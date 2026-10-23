#!/usr/bin/env Rscript

# collapse2Haplotypes (v3)
# P.D. Heintzman (with assistance from chatGPT)
# 20261023

# Collapses compatible haplotypes.
#
# If substringCollapse == TRUE:
#   - Sequences are processed in descending order of unambiguous nucleotide count (A, C, G, T).
#   - IUPAC ambiguity codes are interpreted as sets of possible bases.
#   - Two sequences are compatible if their possible nucleotide sets overlap at every position where both have a recognized IUPAC code.
#   - Gaps and unrecognized characters are ignored.
#   - Compatible sequences are assigned to an existing representative.
#
# If substringCollapse == FALSE:
#   - Only exactly identical aligned sequences are collapsed.
#
# Original FASTA identifiers are preserved in the mapping file.
# Duplicate identifiers are made unique internally.

# Started with https://github.com/legalLab/protocols-scripts/blob/master/scripts/hapCollapse.R
# Completely re-written by P.D. Heintzman to:
#	1. Removed the cleaning option 
#	2. Record which sequences are considered duplicates and which sequence they collapse to

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: collapse2Haplotypes.R <input_fasta> [--substringCollapse]")
}

input_fasta <- args[1]
substringCollapse <- "--substringCollapse" %in% args

library(ape)

collapse2Haplotypes <- function(
    data,
    substringCollapse = FALSE,
    original_names = rownames(data)
) {

    if (!inherits(data, "DNAbin")) {
        stop("Input must be a DNAbin object")
    }

    seq_names <- rownames(data)
    seq_strings <- apply(as.character(data), 1, paste0, collapse = "")

    if (substringCollapse) {

        # Only A, C, G and T count as unambiguous bases when selecting representatives.
        # IUPAC ambiguity codes are interpreted explicitly during compatibility testing.
        # Gaps and unrecognized characters are ignored during compatibility testing.
        informative_bases <- c("A", "C", "G", "T")

        seq_strings <- toupper(seq_strings)

        # Count informative bases in each sequence.
        informative_counts <- vapply(
            strsplit(seq_strings, "", fixed = TRUE),
            function(s) sum(s %in% informative_bases),
            integer(1)
        )

        # Process the most informative sequences first.
        ord <- order(informative_counts, decreasing = TRUE)

        # Define the nucleotides represented by each DNA IUPAC code.
        iupac_bases <- list(
            A = c("A"),
            C = c("C"),
            G = c("G"),
            T = c("T"),
            R = c("A", "G"),
            Y = c("C", "T"),
            S = c("G", "C"),
            W = c("A", "T"),
            K = c("G", "T"),
            M = c("A", "C"),
            B = c("C", "G", "T"),
            D = c("A", "G", "T"),
            H = c("A", "C", "T"),
            V = c("A", "C", "G"),
            N = c("A", "C", "G", "T")
        )
        
        # Two sequences are compatible if their possible nucleotide sets
        # overlap at every position where both have a recognised IUPAC code.
        # Gaps and other unrecognised characters are treated as missing.
        compatible <- function(s1, s2) {
        
            b1 <- strsplit(toupper(s1), "", fixed = TRUE)[[1]]
            b2 <- strsplit(toupper(s2), "", fixed = TRUE)[[1]]
        
            for (k in seq_along(b1)) {
        
                bases1 <- iupac_bases[[b1[k]]]
                bases2 <- iupac_bases[[b2[k]]]
        
                # Ignore positions containing gaps or unrecognised characters.
                if (is.null(bases1) || is.null(bases2)) {
                    next
                }
        
                # If the possible nucleotide sets do not overlap,
                # the sequences are incompatible.
                if (length(intersect(bases1, bases2)) == 0) {
                    return(FALSE)
                }
            }
        
            TRUE
        }

        # Indices of representative sequences, in the original data.
        representatives <- integer(0)

        # Original sequence index -> representative sequence index.
        representative_for <- integer(length(seq_strings))

        for (i in ord) {

            if (length(representatives) > 0) {

                matches <- representatives[
                    vapply(
                        representatives,
                        function(j) compatible(
                            seq_strings[i],
                            seq_strings[j]
                        ),
                        logical(1)
                    )
                ]

            } else {
                matches <- integer(0)
            }

            if (length(matches) > 0) {

                # Collapse into the first compatible representative.
                representative_for[i] <- matches[1]

            } else {

                # No compatible representative: retain this sequence.
                representatives <- c(representatives, i)
                representative_for[i] <- i
            }
        }

        # Keep representative sequences in their original input order.
        representatives <- sort(representatives)

        collapsed <- data[representatives, , drop = FALSE]
        rownames(collapsed) <- seq_names[representatives]

        # Map each input record to its representative.
        mapping <- data.frame(
            original_index = seq_along(seq_names),
            original = original_names,
            unique_id = seq_names,
            representative_input_index = representative_for,
            representative = seq_names[representative_for],
            stringsAsFactors = FALSE
        )

    } else {

        keep <- !duplicated(seq_strings)

        collapsed <- data[keep, , drop = FALSE]
        rownames(collapsed) <- seq_names[keep]

        mapping <- data.frame(
            original_index = seq_along(seq_names),
            original = original_names,
            unique_id = seq_names,
            representative_input_index = match(seq_strings, seq_strings[keep]),
            representative = seq_names[
                match(seq_strings, seq_strings[keep])
            ],
            stringsAsFactors = FALSE
        )
    }

    list(collapsed = collapsed, mapping = mapping)
}

# Run

# Read FASTA records as plain text first, so unequal sequence lengths can be detected before ape::read.dna() attempts to parse the alignment.

fasta_lines <- readLines(input_fasta, warn = FALSE)

header_positions <- which(startsWith(fasta_lines, ">"))

if (length(header_positions) == 0L) {
    stop("Input FASTA contains no sequence headers")
}

if (header_positions[1] != 1L) {
    stop("Invalid FASTA: content appears before the first sequence header")
}

# Extract identifiers from FASTA headers.
original_names <- sub(
    "^>(\\S+).*$",
    "\\1",
    fasta_lines[header_positions]
)

# Determine the sequence lines belonging to each record.
sequence_lengths <- integer(length(header_positions))

for (i in seq_along(header_positions)) {

    start <- header_positions[i] + 1L

    end <- if (i < length(header_positions)) {
        header_positions[i + 1L] - 1L
    } else {
        length(fasta_lines)
    }

    if (start <= end) {
        sequence_lines <- fasta_lines[start:end]
        sequence_lines <- sequence_lines[nzchar(trimws(sequence_lines))]
        sequence_lengths[i] <- sum(nchar(gsub("[[:space:]]", "", sequence_lines)))
    } else {
        sequence_lengths[i] <- 0L
    }
}

if (any(sequence_lengths == 0L)) {
    stop("Input FASTA contains an empty sequence")
}

if (length(unique(sequence_lengths)) != 1L) {
    stop("Input sequences have unequal lengths; provide an aligned FASTA file")
}

# Read the validated alignment.
dna <- read.dna(input_fasta, format = "fasta")

if (is.null(dim(dna)) || nrow(dna) != length(original_names)) {
    stop("Could not read all FASTA sequences correctly")
}

# Make identifiers unique internally.
unique_names <- make.unique(original_names, sep = "_")
rownames(dna) <- unique_names

res <- collapse2Haplotypes(
    dna,
    substringCollapse = substringCollapse,
    original_names = original_names
)

base <- sub("\\.fasta$", "", input_fasta)

write.dna(
    res$collapsed,
    file = paste0(base, ".haplotypes.fasta"),
    format = "fasta",
    nbcol = -1,
    colsep = ""
)

write.table(
    res$mapping,
    file = paste0(base, ".haplotype_mapping.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
