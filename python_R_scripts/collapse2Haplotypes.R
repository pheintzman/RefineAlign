#!/usr/bin/env Rscript

# collapse2Haplotypes (v2)
# P.D. Heintzman (with assistance from chatGPT)
# 20260118

# script collapses haplotypes (i.e. removes duplicate sequences)
# if `substringCollapse==TRUE` the function will consider shorter but identical sequences as the same haplotype and collapse them, returning the longest sequence
# if `substringCollapse==FALSE` the function will consider shorter but identical sequences as different haplotypes and will keep them.

# Started with https://github.com/legalLab/protocols-scripts/blob/master/scripts/hapCollapse.R
# Completely re-written by P.D. Heintzman to:
#	1. Removed the cleaning option 
#	2. Record which sequences are considered duplicates and which sequence they collapse to
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: collapse2Haplotypes.R <input_fasta> [--substringCollapse]")
}

input_fasta <- args[1]
substringCollapse <- "--substringCollapse" %in% args

library(ape)
library(stringr)

collapse2Haplotypes <- function(data, substringCollapse = FALSE) {

    if (!inherits(data, "DNAbin")) {
        stop("Input must be a DNAbin object")
    }

    seq_names <- rownames(data)
    seq_strings <- apply(as.character(data), 1, paste0, collapse = "")
    seq_lengths <- nchar(seq_strings)

    if (substringCollapse) {

        ord <- order(seq_lengths, decreasing = TRUE)
        seq_strings_ord <- seq_strings[ord]
        data_ord <- data[ord, , drop = FALSE]
        names_ord <- seq_names[ord]

        rep_idx_all <- vapply(
            seq_strings_ord,
            function(s) {
                which(str_detect(seq_strings_ord, fixed(s)))[1]
            },
            integer(1)
        )

        rep_idx <- unique(rep_idx_all)

        collapsed <- data_ord[rep_idx, , drop = FALSE]
        rownames(collapsed) <- names_ord[rep_idx]

        mapping <- data.frame(
            original = names_ord,
            representative = names_ord[rep_idx_all],
            stringsAsFactors = FALSE
        )

    } else {

        keep <- !duplicated(seq_strings)

        collapsed <- data[keep, , drop = FALSE]
        rownames(collapsed) <- seq_names[keep]

        mapping <- data.frame(
            original = seq_names,
            representative = seq_names[match(seq_strings, seq_strings[keep])],
            stringsAsFactors = FALSE
        )
    }

    list(collapsed = collapsed, mapping = mapping)
}

# Run

dna <- read.dna(input_fasta, format = "fasta")

res <- collapse2Haplotypes(
    dna,
    substringCollapse = substringCollapse
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
