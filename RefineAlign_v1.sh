#!/bin/bash -l
set -e

#############################
#                           #
#       RefineAlign         #
#       P.D. Heintzman      #
#       20260418            #
#       version 1.0         #
#                           #
#############################


# Requires:
	# A concatenated fasta file of the sequences to be aligned (INPUT_FASTA)
	# A fasta file of the sequence that the alignment coordinates are to be based on (REFERENCE_ID)
	# The REFERENCE_ID fasta needs to be in INPUT_FASTA with a header name identical to the REFERENCE_ID file name

# Generates:
	# A final alignment of all sequences
	# A collapsed alignment of haplotypes
	# A sequence to haplotype lookup

# Before starting:
	# Dependencies: mafft, python3, and R (tested with R v4.4.2)
	# Update the ml commands below for your system
	# Update the SCRIPT_DIRECTORY path to the locations of the python and R scripts
	

# Variables
INPUT_FASTA=$1
	# File name including path (with .fasta extension)
REFERENCE_ID=$2
	# Just the file name (without .fasta extension)


# Scripts and tools
SCRIPT_DIRECTORY=tools/RefineAlign/python_R_scripts
	# For python and R scripts


# Generating an aligned reference panel of haplotypes from a FASTA file
INPUT_FASTA=${INPUT_FASTA%.fasta}


echo ">>> Running Reference panel generation..."
# Create an initial alignment with mafft
# Also checks for and aligns reverse complemented sequences
module load mafft/7.526
mafft \
	--auto \
	--adjustdirection \
	${INPUT_FASTA}.fasta \
	> ${INPUT_FASTA}.aln.fasta


# Place everything on the same alignment coordinates as REFERENCE_ID

# If upstream gaps in REFERENCE_ID, then alignment in this region is moved to the end
# Then all gaps are removed from the alignment
module purge
ml bioinfo-tools biopython/1.68-py3 PDCOLD/23.12 ete/3.1.3-cpeGNU-23.12
python3 ${SCRIPT_DIRECTORY}/rotate_start_gaps.py \
	${INPUT_FASTA}.aln.fasta \
	${REFERENCE_ID} \
	${INPUT_FASTA}.aln.start_rotated.fasta

# Create an intermediate alignment with mafft
module purge
module load mafft/7.526
mafft \
	--auto \
	${INPUT_FASTA}.aln.start_rotated.fasta \
	> ${INPUT_FASTA}.aln.start_rotated.aln.fasta

# If downstream gaps in REFERENCE_ID, then alignment in this region is moved to the start
# Then all gaps are removed from the alignment
module purge
ml bioinfo-tools biopython/1.68-py3 PDCOLD/23.12 ete/3.1.3-cpeGNU-23.12
python3 ${SCRIPT_DIRECTORY}/rotate_end_gaps.py \
	${INPUT_FASTA}.aln.start_rotated.aln.fasta \
	${REFERENCE_ID} \
	${INPUT_FASTA}.aln.end_rotated.fasta

# Create an intermediate alignment with mafft
module purge
module load mafft/7.526
mafft \
	--auto \
	${INPUT_FASTA}.aln.end_rotated.fasta \
	> ${INPUT_FASTA}.aln.end_rotated.aln.fasta

# Final iteration to ensure alignment starts with REFERENCE_ID
# Then all gaps are removed from the alignment
module purge
ml bioinfo-tools biopython/1.68-py3 PDCOLD/23.12 ete/3.1.3-cpeGNU-23.12
python3 ${SCRIPT_DIRECTORY}/rotate_start_gaps.py \
	${INPUT_FASTA}.aln.end_rotated.aln.fasta \
	${REFERENCE_ID} \
	${INPUT_FASTA}.aln.coordinates_fixed.fasta

# Create the final alignment with mafft
module purge
module load mafft/7.526
mafft \
	--auto \
	${INPUT_FASTA}.aln.coordinates_fixed.fasta \
	> ${INPUT_FASTA}.aln.coordinates_fixed.aln.fasta


# Remove positions from the alignment that only consist of gaps and Ns
# Will not remove positions with an N call in the REFERENCE_ID
# Optional: remove all positions in the alignment where the REFERENCE_ID has a gap
	# Add --remove-reference-gaps to enable this option
module purge
ml bioinfo-tools biopython/1.68-py3 PDCOLD/23.12 ete/3.1.3-cpeGNU-23.12
python3 ${SCRIPT_DIRECTORY}/remove_all_gap_and_N_columns.py \
	${INPUT_FASTA}.aln.coordinates_fixed.aln.fasta \
	${REFERENCE_ID} \
	${INPUT_FASTA}.aln.curated.fasta \
	--remove-reference-gaps


# Collapse identical sequences to haplotypes
# Flag to collapse substrings. This means that if two sequences are identical but differ only by Ns or gaps then they are collapsed (if flag is not specified then defaults to FALSE)
module purge
module load PDCOLD/24.11 R/4.4.2-cpeGNU-24.11
${SCRIPT_DIRECTORY}/collapse2Haplotypes.R \
	${INPUT_FASTA}.aln.curated.fasta \
	--substringCollapse
echo ">>>  Reference panel generation complete"


# Remove intermediate files
rm ${INPUT_FASTA}.aln.fasta
rm ${INPUT_FASTA}.aln.start_rotated*.fasta
rm ${INPUT_FASTA}.aln.end_rotated*.fasta
rm ${INPUT_FASTA}.aln.coordinates_fixed.fasta
