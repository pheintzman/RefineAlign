#!/bin/bash -l
set -e

#############################
#                           #
#       RefineAlign         #
#       P.D. Heintzman      #
#       20260809            #
#       version 2.1         #
#                           #
#############################


# Requires:
	# A concatenated fasta file of the sequences to be aligned (INPUT_FASTA)
	# A fasta file of the sequence that the alignment coordinates are to be based on (REFERENCE_ID)
	# The REFERENCE_ID fasta needs to be in INPUT_FASTA with a header name identical to the REFERENCE_ID file name

# Generates:
	# A final alignment of all sequences (INPUT_FASTA.aln.coordinates_fixed.aln.fasta)
	# A final alignment of all sequences with superfluous gaps and N positions removed (INPUT_FASTA.aln.curated.fasta)
	# A collapsed alignment of haplotypes
	# A sequence to haplotype lookup

# Before starting:
	# Dependencies: python, mafft, python3, and R (tested with R v4.4.2)
	# Update the ml commands below for your system (works on Dardel 2026-08)
	# Update the SCRIPT_DIRECTORY path to the locations of the python and R scripts
	

# Variables
INPUT_FASTA=$1
	# File name including path (with .fasta extension)
MISSING_FILTER_PROPORTION=$2
	# Just the file name (without .fasta extension)
REFERENCE_ID=$3
	# Just the file name (without .fasta extension)


# Scripts and tools
SCRIPT_DIRECTORY=tools/RefineAlign/python_R_scripts
	# For python and R scripts


# Generating an aligned reference panel of haplotypes from a FASTA file
INPUT_FASTA=${INPUT_FASTA%.fasta}


# Run an initial missingness filter to remove sequences with a maximum proportion of Ns given in $2
# Script by Pontus Skoglund (https://github.com/pontussk/fasta_nomissing.py) distributed under a GNU General Public License v3.0
module purge
python ${SCRIPT_DIRECTORY}/fasta_nomissing_v2.py \
	--fastafile=${INPUT_FASTA}.fasta \
	--maxmissing_ind=${MISSING_FILTER_PROPORTION} \
	> ${INPUT_FASTA}.missingfiltered.fasta


echo ">>> Running Reference panel generation..."
# Create an initial alignment with mafft
# Also checks for and aligns reverse complemented sequences
module load mafft/7.526
mafft \
	--auto \
	--adjustdirection \
	${INPUT_FASTA}.missingfiltered.fasta \
	> ${INPUT_FASTA}.missingfiltered.aln.fasta


# Place everything on the same alignment coordinates as REFERENCE_ID

# If upstream gaps in REFERENCE_ID, then alignment in this region is moved to the end
# Then all gaps are removed from the alignment
module purge
ml bioinfo-tools biopython/1.68-py3 PDCOLD/23.12 ete/3.1.3-cpeGNU-23.12
python3 ${SCRIPT_DIRECTORY}/rotate_start_gaps.py \
	${INPUT_FASTA}.missingfiltered.aln.fasta \
	${REFERENCE_ID} \
	${INPUT_FASTA}.missingfiltered.aln.start_rotated.fasta

# Create an intermediate alignment with mafft
module purge
module load mafft/7.526
mafft \
	--auto \
	${INPUT_FASTA}.missingfiltered.aln.start_rotated.fasta \
	> ${INPUT_FASTA}.missingfiltered.aln.start_rotated.aln.fasta

# If downstream gaps in REFERENCE_ID, then alignment in this region is moved to the start
# Then all gaps are removed from the alignment
module purge
ml bioinfo-tools biopython/1.68-py3 PDCOLD/23.12 ete/3.1.3-cpeGNU-23.12
python3 ${SCRIPT_DIRECTORY}/rotate_end_gaps.py \
	${INPUT_FASTA}.missingfiltered.aln.start_rotated.aln.fasta \
	${REFERENCE_ID} \
	${INPUT_FASTA}.missingfiltered.aln.end_rotated.fasta

# Create an intermediate alignment with mafft
module purge
module load mafft/7.526
mafft \
	--auto \
	${INPUT_FASTA}.missingfiltered.aln.end_rotated.fasta \
	> ${INPUT_FASTA}.missingfiltered.aln.end_rotated.aln.fasta

# Final iteration to ensure alignment starts with REFERENCE_ID
# Then all gaps are removed from the alignment
module purge
ml bioinfo-tools biopython/1.68-py3 PDCOLD/23.12 ete/3.1.3-cpeGNU-23.12
python3 ${SCRIPT_DIRECTORY}/rotate_start_gaps.py \
	${INPUT_FASTA}.missingfiltered.aln.end_rotated.aln.fasta \
	${REFERENCE_ID} \
	${INPUT_FASTA}.missingfiltered.aln.coordinates_fixed.fasta

# Create the final alignment with mafft
module purge
module load mafft/7.526
mafft \
	--auto \
	${INPUT_FASTA}.missingfiltered.aln.coordinates_fixed.fasta \
	> ${INPUT_FASTA}.missingfiltered.aln.coordinates_fixed.aln.fasta


# Refine the final alignment with muscle
module purge
ml bioinfo-tools muscle/3.8.1551
muscle \
	-refine \
	-in ${INPUT_FASTA}.missingfiltered.aln.coordinates_fixed.aln.fasta \
	-out ${INPUT_FASTA}.missingfiltered.aln.coordinates_fixed.aln.refined.fasta


# Remove positions from the alignment that only consist of gaps and Ns
# Will not remove positions with an N call in the REFERENCE_ID
# Optional: remove all positions in the alignment where the REFERENCE_ID has a gap
	# Add --remove-reference-gaps to enable this option
module purge
ml bioinfo-tools biopython/1.68-py3 PDCOLD/23.12 ete/3.1.3-cpeGNU-23.12
python3 ${SCRIPT_DIRECTORY}/remove_all_gap_and_N_columns.py \
	${INPUT_FASTA}.missingfiltered.aln.coordinates_fixed.aln.refined.fasta \
	${REFERENCE_ID} \
	${INPUT_FASTA}.aln.curated.fasta \
	--remove-reference-gaps


# Collapse identical sequences to haplotypes
# Flag to collapse substrings. This means that if two sequences are identical but 
	# differ only by Ns or gaps then they are collapsed (if flag is not specified then defaults to FALSE)
module purge
module load PDCOLD/24.11 R/4.4.2-cpeGNU-24.11
${SCRIPT_DIRECTORY}/collapse2Haplotypes.R \
	${INPUT_FASTA}.aln.curated.fasta \
	--substringCollapse
echo ">>>  Reference panel generation complete"


# Remove intermediate files
rm ${INPUT_FASTA}.missingfiltered.aln.fasta
rm ${INPUT_FASTA}.missingfiltered.aln.start_rotated*.fasta
rm ${INPUT_FASTA}.missingfiltered.aln.end_rotated*.fasta
rm ${INPUT_FASTA}.missingfiltered.aln.coordinates_fixed.fasta
rm ${INPUT_FASTA}.missingfiltered.aln.coordinates_fixed.aln.fasta
