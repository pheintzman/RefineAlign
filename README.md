# RefineAlign v2.4
### A shell script for generating an alignment of mitogenome sequences from a fasta file.

RefineAlign:
- Corrects reverse complement issues using mafft
- Ensures that all sequences are on the same alignment coordinates
- Provides a missingness filter to remove sequences with Ns above a set proportion
- Refines the mafft alignment using muscle
- Removes superfluous gaps and N positions
- Collapses identical sequences to haplotypes

Made for use on the Swedish NAISS Dardel cluster but easily modifiable for other local or HPC systems.

## Download and setup:
	git clone https://github.com/pheintzman/RefineAlign.git
	cd RefineAlign
	chmod 755 RefineAlign.sh
	chmod 755 ./python_R_scripts/*

## Before starting:
- Dependencies: python, mafft, python3, and R (tested with R v4.4.2)
- In the RefineAlign.sh script:
	- Update the ml commands below for your system (works on Dardel 2026-10)
	- Update the SCRIPT_DIRECTORY path to the location of /python_R_scripts

## Requires:
- A fasta file of the sequences to be aligned (INPUT_FASTA)

## To run:
	./RefineAlign.sh INPUT_FASTA REFERENCE_ID MISSING_FILTER_PROPORTION [Two optional flags: --remove-reference-gaps --substringCollapse]
	# INPUT_FASTA: File name including full path (with .fasta extension). Needs to include the REFERENCE_ID sequence
	# REFERENCE_ID: Header name of the sequence in INPUT_FASTA for which the alignment coordinates are to be based on
	# MISSING_FILTER_PROPORTION: A value between 0.0 (removes sequences with any missing data) and 1.0 (keeps sequences with all missing data). A value of 0.6 will retain sequences with <60% missing data (i.e. >40% completeness)
	# --remove-reference-gaps: Remove all positions in the alignment where the REFERENCE_ID has a gap. This is useful for co-ordinate based analyses, but should not be applied if the alignment is to be used for phylogenetic or network analyses (only applied if flag is specified)
	# --substringCollapse: Flag to collapse substrings. This means that if two sequences are identical but differ only by Ns or gaps then they are collapsed (only applied if flag is specified)

## Generates:
- A final alignment of all sequences (INPUT_FASTA.aln.coordinates_fixed.aln.missingfilter.refined.fasta)
- A final alignment of all sequences with superfluous gaps and N positions removed (INPUT_FASTA.aln.curated.fasta)
- A collapsed alignment of haplotypes (INPUT_FASTA.aln.curated.haplotypes.fasta)
- A sequence to haplotype lookup (INPUT_FASTA.aln.curated.haplotype_mapping.tsv)
	
## Citation:
RefineAlign first appears in Sharif et al. 2026, so for now please cite:
	
	Sharif, B. et al. (2026) DNAharvester: A Nextflow Pipeline for Analysing Highly Degraded DNA from Ancient and Historical Specimens. Bioinformatics btag782
	https://doi.org/10.1093/bioinformatics/btag782

The missingness filter script (fasta_nomissing_v2.py) is modified from the original by Pontus Skoglund (https://github.com/pontussk/fasta_nomissing.py) that is distributed under a GNU General Public License v3.0.
In the original version, only Ns were recognised as missing data whereas the modification also recognises IUPAC ambiguity codes and asterisks as missing data.
This version is also designed to be used on an alignment by recognising the gap character (-). Gaps are not considered in missingness calculations.
If sequences are of different lengths, then the script will exit with an error. This is in contrast to the original script, whereby sequences were silently truncated to the shortest. Applying the script to a multiple sequence alignment should not encounter this issue.
