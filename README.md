# RefineAlign v2.1
A shell script for generating  alignments of mitogenome sequences from a fasta file.

RefineAlign:
- Provides a missingness filter to remove sequences with Ns above a set proportion
- Corrects reverse complement issues using mafft
- Ensures that all sequences are on the same alignment coordinates
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
	- Update the ml commands below for your system (works on Dardel 2026-08)
	- Update the SCRIPT_DIRECTORY path to the location of /python_R_scripts

## Requires:
- A fasta file of the sequences to be aligned (INPUT_FASTA)

## To run:
	./RefineAlign.sh INPUT_FASTA REFERENCE_ID MISSING_FILTER_PROPORTION
	# INPUT_FASTA: File name including full path (with .fasta extension). Needs to include the REFERENCE_ID sequence
	# REFERENCE_ID: Header name of the sequence in INPUT_FASTA for which the alignment coordinates are to be based on
	# MISSING_FILTER_PROPORTION: A value between 0.0 (removes sequences with any Ns) and 1.0 (keeps sequence with all Ns). A value of 0.6 will remove sequences with >60% Ns. (Ns = missingness)

## Generates:
- A final alignment of all sequences (INPUT_FASTA.aln.coordinates_fixed.aln.missingfilter.refined.fasta)
- A final alignment of all sequences with superfluous gaps and N positions removed (INPUT_FASTA.aln.curated.fasta)
- A collapsed alignment of haplotypes (INPUT_FASTA.aln.curated.haplotypes.fasta)
- A sequence to haplotype lookup (INPUT_FASTA.aln.curated.haplotype_mapping.tsv)
	
## Citation:
RefineAlign first appears in Sharif et al. 2026, so for now please cite:
	
	Sharif, B. et al. (2026) DNAharvester: A Nextflow Pipeline for Analysing Highly Degraded DNA from Ancient and Historical Specimens. biorXiv 2026.04.20.719564
	https://doi.org/10.64898/2026.04.20.719564

The missingness filter script (fasta_nomissing.py) is modified from the priginal by Pontus Skoglund (https://github.com/pontussk/fasta_nomissing.py) that is distributed under a GNU General Public License v3.0. The modification treats IUPAC ambiguity codes and asterisks as missing data
