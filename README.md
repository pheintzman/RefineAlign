# RefineAlign
A shell script for generating  alignments of mitogenome sequences from a concatenated fasta file.
Made for use on the Swedish NAISS Dardel cluster.

RefineAlign:
- Corrects reverse complement issues using mafft
- Ensures that all sequences are on the same alignment coordinates
- Removes superfluous gaps and N positions
- Collapses identical sequences to haplotypes

## Download and setup:
	git clone https://github.com/pheintzman/RefineAlign.git
	cd RefineAlign
	chmod 755 RefineAlign_v1.sh
	chmod 755 ./python_R_scripts/*

## Before starting:
- Dependencies: mafft, python3, and R (tested with R v4.4.2)
- In the RefineAlign_v1.sh script:
	- Update the ml commands for your system
	- Update the SCRIPT_DIRECTORY path to the location of /python_R_scripts

## Requires:
- A concatenated fasta file of the sequences to be aligned (INPUT_FASTA)

## To run:
	./RefineAlign_v1.sh INPUT_FASTA REFERENCE_ID
	# INPUT_FASTA: File name including full path (with .fasta extension)
	# REFERENCE_ID: Header name of the sequence for which the alignment coordinates are to be based on

## Generates:
- A final alignment of all sequences (INPUT_FASTA.aln.coordinates_fixed.aln.fasta)
- A final alignment of all sequences with superfluous gaps and N positions removed (INPUT_FASTA.aln.curated.fasta)
- A collapsed alignment of haplotypes (INPUT_FASTA.aln.curated.haplotypes.fasta)
- A sequence to haplotype lookup (INPUT_FASTA.aln.curated.haplotype_mapping.tsv)
	
## Citation:
RefineAlign first appears in Sharif et al. 2026, so for now please cite:
	
	Sharif, B. et al. (2026) DNAharvester: A Nextflow Pipeline for Analysing Highly Degraded DNA from Ancient and Historical Specimens. biorXiv
