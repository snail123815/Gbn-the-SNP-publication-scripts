# Gbn-the-SNP-publication-scripts

Draft now.

## Dependencies

If you do not want to run the code, please ignore this section.

Recommanded environment:

Linux, WSL2

Partially compatible with windows (lack bedtools) and Mac (`lfcShrink`).

### Oleveler

v1.0 - https://gitlab.services.universiteitleiden.nl/ibl-bioinformatic/oleveler

Commit used:

https://gitlab.services.universiteitleiden.nl/ibl-bioinformatic/oleveler/-/tree/8b848a6e9406d9e82a896bca5934661d71629dc6

### `ChIP_Expression.py`

- bedtools
- pip:
	- pybedtools
	- biopython
	- bcbio-gff
	- matplotlib
	- pandas
	- scipy
	- scikit-learn
	- jupyter lab
	- jupyterthemes

## Notes

Re-analysis of transcriptomics data on a Mac system is not possile for unknown incompatibility  
of `lfcShrink` function in DESeq2 package. You may try older versions of R instead (<4.1.2).
