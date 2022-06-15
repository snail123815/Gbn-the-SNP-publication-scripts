# Gbn-the-SNP-publication-scripts

Data analysis scripts in the article:

[**Du, C.**](https://orcid.org/0000-0003-3447-5293), [Willemse, J.](https://www.universiteitleiden.nl/en/staffmembers/joost-willemse#tab-2), [Erkelens, A.M.](https://orcid.org/0000-0001-9369-5814), [Carrion, V.J.](https://orcid.org/0000-0002-4093-0355), [Dame, R.T.](https://orcid.org/0000-0001-9863-1692), and [Wezel, G.P.v.](https://orcid.org/0000-0003-0341-1561) (2022) System-wide analysis of the GATC-binding nucleoid-associated protein Gbn and its impact on Streptomyces development. mSystems: e00061-00022.

doi: [10.1128/msystems.00061-22](https://doi.org/10.1128/msystems.00061-22)

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
