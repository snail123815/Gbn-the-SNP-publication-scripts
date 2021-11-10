# ChIPSeq data processing

Done in a linux server. Other environems should also work if dependencies are met.

Download raw data and put the four `*.fq.gz` files in `cleanreads` dir.

Rename these files as below:

```tree
ChIP_dataAnalysis/cleanreads/
├── ChIP_25h.fq.gz
├── ChIP_48h.fq.gz
├── gDNA_25h.fq.gz
└── gDNA_48h.fq.gz
```

## 1. Align reads to genome

Create a conda environment using following `.yml` specs:

```yml
name: shortReads
channels:
        - conda-forge
        - bioconda
dependencies:
        - biopython
        - bowtie2
        - cairo
        - ipython
        - jedi
        - macs2
        - meme
        - numpy
        - openpyxl
        - openssl
        - pandas
        - pip
        - python
        - samtools=1.9
        - subread
        - spades=3.13.0
        - blast
        - pilon
        - unicycler=0.4.8
```

Adjust `CPU=10`, `BAMSORTMEM=10`, and line 21-22 for environment activation in script `run_align_index_calCoverage.sh` to fit your condition.  
Run alignment:

```sh
./run_align_index_calCoverage.sh
```

## Call peaks using MACS2

Run:

```sh
python macs2_peak_calling.py -b alignmentBAM/ -o peakCalling -c compare.tsv
```