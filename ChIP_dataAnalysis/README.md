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
        - meme=5.1
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

## Get sequences of peaks

Run:

```sh
python pull_peak_sequences.py -g M145.fa -f 'peakCalling/t25/t25_peaks.xls' 'peakCalling/t48/t48_peaks.xls' --filter summit -o ./
```

## Run meme-chip

Dump files `t25_peaks.fa`, `t48_peaks.fa` to [MEME Suit - MEME_ChIP](https://meme-suite.org/meme/tools/meme-chip) is the best way.  
If you really want to run it locally: (No gaurentee)

Database in this folder `collectf.meme` was from MEME version 4

Run this the first time to prepare for meme run:

```sh
python meme-chip.py
```

The program will end without running meme.

Run following command manually because meme-chip looked for `meme` in a wrong location.

```sh
meme MEME_CHIP_Motif/t25_peaks_anr/seqs-centered -oc MEME_CHIP_Motif/t25_peaks_anr/meme_out -mod anr -nmotifs 3 -minw 6 -maxw 30 -bfile MEME_CHIP_Motif/t25_peaks_anr/background -dna -p 4 -revcomp
meme MEME_CHIP_Motif/t48_peaks_anr/seqs-centered -oc MEME_CHIP_Motif/t48_peaks_anr/meme_out -mod anr -nmotifs 3 -minw 6 -maxw 30 -bfile MEME_CHIP_Motif/t48_peaks_anr/background -dna -p 4 -revcomp
```

Then run this again:

```sh
python meme-chip.py
```

You should get more completed result by now.