#!/bin/bash

NCPU=10
BAMSORTMEM=10
RAWREADS=cleanreads # dir in source dir
GENOME=M145.fa # file in source dir
ISPE=""

OUTDIRSAM=alignmentSAM
OUTDIRBAM=alignmentBAM
# ===========================================================================================
echo "[$SHELL] #### Starting Job"

# Get start time
START_TIME=$(date +%s)
echo "[$SHELL] Started at: $(date)"

# Enter the designed environment
echo "[$SHELL] Activate shortReads conda env"

eval "$($MAMBA_EXE shell hook -s bash)"
eval "micromamba activate shortReads"

echo "[$SHELL] run alignment"
CMD="python align_to_genome.py -r $RAWREADS -g $GENOME -o $OUTDIRSAM -p $NCPU $ISPE"
echo "[$SHELL] $CMD"
eval "$CMD"

echo "[$SHELL] Done alignment"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"

echo "[$SHELL] run indexing"
CMD="python processSam2Bam.py -p $OUTDIRSAM -o $OUTDIRBAM -t $NCPU -m $BAMSORTMEM"
echo "[$SHELL] $CMD"
eval "$CMD"

echo "[$SHELL] Done indexing"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"

echo "[$SHELL] Calculating coverage"
CMD="python calculate_coverage.py -p $OUTDIRBAM -o coverage.tsv"
echo "[$SHELL] $CMD"
eval "$CMD"

echo "[$SHELL] Done producing coverage"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"

