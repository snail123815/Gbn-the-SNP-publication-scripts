#!/bin/bash
NCPU=16
RAWREADS=cleanreads # dir in source dir
GENOME=M145.fa # file in source dir
ISPE=""

OUTDIRSAM=alignmentSAM
OUTDIRBAM=alignmentBAM # this will be actually copied back to SOURCEDIR. 

# ===========================================================================================
echo "[$SHELL] #### Starting Job"

# Get start time
START_TIME=$(date +%s)
echo "[$SHELL] Started at: $(date)"

# Enter the designed environment
echo "[$SHELL] Activate shortReads conda env"
source ~/.bashrc
micromamba activate shortReads

echo "[$SHELL] run alignment"
CMD="python align_to_genome.py -r $RAWREADS -g $GENOME -o $OUTDIRSAM -p $NCPU $ISPE"
echo "[$SHELL] $CMD"
eval "$CMD"

echo "[$SHELL] Done alignment"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"

echo "[$SHELL] run indexing"
CMD="python processSam2Bam.py -p $OUTDIRSAM -o $OUTDIRBAM -t $NCPU"
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

