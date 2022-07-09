#!/bin/bash
#PBS -d .
#PBS -l walltime=500:00:00,mem=100gb,nodes=node34,ncpus=32

# mem=80gb,nodes=node34,ncpus=33

SCRIPT1=./pipeline/1.bwa.sh
SCRIPT2=./pipeline/2.longest-alignments.pl
SCRIPT3=./pipeline/3.sam2fasta.py
SCRIPT4=./pipeline/4.fasttree.sh
# SCRIPT5=./pipeline/resolve_polytomies_in_ete3.py
# SCRIPT6=./pipeline/6.anceors-reconstruction.sh

pathToFasttreemp=/mnt/lustre/genkvg/fasttree/FastTreeMP
pathToPrank=/mnt/lustre/genkvg/ALIGNMENT/prank/prank
BWAPATH=./tools/bwa/bwa
PATH_TO_SAMTOOLS=./tools/samtools-1.5/bin/samtools
THREADS=32

REF=./data/reference/covid_ref.fasta
WORKING_PATH=./data/
SAMPLE=gisaid_hcov-19_2020-11-22.fasta

START_FASTA="$WORKING_PATH$SAMPLE"

SAMFILE="$START_FASTA.sam"
PURSAMFILE="$START_FASTA.purified.sam"
MULTALIGH="${WORKING_PATH}mulal_$SAMPLE"
# binarySimplifiedTree=./data/multiple_alignment_gisaid.fasta.tre-simple.pruned.resolved


source /home/kpotoh/tools/python_env/bin/activate
echo "$(date) INFO start calculation" | telegram-send --stdin

$SCRIPT1 $BWAPATH $PATH_TO_SAMTOOLS $REF $START_FASTA
echo "$(date) INFO samfile created" | telegram-send --stdin

$SCRIPT2 $SAMFILE > $PURSAMFILE
echo "$(date) INFO samfile purified" | telegram-send --stdin

$SCRIPT3 $REF $PURSAMFILE $MULTALIGH
echo "$(date) INFO mulal done" | telegram-send --stdin

$SCRIPT4 $pathToFasttreemp $MULTALIGH
echo "$(date) INFO tree builded" | telegram-send --stdin

# $SCRIPT5 "$MULTALIGH.tre-simple"
# echo "$(date) INFO polytomies resolving is done" | telegram-send --stdin