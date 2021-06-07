#!/bin/bash
#PBS -d .
#PBS -l walltime=200:00:00

BWAPATH=./tools/bwa/bwa
REF=./data/covid_ref.fasta
THREADS=32
SAMPLE=./data/filtered_gisaid_hcov-19_2020_11_29_high_qualily.fasta

$BWAPATH index $REF
$BWAPATH bwasw -t $THREADS $REF $SAMPLE > "$SAMPLE.sam.raw"