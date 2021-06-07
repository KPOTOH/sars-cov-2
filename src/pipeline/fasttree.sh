#!/bin/bash
#PBS -q eternity ## big eternity mem1t mem512g
#PBS -d .
#PBS -l walltime=200:00:00,mem=47G,nodes=1:ppn=24

date
~/fasttree/FastTreeMP -nt -fastest gisaid_hcov-19_2020_05_20_18.bwasw.sam.fasta > gisaid_hcov-19_2020_05_20_18.bwasw.sam.tre
date
