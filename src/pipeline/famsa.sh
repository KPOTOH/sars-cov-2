#!/bin/bash
#PBS -q eternity ## big eternity mem1t mem512g
#PBS -d .
#PBS -l walltime=200:00:00,mem=47G,nodes=1:ppn=24

date
module load gcc/gcc-6.3 ;
~/ALIGNMENT/FAMSA/famsa -v -t 24 gisaid_hcov-19_2020_05_20_18_andRef.fasta gisaid_hcov-19_2020_05_20_18_andRef.famsa.fasta
date
