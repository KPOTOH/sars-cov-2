#!/bin/bash
#PBS -q eternity ## big eternity mem1t mem512g
#PBS -d .
#PBS -l walltime=200:00:00,mem=47G,nodes=1:ppn=24

date
#~/ALIGNMENT/mash sketch -p 24 -l gisaid_hcov-19_2020_05_20_18_andRef.lst -o gisaid_hcov-19_2020_05_20_18_andRef
~/ALIGNMENT/mash dist -p 24 SARS_CoV_2RefSeq.msh gisaid_hcov-19_2020_05_20_18_andRef.msh > gisaid_hcov-19_2020_05_20_18_andRef.distances.tab
date
