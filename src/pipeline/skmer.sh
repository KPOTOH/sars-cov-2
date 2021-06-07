#!/bin/bash
#PBS -q eternity ## big eternity mem1t mem512g
#PBS -d .
#PBS -l walltime=200:00:00,mem=47G,nodes=1:ppn=24

date
~/anaconda2/bin/skmer reference -p 12 sequences
date
