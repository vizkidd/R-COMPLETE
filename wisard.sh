#!/bin/bash
#$ -N "wisard"
#$ -cwd
#$ -V
#$ -l data
#$ -o "logs/wisard.o" 
#$ -e "logs/wisard.e"

Rscript wisard.R files/oneway/set.tmp files/all2all_final/
