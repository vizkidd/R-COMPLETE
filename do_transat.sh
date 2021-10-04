#!/bin/bash
#$ -cwd
#$ -V
#$ -l data
#$ -l h_vmem=128G
#$ -l m_mem_free=128G
#$ -o "logs/transat.o" 
#$ -e "logs/transat.e"

TRANSAT_PATH=$(grep -i -w "transat_path" parameters.txt | awk -F'=' '{print $2}') 

$TRANSAT_PATH -v -d --mantissa --indep_pvalues --subhelices -a $1 -t $2 -c $3 -o $4
