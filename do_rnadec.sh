#!/bin/bash
#$ -cwd
#$ -V
#$ -l data
#$ -l h_vmem=256G
#$ -l m_mem_free=256G
#$ -o "logs/rnadecoder.o" 
#$ -e "logs/rnadecoder.e"

RNADEC_PATH=$(grep -i -w "rnadecoder_path" parameters.txt | awk -F'=' '{print $2}')

$RNADEC_PATH/RNA-decoder $1