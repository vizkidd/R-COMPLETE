#!/bin/bash
#$ -cwd
#$ -V
#$ -l data
#$ -l m_mem_free 8G

# $1 = FASTA PATH (Optional, if not given, FASTA PATH taken from parameters.txt)

source functions.sh

#Check if scripts exist
if [[ ! -s parameters.txt ]] ; then
  >&2 color_FG_Bold $Red "Missing parameters.txt!"
  exit 1
fi

if [ $# -eq 0 ]; then
	FASTA_PATH=$(grep -i -w "fasta_path" parameters.txt | check_param) 
else
	FASTA_PATH=$1
fi

GROUPS_PATH=$(grep -i -w "groups_path" parameters.txt | check_param) 
n_threads=$(grep -i -w "max_concurrent_jobs" parameters.txt | check_param)
seqID_delimiter=$(grep -i -w "seqID_delimiter" parameters.txt | check_param) 

for f_org in $FASTA_PATH/*; do 
	f_org_name=$(basename $f_org)
	grep -r -h ">" $f_org  | awk -F"$seqID_delimiter" '{print $NF}' | awk '{split($0,a,","); for(key in a) print a[key];}' | sort -u > files/genes/$f_org_name/ORG_CLUSTERS
	parallel --max-procs $n_threads " printf '%s\t%s\n' {1} {2}" :::: <(grep -H -f files/genes/$f_org_name/ORG_CLUSTERS -r $FASTA_PATH/$f_org_name/ | awk -F'[:>]' -v s_delim="$seqID_delimiter" '{split($2,a,s_delim); n=split($1,b,"."); print $1"\t"$2"\t"a[5]"\t"b[n]'}) | parallel  --max-procs 1 --colsep '\t' --recend '\n'  "if [[ -s {1} && ! -z {2} && ! -z {1} && ! -z {3} ]] ; then samtools faidx {1}  {2} >> $GROUPS_PATH/{3}.{4} ; fi" 
done

Coerce gtf_stats.csv of all organisms
find files/genes -iname "gtf_stats.csv" -exec sed 1d {} \; > files/gtf_stats.csv
##ID Alignments - GENERATE numeric ids for FASTA IDS (because they are long and downstream analysis have difficulty taking long names) 
#Only indexing the IDs for now because find_orthologs.sh depends on the long FASTA IDs and cannot be shorted until orthologous transcripts are obtained
#!!!!!!!#CHANGE FASTA IDs to numeric IDs (because some programs dont work well with long FASTA IDs) ONLY BEFORE ALIGNMENT!!!!!!!!
index_fastaIDs files/rna_ids.txt $FASTA_PATH