#!/bin/bash
#1 - org name
#2 - gene name
#3 - Output FASTA file
#4 - Input File
#5 - OPTIONAL Ortholog Clusters from OrthoDB(eg, files/genes/xenopus_tropicalis/odb.final_map)

###ENTRYPOINT 

#Label Transcripts using the format - Transcript_id $seqID_delimiter Gene_name $seqID_delimiter Organism $seqID_delimiter ODB_cluster

if [[ ! -s $4 ]]; then
	exit 1
fi

source functions.sh

rm -f $3

transcript_delimiter=$(grep -i  -w "transcript_delimiter" parameters.txt | check_param)
seqID_delimiter=$(grep -i -w "seqID_delimiter" parameters.txt | check_param)

transcript_id=""
seq=""
all_seqs=()

IFS=">" read -ra all_rec <<< $(cat $4)

for rec in "${all_rec[@]}"; do
	read -r transcript_id seq <<< $(echo "$rec")
	transcript=$(echo $transcript_id | awk -F"[$(echo $transcript_delimiter)]" '{print $1}' | sed 's/>//g')

	if [[ -s $5 && ! -z $5 ]] ; then
		ortho_cluster=$(grep -w $2 $5 | awk -F'\t' '{if (length(c) == 0){c=$1;}else{c=c","$1;}}END{print c}')
	fi
	if [[ -z "$ortho_cluster" ]]; then
		ortho_cluster="ungrouped"
	fi

	if [[ ! -z "$seq" && ! -z "$transcript_id" ]] ; then
		seq_ID=$(printf "%s%s%s%s%s%s%s" $transcript_id $seqID_delimiter $2 $seqID_delimiter $1 $seqID_delimiter $ortho_cluster )
		printf ">%s\n%s\n" $seq_ID $seq >> $3
		#printf "%s\t%s\t%s\n" $1 $2 $seq_ID >> files/genes/$1/transcripts.metadata
	fi

done

rm -f $4

exit 0
