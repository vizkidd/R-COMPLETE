#!/bin/bash
#1 - org name
#2 - gene name
#3 - GTF file path (for extracting the 'real' gene annotation)
#4 - Output FASTA file
#5 - Input File
#6 - OPTIONAL Ortholog Clusters from OrthoDB(eg, files/genes/xenopus_tropicalis/odb.final_map)

rm $4

transcript_delimiter=$(grep -i  -w "transcript_delimiter" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH=$(grep -i -w "python3_path" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH="${PY3_PATH/#\~/$HOME}"
seqID_delimiter=$(grep -i -w "seqID_delimiter" parameters.txt | awk -F'=' '{print $2}')

transcript_id=""
anno_name=""
seq=""
all_seqs=()

##read -r line seq <&0; #<<< $(cat -)
#IFS=">" read -ra all_rec <<< $(cat -)
IFS=">" read -ra all_rec <<< $(cat $5)
for rec in "${all_rec[@]}"; do
  read -r transcript_id seq <<< $(echo "$rec")

#if $(echo $line | grep -q ">") ; then
	
	#transcript=$(echo $transcript_id | awk -F'[||]' '{print $1}' | sed 's/>//g')
	transcript=$(echo $transcript_id | awk -F"[$(echo $transcript_delimiter)]" '{print $1}' | sed 's/>//g')
	#echo $(grep -w $transcript $3 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq)
	anno_name=$(echo $(grep -w $transcript $3 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq))
	#echo $transcript_id $transcript $anno_name
#fi

if [[ -s $6 && ! -z $6 ]] ; then
ortho_cluster=$(grep -w $2 $6 | awk -v OFS="," -F'\t' '{print $1}')
fi
if [[ -z "$anno_name" ]] ; then
	anno_name=$(echo "$2"--)
fi
if [[ -z "$ortho_cluster" ]]; then
	ortho_cluster="ungrouped"
	#ortho_cluster=$(grep $(echo $2 | awk -F'_' '{print $1}' ) $6 | awk -v OFS="\t" -F'\t' '{print $1}')
fi

#echo $seq
#echo $(grep -w $transcript $3 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g')
#echo $transcript_id $2 $anno_name
if [[ ! -z "$seq" && ! -z "$transcript_id" ]] ; then
	#all_seqs+=$(printf ">%s\n%s\n" $seq_ID $seq) 
	seq_ID=$(printf "%s%s%s%s%s%s%s" $transcript_id $seqID_delimiter $2 $seqID_delimiter $1 $seqID_delimiter $anno_name $seqID_delimiter $ortho_cluster )
	printf ">%s\n%s\n" $seq_ID $seq >> $4
	echo $1 $2 $seq_ID >> files/transcripts.metadata
fi

#echo "----------------"
done
#echo ${all_seqs[@]} #| $PY3_PATH filter_sequenceIDs.py $1 $4 $2 $(echo $anno_name | awk '{ gsub(/[[:punct:]]/, "_", $0); print;} ;') ##$org $file #$filename
rm $5
exit 0
