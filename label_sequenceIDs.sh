#!/bin/bash
#1 - org name
#2 - gene name
#3 - GTF file path (for extracting the 'real' gene annotation)
#4 - Output FASTA file
#5 - Input File
#6 - OPTIONAL Ortholog Clusters from OrthoDB(eg, files/genes/xenopus_tropicalis/odb.final_map)

function check_param() {
  local STDIN=$(cat)
  local tmp_name=$(echo $STDIN | awk -F'==' '{print $1}')
  local tmp_var=$(echo $STDIN | awk -F'==' '{print $2}')
  if [[ -z $tmp_var ]]; then
    echo "Parameter Error : $tmp_name is empty"
    exit 1
  fi
  echo $tmp_var
}

export -f check_param

if [[ ! -s $5 ]]; then
	exit 1
fi

rm -f $4

transcript_delimiter=$(grep -i  -w "transcript_delimiter" parameters.txt | check_param)
seqID_delimiter=$(grep -i -w "seqID_delimiter" parameters.txt | check_param)

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
	anno_name=$2 #$(echo $(grep -w $transcript $3 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq))
	#echo $transcript_id $transcript $anno_name
#fi

if [[ -s $6 && ! -z $6 ]] ; then
	ortho_cluster=$(grep -w $2 $6 | awk -F'\t' '{if (length(c) == 0){c=$1;}else{c=c","$1;}}END{print c}')
fi
#if [[ -z "$anno_name" ]] ; then
#	anno_name=$(echo "$2"--)
#fi
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
rm -f $5
exit 0
