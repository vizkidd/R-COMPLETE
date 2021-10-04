#!/bin/bash

transcript_delimiter=$(grep -i  -w "transcript_delimiter" parameters.txt | awk -F'=' '{print $2}')
PY2_PATH=$(grep -i -w "python2_path" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH=$(grep -i -w "python3_path" parameters.txt | awk -F'=' '{print $2}')

transcript_id=""
anno_name=""
seq=""

#read -r line seq <&0; #<<< $(cat -)
IFS=">" read -ra all_rec <<< $(cat -)
for rec in "${all_rec[@]}"; do
  read -r transcript_id seq <<< $(echo "$rec")

#if $(echo $line | grep -q ">") ; then
	
	#transcript=$(echo $transcript_id | awk -F'[||]' '{print $1}' | sed 's/>//g')
	transcript=$(echo $transcript_id | awk -F"[$(echo $transcript_delimiter)]" '{print $1}' | sed 's/>//g')
	#echo $(grep -w $transcript $3 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq)
	anno_name=$(echo $(grep -w $transcript $3 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq))
	#echo $transcript_id $transcript $anno_name
#fi

if [[ "$anno_name" == "" ]] ; then
	anno_name=$(echo "$2"--)
fi

#echo $seq
#echo $(grep -w $transcript $3 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g')
#echo $transcript_id $2 $anno_name
if [[ "$seq" != "" && "$transcript_id" != "" ]] ; then
	echo $1 $2 $anno_name $transcript_id >> files/transcripts.metadata
	printf ">%s\n%s" $transcript_id $seq | $PY2_PATH filter_sequenceIDs.py $1 $5 $2 $(echo $anno_name | awk '{ gsub(/[[:punct:]]/, "_", $0); print;} ;') $transcript_id ##$org $file #$filename
fi

#echo "----------------"
done

exit 0
