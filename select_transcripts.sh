#!/bin/bash
#$ -cwd
#$ -V
#$ -l data
#$ -o "logs/select_transcripts.o"
#$ -e "logs/select_transcripts.e"

#echo $1
FASTA_PATH=$2 #$2 - files/fasta
#$1 - selected_ENSDART.txt

find $FASTA_PATH -type f -name ".*" -exec rm -f {} + #Deleting residues 
find $FASTA_PATH -type f -name "*.fai" -exec rm -f {} +


#sleep 20

while IFS=$'>' read -r line seq_ID
do
	file_path=$(dirname $(echo "${line%?}" | awk -F':>' '{print $1}'))
	file=$(basename $(echo "${line%?}" | awk -F':>' '{print $1}'))
	#seq_ID=$(awk -F':>' '{print $2}') ##DO NOT UNCOMMENT... read command reads into this variable
	#echo "$file" "$file_path"
	echo "$line" "$seq_ID"
	#echo "$file_path/.$file" 
	samtools faidx "$file_path/$file"
	samtools faidx "$file_path/$file" "$seq_ID" >> "$file_path/.$file" ##Save selected sequences to hidden files
done < $1

find $FASTA_PATH -name ".*" -exec rename -v 's|/\.|/|' {} +  # rename -v -f -n 's|/\.|/|' {} + ##Rename hidden files to replace the original files
find $FASTA_PATH -type f -name "*.fai" -exec rm -f {} + #Deleting index files
find $FASTA_PATH -type f -name ".*" -exec rm -f {} + #Deleting residues 

exit 0