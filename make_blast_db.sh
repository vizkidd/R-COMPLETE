#!/bin/bash

make_db() {
echo "$1"
if [ -s "$1" ]; then
	makeblastdb -in "$1" -dbtype nucl  -hash_index || true # -parse_seqids -out $2
else
	echo "$1 empty..."
fi
}

collate_fasta() {
	path=$1
	gene_name=$2
	safe_gene_name=$(echo $gene_name | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')
	reg=$3
	out_dir=$4
	if find $path/*/ -name "$safe_gene_name*" -type f ! -size 0 | grep -w -i "$reg"; then
		cat $(find $path/*/ -name "$safe_gene_name*" -type f ! -size 0 | grep -w -i "$reg" ) > "$out_dir/blastdb/$safe_gene_name.$reg"
		make_db "$out_dir/blastdb/$safe_gene_name.$reg"
	fi
}

#$1 - FASTA PATH (files/fasta)
#$2 - genelist (files/genelist.txt)
#$3 - blast DB output directory (files/)

rm -rf $3/blastdb 
mkdir $3/blastdb

while IFS= read -r gene
do
echo $gene
collate_fasta $1 $gene exons $3
collate_fasta $1 $gene 3utr $3
collate_fasta $1 $gene 5utr $3
collate_fasta $1 $gene cds $3
collate_fasta $1 $gene noncodingexons $3
done < "$2"
