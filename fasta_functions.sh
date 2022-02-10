#!/bin/bash

index_fastaIDs() {
	#Called using : index_fastaIDs files/rna_ids.txt files/fasta/*.aln
	if [ $# -eq 0 ]
  	then
    	echo "Give output filename, (FASTA folder/Alignment files) as input(exiting)."
    	echo "eg: index_fastaIDs files/rna_ids.txt files/fasta/*.aln"
    	return 1
	fi
	
	arg_arr=($@)
	x=10000
	
	grep ">" -H -r "${arg_arr[@]:1}" | awk -F'>' -v x=$x '{ x = ++x; print substr($1, 1, length($1)-1) "\t" $2 "\t" x; }' > $1
}

fastaID_to_number() {
	#SLOWER implementation
	#while IFS= read -r rna_id
	#do
	#	rna_file=$(echo $rna_id | awk '{print $1}')
	#	rna_name=$(echo $rna_id | awk '{print $2}')
	#	rna_num=$(echo $rna_id | awk '{print $3}')
	#	echo $rna_file $rna_name $rna_num
	#	sed --in-place "s/$rna_name/$rna_num/g" $rna_file 
	#done < "files/rna_ids.txt"
	arg_arr=($@)

	if [ $# -eq 0 ]
  	then
    	echo "Executed with fastaID_to_number transcript_metadata (and optional list of files to replace fasta IDs)"
    	echo "Transcript metadata has (filename, FASTA ID, numeric ID), can be generated with index_fastaIDs"
    	return 1
	fi
	if [ $# == 1 ]
	then
		file_ids=$(parallel "basename {}" ::: $(awk '{print $1}' $1 | sort | uniq))
		parallel -j ${#file_ids} "fID_to_num_parallel $1 {}" ::: ${file_ids[@]}
	else
		parallel -j $# "fID_to_num_multi $1 {}" ::: "${arg_arr[@]:1}"
	fi
}

fID_to_num_multi() {
	#echo $@
	if [ $# -lt 2 ]
  	then
    	echo "Please use fastaID_to_number, this is a parallel(multi file) implementation"
    	return 1
	fi

	f_name=$(basename $2)

	##GETS FASTA IDs which are same between metadata and user provided file
	grep -w -z -o "$(echo $(grep "${f_name%%.*}" $1 | awk '{print $2}'))" $2 | while read -r line; do 
		rna_num=$(grep -w "$line" $1 | awk '{print $3}')
		printf "%s %s -> %s\n" $2 $line $rna_num
		sed --in-place "s/\<$line\>/$rna_num/g" $2
	done
}

fID_to_num_parallel() {
	#PARALLEL implementation
	if [ $# -lt 2 ]
  	then
    	echo "Please use fastaID_to_number, this is a parallel implementation"
    	return 1
	fi

	grep -w $2 $1 | while read -r line; do 
		rna_file=$(echo $line | awk '{print $1}')
		rna_name=$(echo $line | awk '{print $2}')
		rna_num=$(echo $line | awk '{print $3}')
		printf "%s %s -> %s\n" $rna_file $rna_name $rna_num
		sed --in-place "s/\<$rna_name\>/$rna_num/g" $rna_file 
	done
}

number_to_fastaID() {
	#SLOWER implementation
	#while IFS= read -r rna_id
	#do
	#	rna_file=$(echo $rna_id | awk '{print $1}')
	#	rna_name=$(echo $rna_id | awk '{print $2}')
	#	rna_num=$(echo $rna_id | awk '{print $3}')
	#	echo $rna_file $rna_name $rna_num
	#	sed --in-place "s/$rna_num/$rna_name/g" $rna_file 
	#done < "files/rna_ids.txt"

	#echo "$# files"
	arg_arr=($@)

	if [ $# -eq 0 ]
  	then
    	echo "Executed with number_to_fastaID transcript_metadata (and optional list of files to replace fasta IDs)"
    	echo "Transcript metadata has (filename, FASTA ID, numeric ID), can be generated with index_fastaIDs"
    	return 1
	fi
	if [ $# == 1 ]
	then
		file_ids=$(parallel "basename {}" ::: $(awk '{print $1}' $1 | sort | uniq))
		parallel -j ${#file_ids} "num_to_fID_parallel $1 {}" ::: ${file_ids[@]}
	else
		parallel -j $# "num_to_fID_multi $1 {}" ::: "${arg_arr[@]:1}"
	fi
}

num_to_fID_multi() {
	#echo $@
	if [ $# -lt 2 ]
  	then
    	echo "Please use number_to_fastaID, this is a parallel(multi file) implementation"
    	return 1
	fi
	##GETS numeric IDs which are same between metadata and user provided file
	f_name=$(basename $2)
	
	#grep -w -z -o "$(seq $(awk 'BEGIN{a=   0}{if ($3>0+a) a=$3} END{print a}' $1) )" $2 | while read -r line; do 
	grep -w -z -o "$(echo $(grep "${f_name%%.*}" $1 | awk '{print $3}'))" $2 | while read -r line; do 
		rna_name=$(grep -w "$line" $1 | awk '{print $2}')
		printf "%s %s -> %s\n" $2 $line $rna_name
		sed --in-place "s/\<$line\>/$rna_name/g" $2
	done
}

num_to_fID_parallel() {
	#PARALLEL implementation
	#id_sub_list=($(grep $2 $1))
	#echo ${id_sub_list[*]}
	if [ $# -lt 2 ]
  	then
    	echo "Please use number_to_fastaID, this is a parallel implementation"
    	return 1
	fi

	grep -w $2 $1 | while read -r line; do 
		rna_file=$(echo $line | awk '{print $1}')
		rna_name=$(echo $line | awk '{print $2}')
		rna_num=$(echo $line | awk '{print $3}')
		printf "%s %s -> %s\n" $rna_file $rna_name $rna_num
		sed --in-place "s/\<$rna_num\>/$rna_name/g" $rna_file 
	done
}

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
		if [[ ! -s $(blastdb_path -db "$out_dir/$safe_gene_name.$reg") ]]; then
			touch "$out_dir/$safe_gene_name.$reg"
			cat $(find $path/*/ -name "$safe_gene_name*" -type f ! -size 0 | grep -w -i "$reg" ) > "$out_dir/$safe_gene_name.$reg"
			make_db "$out_dir/$safe_gene_name.$reg"
		fi
	fi
}

select_transcripts() { 
# $1 - transcript list file
# $2 - intermediate output (transcripts from CURRENT organism)
# $3 - final output (transcripts from ALL organisms)
# $4 - organism fasta path
## $5 - QSUB process IDs --OBSELETE--

OFS=$'\n' grep -I -i -r -f $1 $4 > $2 ##Get the sequence IDs and file paths of all the valid transcripts
cat $2 >> $3
#qsub -V -N $2 -hold_jid $5 ./select_transcripts.sh $2 $4
./select_transcripts.sh $2 $4
}


delete_empty_orgs() {
	#For each org
	#files/fasta - $1
	#files/genelist.txt - $2
	#files/UNAVAILABLE_ORGS -$3
	rm $3
	for org in $1/*
	do
		if [[ $(grep -c -v -i -f $org/MISSING_GENES $2) == 0 ]]; 
		then
		 	echo $org "is empty!...Removing"
			echo $org >> $3
		 	rm -rf $org
		fi 
	done
}

count_seqs4genes() {
	grep -r -c ">" $1 | awk -F'/' '{print $NF}' | awk -F':' '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }'  | sort | grep "$2" > $3
}

count_genes4orgs() {
	#For each org
	#files/fasta - $1
	#files/genelist.txt - $2
	rm $3
	for org in $1/*
	do
		genecount=$(find $org/* | grep -f $2 | awk -F'/' '{print $(NF-1),$(NF)}' | awk -F'.' '{print $1}' | sort | uniq | wc -l)
		org_name=$(echo $org | awk -F'/' '{print $NF}')
		echo $org_name,$genecount >> $3
	done
}

# mask_stops_3utr() {
# 	#python mask_motifs.py -f $1 -s 3 --pos 1 --mask "N"  -r TRUE --add TRUE -cm "TGA,TAA,TAG" -cmm - -o $1 ##Compare mask is negative, so this will add NNN if stop codon doesnt exist
# 	$PY2_PATH mask_motifs.py -f $1 -s 3 -m "N" -p 1 -cm "TGA,TAA,TAG" -r True -rf 1 -o $1
# }

# mask_stops_cds() {
# 	$PY2_PATH mask_motifs.py -f $1 -s 3 -m "N" -p 0 -cm "TGA,TAA,TAG" -r True -rf 1 -o $1
# 	$PY2_PATH mask_motifs.py -f $1 -s 3 -m "N" -p 0 -cm "TGA,TAA,TAG" -r True -rf 2 -o $1
# 	$PY2_PATH mask_motifs.py -f $1 -s 3 -m "N" -p 0 -cm "TGA,TAA,TAG" -r True -rf 3 -o $1
# }

export -f number_to_fastaID
export -f fastaID_to_number
export -f index_fastaIDs
export -f make_db
export -f collate_fasta
export -f select_transcripts
export -f delete_empty_orgs
export -f count_genes4orgs
export -f count_seqs4genes

export -f fID_to_num_parallel
export -f fID_to_num_multi
export -f num_to_fID_parallel
export -f num_to_fID_multi