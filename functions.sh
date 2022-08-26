#!/bin/bash

##TERMINAL COLORS
# Reset
Color_Off='\033[0m'       # Text Reset

# High intensity FG Colors

Default=39
Black=90
Red=91
Green=92
Yellow=93
Blue=94
Purple=95
Cyan=96
White=97

# BG Colors

BG_Default=49
BG_Black=100
BG_Red=101
BG_Green=102
BG_Yellow=103
BG_Blue=104
BG_Purple=105
BG_Cyan=106
BG_White=107

Bold_Style=1
Reset_Style=0

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

function index_fastaIDs() {
	#Called using : index_fastaIDs files/rna_ids.txt files/fasta/*.aln
	if [ $# -eq 0 ]
  	then
    	echo "Give output filename, (FASTA folder/Alignment files) as input(exiting)."
    	echo "eg: index_fastaIDs files/rna_ids.txt files/fasta/*.aln"
    	return 1
	fi
	
	local arg_arr=($@)
	local x=10000
	
	grep ">" -H -r "${arg_arr[@]:1}" | awk -F'>' -v x=$x '{ x = ++x; print substr($1, 1, length($1)-1) "\t" $2 "\t" "i"x; }' > $1
}

function fastaID_to_number() {
	#SLOWER implementation
	#while IFS= read -r rna_id
	#do
	#	rna_file=$(echo $rna_id | awk '{print $1}')
	#	rna_name=$(echo $rna_id | awk '{print $2}')
	#	rna_num=$(echo $rna_id | awk '{print $3}')
	#	echo $rna_file $rna_name $rna_num
	#	sed --in-place "s/$rna_name/$rna_num/g" $rna_file 
	#done < "files/rna_ids.txt"
	local arg_arr=($@)

	if [ $# -eq 0 ]
  	then
    	echo "Executed with fastaID_to_number transcript_metadata (and optional list of files to replace fasta IDs)"
    	echo "If optional files are given then the IDs in files from transcript_metadata are not replaced "
    	echo "Transcript metadata has (filename, FASTA ID, numeric ID), can be generated with index_fastaIDs"
    	return 1
	fi
	if [ $# == 1 ]
	then
		local file_ids=$(parallel "basename {}" ::: $(awk '{print $1}' $1 | sort | uniq))
		parallel -j ${#file_ids} "fID_to_num_parallel $1 {}" ::: ${file_ids[@]}
	else
		parallel -j $# "fID_to_num_multi $1 {}" ::: "${arg_arr[@]:1}"
	fi
}

function fID_to_num_multi() {
	#echo $@
	if [ $# -lt 2 ]
  	then
    	echo "Please use fastaID_to_number, this is a parallel(multi file) implementation"
    	return 1
	fi

	local f_name=$(basename $2)

	##GETS FASTA IDs which are same between metadata and user provided file
	grep -w -z -o "$(grep "${f_name%%.*}" $1 | awk '{print $2}')" $2 | while read -r line; do 
		local rna_num=$(grep -w "$line" $1 | awk '{print $3}')
		printf "%s %s -> %s\n" $2 $line $rna_num
		sed --in-place "s/\<$line\>/$rna_num/g" $2
	done
}

function fID_to_num_parallel() {
	#PARALLEL implementation
	if [ $# -lt 2 ]
  	then
    	echo "Please use fastaID_to_number, this is a parallel implementation"
    	return 1
	fi

	grep -w $2 $1 | while read -r line; do 
		local rna_file=$(echo $line | awk '{print $1}')
		local rna_name=$(echo $line | awk '{print $2}')
		local rna_num=$(echo $line | awk '{print $3}')
		printf "%s %s -> %s\n" $rna_file $rna_name $rna_num
		sed --in-place "s/\<$rna_name\>/$rna_num/g" $rna_file 
	done
}

function number_to_fastaID() {
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
	local arg_arr=($@)

	if [ $# -eq 0 ]
  	then
    	echo "Executed with number_to_fastaID transcript_metadata (and optional list of files to replace fasta IDs)"
    	echo "If optional files are given then the IDs in files from transcript_metadata are not replaced "
    	echo "Transcript metadata has (filename, FASTA ID, numeric ID), can be generated with index_fastaIDs"
    	return 1
	fi
	if [ $# == 1 ]
	then
		local file_ids=$(parallel "basename {}" ::: $(awk '{print $1}' $1 | sort | uniq))
		parallel -j ${#file_ids} "num_to_fID_parallel $1 {}" ::: ${file_ids[@]}
	else
		parallel -j $# "num_to_fID_multi $1 {}" ::: "${arg_arr[@]:1}"
	fi
}

function num_to_fID_multi() {
	#echo $@
	if [ $# -lt 2 ]
  	then
    	echo "Please use number_to_fastaID, this is a parallel(multi file) implementation"
    	return 1
	fi
	##GETS numeric IDs which are same between metadata and user provided file
	local f_name=$(basename $2)
	echo $@
	echo $f_name
	#grep -w -z -o "$(seq $(awk 'BEGIN{a=   0}{if ($3>0+a) a=$3} END{print a}' $1) )" $2 | while read -r line; do 
	grep -w -z -o "$(grep "${f_name%%.*}" $1 | awk '{print $3}')" $2 | while read -r line; do 
		local rna_name=$(grep -w "$line" $1 | awk '{print $2}')
		printf "%s %s -> %s\n" $2 $line $rna_name
		sed --in-place "s/\<$line\>/$rna_name/g" $2
	done
}

function num_to_fID_parallel() {
	#PARALLEL implementation
	#id_sub_list=($(grep $2 $1))
	#echo ${id_sub_list[*]}
	if [ $# -lt 2 ]
  	then
    	echo "Please use number_to_fastaID, this is a parallel implementation"
    	return 1
	fi

	grep -w $2 $1 | while read -r line; do 
		local rna_file=$(echo $line | awk '{print $1}')
		local	rna_name=$(echo $line | awk '{print $2}')
		local	rna_num=$(echo $line | awk '{print $3}')
		printf "%s %s -> %s\n" $rna_file $rna_name $rna_num
		sed --in-place "s/\<$rna_num\>/$rna_name/g" $rna_file 
	done
}

function make_db() {
echo "$1"
if [ -s "$1" ]; then
	makeblastdb -in "$1" -dbtype nucl  -hash_index || true # -parse_seqids -out $2
else
	echo "$1 empty..."
fi
}

function collate_fasta() {
	local path=$1
	local gene_name=$2
	local safe_gene_name=$(echo $gene_name | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')
	local reg=$3
	local out_dir=$4
	if find $path/*/ -name "$safe_gene_name*" -type f ! -size 0 | grep -w -i "$reg"; then
		if [[ ! -s $(blastdb_path -db "$out_dir/$safe_gene_name.$reg") ]]; then
			touch "$out_dir/$safe_gene_name.$reg"
			cat $(find $path/*/ -name "$safe_gene_name*" -type f ! -size 0 | grep -w -i "$reg" ) > "$out_dir/$safe_gene_name.$reg"
			make_db "$out_dir/$safe_gene_name.$reg"
		fi
	fi
}

function select_transcripts() { 
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


function delete_empty_orgs() {
	#For each org
	#files/fasta - $1
	#files/genelist.txt - $2
	#files/UNAVAILABLE_ORGS -$3
	rm $3
	#for org in $1/*
	#do
	#	if [[ $(grep -c -v -i -f $org/MISSING_GENES $2) == 0 ]]; 
	#	then
	#	 	echo $org "is empty!...Removing"
	#		echo $org >> $3
	#	 	rm -rf $org
	#	fi 
	#done
	for org in $(cat $2 | grep -c -v -i -z $(find $1 -name MISSING_GENES))
	do
		if [[ $(echo $org | awk -F':' '{print $2}') == 0 ]];
		then
			echo $org "is empty!...Removing"
			echo $org >> $3
			rm -rf $org
		fi
	done
}

function count_seqs4genes() {
	if [ $# -eq 0 ]
  	then
    	echo "Give path to sequence files/alignment files, region/' ', out_file to get the distribution of sequences across genes"
    	echo "eg, count_seqs4genes files/fasta/ 3utr files/3utr_seqs_meta.txt"
    	echo "eg, count_seqs4genes files/alns/*.aln ' ' files/aln_seqs_meta.txt"
    	return 1
	fi
	grep -r -c ">" $1 | awk -F'/' '{print $NF}' | awk -F':' '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }'  | sort | grep "$2" > $3
}

function count_genes4orgs() {
	#For each org
	#files/fasta - $1
	#files/genelist.txt - $2
	if [ $# -eq 0 ]
  	then
    	echo "Give path to organism folders, genelist, out_file to get the distribution of genes across organisms"
    	echo "eg, count_genes4orgs files/fasta files/genelist.txt files/gene_counts.txt"
    	return 1
	fi
	rm $3
	for org in $1/*
	do
		local genecount=$(find $org/* | grep -f $2 | awk -F'/' '{print $(NF-1),$(NF)}' | awk -F'.' '{print $1}' | sort | uniq | wc -l)
		local org_name=$(echo $org | awk -F'/' '{print $NF}')
		echo $org_name,$genecount >> $3
	done
}

function get_length_dist() {
 	if [ $# -eq 0 ]
  	then
    	echo "Give path to sequence files to get the length distribution, use -v to print file names"
    	return 1
	fi
  	local f_arr=($@); 
  	#print_fn=$(echo $@ | grep -w -q "\-v");
  	local seq_files=( “${f_arr[@]/\-v}” ); 
  	#echo ${seq_files[@]};
  	for seq_f in ${seq_files[@]};
  	do
	  	if [[ -s $seq_f ]];
  		then
	  		if $(echo $@ | grep -q -w "\-v") ;
	  		then
	  			echo $seq_f;
  			fi;
  			awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $seq_f | grep '>' -v | sort | uniq ;
  		fi
  	done; 
}

function get_count_dist() {
 	if [ $# -eq 0 ]
  	then
    	echo "Give path to sequence files to get the count distribution(number of sequences in each file), use -v to print file names"
    	return 1
	fi
  	local f_arr=($@); 
  	#print_fn=$(echo $@ | grep -w -q "\-v");
  	local seq_files=( “${f_arr[@]/\-v}” ); 
  	#echo ${seq_files[@]};
  	for seq_f in ${seq_files[@]};
  	do
	  	if [[ -s $seq_f ]];
  		then
	  		if $(echo $@ | grep -q -w "\-v") ;
	  		then
	  			echo $seq_f;
  			fi;
  			grep '>' $seq_f | wc -l;
  		fi
  	done; 
}

function get_pairwise_pid() {
	#INPUT: fasta/aln files
	#For each sequence in the alignment get its pid's from BLAST hits and in the end average them
	local arg_arr=($@);
	echo ${arg_arr[0]}
	return 0
	for seq_file in "${arg_arr[@]}-1";
	do
		echo $seq_file
	done;
}

function get_fasta() {
#1 - gene name
#2 - transcript region (cds/5utr/3utr/exons/noncodingexons)
#3 - basename for bedfiles
#4 - genome fasta file
#5 - org name
#6 - FASTA OUPUT FOLDER
#7 - Formatted file name for gene (without punctuation)
#8 - GTF file path (for extracting the 'real' gene annotation)
#9 - TEMP PATH

local gene_full=$1
local s_name=$(echo $gene_full | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')
local reg=$2
local base_bed=$3
local genome_fa=$4
local org_name=$5
local FASTA_PATH=$6
local GTF_PATH=$7
local TEMP_PATH=$8
local LABEL_FASTA=$9

if [ ! -s $FASTA_PATH/$s_name.$reg ]; then #$FASTA_PATH/$file_out.cds
  >&2 echo $s_name $reg
  
  if [[ -s "$base_bed"_"$reg.bed" ]]; then
  	grep -i -w $gene_full files/genes/$org_name/gtf_stats.csv | awk -F',' '{print $3}' | grep -w -f - "$base_bed"_"$reg.bed" > $TEMP_PATH/"$s_name"_"$reg.bed"
  fi

  ##TO get flanks 
  #grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3"

  ##TO find MEDIAN values
  #grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | sort -n | awk '{arr[NR]=$1} END {if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}'

  ##TO find MEAN
  #grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'

  if [[ -s $TEMP_PATH/"$s_name"_"$reg.bed" ]]; then 
    ln -r -f -s $TEMP_PATH/"$s_name"_"$reg.bed" $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed
  elif [[ -s $TEMP_PATH/"$s_name"_cds.bed && $reg!="cds" ]]; then
    if [[ "$reg" == "3utr" ]]; then
      local flank_len=$(grep -w -i "$gene_full" files/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($9 > 3) sum += $9; n++ } END { if (n > 0) print sum / n; }') ##3' UTR length must be > 3 (because of existence of stop codons)
      if [[ $flank_len == "" ]]; then
        return 2
      fi
      bedtools flank -s -l $flank_len -r 0 -i $TEMP_PATH/"$s_name"_cds.bed -g $genome_fa.fai > $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed 
    elif [[ "$reg" == "5utr" ]]; then
      local flank_len=$(grep -w -i "$gene_full" files/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($8 > 0) sum += $8; n++ } END { if (n > 0) print sum / n; }') ##5' UTR length must be > 0
      if [[ $flank_len == "" ]]; then
        return 2
      fi
      bedtools flank -s -r $flank_len -l 0 -i $TEMP_PATH/"$s_name"_cds.bed -g $genome_fa.fai > $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed 
    fi
    if [[ -s $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed ]]; then
      sed -i "s/cds/$(echo $reg)_FLANK/g" $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed
    else
      >&2 echo "$reg FLANKS not found for $s_name"
      return 1
    fi
  fi

  if [[ -s $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed || -L $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed ]]; then
    if [[ $LABEL_FASTA ==  "TRUE" ]] ; then
      bedtools getfasta -s -split -fi $genome_fa -bed $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed -nameOnly -fullHeader > "$FASTA_PATH/$s_name.$reg.tmp" #NOTUSING name+ because it also gives coordinates
    else
      bedtools getfasta -s -split -fi $genome_fa -bed $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed -nameOnly -fullHeader > "$FASTA_PATH/$s_name.$reg" #NOTUSING name+ because it also gives coordinates
    fi
  fi
fi

return 0
}

# function index_genome(){
# 	# 1 - GENOME FILE
# 	# 2 - FIFO pipe name
# 	local GENOME_FILE=$1
# 	local FIFO_FILE=$2
# 	mkfifo $FIFO_FILE
# 	if [[ ${GENOME_FILE##*.} == "gz" ]] ; then
#   	local gfile_name=${GENOME_FILE%.*}
# 	  zcat -f $GENOME_FILE | tee $FIFO_FILE $gfile_name &
# 	else
# 		local gfile_name=$GENOME_FILE
# 		zcat -f $GENOME_FILE > $FIFO_FILE &
# 	fi

#   time samtools faidx --fai-idx $gfile_name.fai $FIFO_FILE
# }

function color_FG(){
	# 1 - FG color
	# 2 - Message
	echo -e '\033[0;'$1';49m '$2 ${Color_Off}
}

function color_FG_Bold(){
	# 1 - FG color
	# 2 - Message
	echo -e '\033[1;'$1';49m '$2 ${Color_Off}
}

function color_BG(){
	# 1 - BG color
	# 2 - Message
	echo -e '\033[0;39;'$1'm '$2 ${Color_Off}
}

function color_FG_BG(){
	# 1 - FG color
	# 2 - BG color
	# 3 - Message
	echo -e '\033[0;'$1';'$2'm '$3 ${Color_Off}
}

function color_FG_BG_Bold(){
	# 1 - FG color
	# 2 - BG color
	# 3 - Message
	echo -e '\033[1;'$1';'$2'm '$3 ${Color_Off}
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
#export -f index_genome
export -f get_fasta
export -f get_length_dist
export -f get_count_dist
export -f check_param
#export -f checkForVariable

export -f fID_to_num_parallel
export -f fID_to_num_multi
export -f num_to_fID_parallel
export -f num_to_fID_multi

export -f color_FG
export -f color_FG_Bold
export -f color_BG
export -f color_FG_BG
export -f color_FG_BG_Bold
