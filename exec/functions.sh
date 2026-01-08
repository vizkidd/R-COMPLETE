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
  
  #if [[ -z $tmp_var ]]; then
  #  echo "Parameter Error : $tmp_name is empty"
  #  exit 1
  #fi
  
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
	
	grep ">" -H -r "${arg_arr[@]:1}" | awk -F'>' -v x=$x '{ x = ++x; print substr($1, 1, length($1)-1) "\t" $2 "\t" "i"x; }' > "${arg_arr[0]}" #$1
}

function shorten_fastaIDs() {
	local arg_arr=($@)
	export GNU_PARALLEL=${arg_arr[1]}
	local transcript_meta=${arg_arr[0]}
	local fasta_files=(${arg_arr[@]:2})
	echo "${fasta_files[@]}"
	# if [ $# -eq 0 ]
  # 	then
  #   	echo "Executed with shorten_fastaIDs transcript_metadata (or optional list of files to replace fasta IDs)"
  #   	echo "If optional files are given then the IDs in files from transcript_metadata are not replaced "
  #   	echo "Transcript metadata has (filename, FASTA ID, numeric ID), can be generated with index_fastaIDs"
  #   	return 1
	# fi
if [[ ${#fasta_files} > 0 ]]; then
	$GNU_PARALLEL -j0 "sed -i -f <(grep -w {} $transcript_meta | awk -F'\t' '{print \"s/\"\$2\"/\"\$3\"/g\"}') {}" ::: $(grep -f <(echo "${fasta_files[@]}") $transcript_meta | awk -F'\t' '{print $1}' | sort -u)
else
	$GNU_PARALLEL -j0 "sed -i -f <(grep -w {} $transcript_meta | awk -F'\t' '{print \"s/\"\$2\"/\"\$3\"/g\"}') {}" ::: "$(awk -F'\t' '{print $1}' $transcript_meta | sort -u)"
fi

}

function lengthen_fastaIDs() {
	local arg_arr=($@)
	export GNU_PARALLEL=${arg_arr[1]}
	local transcript_meta=${arg_arr[0]}
	local fasta_files=(${arg_arr[@]:2})

	# if [ $# -eq 0 ]
  # 	then
  #   	echo "Executed with lengthen_fastaIDs transcript_metadata (and optional list of files to replace fasta IDs)"
  #   	echo "If optional files are given then the IDs in files from transcript_metadata are not replaced "
  #   	echo "Transcript metadata has (filename, FASTA ID, numeric ID), can be generated with index_fastaIDs"
  #   	return 1
	# fi
	
if [[ ${#fasta_files} > 0 ]]; then
	$GNU_PARALLEL -j0 "sed -i -f <(grep -w {} $transcript_meta | awk -F'\t' '{print \"s/\"\$3\"/\"\$2\"/g\"}') {}" ::: $(grep -f <(echo "${fasta_files[@]}") $transcript_meta | awk -F'\t' '{print $1}' | sort -u)
else
	$GNU_PARALLEL -j0 "sed -i -f <(grep -w {} $transcript_meta | awk -F'\t' '{print \"s/\"\$3\"/\"\$2\"/g\"}') {}" ::: "$(awk -F'\t' '{print $1}' $transcript_meta | sort -u)"
fi

}

function check_DB(){
	local script_args=($(echo $@))
	local fasta_file=${script_args[0]}
	local blast_bin=${script_args[1]}
	$blast_bin/blastdb_path -db $fasta_file
	return 0
}

function make_BLAST_DB() {
	# $1 - FASTA file
	# $2 - path to BLAST bin
	local script_args=($(echo $@))
	local fasta_file=${script_args[0]}
	local blast_bin=${script_args[1]}
	#echo "$1"
	if [[ ! -s "$($blast_bin/blastdb_path -db "$fasta_file")" ]]; then #if [ -s "$1" ]; then
		>&2 $blast_bin/makeblastdb -in "$fasta_file" -dbtype nucl -hash_index #-parse_seqids #|| true # -parse_seqids -out $2
	elif [[ ! -s $fasta_file ]]; then
		>&2 echo "$fasta_file empty..."
		exit 255
	else
		>&2 echo "BLAST DB for $fasta_file exists!"
		exit 0
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
			make_BLAST_DB "$out_dir/$safe_gene_name.$reg"
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

function get_FASTA() {
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
	local OUT_PATH=$9
	#local LABEL_FASTA=$9

	if [ ! -s $FASTA_PATH/$s_name.$reg ]; then #$FASTA_PATH/$file_out.cds
	  #>&2 echo $s_name $reg
	  
	  if [[ -s "$base_bed"_"$reg.bed" ]]; then
	  	grep -i -w $gene_full $OUT_PATH/genes/$org_name/gtf_stats.csv | awk -F',' '{print $3}' | grep -w -f - "$base_bed"_"$reg.bed" > $TEMP_PATH/"$s_name"_"$reg.bed"
	  fi

	  ##TO get flanks 
	  #grep -i -f $OUT_PATH/some_org/cat1.rna_list $OUT_PATH/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3"

	  ##TO find MEDIAN values
	  #grep -i -f $OUT_PATH/some_org/cat1.rna_list $OUT_PATH/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | sort -n | awk '{arr[NR]=$1} END {if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}'

	  ##TO find MEAN
	  #grep -i -f $OUT_PATH/some_org/cat1.rna_list $OUT_PATH/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'

	  if [[ -s $TEMP_PATH/"$s_name"_"$reg.bed" ]]; then 
	    ln -r -f -s $TEMP_PATH/"$s_name"_"$reg.bed" $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed
	  elif [[ -s $TEMP_PATH/"$s_name"_cds.bed && $reg!="cds" ]]; then
	    if [[ "$reg" == "3utr" ]]; then
	      local flank_len=$(grep -w -i "$gene_full" $OUT_PATH/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($9 > 3) sum += $9; n++ } END { if (n > 0) print sum / n; }') ##3' UTR length must be > 3 (because of existence of stop codons)
	      if [[ $flank_len == "" ]]; then
	        return 2
	      fi
	      bedtools flank -s -l $flank_len -r 0 -i $TEMP_PATH/"$s_name"_cds.bed -g $genome_fa.fai > $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed 
	    elif [[ "$reg" == "5utr" ]]; then
	      local flank_len=$(grep -w -i "$gene_full" $OUT_PATH/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($8 > 0) sum += $8; n++ } END { if (n > 0) print sum / n; }') ##5' UTR length must be > 0
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
	    #if [[ $LABEL_FASTA ==  "TRUE" ]] ; then
	      bedtools getfasta -s -split -fi $genome_fa -bed $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed -name -fullHeader > "$FASTA_PATH/$s_name.$reg" #-split #NOTUSING name+ because it also gives coordinates
	    #else
	    #  bedtools getfasta -s -split -fi $genome_fa -bed $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed -nameOnly -fullHeader > "$FASTA_PATH/$s_name.$reg" #NOTUSING name+ because it also gives coordinates
	    #fi
	  fi
	fi

	return 0
}

function accumulate_clusters(){
	# $1 - ORGANISM FASTA PATH (../files/fasta/danio_rerio/)
	# $2 - Transcript sequence ID delimiter (::)
	# $3 - Output file ($OUT_PATH/genes/danio_rerio/ALL_CLUSTERS)
	grep -h ">" "$1"/* | awk -F"$2" '{print $NF}' | awk '{split($0,a,","); for(key in a) print a[key];}' | sort -u > $3
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

function index_genome(){
	# This function also extracts genome because bedtools getfasta requires it to be
	# 1 - Genome.fa File path

	local GENOME_FILE=$1
	#FIFO_FILE="$TEMP_PATH/$f_org_name/genome_pipe"
	#rm -f $FIFO_FILE
	#mkfifo $FIFO_FILE
	if [[ ${GENOME_FILE##*.} == "gz" ]] ; then
  	gfile_name=${GENOME_FILE%.*}
 	 zcat -f $GENOME_FILE > $gfile_name #| tee $gfile_name > $FIFO_FILE  &
  	#genome_ext_proc=$(echo $!)
	else
  	gfile_name=$GENOME_FILE
  	#zcat -f $GENOME_FILE > $FIFO_FILE &
	fi

	#if [[ ! -s $gfile_name.fai ]]; then
  #samtools faidx --fai-idx $gfile_name.fai $gfile_name #&
  samtools faidx $gfile_name #--fai-idx $gfile_name.fai 
  	#genome_index_proc=$(echo $!)
	#else
	#  rm -f $FIFO_FILE
	#fi
}

# function group_FASTA_clusters(){
# 	local script_args=($(echo $@))
# 	local PARALLEL_PATH=${script_args[0]}) 
# 	local FASTA_PATH=${script_args[1]}) 
# 	local GROUPS_PATH=${script_args[2]}) 
# 	local seqID_delimiter=${script_args[3]} 
# 	local n_threads=${script_args[4]} 
# 	local OUT_PATH=${script_args[5]} 

# 	for f_org in $FASTA_PATH/*; do 
# 	f_org_name=$(basename $f_org)
# 		grep -r -h ">" $f_org  | awk -F"$seqID_delimiter" '{print $NF}' | awk '{split($0,a,","); for(key in a) print a[key];}' | sort -u > $OUT_PATH/genes/$f_org_name/ORG_CLUSTERS
# 		$PARALLEL_PATH --max-procs $n_threads " printf '%s\t%s\n' {1} {2}" :::: <(grep -H -f $OUT_PATH/genes/$f_org_name/ORG_CLUSTERS -r $FASTA_PATH/$f_org_name/ | awk -F'[:>]' -v s_delim="$seqID_delimiter" '{split($2,a,s_delim); n=split($1,b,"."); print $1"\t"$2"\t"a[5]"\t"b[n]'}) | parallel  --max-procs 1 --colsep '\t' --recend '\n'  "if [[ -s {1} && ! -z {2} && ! -z {1} && ! -z {3} ]] ; then samtools faidx {1} {2} >> $GROUPS_PATH/{3}.{4} ; fi" 
# 	done
# }

function concatenate_FASTA_groups(){
	local script_args=($(echo $@))
	local GROUPS_PATH=${script_args[0]}
	local PARALLEL_PATH=${script_args[1]} 
	if [[ ! -z $GROUPS_PATH && ! -z $PARALLEL_PATH ]]; then
		group_files=($(awk -F'.' '{print $1"."$2}' <(ls -1 $GROUPS_PATH) | sort -u))
		$PARALLEL_PATH -j0 "cat $GROUPS_PATH/{}* | awk 'NF' > $GROUPS_PATH/{}" ::: $(echo ${group_files[@]})
		rm -f $GROUPS_PATH/*.*.*
	else
		>&2 echo "concatenate_FASTA_groups() - Requires GROUPS_PATH and PARALLEL_PATH."
	fi
}

function group_FASTA_seqs(){
	local script_args=($(echo $@))
	local PARALLEL_PATH=${script_args[0]} 
	local FASTA_FILE=${script_args[1]} 
	local GROUPS_PATH=${script_args[2]} 
	local seqID_delimiter=${script_args[3]} 
	local n_threads=${script_args[4]} 
	local OUT_PATH=${script_args[5]} 
	local run_mode=${script_args[6]}
	local FASTA_PATH=${script_args[7]}

	local f_org_name=$(basename $(dirname $FASTA_FILE))
	grep ">" $FASTA_FILE | sed "s/>/\n>/g" | awk -F"$seqID_delimiter" '{print $NF}' | awk '{split($0,a,","); for(key in a) print a[key];}' | sort -u > $OUT_PATH/genes/$f_org_name/ORG_CLUSTERS.$run_mode
	#grep ">" $FASTA_FILE | sed "s/>/\n>/g" | awk -F">" -v s_delim="$seqID_delimiter" -v run_mode=$run_mode -v f_file=$FASTA_FILE '{split($2,a,s_delim); n=split(f_file,s,"."); split(a[run_mode],b,","); for (key in b) print f_file"\t"$2"\t"a[run_mode]"\t"b[key]"\t"s[n] ;}' | $PARALLEL_PATH  --max-procs $n_threads --colsep '\t' --recend '\n' "samtools faidx {1} {2} | sed \"s/>/\n>/g\" >> $GROUPS_PATH/\"{4}\".\"{5}\".$RANDOM"
	grep ">" $FASTA_FILE | sed "s/>/\n>/g" | awk -F">" -v s_delim="$seqID_delimiter" -v run_mode=$run_mode -v f_file=$FASTA_FILE '{split($2,a,s_delim); n=split(f_file,s,"."); split(a[run_mode],b,","); for (key in b) print f_file"\t"$2"\t"a[run_mode]"\t"b[key]"\t"s[n] ;}' | $PARALLEL_PATH --max-procs $n_threads --colsep '\t' --recend '\n' "samtools faidx {1} {2} | sed \"s/>/\n>/g\" >> $GROUPS_PATH/\"{4}\".\"{5}\".$RANDOM"
	rm -f $FASTA_FILE.fai

	return 0

}

function do_BLAST() {
	# $1 - Name of the BLAST run
	# $2 - Query FASTA
	# $3 - Subject BLAST DB/FASTA
	# $4 - BLAST output file
	# $5 - BLAST program
	# $6 - BLAST options

	#echo $1 $2 $3 $4 $5 $6 $7
	#>&2 echo $@
	local script_args=($(echo $@))
	local parallel_path=${script_args[0]}
	local run_name=${script_args[1]}
	local query=${script_args[2]}
	local DB=${script_args[3]}
	local BLAST_output=${script_args[4]}
	local prog=${script_args[5]}
	local n_threads=${script_args[6]}
	local prog_path=$(dirname $prog)
	local blast_options=$(echo "${script_args[@]: 7:${#script_args[@]}}")

		if [[ -s "$BLAST_output" ]]; then
			rm $BLAST_output
		fi
		if [[ ! -s $($prog_path/blastdb_path -db $DB) ]]; then
			>&2 $prog_path/makeblastdb -in "$DB" -dbtype nucl  -hash_index #|| true
		fi

		#if [[ -z $parallel_path ]]; then
		#	local parallel_path=$(which parallel)
		#fi

		#>&2 echo "$run_name Started..."
		#>&2 echo "$prog -db $DB $blast_options -out - "
		if [[ -s "$query" && -s "$DB" ]]; then
 			rm -f "$BLAST_output"
 			mkfifo "$BLAST_output" #_pipe
 	if [[ -z $parallel_path ]]; then
 		$prog -query $query -db $DB $blast_options -out - 1> $BLAST_output & #_pipe & #  #-outfmt 11 -out $BLAST_output 
 		proc=$!
 	else
 		#$parallel_path -j1 --joblog $(dirname $DB)/parallel_JOBLOG.txt --compress --pipepart -a "$query" --recstart '>' --block -1 "$prog -db $DB -outfmt 11 $blast_options -out $BLAST_output " #-word_size 5 -evalue 1e-25
 		$parallel_path -j$n_threads --joblog $(dirname $DB)/parallel_JOBLOG.txt --compress --pipepart -a "$query" --recstart '>' --line-buffer --block -1 "$prog -db $DB $blast_options -out - "  1> "$BLAST_output" & #-j 1 #-outfmt 11 -out $BLAST_output #--block -1 #--group #--keep-order #$n_threads
 		proc=$!
 	fi
		
		wait $proc
		rm -f "$BLAST_output" 

		#>&2 echo "$run_name is done"
		else
			>&2 echo "($run_name) Error: Either ($query)/($DB) is empty/not found"
			return 255
		fi
		return 0
}

function all2allblast() {
	#all2allblast $reference_ORGS $fasta_path files/genelist.txt files/all2all $blastdb_path cds cds 0 tblastx
	#1 - list of species
	#2 - FASTA path
	#3 - gene list
	#5 - blastdb path
	#4 - output dir
	#6 - query region (cds/3utr/5utr)
	#7 - subject region (cds/3utr/5utr)
	#8 - BLAST program
	proc_list=()
	mkdir $4
	while IFS= read -r subject
	do
	#cat $2/$subject/*.$7 > $5/$subject.$7
	if [[ ! -s $(blastdb_path -db $5/$subject.$7) ]]; then
		touch $5/$subject.$7
		#cat $(grep -r -i -l -f $3 $2/$subject | grep ".$7") > $5/$subject.$7
		cat $(find $2/$subject -name "*.*" | grep -i -f $3 | grep ".$7" ) > $5/$subject.$7
		makeblastdb -in "$5/$subject.$7" -dbtype nucl  -hash_index || true 
	fi
	while IFS= read -r query
	do
	if [[ ! -s $4/all2all.$query.$subject ]]; then
		if [[ "$query" != "$subject" ]]; then
			#cat $2/$query/*.$6 > $5/$query.$6
			if [[ ! -s $(blastdb_path -db $5/$query.$6) ]]; then
				touch $5/$query.$6
				#cat $(grep -r -i -l -f $3 $2/$query | grep ".$6") > $5/$query.$6
				cat $(find $2/$query -name "*.*" | grep -i -f $3 | grep ".$6" ) > $5/$query.$6
				makeblastdb -in "$5/$query.$6" -dbtype nucl  -hash_index || true 
			fi
			printf "$query\n$subject" > files/pairwise_ORGS.txt
			nohup ./jobhold.sh "$query-$subject do_blast.sh $4 $5/$query.$6 $5/$subject.$7 $8 all2all $query $subject" &>> logs/job_status.o& #n_proc=2 #/fast/gridengine/uge-8.6.14/bin/lx-amd64/
			proc_list+=("$!")
		fi
	fi
	done < "$1"
	done < "$1"

	for proc_id in "${proc_list[@]}"
	do
		echo $proc_id
		wait $proc_id #|| true
	done
	echo "Submitted jobs complete(all-all $8)"
}

function onewayblast() {
	#onewayblast files/selected_ORGS.txt $fasta_path $3 files/tables $blastdb_path cds cds 6 tblastx
	#1 - list of species
	#2 - FASTA path
	#3 - gene list
	#5 - blastdb path
	#4 - output dir
	#6 - query region (cds/3utr/5utr)
	#7 - subject region (cds/3utr/5utr)
	#8 - BLAST program
	proc_list=()
	mkdir $4
	##FOR EACH GENE DO BLAST 
	while IFS= read -r line; do
		collate_fasta $2 $line $7 $5
		while IFS= read -r species; do
	if [[ ! -s $4/$line.$species.all ]]; then
	  #gene_list+=("$line")
	  #j_name=$(echo $line | sed "s/://g")
	  j_name=$(echo $line | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')
	  ##THIS is how you submit qsub jobs from within a script
	  nohup ./jobhold.sh "$j_name do_blast.sh $4 $2/$species/$j_name.$6 $5/$j_name.$7 $8 $j_name $species all" &>> logs/job_status.o& #n_proc=2 #/fast/gridengine/uge-8.6.14/bin/lx-amd64/
	  ##proc_id=$("qsub -V -cwd -l data -N "$j_name" do_blast.sh $1 $2 "$line" $4 $5 $6 $7 $n_proc") #n_proc=2
	  #id=`echo $proc_id | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'` ##GET QSUB IDs
	  ##proc_list+=("$id",)
	  proc_list+=("$!")
	fi
	done < $1
	done < $3

	#parallel -j+$(nproc) do_blast $1 $2 {} $4 $5 $6 $7 $n_proc ::: "${gene_list[@]}"
	#echo "${proc_list[@]}"
	for proc_id in "${proc_list[@]}"
	do
		echo $proc_id
		wait $proc_id #|| true
	done
	echo "Submitted jobs complete(one-way $8)"
}

function all2all_refblast() {
	#all2all_refblast $reference_ORGS $fasta_path files/all2all/all2all.genelist files/all2all_final $blastdb_path cds cds tblastx files/oneway/SET1
	#all2all_refblast $reference_ORGS $fasta_path files/all2all/all2all.genelist files/all2all_final $blastdb_path $region $region tblastx files/oneway/set.tmp
	#1 - list of ref species
	#2 - FASTA path
	#3 - gene list
	#5 - blastdb path
	#4 - output dir
	#6 - query region (cds/3utr/5utr)
	#7 - subject region (cds/3utr/5utr)
	#8 - BLAST program
	#9 - Organism SET
	proc_list=()
	mkdir $4
	while IFS= read -r subject
	do
	#cat $2/$subject/*.$7 > $5/$subject.$7
	if [[ ! -s $(blastdb_path -db $5/$subject.$7) ]]; then
		touch $5/$subject.$7
		#cat $(grep -r -i -l -f $3 $2/$subject | grep ".$7") > $5/$subject.$7
		cat $(find $2/$subject -name "*.*" | grep -i -f $3 | grep ".$7" ) > $5/$subject.$7
		makeblastdb -in "$5/$subject.$7" -dbtype nucl  -hash_index || true 
	fi
	while IFS= read -r query
	do
	if [[ ! -s $4/all2all.$query.$subject ]] || [[ ! -s $4/all2all.$subject.$query ]]; then
		if [[ "$query" != "$subject" ]]; then
			#cat $2/$query/*.$6 > $5/$query.$6
			if [[ ! -s $(blastdb_path -db $5/$query.$6) ]]; then
				touch $5/$query.$6
				#cat $(grep -r -i -l -f $3 $2/$query | grep ".$6") > $5/$query.$6
				cat $(find $2/$query -name "*.*" | grep -i -f $3 | grep ".$6" ) > $5/$query.$6
				makeblastdb -in "$5/$query.$6" -dbtype nucl  -hash_index || true 
			fi
			printf "$query\n$subject" > files/pairwise_ORGS.txt
			nohup ./jobhold.sh "$query-$subject do_blast.sh $4 $5/$query.$6 $5/$subject.$7 $8 all2all $query $subject" &>> logs/job_status.o& #n_proc=2 #/fast/gridengine/uge-8.6.14/bin/lx-amd64/
			proc_list+=("$!")
			nohup ./jobhold.sh "$subject-$query do_blast.sh $4 $5/$subject.$7 $5/$query.$6 $8 all2all $subject $query" &>> logs/job_status.o& #n_proc=2 #/fast/gridengine/uge-8.6.14/bin/lx-amd64/
			proc_list+=("$!")
		fi
	fi
	done < "$9"
	done < "$1"

	for proc_id in "${proc_list[@]}"
	do
		echo $proc_id
		wait $proc_id #|| true
	done
	echo "Submitted jobs complete(all-all(ref) $8)"
}

function oneway_RBH() {
	# $1 - FOLDER PATH
	# $2 - Query organism
	# $3 - Subject organism
	# $4 $5 $6 - Python path, RBH script and SAME_GENE flage for RBH script
	path=$1
	query=$2
	subject=$3
	rbh_script=$(echo "$4 $5")
	echo $rbh_script $6
	proc_list=()
	if [[ "$query" != "$subject" ]]; then
		if [[ ! -s "$path/$query-$subject.orths" ]]; then
			nohup blast_formatter -archive $path/all2all.$query.$subject -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out $path/$query-$subject.out &>> /dev/null&
			proc_list+=("$!")
			nohup blast_formatter -archive $path/all2all.$subject.$query -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out $path/$subject-$query.out &>> /dev/null&
			proc_list+=("$!")
			for proc_id in "${proc_list[@]}"
			do
				echo $proc_id
				wait $proc_id #|| true
			done
			#/data/meyer/viz/tools/miniconda3/envs/local_root/bin/python transcriptologs.py -i1 $query-$subject.out -i2 $subject-$query.out -o $query-$subject.orths
			nohup $rbh_script -i1 $path/$query-$subject.out -i2 $path/$subject-$query.out -o $path/$query-$subject.orths $6 &>> logs/oneway_RBH.o&
			wait "$!"
		fi
	fi
}

function twoway_RBH() {
	# $1 - FOLDER PATH
	# $2 - Query organism
	# $3 - Subject organism
	# $4 $5 $6 - Python path, RBH script and SAME_GENE flage for RBH script
	path=$1
	query=$2
	subject=$3
	rbh_script=$(echo "$4 $5")
	proc_list=()
	if [[ "$query" != "$subject" ]]; then
		if [[ ! -s "$path/$query-$subject.orths" ]] || [[ ! -s "$path/$subject-$query.orths" ]]; then
			nohup blast_formatter -archive $path/all2all.$query.$subject -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out $path/$query-$subject.out &>> /dev/null&
			proc_list+=("$!")
			nohup blast_formatter -archive $path/all2all.$subject.$query -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out $path/$subject-$query.out &>> /dev/null&
			proc_list+=("$!")
			for proc_id in "${proc_list[@]}"
			do
				echo $proc_id
				wait $proc_id #|| true
			done
			proc_list=()
			#/data/meyer/viz/tools/miniconda3/envs/local_root/bin/python transcriptologs.py -i1 $query-$subject.out -i2 $subject-$query.out -o $query-$subject.orths
			nohup $rbh_script -i1 $path/$query-$subject.out -i2 $path/$subject-$query.out -o $path/$query-$subject.orths $6 &>> logs/twoway_RBH.o&
			proc_list+=("$!")
			nohup $rbh_script -i1 $path/$subject-$query.out -i2 $path/$query-$subject.out -o $path/$subject-$query.orths $6 &>> logs/twoway_RBH.o&
			proc_list+=("$!")
			for proc_id in "${proc_list[@]}"
			do
				echo $proc_id
				wait $proc_id #|| true
			done
		fi
	fi
}

function merge_OG2genes_OrthoDB(){
	# set -eo pipefail
	local script_args=($(echo $@))
	local ORTHODB_PATH_PREFIX=$(realpath ${script_args[0]}) #$1 #$(grep -i -w "orthodb_path_prefix" $1 | check_param) 
	local ORTHODB_PATH=$(dirname "$ORTHODB_PATH_PREFIX")
	local CLEAN_EXTRACT=${script_args[1]} #$2 #$(grep -i -w "clean_extract" $1 | check_param) 
	local n_threads=${script_args[2]} #$3 #$(grep -i -w "max_concurrent_jobs" $1 | check_param)
	local gene_list=${script_args[3]} #$4
	
	# echo "pwd:"$(pwd)
	if [[ $(ls -1 "$ORTHODB_PATH_PREFIX"*.gz | awk "END{print NR}") == 0 ]] ; then
		>&2 echo $(color_FG_BG_Bold $Red $BG_White "Error : ODB Path not found!") #| tee >(cat >&2)
	  exit 1
	fi

	# if [[ $CLEAN_EXTRACT == "TRUE" ]]; then
	# #	rm -f "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz
	# 	>&1  echo $(color_FG_BG_Bold $Yellow $BG_White "Removing "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz! due to CLEAN_EXTRACT option in parameters")
	# 	rm -f "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
	# fi

	if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz &> /dev/null ; then
		time join -t $'\t' -1 2 -2 1 <(zcat "$ORTHODB_PATH_PREFIX"_OG2genes.tab.gz | sort --parallel=$n_threads -k 2) <(zcat "$ORTHODB_PATH_PREFIX"_genes.tab.gz | sort --parallel=$n_threads -k 1) | awk -F "\t" '{if($2 in a)a[$2]=a[$2]","$1"||"$5;else a[$2]=$1"||"$5 ;}END{for(key in a)print key"\t"a[key];}' |  awk '{split($2,a,","); delete c; for(key in a){split(a[key],b,/\|\|/); if(b[2] in c==0){c[b[2]]=0;} } e=""; for(gene in c){ if(length(e)==0){e=gene;}else{e=gene","e;} } print $1"\t"$2"\t"e }' | gzip -c > "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz 
	else
		>&1 echo $(color_FG $Green "File exists!")	
	fi
	>&1 echo $(color_FG $Green "Fixed ODB file stored in : ")$(color_FG_BG_Bold $White $BG_Purple "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz)

	if [ ${#script_args[@]} -eq 4 ]; then
		#   echo "here 2"
		if [[ -s "$ORTHODB_PATH"/odb.gene_list ]]; then
			#if old list is different from new list, recreate "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
			if [[ ! -z $(diff -q files/genelist.txt ./files/reference_ORGS.txt) ]]; then
				>&1 echo $(color_FG $Green "New gene list detected: Recreating ")$(color_FG_BG_Bold $White $BG_Purple "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz)
				#gene lists are different, recalculate
				zcat -f $gene_list | sort | uniq | grep -v -w -i "gene" | grep -v '^$' > "$ORTHODB_PATH"/odb.gene_list
				cp "$ORTHODB_PATH"/odb.gene_list $gene_list
				time zgrep -f "$ORTHODB_PATH"/odb.gene_list "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz | gzip -c > "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
			fi
			>&1 echo $(color_FG $Green "(User gene list) Fixed ODB file found in : ")$(color_FG_BG_Bold $White $BG_Purple "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz)
		else
			if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz &> /dev/null ; then
				# echo "here 2.1"
				zcat -f $gene_list | sort | uniq | grep -v -w -i "gene" | grep -v '^$' > "$ORTHODB_PATH"/odb.gene_list	
				cp "$ORTHODB_PATH"/odb.gene_list $gene_list
				time zgrep -f "$ORTHODB_PATH"/odb.gene_list "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz | gzip -c > "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
			else
				>&1  echo $(color_FG $Green "File exists!")
			fi
			>&1 echo $(color_FG $Green "(User gene list) Fixed ODB file stored in : ")$(color_FG_BG_Bold $White $BG_Purple "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz)
		fi
	fi
}

function extract_transcript_regions(){
	# $1 = genome fasta file
	# $2 = formatted gtf/gff3(UNTESTED) file
	# $3 = input gene list
	# $4 = org name(eg x_laevis)
	# $5 = Path to parameters file

	local script_args=($(echo $@))
	local GENOME_FILE=${script_args[0]} #$1
	local ANNO_FILE=${script_args[1]} #$2
	local GENE_LIST=${script_args[2]} #$3
	local f_org_name=${script_args[3]} #$4
	local param_file=${script_args[4]} #$5
	local PARALLEL_PATH=${script_args[5]} #$5
	local org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0); print $1" "$2; } ;')

	#chmod a+x $GENOME_FILE
	#chmod a+x $ANNO_FILE
	if [[ ! -s $GENOME_FILE || ! -s $ANNO_FILE ]]; then
		echo $(color_FG_BG_Bold $Red $BG_White "Error : Genome or Annotation missing!.") #| tee >(cat >&2)
		exit 1
	fi


	local GENOMES_PATH=$(grep -i -w "genomes_path" $param_file | check_param)
	local ANNOS_PATH=$(grep -i -w "annos_path" $param_file | check_param)
	local BED_PATH=$(grep -i -w "bed_path" $param_file | check_param)
	local FASTA_PATH=$(grep -i -w "fasta_path" $param_file | check_param)
	local OUT_PATH=$(grep -i -w "out_path" $param_file | check_param)
	local TEMP_PATH=$(grep -i -w "temp_path" $param_file | check_param)
	local CLEAN_EXTRACT=$(grep -i -w "clean_extract" $param_file | check_param) 
	local TRANSCRIPT_REGIONS=($(grep -i -w "transcript_regions" $param_file | check_param | sed "s/,/\n/g" ))
	local GENE_SEARCH_MODE=$(grep -i -w "gene_search_mode" $param_file | awk -F"==" '{print $2}') # | check_param)
	local ORTHODB_PATH_PREFIX=$(grep -i -w "orthodb_path_prefix" $param_file | check_param)
	local REF_ORGS=$(grep -i -w "ref_orgs" $param_file | check_param)
	local seqID_delimiter=$(grep -i -w "seqID_delimiter" $param_file | check_param) 
	#LABEL_SEQS=$(grep -i -w "label_sequence_IDs" $5 | check_param) 
	#n_threads=$(nproc --all)
	local n_threads=$(grep -i -w "max_concurrent_jobs" $param_file | check_param)
	if [[ -z $n_threads || $n_threads == 0 || $n_threads == " " ]]; then
		n_threads=$(nproc)
	fi
	local bed_prefix="$BED_PATH/$f_org_name"
	local MODE=""
	if [[ $GENE_SEARCH_MODE=="HARD" || -z $GENE_SEARCH_MODE || $GENE_SEARCH_MODE == " "  ]]; then
	  local MODE="-w"
	fi

	echo "Command : $0 extract_transcript_regions $@"

	if [[ $CLEAN_EXTRACT ==  "TRUE" ]] ; then
	  rm -rf $FASTA_PATH/$f_org_name
	  rm -rf $OUT_PATH/genes/$f_org_name
	  rm -rf $TEMP_PATH/$f_org_name
	  rm -f $GENOME_FILE.fai
	  rm -f $OUT_PATH/genes/$f_org_name/*
	  rm -f $bed_prefix/*
	fi

	mkdir -p $FASTA_PATH/$f_org_name
	mkdir -p $OUT_PATH/genes/$f_org_name/
	mkdir -p $TEMP_PATH/$f_org_name
	mkdir -p $bed_prefix

	#if [[ $(echo $ANNO_FILE | grep -q -i "gtf") != 0 ]] ; then
	if ! basename $ANNO_FILE | grep -q -i "gtf" ; then 
	  local file_name=${ANNO_FILE%.*}
	  #zcat -f $ANNO_FILE | gffread - -T -O -E -o - | gzip -c > $ANNOS_PATH/"$f_org_name".gtf.gz & ## -O ONLY FOR GFF3
	  zcat -f $ANNO_FILE | gffread - -T -E -o - | gzip -c > $ANNOS_PATH/"$f_org_name".gtf.gz &
	  local anno_proc_id=$(echo $!)
	  local ANNO_FILE=$ANNOS_PATH/"$f_org_name".gtf.gz
	fi

	# while [[ $(lsof -R -p $$ | wc -l) -gt $(($(ulimit -Sn)-200)) ]]; #limit file descriptors for child & parent process to 1024-200 
	# do
	#   printf %s"\t"%s"\n" "Waiting for file handles to close.." $(lsof -R -p $$ | wc -l)
	#   sleep 5
	# done

  # Get the soft limit once to avoid repeating the subshell
  limit=$(ulimit -Sn)
  threshold=$((limit - 200))
  
  # Use /proc/$$/fd to count open handles
  while [[ $(ls /proc/$$/fd | wc -l) -gt $threshold ]]; 
  do
      current_count=$(ls /proc/$$/fd | wc -l)
      printf "Waiting for file handles to close... Current: %s / Max: %s\n" "$current_count" "$limit"
      sleep 5
  done

	>&1 color_FG_BG_Bold $Black $BG_Yellow "Extracting Genome & Building Index (Samtools faidx)..."

	if [[ ${GENOME_FILE##*.} == "gz" ]] ; then
	  local gfile_name=${GENOME_FILE%.*}
	else
	  local gfile_name=$GENOME_FILE
	fi

	index_genome $GENOME_FILE &
	local genome_index_proc=$(echo $!)

	>&1 color_FG $Yellow "Genome : $gfile_name\nAnnotation : $ANNO_FILE"

	if [[ $(zgrep -m 1 -hPo 'gene_name "\K[^"]+' $ANNO_FILE | awk '{print NR}') == 0 ]] ; then
		echo $(color_FG_BG_Bold $Red $BG_White "Error : GTF does not contain gene_name attribute") #| tee >(cat >&2)
		exit 1
	fi

	###################################################################################################

	local gene_list=($(cat $GENE_LIST | sort | uniq | grep -v -w -i "gene" | grep -v '^$'))

	>&1 color_FG_BG_Bold $Black $BG_Yellow "1. Checking Gene Names & Splitting GTFs for parallel processing..."

	touch $OUT_PATH/genes/$f_org_name/1.list $OUT_PATH/genes/$f_org_name/2.list $OUT_PATH/genes/$f_org_name/odb.list

	if [[ ! -z $anno_proc_id ]]; then
	  wait $anno_proc_id
	fi

	if [[ ! -s $OUT_PATH/genes/$f_org_name/1.list || ! -s $OUT_PATH/genes/$f_org_name/2.list || ! -s $OUT_PATH/genes/$f_org_name/gtf_stats.csv || ! -s $BED_PATH/$f_org_name/$f_org_name.bed ]] ; then
	  local eexp_gene=$(printf -- '%s\n' "${gene_list[@]}")
	  
	  if [[ ${#eexp_gene[@]} > 0 ]]; then
	  	#Splitting GTF into multiple parts based on grep output for downstream parallel processing in extract_gtf_info.R
	  time zgrep -i $MODE -A 0 --group-separator='>' -f <(echo "${eexp_gene[@]}") $ANNO_FILE | csplit --quiet -z --suffix-format="%0d.gtf_slice" --prefix="$TEMP_PATH/$f_org_name/1." --suppress-matched - '/>/' '{*}' #> $TEMP_PATH/$f_org_name/1.gtf_slice 
	  fi
	  
	  zgrep -hPo 'gene_name "\K[^"]+' $TEMP_PATH/$f_org_name/*.gtf_slice | sort | uniq | awk 'NF' | awk '{print tolower($0)}' > $OUT_PATH/genes/$f_org_name/1.list
	  printf -- "%s\n" ${gene_list[@]/($(cat $OUT_PATH/genes/$f_org_name/1.list)))} | awk 'NF' | awk '{print tolower($0)}' > $OUT_PATH/genes/$f_org_name/2.list
	fi

	if [[ -s $OUT_PATH/genes/$f_org_name/1.list || -s $OUT_PATH/genes/$f_org_name/2.list ]]; then
	  >&1 echo $(color_FG $Green "1. DONE : Available Genes : ")$(color_FG_BG_Bold $White $BG_Purple "$OUT_PATH/genes/$f_org_name/1.list")$(color_FG $Green ", Genes not found in annotation: ")$(color_FG_BG_Bold $White $BG_Purple "$OUT_PATH/genes/$f_org_name/2.list ")
	else
	  echo $(color_FG_BG_Bold $Red $BG_White "1. Error : Step 1 Failed (Check if the GTF file exists and if the size is right)") #| tee >(cat >&2)
	  exit 1
	fi

	#######################################################################################################

	>&1 color_FG_BG_Bold $Black $BG_Yellow "2. Checking OrthoDB for missing & orthologous genes..."

	if [[ ! -s $OUT_PATH/genes/$f_org_name/odb.list || ! -s $OUT_PATH/genes/$f_org_name/final.list || ! -s $OUT_PATH/genes/$f_org_name/odb.final_map ]] ; then
	  time check_OrthoDB $f_org_name $GENE_LIST $OUT_PATH/genes/$f_org_name/odb.list $OUT_PATH/genes/$f_org_name/odb.final_map $param_file
	fi

	if [[ $(awk 'END{print NR;}' "$OUT_PATH/genes/$f_org_name/odb.list" | awk '{print $1}') !=  0 ]] ; then
	    local odb_gene_list=($(cat "$OUT_PATH/genes/$f_org_name/odb.list" | grep -v -w -i "gene" | grep -v '^$' | awk '{print tolower($0)}'))
	    local existing_list=($(cat $OUT_PATH/genes/$f_org_name/1.list | sort | uniq | awk '{print tolower($0)}')) 
	    local short_list=($(echo ${odb_gene_list[@]/${existing_list[@]}}))

	  if [[ "${#short_list[@]}" > 0 ]] ; then
	    local eexp_gene=$(printf -- '%s\n' "${short_list[@]}")
	    time echo "${eexp_gene[@]}" | zgrep -i $MODE -A 0 --group-separator='>' -f - $ANNO_FILE | csplit --quiet -z --suffix-format="%0d.gtf_slice" --prefix="$TEMP_PATH/$f_org_name/2." --suppress-matched - '/>/' '{*}' #> $TEMP_PATH/$f_org_name/2.gtf_slice 
	  fi

	  ##REFRESH geene list based on file names
	  cat $OUT_PATH/genes/$f_org_name/1.list $OUT_PATH/genes/$f_org_name/2.list $OUT_PATH/genes/$f_org_name/odb.list | awk 'NF' | awk '{print tolower($0)}' > $OUT_PATH/genes/$f_org_name/full.list

	  >&1 echo $(color_FG $Green "2. DONE : Full List : ")$(color_FG_BG_Bold $White $BG_Purple "$OUT_PATH/genes/$f_org_name/full.list")$(color_FG $Green ", List from ODB : ")$(color_FG_BG_Bold $White $BG_Purple "$OUT_PATH/genes/$f_org_name/odb.list")$(color_FG $Green ", ODB Cluster to Genes Map : ")$(color_FG_BG_Bold $White $BG_Purple "$OUT_PATH/genes/$f_org_name/odb.final_map")

	else
	  echo $(color_FG_Bold $Red "2. Warning : ")$(color_FG_BG_Bold $White $BG_Red "$OUT_PATH/genes/$f_org_name/odb/odb.list")$(color_FG_Bold $Red " missing, Possibly orthologous genes were not found ") #| tee >(cat >&2)
	  color_FG_Bold $Red "2. If unsure, re-run command: $0 check_OrthoDB $f_org_name $GENE_LIST $ANNO_FILE" #| tee >(cat >&2)
	  cat $OUT_PATH/genes/$f_org_name/1.list $OUT_PATH/genes/$f_org_name/2.list | awk 'NF' > $OUT_PATH/genes/$f_org_name/full.list
	fi

	#######################################################################################################

	>&1 color_FG_BG_Bold $Black $BG_Yellow "3. Extracting Transcript Stats from GTF_Slices..." #(log:$TEMP_PATH/$f_org_name/get_GTF_info.[o/e])

	time Rscript --vanilla --verbose $(echo $(dirname $0))/extract_gtf_info.R $TEMP_PATH/$f_org_name/ $f_org_name $OUT_PATH/genes/$f_org_name/gtf_stats.csv $param_file #1> $TEMP_PATH/$f_org_name/get_GTF_info.o 2> $TEMP_PATH/$f_org_name/get_GTF_info.e
	r_exit_code="$?"

	if [[ ! -s $OUT_PATH/genes/$f_org_name/gtf_stats.csv || $r_exit_code != 0 ]] ; then #|| ! -s $OUT_PATH/genes/$f_org_name/final.list
	  >&2 color_FG_Bold $Red "3. ERROR: Extraction of transcript stats failed... Possibly no genes were found. Check if GTF file has gene_name attribute"
	  #>&2 color_FG_Bold $Red "3. Check $TEMP_PATH/$f_org_name/get_GTF_info.[o/e]"
	  >&2 color_FG_Bold $Red "3. (Possible Fix) : Remove $TEMP_PATH/$f_org_name/ & $OUT_PATH/genes/$f_org_name/gtf_stats.csv and re-run the pipeline"
	  exit 255
	fi

	sed 1d $OUT_PATH/genes/$f_org_name/gtf_stats.csv | awk -F',' '{print $1"\n"}' | sort | uniq | awk 'NF' | awk '{print tolower($0)}' > $OUT_PATH/genes/$f_org_name/final.list

	if [[ -s $OUT_PATH/genes/$f_org_name/gtf_stats.csv && -s $OUT_PATH/genes/$f_org_name/final.list && $r_exit_code == 0 ]] ; then
	>&1 echo $(color_FG $Green "3. DONE : Final List : ")$(color_FG_BG_Bold $White $BG_Purple "$OUT_PATH/genes/$f_org_name/final.list")$(color_FG $Green ", GTF stats : ")$(color_FG_BG_Bold $White $BG_Purple "$OUT_PATH/genes/$f_org_name/gtf_stats.csv")
	else
	  echo $(color_FG_BG_Bold $Red $BG_White "3. Error : Step 3 Failed") #| tee >(cat >&2)
	  exit 1
	fi

	#######################################################################################################

	local gene_list=($(cat $OUT_PATH/genes/$f_org_name/final.list | awk 'NF' | awk '{print tolower($0)}'))
	local s_names=($(awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo ${gene_list[@]}) | awk '{print tolower($0)}'))

	if [[ ! -z $genome_index_proc ]]; then
	  wait $genome_index_proc
	fi

	>&1 color_FG_BG_Bold $Black $BG_Yellow "4. Fetching sequences from Genome..."
	>&1 color_FG $Yellow "Transcript Regions : $(IFS=","; echo ${TRANSCRIPT_REGIONS[*]}) "

	if [[ -s $gfile_name ]]; then

	  time $PARALLEL_PATH --max-procs $n_threads "get_FASTA {1} {2} $bed_prefix/$f_org_name $gfile_name $f_org_name $FASTA_PATH/$f_org_name $ANNO_FILE $TEMP_PATH/$f_org_name $OUT_PATH ;" ::: ${gene_list[@]} ::: ${TRANSCRIPT_REGIONS[@]}

	  >&1 echo $(color_FG $Green "4. DONE : FASTA PATH : ")$(color_FG_BG_Bold $White $BG_Purple "$FASTA_PATH/$f_org_name")
	else
	  echo $(color_FG_BG_Bold $Red $BG_White "4. Error : Step 4 Failed, Genome not found!") #| tee >(cat >&2)
	  exit 1
	fi

	#############################################################################################################

	#if [[ $LABEL_SEQS ==  "TRUE" && -s $OUT_PATH/genes/$f_org_name/odb.final_map ]] ; then
	  >&1 color_FG_BG_Bold $Black $BG_Yellow "4.1 Labelling sequences..."
	  #local tmp_names=($(parallel --link --max-procs $n_threads "echo {1},{2}" ::: ${gene_list[@]} ::: ${s_names[@]}))
	  #time parallel --max-procs $n_threads "printf -- %s,%s\\\n {1} {2}" ::: ${tmp_names[@]} ::: ${TRANSCRIPT_REGIONS[@]} | parallel --colsep "," --max-procs $n_threads "label_sequenceIDs $f_org_name {1} $FASTA_PATH/$f_org_name/{2}.{3} $FASTA_PATH/$f_org_name/{2}.{3}.tmp $param_file $OUT_PATH/genes/$f_org_name/odb.final_map" 
	  time Rscript --vanilla --verbose $(echo $(dirname $0))/label_FASTA_files.R $FASTA_PATH/$f_org_name/ $f_org_name $OUT_PATH/genes/$f_org_name/final.list $param_file $OUT_PATH/genes/$f_org_name/odb.final_map
	#fi

	#######################################################################################################

	>&1 color_FG_BG_Bold $Black $BG_Yellow "5. Generating Metadata and Cleaning up..."

	find $FASTA_PATH/$f_org_name/ -type f  -empty -delete

	#sed 1d $OUT_PATH/genes/$f_org_name/gtf_stats.csv | awk -F',' '{print $1}' | sort | uniq >  $OUT_PATH/genes/$f_org_name/AVAILABLE_GENES 
	ls -1 $FASTA_PATH/$f_org_name/ | sed -e 's/\.[^.]*$//' | sort -u | awk '{print tolower($0)}' >  $OUT_PATH/genes/$f_org_name/AVAILABLE_GENES 
	grep -v -i -f $OUT_PATH/genes/$f_org_name/AVAILABLE_GENES $GENE_LIST | sort | uniq | awk '{print tolower($0)}' > $OUT_PATH/genes/$f_org_name/MISSING_GENES

	find $FASTA_PATH/$f_org_name/ -type f -name "*.fai" -exec rm -f {} +

	if [[ ${GENOME_FILE##*.} == "gz" ]]; then #&& ! -z $gfile_name
	  rm -f $gfile_name
	fi

	if [[ ! -s $GENOMES_PATH/$f_org_name.fa.gz ]]; then
	     zcat -f $GENOME_FILE | gzip -c > $GENOMES_PATH/$f_org_name.fa.gz
	     >&1 echo $(color_FG $Green "Genome saved to : ")$(color_FG_BG_Bold $White $BG_Purple "$GENOMES_PATH/$f_org_name.fa.gz")
	     rm -f $GENOME_FILE
	fi
	if [[ ! -s $ANNOS_PATH/$f_org_name.gtf.gz ]]; then
	     zcat -f $ANNO_FILE | gzip -c > $ANNOS_PATH/$f_org_name.gtf.gz
	     >&1 echo $(color_FG $Green "Annotation saved to : ")$(color_FG_BG_Bold $White $BG_Purple "$ANNOS_PATH/$f_org_name.gtf.gz")
	     rm -f $ANNO_FILE
	fi

	rm -rf $TEMP_PATH/$f_org_name/
	rm -f $bed_prefix/"$f_org_name"_*

	>&1 color_FG_BG_Bold $Purple $BG_White "Extraction DONE for organism : $f_org_name"

	exit 0
}

function check_OrthoDB(){
	# # INTERNAL SCRIPT FOR finding orthologous genes
	# $1 - reference organism
	# $2 - gene list (to look for orthologous genes)
	# $3 - output file for gene list from OrthoDB
	# $4 - output file for gene clusters from OrthoDB
	# $5 - Path to parameters.txt
	# $6 - TRUE/FALSE. Select all genes from each cluster for the organism? (otherwise only the user specified genes are selected)

	##ENTRYPOINT
	# set -xeo pipefail
	local script_args=($(echo $@))
	local f_org_name=${script_args[0]}  #$1
	#readarray gene_list < $2
	local gene_list=($(cat ${script_args[1]} | grep -v -w -i "gene" | grep -v '^$')) #$2
	local out_gene_list=${script_args[2]} #$3
	local out_gene_clusters=${script_args[3]} #$4 #$5 - is parameters.txt
	local param_file=${script_args[4]}
	local select_all_genes=${script_args[5]}

	local TEMP_PATH=$(grep -i -w "temp_path" $param_file | check_param)
	local ORTHODB_PATH_PREFIX=$(realpath $(grep -i -w "orthodb_path_prefix" $param_file | check_param))
	local REF_ORGS=$(grep -i -w "ref_orgs" $param_file | check_param)
	local GENE_SEARCH_MODE=$(grep -i -w "gene_search_mode" $param_file | awk -F"==" '{print $2}' ) #| check_param)
	local n_threads=$(grep -i -w "max_concurrent_jobs" $param_file | check_param)
	
	# echo $ORTHODB_PATH_PREFIX

	if [[ -z $ORTHODB_PATH_PREFIX ]]; then
		>&2 color_FG_Bold $Yellow "Not using OrthoDB!"
		return 0
	fi
	
	if [[ -z $n_threads || $n_threads == 0 || $n_threads == " " ]]; then
		n_threads=$(nproc)
	fi

	local MODE=""
	if [[ $GENE_SEARCH_MODE=="HARD" || -z $GENE_SEARCH_MODE || $GENE_SEARCH_MODE == " " ]]; then
	  local MODE="-w"
	fi

	mkdir -p $TEMP_PATH/$f_org_name/

	#Check if ODB files exist
	if [[ ( ! -s "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz || ! -s "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz ) && ! -s "$ORTHODB_PATH_PREFIX"_species.tab.gz ]] ; then
	  >&2 color_FG_Bold $Red "2. OrthoDB files missing/corrupt"
	  >&2 color_FG_Bold $Red "2. Require: "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz, "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz, "$ORTHODB_PATH_PREFIX"_species.tab.gz"
	  exit -1
	fi

	# if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz ; then
	#   if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz ; then
	#     >&2 color_FG_Bold $Red "2. Error in ODB files! Missing "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz : Rerun merge_OG2genes_OrthoDB($ORTHODB_PATH_PREFIX FALSE $n_threads ${script_args[1]})"
	# 	## merge_OG2genes_OrthoDB $ORTHODB_PATH_PREFIX FALSE $n_threads ${script_args[1]}
	# # 	exit -1
	# #   else
	# #     local ODB_FILE="$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz
	# #     >&1 echo $(color_FG $Yellow "2. Selected Fixed ODB file : ")$(color_FG_BG_Bold $White $BG_Purple "$ODB_FILE")
	#   fi
	#   >&2 color_FG_Bold $Red "2. Error in ODB files! Corrupt "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz : Rerun merge_OG2genes_OrthoDB($ORTHODB_PATH_PREFIX FALSE $n_threads ${script_args[1]})"
	#   exit -1
	# else
	#   local ODB_FILE="$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
	#   >&1 echo $(color_FG $Yellow "2. Selected Fixed ODB (User) file : ")$(color_FG_BG_Bold $White $BG_Purple "$ODB_FILE")
	# fi

	if [ -s "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz ]; then
		local ODB_FILE="$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
		>&1 echo $(color_FG $Yellow "2. Selected Fixed ODB (User) file : ")$(color_FG_BG_Bold $White $BG_Purple "$ODB_FILE")
	else
		local ODB_FILE="$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz
		>&1 echo $(color_FG $Yellow "2. Selected Fixed ODB file : ")$(color_FG_BG_Bold $White $BG_Purple "$ODB_FILE")
	fi
	
	readarray refs < $REF_ORGS

	local s_names=($(awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo "${gene_list[@]}") | awk '{print tolower($0)}'))
	local genes_strip=($(awk -F'_' '{print $1;}' <(echo "${s_names[@]}")))
	local lookup_genes=($(echo "${genes_strip[@]}" "${gene_list[@]}" | sort | uniq | awk 'NF' | awk '{print tolower($0)}' ))
	local org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;' | awk -F" " '{print $1" "$2}')

	#Find ODB ORGANISM ID for the organism
	local org_ID=$(zgrep -i -P "\t\b($org_name)\b\t" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
	if [[ -z $org_ID ]] ; then
	  local org_ID=$(zgrep -i -P "($org_name)" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
	fi

	# echo "ODB_FILE:"$ODB_FILE "out file:"$out_gene_list "GENES:" ${gene_list[@]}

	if [[ ! -z $org_ID ]] ; then

	  #Find ODB ORGANISM ID for the reference organism(s)
	  for ref in "${refs[@]}"; do 
	    local ref_org_name=$(echo $ref | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
	    local ref_org_ID=""
	    local ref_org_ID=$(zgrep -i -P "\t\b($ref_org_name)\b\t" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
	    if [[ -z $ref_org_ID ]] ; then
	      local ref_org_ID=$(zgrep -i -P "($ref_org_name)" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
	    fi
	    if [[ ! -z $ref_org_ID && $ref_org_ID != $f_org_name ]] ; then
	      ref_org_IDs+=($ref_org_ID)
	    fi
	  done

	#Get only the ODB GENE IDs, genes and ODB clusters of the reference organisms from the reduced ODB FILE which is generated by merge_OG2genes_OrthoDB.sh
	zgrep -f <(printf -- '%s\n' "${ref_org_IDs[@]}") $ODB_FILE | gzip -c > $TEMP_PATH/$f_org_name/odb.clusters.gz

	if [[ $select_all_genes == "TRUE" ]]; then
		#ODB GENE IDs and GENE NAMES are delimited with || which I am splitting with awk regexp \|/|, and selecting the elements which match the organisms and genes
		zcat -f $TEMP_PATH/$f_org_name/odb.clusters.gz | awk '{print $2}' | sed 's/,/\n/g' | grep -w "$org_ID" - | awk '{split($0,a,/\|\|/); print a[2]}' | sort -u > $out_gene_list
	else
		#ODB GENE IDs and GENE NAMES are delimited with || which I am splitting with awk regexp \|/|, and selecting the elements which match the organisms and genes
		zcat -f $TEMP_PATH/$f_org_name/odb.clusters.gz | awk '{print $2}' | sed 's/,/\n/g' | grep -w "$org_ID" - | grep -i $MODE -f <(printf -- '%s\n' "${lookup_genes[@]}") - | awk '{split($0,a,/\|\|/); print a[2]}' | sort -u | awk '{print tolower($0)}' > $out_gene_list
	fi

	if [[ -s $out_gene_list ]]; then
	  >&1 echo $(color_FG $Green "2. DONE: Found Orthologous genes : ")$(color_FG_BG_Bold $White $BG_Purple "$out_gene_list")
	else
	  >&2 echo $(color_FG $Red "2. Error: Orthologous genes couldn't be found/saved to $out_gene_list")
	  exit 1
	fi

	#Save selected clusters and the participating genes(delimiter=",") in odb.final_map
	zgrep "$org_ID" $TEMP_PATH/$f_org_name/odb.clusters.gz | awk '{split($3,a,","); for(key in a){print $1"\t"a[key]}}' | grep -i $MODE -f <(cat ${script_args[1]} $out_gene_list | grep -v -w -i "gene" | grep -v '^$') | sort -u | awk '{if($1 in a){a[$1]=a[$1]","$2}else{a[$1]=$2;} } END{for(key in a){print key"\t"a[key]}}' > $out_gene_clusters #> $OUT_PATH/genes/$f_org_name/odb.final_map

	if [[ -s $out_gene_clusters ]]; then #$OUT_PATH/genes/$f_org_name/odb.final_map
	  >&1 echo $(color_FG $Green "2. DONE: Mapped gene names to clusters : ")$(color_FG_BG_Bold $White $BG_Purple "$out_gene_clusters")
	else
	  >&2 echo $(color_FG $Red "2. Error: OG Cluster to gene mappings couldn't be saved to $out_gene_clusters")
	  exit 1
	fi


	rm -f $TEMP_PATH/$f_org_name/odb.clusters.gz

	else
	  >&2 color_FG_Bold $Red "2. $org_name not found in OrthoDB files"
	fi

	#exit 0
}

function cat_files(){
	local script_args=($(echo $@))
	local output_file=${script_args[0]}
	local infiles=($(echo "${script_args[@]: 1:${#script_args[@]}}"))

	cat "${infiles[@]}" > $output_file
return 0
}

function sed_replace(){
	local script_args=($(echo $@))
	local input_file=${script_args[0]}
	local old_name=${script_args[1]}
	local new_name=${script_args[2]}

	sed -i "s/$old_name/$new_name/g" $input_file
}

function convert_BLAST_format(){
	local script_args=($(echo $@))
	local blast_bin_path=${script_args[0]}
	local blast_formatter_path=$(echo "$blast_bin_path/blast_formatter")
	local blast_archive=${script_args[1]}  
	local out_file=${script_args[2]} 
	local out_fmt=${script_args[3]}  
	local blast_cols=$(echo "${script_args[@]: 4:${#script_args[@]}}")

	if [[ -s $blast_archive && -s $blast_formatter_path ]]; then
		$blast_formatter_path -archive $blast_archive -outfmt "$out_fmt $blast_cols" -out $out_file #&> /dev/null
	else
		>&2 echo "blast_formatter not found in $blast_formatter_path"
		exit 1
	fi
}

function select_ref_org_groups(){
	local script_args=($(echo $@))
	local ref_orgs_file=${script_args[0]}
	local groups_path=${script_args[1]}  
	local parallel_path=${script_args[2]}  
	local any_or_all_refs=${script_args[3]}  
	local n_threads=${script_args[4]}

if [[ $any_or_all_refs == "ANY" ]]; then
	rm -f $(grep -r -L -i -f $ref_orgs_file $groups_path)
elif [[ $any_or_all_refs == "ALL" ]]; then
	readarray refs < $ref_orgs_file	
	$parallel_path --max-procs $n_threads "grep -r -l -i {1} $groups_path | sort > $groups_path/{1}.tmp" ::: "${refs[@]}" 
	#rm -f $(grep -r -L -i -f $(cat $groups_path/*.tmp | sort | uniq -d))
	rm -f $(cat $groups_path/*.tmp | grep -v -i -f <(cat $groups_path/*.tmp | sort | uniq -d))
	rm -f $groups_path/*.tmp
fi
}

function get_all_odb_genes(){
	local script_args=($(echo $@))
	local gene_list=${script_args[0]}
	local odb_files_path=${script_args[1]} #($(echo "${script_args[@]: 1:${#script_args[@]}}"))
	#parallel --max-procs 8 --keep-order "grep -o -h -i -f $gene_list {}" ::: "${odb_map_path[@]}"
	grep -h -i -f $gene_list <(cat "$odb_files_path") | awk '{print $2}' | sed 's/,/\n/g'
	return 0
}

function install_parallel(){
	if [[ -z $(which parallel) || $(which parallel | awk '{print NR}') == 0 ]] ; then
		#(wget -O - pi.dk/3 || lynx -source pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3 ) > $(dirname $0)/install.sh
		if [[ ! -z $SHELL ]]; then
			#$SHELL $(dirname $0)/install.sh
			(wget -O - pi.dk/3 || lynx -source pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3 ) | $SHELL
		else
			echo "Could not install GNU parallel..SHELL not set in path or SHELL not bash"
			return 255
		fi
	fi
	return 0
}
# function test_print(){
# 	echo $1
# 	echo $(pwd)
# 	echo $0
# 	echo $(dirname $0)
# }

# mask_stops_3utr() {
# 	#python mask_motifs.py -f $1 -s 3 --pos 1 --mask "N"  -r TRUE --add TRUE -cm "TGA,TAA,TAG" -cmm - -o $1 ##Compare mask is negative, so this will add NNN if stop codon doesnt exist
# 	$PY2_PATH mask_motifs.py -f $1 -s 3 -m "N" -p 1 -cm "TGA,TAA,TAG" -r True -rf 1 -o $1
# }

# mask_stops_cds() {
# 	$PY2_PATH mask_motifs.py -f $1 -s 3 -m "N" -p 0 -cm "TGA,TAA,TAG" -r True -rf 1 -o $1
# 	$PY2_PATH mask_motifs.py -f $1 -s 3 -m "N" -p 0 -cm "TGA,TAA,TAG" -r True -rf 2 -o $1
# 	$PY2_PATH mask_motifs.py -f $1 -s 3 -m "N" -p 0 -cm "TGA,TAA,TAG" -r True -rf 3 -o $1
# }

export -f install_parallel
export -f get_all_odb_genes
export -f cat_files
export -f sed_replace
##DATA EXTRACTION FUNCTIONS
#export -f group_FASTA_clusters
export -f concatenate_FASTA_groups
export -f group_FASTA_seqs
export -f lengthen_fastaIDs
export -f shorten_fastaIDs
export -f index_fastaIDs
export -f index_genome
export -f collate_fasta
export -f select_transcripts
export -f delete_empty_orgs
export -f count_genes4orgs
export -f count_seqs4genes
#export -f index_genome
export -f get_FASTA
export -f get_length_dist
export -f get_count_dist
export -f check_param
export -f merge_OG2genes_OrthoDB
export -f extract_transcript_regions
export -f check_OrthoDB
#export -f label_sequenceIDs
#export -f test_print

##FIND ORTHOLOGS FUNCTIONS
export -f all2allblast
export -f onewayblast
export -f all2all_refblast
export -f oneway_RBH
export -f twoway_RBH
export -f do_BLAST
export -f make_BLAST_DB
export -f convert_BLAST_format
export -f select_ref_org_groups

##FUNCTIONS FOR COLORING CONSOLE OUTPUT
export -f color_FG
export -f color_FG_Bold
export -f color_BG
export -f color_FG_BG
export -f color_FG_BG_Bold

export BLASTDB_LMDB_MAP_SIZE=3000000

#export GNU_PARALLEL=$(which parallel)

if [ $# -gt 0 ] ; then
	script_args=($(echo $@))
	#>&2 echo $@
	#echo "${script_args[0]}" "$(echo ${script_args[@]:1:${#script_args[@]}})" 
	#export -f "${script_args[0]}"
	time "${script_args[0]}" "$(echo ${script_args[@]:1:${#script_args[@]}})" 
	exit 0
fi
