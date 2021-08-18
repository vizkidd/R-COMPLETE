#!/bin/bash
#$ -cwd
#$ -V
#$ -l data
# $ -pe smp 1
#$ -l h_vmem=300G
#$ -l m_mem_free=32G
# $ -l full_node=TRUE
#$ -o "logs/blast.o" 
#$ -e "logs/blast.e"

# $ -l h_vmem=22G
# $ -l mem_free=8G
# $ -pe smp 4

# $1 - table output dir (files/tables)
# $2 - Query FASTA
# $3 - Subject BLAST DB
# $4 - BLAST program
# $5 - Gene name
# $6 - Query name
# $7 - Subject name

(source ./parallel_script.sh || true)

#echo "do_blast.sh"
#echo $1 $2 $3 $4 $5 $6 $7
gene=$(echo $5 | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')
blast_options=$(grep -i -w "blast_options" parameters.txt | awk -F'=' '{print $2}' )
tables_path=$1
query=$2
DB=$3
prog=$4
query_name=$6
subject_name=$7
n_proc=$((1+$(nproc)/8))

#alias tempfile='mktemp'
#export -f tempfile

	if [ -s "$tables_path/$gene.$query_name.$subject_name" ]; then
	rm $tables_path/$gene.$query_name.$subject_name
	fi
	if [ ! -s $(blastdb_path -db $DB) ]; then
	makeblastdb -in "$DB" -dbtype nucl  -hash_index || true
	fi
#while IFS='' read -r line # $line is the org
#do 
#	echo $line - $gene
	#echo "($gene) $query_name-$subject_name Started..."
	##tblastx -query "$fasta_path/$line/$gene.$region" -db "$DB/$gene.$region" -outfmt 6 -word_size 5 > $tables_path/$gene.$line
	##cat "$fasta_path/$line/$gene.$region" | parallel --round-robin --pipe --recstart '>' 'tblastx -db "$DB/$gene.$region" -outfmt 6 -word_size 5 > $tables_path/$gene.$line' > $tables_path/$gene.$line
	##nohup parallel -a files/fasta/danio_rerio/sh2d.cds --round-robin --pipepart --recstart '>' 'blastn -db files/blastdb/sh2d.cds -word_size 9 -outfmt 6'  > sh2d.out & >nohup.out &
	#parallel -a "$fasta_path/$line/$gene.$region" -j1 --load 95% --memfree 4G --retries 3 --compress --tmpdir "files/temp" --pipepart --recstart '>' tblastx -db "$DB/$gene.$region" -word_size 5 -outfmt 6  > $tables_path/$gene.$line  # -j0 -j$n_proc --memfree 4G
	#tblastx -query  "$fasta_path/$line/$gene.$region" -db "$DB/$gene.$region" -word_size 5 -outfmt 6  > $tables_path/$gene.$line  # -j0 -j$n_proc --memfree 4G
	>&2 echo "echo ($gene) $query_name-$subject_name Started..."
#	tmpout=$(mktemp -u $tables_path/$gene.$line.XXXXXXX)
	#input=$(echo `find $fasta_path/$line/ -iname "$gene*"  -type f ! -size 0 | grep -i "$region"`)
#	input=$(echo `find $fasta_path/$line/ -iname "$gene*"  -type f ! -size 0 | grep -w -i "$gene" | grep -i "$region"`)
#	echo $input
	#blast_cmd=(time parallel -a "$input" -j1 --joblog files/temp/$gene --tmpdir "/fast/home/vsuresh/temp" --pipepart --recstart '>' --block $n_proc "$prog -db $DB/$gene.$to_region -num_threads $n_proc -word_size 5 -outfmt 11 -evalue 0.01 -out $tmpout")    # -j0 -j$n_proc --memfree 4G --compress --limit "mem 1" --limit "io 1" --retry-failed --retries 2 
	#blast_cmd=(parallel -a "$input" -j1 --joblog files/JOBLOG --tmpdir "/fast/home/vsuresh/temp" --compress --compress-program "bzip2" --pipepart --recstart '>' --block -1 "$prog -db $DB/$gene.$to_region -num_threads $n_proc -word_size 5 -outfmt 11 -evalue 0.01 -out $tmpout")    # -j0 -j$n_proc --memfree 4G --compress --limit "mem 1" --limit "io 1" --retry-failed --retries 2 --tmpdir "/fast/home/vsuresh/temp" --cat
	#echo "${blast_cmd[@]}"
	#"${blast_cmd[@]}"
	if [ -s "$query" ] && [ -s "$DB" ]; then
	#parallel -a "$input" -j1 --joblog files/JOBLOG --tmpdir "/fast/home/vsuresh/temp" --compress --compress-program "bzip2" --pipepart --recstart '>' --block 1 "$prog -db $DB/$gene.$to_region -word_size 5 -outfmt 11 -evalue 0.01 -out $tmpout" # -num_threads $n_proc
	#parallel -a "$query" -j1 --joblog files/JOBLOG --retry-failed  --compress --pipepart --recstart '>' --block 1 "$prog -db $DB -word_size 5 -outfmt 11 -evalue 0.01 -out $tables_path/$gene.$query_name.$subject_name" # -num_threads $n_proc --tmpdir "/fast/home/vsuresh/temp" --compress-program "bzip2" 
	#cat "$query" | parallel -j1 --joblog files/JOBLOG --compress --pipe --recstart '>' --block 1 "$prog -db $DB -word_size 5 -outfmt 11 -evalue 0.01 -out $tables_path/$gene.$query_name.$subject_name" # -num_threads $n_proc --tmpdir "/fast/home/vsuresh/temp" --compress-program "bzip2" 
	parallel -j1 --joblog files/JOBLOG --compress --pipepart -a "$query" --recstart '>' --block -1 "$prog -db $DB -outfmt 11 $blast_options -out $tables_path/$gene.$query_name.$subject_name " #-word_size 5 -evalue 1e-25
	echo "($gene) $query_name-$subject_name is done"
	else
		echo "($gene) Error: Either ($query)/($DB) is empty/not found"
	fi
#done < "$selected_orgs"

#fi

exit 0
