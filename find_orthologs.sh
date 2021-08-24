#!/bin/bash
#$ -N "orths" 
#$ -cwd 
#$ -o "logs/orths.o" 
#$ -e "logs/orths.e" 
#$ -m ea

##This script also 'selects' the transcripts
##Inputs
# $2 genelist (files/genelist.txt)
##Intermediate inputs
# List of selected organisms from $fasta_path (selected_ORGS.txt)
##Taken from parameters.txt
# Path of FASTA ($fasta_path) 
# Path of BLAST DB ($blastdb_path)
# List of reference organisms (reference_ORGS.txt)
#FROM & TO regions allow for cross-region blasting (old)

#convert_fmt() {
#	blast_formatter -archive $tmpout -outfmt $8 -out $tmpout.dup
#}

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

all2allblast() {
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
	cat $(grep -r -i -l -f $3 $2/$subject | grep ".$7") > $5/$subject.$7
	makeblastdb -in "$5/$subject.$7" -dbtype nucl  -hash_index || true 
fi
while IFS= read -r query
do
if [[ ! -s $4/all2all.$query.$subject ]]; then
	if [[ "$query" != "$subject" ]]; then
		#cat $2/$query/*.$6 > $5/$query.$6
		if [[ ! -s $(blastdb_path -db $5/$query.$6) ]]; then
			touch $5/$query.$6
			cat $(grep -r -i -l -f $3 $2/$query | grep ".$6") > $5/$query.$6
			makeblastdb -in "$5/$query.$6" -dbtype nucl  -hash_index || true 
		fi
		printf "$query\n$subject" > files/pairwise_ORGS.txt
		nohup ./jobhold.sh $query-$subject" do_blast.sh $4 $5/$query.$6 $5/$subject.$7 $8 all2all $query $subject" &>> logs/job_status.o& #n_proc=2 #/fast/gridengine/uge-8.6.14/bin/lx-amd64/
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

onewayblast() {
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
  nohup ./jobhold.sh $j_name" do_blast.sh $4 $2/$species/$j_name.$6 $5/$j_name.$7 $8 $j_name $species all" &>> logs/job_status.o& #n_proc=2 #/fast/gridengine/uge-8.6.14/bin/lx-amd64/
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
	wait $proc_id || true
done
echo "Submitted jobs complete(one-way $8)"

}

all2all_refblast() {
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
	cat $(grep -r -i -l -f $3 $2/$subject | grep ".$7") > $5/$subject.$7
	makeblastdb -in "$5/$subject.$7" -dbtype nucl  -hash_index || true 
fi
while IFS= read -r query
do
if [[ ! -s $4/all2all.$query.$subject ]] || [[ ! -s $4/all2all.$subject.$query ]]; then
	if [[ "$query" != "$subject" ]]; then
		#cat $2/$query/*.$6 > $5/$query.$6
		if [[ ! -s $(blastdb_path -db $5/$query.$6) ]]; then
			touch $5/$query.$6
			cat $(grep -r -i -l -f $3 $2/$query | grep ".$6") > $5/$query.$6
			makeblastdb -in "$5/$query.$6" -dbtype nucl  -hash_index || true 
		fi
		printf "$query\n$subject" > files/pairwise_ORGS.txt
		nohup ./jobhold.sh $query-$subject" do_blast.sh $4 $5/$query.$6 $5/$subject.$7 $8 all2all $query $subject" &>> logs/job_status.o& #n_proc=2 #/fast/gridengine/uge-8.6.14/bin/lx-amd64/
		proc_list+=("$!")
		nohup ./jobhold.sh $subject-$query" do_blast.sh $4 $5/$subject.$7 $5/$query.$6 $8 all2all $subject $query" &>> logs/job_status.o& #n_proc=2 #/fast/gridengine/uge-8.6.14/bin/lx-amd64/
		proc_list+=("$!")
	fi
fi
done < "$9"
done < "$1"

for proc_id in "${proc_list[@]}"
do
	echo $proc_id
	wait $proc_id || true
done
echo "Submitted jobs complete(all-all(ref) $8)"
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
oneway_RBH() {
# $1 - FOLDER PATH
# $2 - Query organism
# $3 - Subject organism
path=$1
query=$2
subject=$3
proc_list=()
if [[ "$query" != "$subject" ]]; then
	if [[ ! -s $path/$query-$subject.orths ]]; then
		nohup blast_formatter -archive $path/all2all.$query.$subject -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out $path/$query-$subject.out &>> /dev/null&
		proc_list+=("$!")
		nohup blast_formatter -archive $path/all2all.$subject.$query -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out $path/$subject-$query.out &>> /dev/null&
		proc_list+=("$!")
		for proc_id in "${proc_list[@]}"
		do
			echo $proc_id
			wait $proc_id || true
		done
		#/data/meyer/viz/tools/miniconda3/envs/local_root/bin/python transcriptologs.py -i1 $query-$subject.out -i2 $subject-$query.out -o $query-$subject.orths
		nohup python RBH-v1.py $path/$query-$subject.out $path/$subject-$query.out $path/$query-$subject.orths &>> logs/oneway_RBH.o&
		wait "$!"
	fi
fi
}
twoway_RBH() {
# $1 - FOLDER PATH
# $2 - Query organism
# $3 - Subject organism
path=$1
query=$2
subject=$3
proc_list=()
if [[ "$query" != "$subject" ]]; then
	if [[ ! -s $path/$query-$subject.orths ]] || [[ ! -s $path/$subject-$query.orths ]]; then
		nohup blast_formatter -archive $path/all2all.$query.$subject -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out $path/$query-$subject.out &>> /dev/null&
		proc_list+=("$!")
		nohup blast_formatter -archive $path/all2all.$subject.$query -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out $path/$subject-$query.out &>> /dev/null&
		proc_list+=("$!")
		for proc_id in "${proc_list[@]}"
		do
			echo $proc_id
			wait $proc_id || true
		done
		proc_list=()
		#/data/meyer/viz/tools/miniconda3/envs/local_root/bin/python transcriptologs.py -i1 $query-$subject.out -i2 $subject-$query.out -o $query-$subject.orths
		nohup python RBH-v1.py $path/$query-$subject.out $path/$subject-$query.out $path/$query-$subject.orths &>> logs/twoway_RBH.o&
		proc_list+=("$!")
		nohup python RBH-v1.py $path/$subject-$query.out $path/$query-$subject.out $path/$subject-$query.orths &>> logs/twoway_RBH.o&
		proc_list+=("$!")
		for proc_id in "${proc_list[@]}"
		do
			echo $proc_id
			wait $proc_id || true
		done
	fi
fi
}

##ENTRYPOINT

#mkdir $4
#mkdir $7
rm files/all_ENSDART.txt
rm files/JOBLOG
touch files/gene_thresholds.txt
#export -f do_blast ##Need for function to work with sem
#export -f select_transcripts

#n_proc=$(($(nproc))) #/2))
#echo $n_proc

#gene_list=()


#orth_path=$5
selected_orgs=$1
fasta_path=$(grep -i "fasta_path" parameters.txt | awk -F'=' '{print $2}')
tables_path="files/oneway"
blastdb_path=$(grep -i "blastdb_path" parameters.txt | awk -F'=' '{print $2}')
region=$(grep -i "blast_region" parameters.txt | awk -F'=' '{print $2}')  #"cds"
reference_ORGS=$(grep -i "ref_orgs" parameters.txt | awk -F'=' '{print $2}')  
clean_download=$(grep -i "clean_download" parameters.txt | awk -F'=' '{print $2}') 
clean_extract=$(grep -i "clean_extract" parameters.txt | awk -F'=' '{print $2}') 
seqID_delimiter=$(grep -i "seqID_delimiter" parameters.txt | awk -F'=' '{print $2}') 
e_value=$(grep -i "e_value" parameters.txt | awk -F'=' '{print $2}')

if [[ $clean_download == "TRUE" ]] || [[ $clean_extract == "TRUE" ]] ; then 
	rm -rf files/oneway
	rm -rf files/all2all
	rm -rf files/all2all_final
	rm files/gene_thresholds.txt
fi

rm -rf $blastdb_path
mkdir $blastdb_path

if [[ $(wc -l $reference_ORGS) > 2 ]]; then

all2allblast $reference_ORGS $fasta_path $2 files/all2all $blastdb_path $region $region tblastx

while IFS= read -r subject
do
while IFS= read -r query
do
		#if [[ ! -s files/all2all/$query-$subject.orths ]]; then
		#	blast_formatter -archive files/all2all/all2all.$query.$subject -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out files/all2all/$query-$subject.out
		#	blast_formatter -archive files/all2all/all2all.$subject.$query -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out files/all2all/$subject-$query.out
		#	#/data/meyer/viz/tools/miniconda3/envs/local_root/bin/python transcriptologs.py -i1 $query-$subject.out -i2 $subject-$query.out -o $query-$subject.orths
		#	python RBH-v1.py files/all2all/$query-$subject.out files/all2all/$subject-$query.out files/all2all/$query-$subject.orths
		#fi
		oneway_RBH files/all2all $query $subject
done < "$reference_ORGS"
done < "$reference_ORGS"

Rscript merge_orths.R $reference_ORGS files/all2all/ files/all2all/final_df.txt $e_value

awk '{for(i=1;i<=2;i++) print $i}' files/all2all/final_df.txt > files/all2all/all2all.selected

select_transcripts files/all2all/all2all.selected files/temp/all2all.selected files/all_ENSDART.txt $fasta_path/ #$(echo "${proc_list[@]}" | sed 's/ //g')

awk -F$seqID_delimiter '{print $2}' files/all2all/all2all.selected | sort | uniq > files/all2all/all2all.genelist

else
	mkdir files/all2all/
	cat $2 > files/all2all/all2all.genelist
fi

onewayblast $reference_ORGS $fasta_path $2 files/oneway $blastdb_path $region $region tblastx

while IFS= read -r gene
do
if grep -q $(echo $gene | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;') "files/all2all/all2all.genelist"; then
	rm $tables_path/$gene
	for f in $(ls $tables_path/$gene.*.*); 
	do 
  	blast_formatter -archive $f -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" >> files/oneway/$gene
  	##rm f
	done;

	
		if [ -s "files/gene_thresholds.txt" ]; then
			#if grep -q "$gene" "files/gene_thresholds.txt"; then
			if [[ $(grep -q "$gene" files/gene_thresholds.txt) ]] || [[ $(grep "$gene" "files/gene_thresholds.txt" | awk -F',' '{print $2}') == "NA" ]]; then
				echo "Score already calculated for $gene"
				#exit 0
			else
				#/gnu/store/71wrbjz4xx53mhd40qfzm5czfzsc89sx-profile/bin/Rscript calculate_gene_conservation.R $tables_path/$gene $gene "files/idents" "files/gene_thresholds.txt" "files/plots" "files/available_orgs.txt" "files/selected_ORGS.txt" "files/orgs" $(grep ">" $blastdb_path/$gene.$region | wc -l)
				Rscript calculate_gene_conservation.R $tables_path/$gene $gene "files/idents" "files/gene_thresholds.txt" "files/available_orgs.txt" "files/orgs" "$(grep ">" $blastdb_path/$gene.$region | wc -l)"
			fi
		else
			Rscript calculate_gene_conservation.R $tables_path/$gene $gene "files/idents" "files/gene_thresholds.txt" "files/available_orgs.txt" "files/orgs" "$(grep ">" $blastdb_path/$gene.$region | wc -l)"
		fi

	thres=$(($(grep -i $gene files/gene_thresholds.txt | awk -F',' '{print $2}'))) #-5)) ##Subtracting threshhold by 5 for padding
	echo "$gene - $thres"
	#/gnu/store/71wrbjz4xx53mhd40qfzm5czfzsc89sx-profile/bin/Rscript select_orthologs.R $gene $tables_path/$gene $orth_path $thres ##/gnu/store/71wrbjz4xx53mhd40qfzm5czfzsc89sx-profile/bin/Rscript
	#awk -F'_' '{print $1}' $orth_path/"$gene".orths > $orth_path/"$gene".temp
	#cat $orth_path/"$gene".temp > $orth_path/"$gene".orths
	#rm $orth_path/"$gene".temp
		#'Select' transcripts, thereby deleting all the irrelavnt (non-orthologous) sequences
	#FOR each gene, we select orthologous transcripts of all organisms. (the tables with the orth transcript list is in files/tables)
	#select_transcripts $7/"$gene".orths files/temp/$gene.selected files/all_ENSDART.txt $fasta_path/ #$(echo "${proc_list[@]}" | sed 's/ //g')
fi
done < $2

Rscript create_sets.R $2 "files/all2all/all2all.genelist"

cat $reference_ORGS files/oneway/SET > files/oneway/set.tmp

all2all_refblast $reference_ORGS $fasta_path files/all2all/all2all.genelist files/all2all_final $blastdb_path $region $region tblastx files/oneway/set.tmp

while IFS= read -r subject
do
#	new_subject=$old_subject
#	old_subject=$subject
while IFS= read -r query
do
		#if [[ ! -s files/all2all_final/$query-$subject.orths ]]; then
		#	blast_formatter -archive files/all2all_final/all2all.$query.$subject -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out files/all2all_final/$query-$subject.out
		#	blast_formatter -archive files/all2all_final/all2all.$subject.$query -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive" -out files/all2all_final/$subject-$query.out
		#	#blast_formatter -archive files/all2all_final/all2all.$query.$subject -outfmt "5" -out files/all2all_final/$query-$subject.xml
		#	#blast_formatter -archive files/all2all_final/all2all.$subject.$query -outfmt "5" -out files/all2all_final/$subject-$query.xml
		#	#/data/meyer/viz/tools/miniconda3/envs/local_root/bin/python transcriptologs.py -i1 $query-$subject.out -i2 $subject-$query.out -o $query-$subject.orths
		#	python RBH-v1.py files/all2all_final/$query-$subject.out files/all2all_final/$subject-$query.out files/all2all_final/$query-$subject.orths
		#	python RBH-v1.py files/all2all_final/$subject-$query.out files/all2all_final/$query-$subject.out files/all2all_final/$subject-$query.orths
		#fi
		twoway_RBH files/all2all_final $query $subject
done < "files/oneway/set.tmp"
#	if [[ "$new_subject" != "-1" ]]; then
#		twoway_RBH files/all2all_final $old_subject $subject
#	fi
done < "$reference_ORGS"

#Rscript merge_orths.R files/oneway/set.tmp files/all2all_final/ files/all2all_final/final_df.txt 1e-05

Rscript wisard.R files/oneway/set.tmp files/all2all_final

exit
