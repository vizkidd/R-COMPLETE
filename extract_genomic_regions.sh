#!/bin/bash
#$1 = genome fasta file
#$2 = formatted gtf/gff3(UNTESTED) file
#$3 = basename for bedfiles
#$4 = input gene list
#$5 = org name(eg x_laevis)

# parallel_checkODB() {
#     #echo $1 $2 $3
#     ref_org_name=$(echo $1 | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
#   ref_org_ID=$(grep -i -P "\t\b($ref_org_name)\b\t" $ORTHODB_PATH/$ORTHODB_PREFIX_species.tab | awk '{print $2}')
#   if [[ $ref_org_ID ==  "" ]] ; then
#     ref_org_ID=$(grep -i -P "($ref_org_name)" "$ORTHODB_PATH/$ORTHODB_PREFIX"_species.tab | awk '{print $2}')
#   fi
#   ref_gene_ID=$(grep -i -w $ref_org_ID $ORTHODB_PATH/$ORTHODB_PREFIX_genes.tab | grep -i -P "\b($2)\b\t" | awk '{print $1}')
#   #echo $ref_org_ID $ref_org_name $ref_gene_ID $ORTHODB_PATH
#   if [[ $ref_gene_ID !=  "" ]] ; then
#   echo $ref_org_name $ref_gene_ID #$ref_org_ID
#   grep "$ref_gene_id" $ORTHODB_PATH/$ORTHODB_PREFIX_OG2genes.tab | awk '{print $1}' > "$TEMP_PATH/$f_org_name.$s_name.clusters"
#   #cat "$TEMP_PATH/$f_org_name.$s_name.clusters"
#   grep -w -f "$TEMP_PATH/$f_org_name.$s_name.clusters" $ORTHODB_PATH/$ORTHODB_PREFIX_OG2genes.tab | grep -w $3 | awk '{print $2}' | sort | uniq > "$TEMP_PATH/$f_org_name.$s_name.genes"
#   grep -w -f "$TEMP_PATH/$f_org_name.$s_name.genes" $ORTHODB_PATH/$ORTHODB_PREFIX_genes.tab | awk '{print $4}' > "$TEMP_PATH/$f_org_name.$s_name.names"
#   cat "$TEMP_PATH/$f_org_name.$s_name.names"
#   grep -i -f "$TEMP_PATH/$f_org_name.$s_name.names" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq >> files/genes/$f_org_name/$file_out.rna_list
# fi
# }

##FUNCTIONS
check_OrthoDB() {
  #1 - gene (name)
 # proc_list=()
rm files/genes/$f_org_name/$file_out.rna_list
readarray refs < $REF_ORGS
s_name=$(echo $1 | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')
org_ID=$(grep -i -P "\t\b($org_name)\b\t" "$ORTHODB_PATH/$ORTHODB_PREFIX"_species.tab | awk '{print $2}')
#gene_ID=$(zcat $ORTHODB_PATH/odb10v1_genes.tab.gz | grep $org_ID | grep -i -P "\b($1)\b\t" | awk '{print $1}')
if [[ $org_ID ==  "" ]] ; then
  org_ID=$(grep -i -P "($org_name)" "$ORTHODB_PATH/$ORTHODB_PREFIX"_species.tab | awk '{print $2}')
fi
if [[ $org_ID !=  "" ]] ; then
  #echo $org_ID
  #parallel -j1 --compress --pipepart -a "$REF_ORGS" --recstart '\n' --block -1 "./check_OrthoDB.sh {} $1 $org_ID $f_org_name $s_name"
  parallel -j ${#refs[@]} "./check_OrthoDB.sh {} $1 $org_ID $f_org_name $s_name $ANNO_FILE" ::: ${refs[@]}
  #echo ${refs[@]} | parallel  -j ${#refs[@]} --recstart ' ' --block -1 "./check_OrthoDB.sh {} $1 $org_ID $f_org_name $s_name " 
  cat $TEMP_PATH/*.$f_org_name.$s_name > files/genes/$f_org_name/$file_out.rna_list
  rm $TEMP_PATH/*.$f_org_name.$s_name*
else
  echo "$org_name not found in OrthoDB files/OrthoDB files don't exist"
  ORG_IN_ODB="FALSE"
fi
#for id in ${proc_list[@]}
#do
#  echo $id
#  wait $id
#done
}

check_gene() {
  gene_full=$1
  f_org_name=$2
  bed_prefix=$3
  ANNO_FILE=$4
  GENOME_FILE=$5
  TEMP_PATH=$6
  FASTA_PATH=$7
  transcript_regions=("cds" "exons" "noncodingexons" "3utr" "5utr")

  if [[ $gene_full !=  "" ]] ; then
  gene=$(echo $gene_full | sed 's/.$//')
  file_out=$(echo $gene_full | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')
  echo $gene $gene_full #$f_org_name $file_out $bed_prefix
  #Saving gene info list in a file
  grep -i "$gene" $ANNO_FILE > files/genes/$f_org_name/$file_out.gtf_slice
  #grep -i "$gene" $2 | awk -F'\t' '{print $9}' | awk -F';' '{print $3}' | awk '{print $2}' | sed 's/"//g' | sort | uniq > files/genes/$5/$gene.rna_list
  if grep -q -w -i "$gene_full" files/genes/$f_org_name/$file_out.gtf_slice ; then
    echo "1 - Gene name has perfect match"
    grep -i -w "$gene_full" files/genes/$f_org_name/$file_out.gtf_slice | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$f_org_name/$file_out.rna_list
  elif grep -q -i "$gene_full" files/genes/$f_org_name/$file_out.gtf_slice ; then
    #grep -i "$gene_full" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
    if [[ $ORG_IN_ODB ==  "TRUE" ]] ; then
      echo "2 - Gene name matches a supergroup or subgroup, checking OrthoDB"
     time check_OrthoDB $gene_full
    fi
  else
    #grep -i "$gene" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
    if [[ $ORG_IN_ODB ==  "TRUE" ]] ; then
      echo "3 - Gene name has a partial/no match, checking OrthoDB"
     time check_OrthoDB $gene_full
    fi
  fi
  if [ ! -s files/genes/$f_org_name/$file_out.rna_list ]; then
    grep -i "$gene_full" files/genes/$f_org_name/$file_out.gtf_slice | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$f_org_name/$file_out.rna_list
  fi
  cat files/genes/$f_org_name/$file_out.rna_list
  if [ -s files/genes/$f_org_name/$file_out.rna_list ]
  then
    #if [ ! -s $FASTA_PATH/$5/$file_out.cds ]; then #$6/$5/$7.$2
    #  get_fasta $gene_full cds $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE ##CDS should come first because the 3utr flanks depend on cds coordinates
    #fi
    #if [ ! -s $FASTA_PATH/$5/$file_out.exons ]; then
    # get_fasta $gene_full exons $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE
    #fi
    #if [ ! -s $FASTA_PATH/$5/$file_out.noncodingexons ]; then
    #  get_fasta $gene_full noncodingexons $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE
    #fi
    #if [ ! -s $FASTA_PATH/$5/$file_out.3utr ]; then
    #  get_fasta $gene_full 3utr $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE
    #fi
    #if [ ! -s $FASTA_PATH/$5/$file_out.5utr ]; then
    #  get_fasta $gene_full 5utr $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE
    #fi

    time parallel -j ${#transcript_regions[@]} "get_fasta $gene_full {} $bed_prefix $GENOME_FILE $f_org_name $FASTA_PATH $file_out 'files/genes/$f_org_name/$file_out.gtf_slice' $TEMP_PATH " ::: ${transcript_regions[@]} #$8
  fi

fi
}

get_fasta() {
#1 - gene name
#2 - transcript region (cds/5utr/3utr/exons/noncodingexons)
#3 - basename for bedfiles
#4 - genome fasta file
#5 - org name
#6 - FASTA OUPUT FOLDER
#7 - Formatted file name for gene
#8 - GTF file path (for extracting the 'real' gene annotation)
#9 - TEMP PATH
#10 - Default flank_len
#echo $1 $2 $3 $4 $5 $6 $7 $8 $9 
#echo $(grep -i -w -f files/genes/$5/$7.rna_list $3_$2.bed)
#echo $1 $2 $3 $4 $5 $6 $7 $8 $9 
#flank_len=${10}
#echo $flank_len
if [ ! -s $6/$5/$7.$2 ]; then #$FASTA_PATH/$5/$file_out.cds
grep -i -w -f files/genes/$5/$7.rna_list $3_$2.bed > $9/$5_$1_$2.bed

##TO get flanks 
#grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3"

##TO find MEDIAN values
#grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | sort -n | awk '{arr[NR]=$1} END {if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}'

##TO find MEAN
#grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'


##EXTEND regions using bedtools slop
if [[ "$2" == "3utr" ]]; then
  ##SHIFT by 3bp to omit stop codons in 3UTRs
  #bedtools shift -s 3 -i $TEMP_PATH/$5_"$1"_$2.bed -g $TEMP_PATH/$5_genomeFile.txt ##WONT WORK, and dont need to because stop codons are part of CDS
  ##GET AVERAGE UTR LENGTH for getting FLANK LEN
  flank_len=$(grep -i -f files/genes/$5/$7.rna_list files/genes/$5/$7.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
  
  if [ -s $9/$5_$1_$2.bed ]; then ## if bed file is empty then we will have to take flanks from CDS
  	grep -i "$2_FLANK" $9/$5_$1_$2.bed > $9/$5_$1_$2_FLANK.tmp
  	##REMOVE the line because we stored flanking regions in a seperate bed file
  	sed -i'' "/$2_FLANK/d" $9/$5_$1_$2.bed
  	bedtools slop -s -r $flank_len -l 0 -i $9/$5_"$1"_$2_FLANK.tmp -g $9/$5_genomeFile.txt > $9/$5_$1_$2_FLANK.bed
  else
  	bedtools flank -s -r $flank_len -l 0 -i $9/$5_"$1"_cds.bed -g $9/$5_genomeFile.txt > $9/$5_"$1"_$2_FLANK.bed
  	sed -i'' "s/cds/3utr_FLANK/g" $9/$5_"$1"_$2_FLANK.bed
  fi
  bedtools getfasta -s -split -fi $4 -bed $9/$5_"$1"_$2_FLANK.bed -nameOnly -fullHeader >> $6/$5/$7.$2.tmp #NOTUSING name+ because BLAST doesn't like long IDs
elif [[ "$2" == "5utr" ]]; then
  ##GET AVERAGE UTR LENGTH for getting FLANK LEN
  flank_len=$(grep -i -f files/genes/$5/$7.rna_list files/genes/$5/$7.gtf_slice | grep -i "utr" | grep -i "five\|5" | awk '{print $5-$4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
  if [ -s $9/$5_"$1"_$2.bed ]; then ## if bed file is empty then we will have to take flanks from CDS
  	grep -i "$2_FLANK" $9/$5_"$1"_$2.bed > $9/$5_"$1"_$2_FLANK.tmp
  	##REMOVE the line because we stored flanking regions in a seperate bed file
  	sed -i'' "/$2_FLANK/d" $9/$5_"$1"_$2.bed
  	bedtools slop -s -l $flank_len -r 0 -i $9/$5_"$1"_$2_FLANK.tmp -g $9/$5_genomeFile.txt > $9/$5_"$1"_$2_FLANK.bed
  else
  	bedtools flank -s -l $flank_len -r 0 -i $9/$5_"$1"_cds.bed -g $9/$5_genomeFile.txt > $9/$5_"$1"_$2_FLANK.bed
  	sed -i'' "s/cds/5utr_FLANK/g" $9/$5_"$1"_$2_FLANK.bed
  fi
  bedtools getfasta -s -split -fi $4 -bed $9/$5_"$1"_$2_FLANK.bed -nameOnly -fullHeader >> $6/$5/$7.$2.tmp #NOTUSING name+ because BLAST doesn't like long IDs
fi

bedtools getfasta -s -split -fi $4 -bed $9/$5_"$1"_$2.bed -nameOnly -fullHeader >> $6/$5/$7.$2.tmp #NOTUSING name+ because BLAST doesn't like long IDs

#python label_sequenceIDs.py $5 $FASTA_PATH/$5/$7.$2.tmp $1 $8 $(grep -i "filter" parameters.txt | awk -F'=' '{print $2}') ##$org $file #$filename

#while IFS='>' read -r rec;
#do
#echo $FASTA_PATH/$5/$7.$2
parallel -j1 --compress --pipepart -a "$6/$5/$7.$2.tmp" --recstart '>' --block -1 "sh label_sequenceIDs.sh $5 $1 $8 $6/$5/$7.$2.tmp $6/$5/$7.$2" #-word_size 5 -evalue 1e-25
#read -r line seq <<< $(echo "$rec")
#printf "%s\n%s" $line $seq | parallel -j1 "sh label_sequenceIDs.sh $5 $1 $8 $FASTA_PATH/$5/$7.$2.tmp $FASTA_PATH/$5/$7.$2" #-word_size 5 -evalue 1e-25
#done< "$FASTA_PATH/$5/$7.$2.tmp"
##sh label_sequenceIDs.sh $5 $1 $8 $FASTA_PATH/$5/$7.$2.tmp $FASTA_PATH/$5/$7.$2 #-word_size 5 -evalue 1e-25
#echo $flank_len
rm $6/$5/$7.$2.tmp
fi
}
####

##ENTRYPOINT
GENOMES_PATH=$(grep -i -w "genomes_path" parameters.txt | awk -F'=' '{print $2}') 
ANNOS_PATH=$(grep -i -w "annos_path" parameters.txt | awk -F'=' '{print $2}') 
FASTA_PATH=$(grep -i -w "fasta_path" parameters.txt | awk -F'=' '{print $2}') 
TEMP_PATH=$(grep -i -w "temp_path" parameters.txt | awk -F'=' '{print $2}') 
CLEAN_EXTRACT=$(grep -i -w "clean_extract" parameters.txt | awk -F'=' '{print $2}') 
ORTHODB_PATH=$(grep -i -w "orthodb_files_path" parameters.txt | awk -F'=' '{print $2}') 
ORTHODB_PREFIX=$(grep -i -w "orthodb_prefix" parameters.txt | awk -F'=' '{print $2}') 
flank_len=$(grep -i -w "flank_len" parameters.txt | awk -F'=' '{print $2}')   #2000
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'=' '{print $2}') 
#PY2_PATH=$(grep -i -w "python2_path" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH=$(grep -i -w "python3_path" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH="${PY3_PATH/#\~/$HOME}"
ORG_IN_ODB="TRUE"
transcript_regions=("cds" "exons" "noncodingexons" "3utr" "5utr")

export -f get_fasta
export -f check_gene
export -f check_OrthoDB
#export -f parallel_checkODB

f_org_name=$5
bed_prefix=$3
#org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0); print $1" "$2; } ;')
echo $org_name

mkdir files
mkdir $FASTA_PATH
mkdir files/bed
mkdir $TEMP_PATH
mkdir files/genes
mkdir files/genes/$5/

GENOME_FILE=$1
ANNO_FILE=$2

if [[ ! -s $GENOME_FILE || ! -s $ANNO_FILE ]] ; then
  echo "$GENOME_FILE (or) $ANNO_FILE doesn't exist"
  exit 1
fi

if [[ ${GENOME_FILE##*.} ==  "gz" ]] ; then
file_name=${GENOME_FILE%.*}
##gunzip $GENOME_FILE
gunzip -f $GENOME_FILE
GENOME_FILE=$file_name
fi

if [[ ${ANNO_FILE##*.} ==  "gz" ]] ; then
file_name=${ANNO_FILE%.*}
#gunzip $ANNO_FILE
gunzip -f $ANNO_FILE
ANNO_FILE=$file_name
fi

if [[ ${ANNO_FILE##*.} !=  "gtf" ]] ; then
file_name=${ANNO_FILE%.*}
#gunzip $ANNO_FILE
gffread $ANNO_FILE -T -O -E -o $ANNOS_PATH/"$5".gtf
ANNO_FILE=$ANNOS_PATH/"$5".gtf
#echo $ANNO_FILE
fi

echo $GENOME_FILE
echo $ANNO_FILE

#if [ ! -f "$GENOME_FILE".fai ]; then
samtools faidx $GENOME_FILE
#fi

awk -v OFS='\t' {'print $1,$2'} $GENOME_FILE.fai > $TEMP_PATH/$5_genomeFile.txt
#rm $3*
rm $bed_prefix*
$PY3_PATH extract_transcript_regions.py -i $ANNO_FILE --gtf  -o $3

#if [[ "$7" == "TRUE" ]]; then
if [[ $CLEAN_EXTRACT ==  "TRUE" ]] ; then
  rm -rf $FASTA_PATH/$5
fi
mkdir $FASTA_PATH/$5

#FOR each gene
#while IFS= read -r gene_full
#do
#if [[ $gene_full !=  "" ]] ; then
#  gene=$(echo $gene_full | sed 's/.$//')
#  echo $gene $gene_full
#  file_out=$(echo $gene_full | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')
#
#  #grep -i "$gene" $2 | awk -F'\t' '{print $9}' | awk -F';' '{print $3}' | awk '{print $2}' | sed 's/"//g' | sort | uniq > files/genes/$5/$gene.rna_list
#  if grep -q -w -i "$gene_full" $ANNO_FILE ; then
#    echo "1 - Gene name has perfect match"
#    grep -i -w "$gene_full" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
#  elif grep -q -i "$gene_full" $ANNO_FILE ; then
#    #grep -i "$gene_full" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
#    if [[ $ORG_IN_ODB ==  "TRUE" ]] ; then
#      echo "2 - Gene name matches a supergroup or subgroup, checking OrthoDB"
#     check_OrthoDB $gene_full
#    fi
#  else
#    #grep -i "$gene" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
#    if [[ $ORG_IN_ODB ==  "TRUE" ]] ; then
#      echo "3 - Gene name has a partial/no match, checking OrthoDB"
#      check_OrthoDB $gene_full
#    fi
#  fi
#  if [ ! -s files/genes/$5/$file_out.rna_list ]; then
#    grep -i "$gene_full" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
#  fi
#  cat files/genes/$5/$file_out.rna_list
#  if [ -s files/genes/$5/$file_out.rna_list ]
#  then
#    #if [ ! -s $FASTA_PATH/$5/$file_out.cds ]; then #$6/$5/$7.$2
#    #  get_fasta $gene_full cds $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE ##CDS should come first because the 3utr flanks depend on cds coordinates
#    #fi
#    #if [ ! -s $FASTA_PATH/$5/$file_out.exons ]; then
#    # get_fasta $gene_full exons $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE
#    #fi
#    #if [ ! -s $FASTA_PATH/$5/$file_out.noncodingexons ]; then
#    #  get_fasta $gene_full noncodingexons $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE
#    #fi
#    #if [ ! -s $FASTA_PATH/$5/$file_out.3utr ]; then
#    #  get_fasta $gene_full 3utr $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE
#    #fi
#    #if [ ! -s $FASTA_PATH/$5/$file_out.5utr ]; then
#    #  get_fasta $gene_full 5utr $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE
#    #fi
#    parallel -j ${#transcript_regions[@]} "get_fasta $gene_full {} $3 $GENOME_FILE $5 $FASTA_PATH $file_out $ANNO_FILE " ::: ${transcript_regions[@]}
#  fi
#fi
#done < "$4"

readarray gene_list < "$4"

time parallel -j ${#gene_list[@]} "check_gene {} $f_org_name $bed_prefix $ANNO_FILE $GENOME_FILE $TEMP_PATH $FASTA_PATH $flank_len" ::: ${gene_list[@]} #get_fasta is called from check_gene function after preliminary checks

find files/genes/$5 -empty | awk -F'/' '{print $NF}' | sed 's/.rna_list//g' > $FASTA_PATH/$5/MISSING_GENES
grep -v -i -f $FASTA_PATH/$5/MISSING_GENES $4 | sort > $FASTA_PATH/$5/AVAILABLE_GENES
find files/genes/$5 -empty -delete
find $FASTA_PATH/$5 -empty -delete
#rm $TEMP_PATH/$5*
exit 0
mv $GENOME_FILE $GENOMES_PATH/$5.fa
mv $ANNO_FILE $ANNOS_PATH/$5.gtf
#bgzip -i $GENOMES_PATH/$5.fa
#bgzip -i $ANNOS_PATH/$5.gtf
gzip -f $GENOMES_PATH/$5.fa
gzip -f $ANNOS_PATH/$5.gtf
#rm $GENOMES_PATH/$5.fa
#rm $ANNOS_PATH/$5.gtf
rm $TEMP_PATH/$f_org_name.*.*
rm files/genes/$5/*.gtf_slice
#rm $2
rm $3*
#rm files/bed/$5*
##rm $GENOMES_PATH/doc_*
##rm $GENOMES_PATH/$1.fai
##rm $ANNOS_PATH/*.tsv
##rm $ANNOS_PATH/*.txt

exit 0
