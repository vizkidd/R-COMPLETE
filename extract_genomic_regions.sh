#!/bin/bash
#$ -cwd
#$ -V
#$ -l data
#$ -l m_mem_free 8G

# $1 = genome fasta file
# $2 = formatted gtf/gff3(UNTESTED) file
# $3 = input gene list
# $4 = org name(eg x_laevis)

##FUNCTIONS

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
  
  grep -i -w $gene_full files/genes/$org_name/gtf_stats.csv | awk -F',' '{print $3}' | grep -w -f - "$base_bed"_"$reg.bed" > $TEMP_PATH/"$s_name"_"$reg.bed"

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

####

##ENTRYPOINT

#Check if scripts exist
if [[ ! -s extract_gtf_info.R || ! -s check_OrthoDB.sh || ! -s label_sequenceIDs.sh || ! -s parameters.txt ]] ; then
  >&2 color_FG_Bold $Red "Missing scripts &/or parameters.txt!"
  >&2 color_FG_Bold $Red  "Need : extract_gtf_info.R, check_OrthoDB.sh, label_sequenceIDs.sh, parameters.txt"
  exit -1
fi

GENOME_FILE=$1
ANNO_FILE=$2
f_org_name=$4
org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0); print $1" "$2; } ;')

export -f get_fasta

source fasta_functions.sh

GENOMES_PATH=$(grep -i -w "genomes_path" parameters.txt | check_param) 
ANNOS_PATH=$(grep -i -w "annos_path" parameters.txt | check_param) 
BED_PATH=$(grep -i -w "bed_path" parameters.txt | check_param) 
FASTA_PATH=$(grep -i -w "fasta_path" parameters.txt | check_param) 
TEMP_PATH=$(grep -i -w "temp_path" parameters.txt | check_param) 
CLEAN_EXTRACT=$(grep -i -w "clean_extract" parameters.txt | check_param) 
TRANSCRIPT_REGIONS=($(grep -i -w "transcript_regions" parameters.txt | check_param | sed "s/,/\n/g" ))
GENE_SEARCH_MODE=$(grep -i -w "gene_search_mode" parameters.txt | check_param)
ORTHODB_PATH_PREFIX=$(grep -i -w "orthodb_path_prefix" parameters.txt | check_param) 
REMOVE_DOWNLOADS=$(grep -i -w "remove_downloads" parameters.txt | check_param) 
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | check_param) 
seqID_delimiter=$(grep -i -w "seqID_delimiter" parameters.txt | check_param) 
LABEL_SEQS=$(grep -i -w "label_sequence_IDs" parameters.txt | check_param) 
#n_threads=$(nproc --all)
n_threads=$(grep -i -w "max_concurrent_jobs" parameters.txt | check_param)

bed_prefix="$BED_PATH/$f_org_name"
MODE=""
if [[ $GENE_SEARCH_MODE=="HARD" ]]; then
  MODE="-w"
fi

echo "Command : $0 $@"

if [[ $CLEAN_EXTRACT ==  "TRUE" ]] ; then
  rm -rf $FASTA_PATH/$5
  rm -rf files/genes/$f_org_name
  rm -rf $TEMP_PATH/$f_org_name
  rm -f $GENOME_FILE.fai
  rm -f files/genes/$f_org_name/*
  rm -f $bed_prefix/*
fi

mkdir -p $FASTA_PATH/$f_org_name
mkdir -p files/genes/$f_org_name/odb/
mkdir -p $TEMP_PATH/$f_org_name
mkdir -p $bed_prefix

genome_pipe="$TEMP_PATH/$f_org_name/genome_pipe"
rm -f $genome_pipe
mkfifo $genome_pipe

gfile_name=$GENOME_FILE
if [[ ${GENOME_FILE##*.} ==  "gz" ]] ; then
  gfile_name=${GENOME_FILE%.*}
    zcat -f $GENOME_FILE > $genome_pipe & 
    cat < $genome_pipe > $gfile_name & 
    genome_proc_id=$(echo $!)
fi

if [[ $(echo $ANNO_FILE | grep -q -i "gtf") ]] ; then
  file_name=${ANNO_FILE%.*}
  nohup zcat -f $ANNO_FILE | gffread - -T -O -E | gzip -c - > $ANNOS_PATH/"$5".gtf.gz &> /dev/null &
  anno_proc_id=$(echo $!)
  ANNO_FILE=$ANNOS_PATH/"$5".gtf.gz
fi

>&1 color_FG_BG_Bold $Black $BG_Yellow "0.1 Building Reference Genome Index (Samtools faidx)..."

if [[ ! -s $GENOME_FILE.fai ]] ; then
  time samtools faidx --fai-idx $GENOME_FILE.fai $genome_pipe & 
  genome_index_proc=$(echo $!)
fi

while [[ $(lsof -R -p $$ | wc -l) -gt $(($(ulimit -Sn)-200)) ]]; #limit file descriptors for child & parent process to 1024-200 
do
  printf %s"\t"%s"\n" "Waiting for file handles to close.." $(lsof -R -p $$ | wc -l)
  sleep 5
done

###################################################################################################

gene_list=($(cat $3 | sort | uniq | grep -v -w -i "gene"))

>&1 color_FG_BG_Bold $Black $BG_Yellow "1. Checking Gene Names & Splitting GTFs for parallel processing..."


if [[ ! -z $anno_proc_id ]]; then
  wait $anno_proc_id
fi

if [[ ! -s files/genes/$f_org_name/1.list || ! -s files/genes/$f_org_name/2.list || ! -s files/genes/$f_org_name/gtf_stats.csv || ! -s $BED_PATH/$f_org_name/$f_org_name.bed ]] ; then
  eexp_gene=$(printf -- '%s\n' "${gene_list[@]}")
  
  #Splitting GTF into multiple parts based on grep output for downstream parallel processing in extract_gtf_info.R
  time zgrep -i $MODE -A 0 --group-separator='>' -f <(echo "${eexp_gene[@]}") $ANNO_FILE | csplit --quiet -z --suffix-format="%0d.gtf_slice" --prefix="$TEMP_PATH/$f_org_name/1." --suppress-matched - '/>/' '{*}'
  
  zgrep -hPo 'gene_name "\K[^"]+' $TEMP_PATH/$f_org_name/*.gtf_slice | sort | uniq | awk 'NF' > files/genes/$f_org_name/1.list
  echo ${gene_list[@]/($(cat files/genes/$f_org_name/1.list)))} | awk 'NF' > files/genes/$f_org_name/2.list
fi

>&1 echo $(color_FG $Green "1. DONE : Available Genes : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/1.list")$(color_FG $Green ", Genes not found : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/2.list ")

#######################################################################################################

>&1 color_FG_BG_Bold $Black $BG_Yellow "2. Checking OrthoDB for missing & orthologous genes..."

if [[ ! -s files/genes/$f_org_name/odb.list || ! -s files/genes/$f_org_name/final.list || ! -s files/genes/$f_org_name/odb.final_map ]] ; then
  time ./check_OrthoDB.sh $f_org_name $4 files/genes/$f_org_name/odb.list $ANNO_FILE 
fi

if [[ $(awk 'END{print NR;}' "files/genes/$f_org_name/odb.list" | awk '{print $1}') !=  0 ]] ; then
    odb_gene_list=($(cat "files/genes/$f_org_name/odb.list" | grep -v -w -i "gene"))
    existing_list=($(cat files/genes/$f_org_name/1.list | sort | uniq)) 
    short_list=($(echo ${odb_gene_list[@]/${existing_list[@]}}))

  if [[ "${#short_list[@]}" > 0 ]] ; then
    eexp_gene=$(printf -- '%s\n' "${short_list[@]}")
    time echo "${eexp_gene[@]}" | zgrep -i $MODE -A 0 --group-separator='>' -f - $ANNO_FILE | csplit --quiet -z --suffix-format="%0d.gtf_slice" --prefix="$TEMP_PATH/$f_org_name/2." --suppress-matched - '/>/' '{*}' 
  fi
else
  >&2 echo $(color_FG_Bold $Red "2. Error : ")$(color_FG_BG_Bold $White $BG_Red "files/genes/$f_org_name/odb/odb.list")$(color_FG_Bold $Red " missing, Possibly orthologous genes were not found ")
  >&2 color_FG_Bold $Red "2. If unsure, re-run command: ./jobhold.sh check_ortho ./check_OrthoDB.sh $f_org_name $4 $ANNO_FILE"
fi

##REFRESH geene list based on file names
cat files/genes/$f_org_name/1.list files/genes/$f_org_name/2.list files/genes/$f_org_name/odb.list | awk 'NF' > files/genes/$f_org_name/full.list

>&1 echo $(color_FG $Green "2. DONE : Full List : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/full.list")$(color_FG $Green ", List from ODB : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/odb.list")$(color_FG $Green ", ODB Cluster to Genes Map : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/odb.final_map")

#######################################################################################################

>&1 color_FG_BG_Bold $Black $BG_Yellow "3. Extracting Transcript Stats from GTF_Slices...(log:$TEMP_PATH/$f_org_name/get_GTF_info.[o/e])"

time Rscript --vanilla --verbose extract_gtf_info.R $TEMP_PATH/$f_org_name/ $f_org_name files/genes/$f_org_name/gtf_stats.csv 1> $TEMP_PATH/$f_org_name/get_GTF_info.o 2> $TEMP_PATH/$f_org_name/get_GTF_info.e
sed 1d files/genes/$f_org_name/gtf_stats.csv | awk -F',' '{print $1"\n"}' | sort | uniq | awk 'NF' > files/genes/$f_org_name/final.list

if [[ ! -s files/genes/$f_org_name/gtf_stats.csv || ! -s files/genes/$f_org_name/final.list ]] ; then
  >2& color_FG_Bold $Red "3. ERROR: Extraction of transcript stats failed..."
  >2& color_FG_Bold $Red "3. Check $TEMP_PATH/$f_org_name/get_GTF_info.[o/e]"
  >2& color_FG_Bold $Red "3. Remove $TEMP_PATH/$f_org_name/ & files/genes/$f_org_name/gtf_stats.csv and re-run the pipeline"
  exit 1
fi

>&1 echo $(color_FG $Green "3. DONE : Final List : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/final.list")$(color_FG $Green ", GTF stats : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/gtf_stats.csv")

#######################################################################################################

gene_list=($(cat files/genes/$f_org_name/final.list | awk 'NF' ))
s_names=($(awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo ${gene_list[@]})))

if [[ ! -z $genome_index_proc ]]; then
  wait $genome_index_proc
fi
if [[ ! -z $genome_proc_id ]]; then
  wait $genome_proc_id
fi

>&1 color_FG_BG_Bold $Black $BG_Yellow "4. Fetching sequences from Genome..."
>&1 color_FG $Yellow "Transcript Regions : $(IFS=","; echo ${TRANSCRIPT_REGIONS[*]}) "

time parallel --max-procs $n_threads "get_fasta {1} {2} $bed_prefix/$f_org_name $gfile_name $f_org_name $FASTA_PATH/$f_org_name $ANNO_FILE $TEMP_PATH/$f_org_name $LABEL_SEQS ;" ::: ${gene_list[@]} ::: ${TRANSCRIPT_REGIONS[@]}

>&1 echo $(color_FG $Green "4. DONE : FASTA PATH : ")$(color_FG_BG_Bold $White $BG_Purple "$FASTA_PATH/$f_org_name")

#############################################################################################################

if [[ $LABEL_SEQS ==  "TRUE" ]] ; then
  >&1 color_FG_BG_Bold $Black $BG_Yellow "4.1 Labelling sequences..."
  tmp_names=($(parallel --link --max-procs $n_threads "echo {1},{2}" ::: ${gene_list[@]} ::: ${s_names[@]}))
  time parallel --max-procs $n_threads "printf -- %s,%s\\\n {1} {2}" ::: ${tmp_names[@]} ::: ${TRANSCRIPT_REGIONS[@]} | parallel --colsep "," --max-procs $n_threads "./label_sequenceIDs.sh $f_org_name {1} $FASTA_PATH/$f_org_name/{2}.{3} $FASTA_PATH/$f_org_name/{2}.{3}.tmp files/genes/$f_org_name/odb.final_map" 
fi

#######################################################################################################

>&1 color_FG_BG_Bold $Black $BG_Yellow "5. Generating Metadata and Cleaning up..."

sed 1d files/genes/$f_org_name/gtf_stats.csv | awk -F',' '{print $2}' | sort | uniq >  files/genes/$f_org_name/AVAILABLE_GENES
grep -v -i -f files/genes/$f_org_name/AVAILABLE_GENES $4 | sort | uniq > files/genes/$f_org_name/MISSING_GENES

find $FASTA_PATH/$f_org_name/ -type f -name "*.fai" -exec rm -f {} +

if [[ ! -z $gfile_name ]]; then
  rm -f $gfile_name
fi

rm -rf $TEMP_PATH/$f_org_name/
rm -f $bed_prefix/$f_org_name_*
if [[ $REMOVE_DOWNLOADS ==  "TRUE" ]] ; then
   rm $GENOME_FILE
   rm $ANNO_FILE
fi

>&1 color_FG_BG_Bold $Purple $BG_White "Extraction DONE for organism : $f_org_name"

exit 0
