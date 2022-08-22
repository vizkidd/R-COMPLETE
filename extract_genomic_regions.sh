#!/bin/bash
#$ -cwd
#$ -V
#$ -l data
# $ -pe smp 128
#$ -l m_mem_free 8G

# $1 = genome fasta file
# $2 = formatted gtf/gff3(UNTESTED) file
# $3 = basename for bedfiles
# $4 = input gene list
# $5 = org name(eg x_laevis)

##FUNCTIONS
function remove_duplicate_genes() {
  #$1 - gene_list
  #$2 - TEMP_PATH
  cat $1 | sort | uniq > $1.tmp
  cat $1.tmp > $1
  rm -f $1.tmp
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

#echo $(grep -i -w -f files/genes/$org_name/$s_name.rna_list $base_bed_$reg.bed)
#echo $gene_full $reg $base_bed $genome_fa $org_name $FASTA_PATH $s_name $GTF_PATH $TEMP_PATH 
#flank_len=${10}
#echo $flank_len

if [ ! -s $FASTA_PATH/$s_name.$reg ]; then #$FASTA_PATH/$file_out.cds
  >&2 echo $s_name $reg
  # if grep -q -w -i "$s_name" files/genes/$org_name/gtf_stats.csv ; then
  #   return 2
  # fi

  #grep -i -w -f $GENES_META/$org_name/$s_name.rna_list "$base_bed"_"$reg.bed" > $TEMP_PATH/"$s_name"_"$reg.bed"
  grep -i -w $gene_full files/genes/$org_name/gtf_stats.csv | awk -F',' '{print $3}' | grep -w -f - "$base_bed"_"$reg.bed" > $TEMP_PATH/"$s_name"_"$reg.bed"

  ##TO get flanks 
  #grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3"

  ##TO find MEDIAN values
  #grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | sort -n | awk '{arr[NR]=$1} END {if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}'

  ##TO find MEAN
  #grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'

  if [[ -s $TEMP_PATH/"$s_name"_"$reg.bed" ]]; then 
    #cp -s -f $TEMP_PATH/"$s_name"_"$reg.bed" $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed
    ln -r -f -s $TEMP_PATH/"$s_name"_"$reg.bed" $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed
  elif [[ -s $TEMP_PATH/"$s_name"_cds.bed && $reg!="cds" ]]; then
    ## if bed file is empty then we will have to take flanks from CDS
    if [[ "$reg" == "3utr" ]]; then
      local flank_len=$(grep -w -i "$gene_full" files/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($9 > 3) sum += $9; n++ } END { if (n > 0) print sum / n; }') ##3' UTR length must be > 3 (because of existence of stop codons)
      if [[ $flank_len == "" ]]; then
        return 2
      fi
      bedtools flank -s -l $flank_len -r 0 -i $TEMP_PATH/"$s_name"_cds.bed -g $genome_fa.fai > $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed #files/genes/"$org_name"/genomeFile.txt
    elif [[ "$reg" == "5utr" ]]; then
      local flank_len=$(grep -w -i "$gene_full" files/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($8 > 0) sum += $8; n++ } END { if (n > 0) print sum / n; }') ##5' UTR length must be > 0
      if [[ $flank_len == "" ]]; then
        return 2
      fi
      bedtools flank -s -r $flank_len -l 0 -i $TEMP_PATH/"$s_name"_cds.bed -g $genome_fa.fai > $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed #files/genes/$org_name/genomeFile.txt
    fi
    if [[ -s $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed ]]; then
      sed -i "s/cds/$(echo $reg)_FLANK/g" $TEMP_PATH/"$s_name"_"$reg"_FETCH.bed
    else
      >&2 echo "$reg FLANKS not found for $s_name"
      return 1
    fi
  #else
  #  >&2 echo "CDS not found for $s_name"
  #  return 1
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
####

##ENTRYPOINT

#Check if scripts exist
if [[ ! -s extract_transcript_regions.py || ! -s extract_gtf_info.R || ! -s check_OrthoDB.sh || ! -s label_sequenceIDs.sh || ! -s parameters.txt ]] ; then
  echo "Missing scripts &/or parameters.txt!"
  echo "Need : extract_transcript_regions.py, extract_gtf_info.R, check_OrthoDB.sh, label_sequenceIDs.sh, parameters.txt"
  exit -1
fi

GENOME_FILE=$1
ANNO_FILE=$2
f_org_name=$5
bed_prefix=$3
org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0); print $1" "$2; } ;')

export -f get_fasta
#export -f check_param
export -f remove_duplicate_genes

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
#flank_len=2000 #$(grep -i -w "flank_len" parameters.txt | check_param)   #2000
REMOVE_DOWNLOADS=$(grep -i -w "remove_downloads" parameters.txt | check_param) 
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | check_param) 
#GENES_META=$(grep -i -w "genes_meta" parameters.txt | check_param) 
seqID_delimiter=$(grep -i -w "seqID_delimiter" parameters.txt | check_param) 
LABEL_SEQS=$(grep -i -w "label_sequence_IDs" parameters.txt | check_param) 
#PY2_PATH=$(grep -i -w "python2_path" parameters.txt | check_param)
PY3_PATH=$(grep -i -w "python3_path" parameters.txt | check_param)
PY3_PATH="${PY3_PATH/#\~/$HOME}"
#n_threads=$(nproc --all)
max_block_size_mb=$(grep -i -w "max_block_size_mb" parameters.txt | check_param) 
n_threads=$(grep -i -w "max_concurrent_jobs" parameters.txt | check_param)

MODE=""
if [[ $GENE_SEARCH_MODE=="HARD" ]]; then
  MODE="-w"
fi

#echo $org_name
echo "Command : $0 $@"
#echo "Threads:$n_threads, File Block size(MB):$block_size_mb"

if [[ $CLEAN_EXTRACT ==  "TRUE" ]] ; then
  rm -rf $FASTA_PATH/$5
  rm -rf files/genes/$f_org_name
  rm -rf $TEMP_PATH/$f_org_name
  rm -f $GENOME_FILE.fai
  rm -f files/genes/$f_org_name/gtf_stats.csv
  rm -f files/genes/$f_org_name/1.list
  rm -f files/genes/$f_org_name/2.list
  rm -f files/genes/$f_org_name/full.list
  rm -f files/genes/$f_org_name/odb.list
  rm -f files/genes/$f_org_name/final.list
fi

mkdir -p $FASTA_PATH/$f_org_name
mkdir -p files/genes/$f_org_name/odb/
mkdir -p $TEMP_PATH/$f_org_name

genome_pipe="$TEMP_PATH/$f_org_name/genome_pipe"
rm -f $genome_pipe
mkfifo $genome_pipe

gfile_name=$GENOME_FILE
if [[ ${GENOME_FILE##*.} ==  "gz" ]] ; then
  gfile_name=${GENOME_FILE%.*}
  #if [[ ! -s $gfile_name || $CLEAN_EXTRACT == "TRUE" ]] ; then
    zcat -f $GENOME_FILE > $genome_pipe & #> /dev/null &
    ##zcat -f $GENOME_FILE > $genome_pipe &
    cat < $genome_pipe > $gfile_name & #&> /dev/null &
    ##nohup zcat -f $GENOME_FILE > $gfile_name & #> /dev/null &
    genome_proc_id=$(echo $!)
  #fi
  #GENOME_FILE=$gfile_name
fi

if [[ $(echo $ANNO_FILE | grep -q -i "gtf") ]] ; then
  file_name=${ANNO_FILE%.*}
  nohup zcat -f $ANNO_FILE | gffread - -T -O -E | gzip -c - > $ANNOS_PATH/"$5".gtf.gz &> /dev/null &
  anno_proc_id=$(echo $!)
  ANNO_FILE=$ANNOS_PATH/"$5".gtf.gz
fi

>&1 color_FG_BG_Bold $Black $BG_Yellow "0.1 Building Reference Genome Index (Samtools faidx)...${Color_Off}"

if [[ ! -s $GENOME_FILE.fai ]] ; then
##    # if [[ ${GENOME_FILE##*.} ==  "gz" ]] ; then
##    #   time samtools faidx --fai-idx $GENOME_FILE.fai <(zcat $GENOME_FILE)
##    # else
##    #   time samtools faidx --fai-idx $GENOME_FILE.fai <(cat $GENOME_FILE)
##    # fi
##    #time samtools faidx --fai-idx $GENOME_FILE.fai < $genome_pipe
  time samtools faidx --fai-idx $GENOME_FILE.fai $genome_pipe & #> /dev/null &
  genome_index_proc=$(echo $!)
fi

while [[ $(lsof -R -p $$ | wc -l) -gt $(($(ulimit -Sn)-200)) ]]; #limit file descriptors for child & parent process to 1024-200 #$(lsof -R -p $$ -u $USER | wc -l) > $(($(ulimit -Sn)-200))
do
  printf %s"\t"%s"\n" "Waiting for file handles to close.." $(lsof -R -p $$ | wc -l)
  sleep 5
done

gene_list=($(cat $4 | sort | uniq | grep -v -w -i "gene"))

>&1 color_FG_BG_Bold $Black $BG_Yellow "1. Checking Gene Names & Splitting GTFs for parallel processing...${Color_Off}"

#rm -f files/genes/$f_org_name/1.list
#rm -f files/genes/$f_org_name/2.list
#rm -f files/genes/$f_org_name/full.list
#rm -f files/genes/$f_org_name/final.list
##find $TEMP_PATH/$f_org_name/ -type f -exec rm -f -- {} +


if [[ ! -z $anno_proc_id ]]; then
  wait $anno_proc_id
fi

if [[ ! -s files/genes/$f_org_name/1.list || ! -s files/genes/$f_org_name/2.list ]] ; then
  eexp_gene=$(printf -- '%s\n' "${gene_list[@]}")
  time zgrep -i $MODE -A 0 --group-separator='>' -f <(echo "${eexp_gene[@]}") $ANNO_FILE | csplit --quiet -z --suffix-format="%0d.gtf_slice" --prefix="$TEMP_PATH/$f_org_name/1." --suppress-matched - '/>/' '{*}'
  zgrep -hPo 'gene_name "\K[^"]+' $TEMP_PATH/$f_org_name/*.gtf_slice | sort | uniq | awk 'NF' > files/genes/$f_org_name/1.list
  echo ${gene_list[@]/($(cat files/genes/$f_org_name/1.list)))} | awk 'NF' > files/genes/$f_org_name/2.list
#else
#  >&1 color_FG $Green "1. DONE : Genes already checked ${Color_Off}"
fi

>&1 echo $(color_FG $Green "1. DONE : Available Genes : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/1.list ${Color_Off}")$(color_FG $Green ", Genes not found : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/2.list ")

#######################################################################################################

>&1 color_FG_BG_Bold $Black $BG_Yellow "2. Checking OrthoDB for missing & orthologous genes..."

if [[ ! -s files/genes/$f_org_name/odb.list || ! -s files/genes/$f_org_name/final.list || ! -s files/genes/$f_org_name/odb.final_map ]] ; then
  time ./check_OrthoDB.sh $f_org_name $4 files/genes/$f_org_name/odb.list $ANNO_FILE #1> $TEMP_PATH/$f_org_name/check_ortho.o 2> $TEMP_PATH/$f_org_name/check_ortho.e 
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

#if [[ ! -s files/genes/$f_org_name/gtf_stats.csv || ! -s files/genes/$f_org_name/final.list ]] ; then
  #echo "With R"
  #time ./jobhold.sh get_GTF_info Rscript --vanilla --verbose extract_gtf_info.R $GENES_META/$f_org_name/ $f_org_name files/genes/$f_org_name/gtf_stats.csv ; #1> $TEMP_PATH/get_GTF_info.o 2> $TEMP_PATH/get_GTF_info.e
  time Rscript --vanilla --verbose extract_gtf_info.R $TEMP_PATH/$f_org_name/ $f_org_name files/genes/$f_org_name/gtf_stats.csv 1> $TEMP_PATH/$f_org_name/get_GTF_info.o 2> $TEMP_PATH/$f_org_name/get_GTF_info.e
  sed 1d files/genes/$f_org_name/gtf_stats.csv | awk -F',' '{print $1"\n"}' | sort | uniq | awk 'NF' > files/genes/$f_org_name/final.list
#fi

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

>&1 color_FG_BG_Bold $Black $BG_Yellow "4. Fetching sequences from Genome..." #(log:$TEMP_PATH/$f_org_name/get_FASTA.[o/e])
>&1 color_FG $Yellow "Transcript Regions : $(IFS=","; echo ${TRANSCRIPT_REGIONS[*]}) "

time parallel --max-procs $n_threads "get_fasta {1} {2} $bed_prefix $gfile_name $f_org_name $FASTA_PATH/$f_org_name $ANNO_FILE $TEMP_PATH/$f_org_name $LABEL_SEQS ;" ::: ${gene_list[@]} ::: ${TRANSCRIPT_REGIONS[@]} #1> $TEMP_PATH/$f_org_name/get_FASTA.o 2> $TEMP_PATH/$f_org_name/get_FASTA.e

>&1 echo $(color_FG $Green "4. DONE : FASTA PATH : ")$(color_FG_BG_Bold $White $BG_Purple "$FASTA_PATH/$f_org_name") #$(color_FG $Green " (Logs : $TEMP_PATH/$f_org_name/get_FASTA.[o/e])")

#############################################################################################################

if [[ $LABEL_SEQS ==  "TRUE" ]] ; then
  >&1 color_FG_BG_Bold $Black $BG_Yellow "4.1 Labelling sequences..."
  tmp_names=($(parallel --link --max-procs $n_threads "echo {1},{2}" ::: ${gene_list[@]} ::: ${s_names[@]}))
  time parallel --max-procs $n_threads "printf -- %s,%s\\\n {1} {2}" ::: ${tmp_names[@]} ::: ${TRANSCRIPT_REGIONS[@]} | parallel --colsep "," --max-procs $n_threads "./label_sequenceIDs.sh $f_org_name {1} $ANNO_FILE $FASTA_PATH/$f_org_name/{2}.{3} $FASTA_PATH/$f_org_name/{2}.{3}.tmp files/genes/$f_org_name/odb.final_map" 
#else
#  time parallel --link --max-procs $n_threads "mv -f $FASTA_PATH/$f_org_name/{1}.{2}.tmp $FASTA_PATH/$f_org_name/{1}.{2}" ::: ${gene_list[@]} ::: ${s_names[@]}
fi

#######################################################################################################

>&1 color_FG_BG_Bold $Black $BG_Yellow "5. Generating Metadata and Cleaning up..."
#find files/genes/$5 -iname "*.rna_list" -empty | awk -F'/' '{print $NF}' | sed 's/.rna_list//g' > $FASTA_PATH/$5/MISSING_GENES
#grep -v -i -f $FASTA_PATH/$5/MISSING_GENES $4 | sort > $FASTA_PATH/$5/AVAILABLE_GENES
sed 1d files/genes/$f_org_name/gtf_stats.csv | awk -F',' '{print $2}' | sort | uniq >  files/genes/$f_org_name/AVAILABLE_GENES
grep -v -i -f files/genes/$f_org_name/AVAILABLE_GENES $4 | sort | uniq > files/genes/$f_org_name/MISSING_GENES

find $FASTA_PATH/$f_org_name/ -type f -name "*.fai" -exec rm -f {} +

if [[ ! -z $gfile_name ]]; then
  rm -f $gfile_name
fi

rm -rf $TEMP_PATH/$f_org_name/
if [[ $REMOVE_DOWNLOADS ==  "TRUE" ]] ; then
   rm $GENOME_FILE
   rm $ANNO_FILE
fi

>&1 color_FG_BG_Bold $Purple $BG_White "Extraction DONE for organism : $f_org_name"

#wait $odb_proc_id

# ##Populate clusters(files/genes/<org>/ALL_CLUSTERS) and group FASTA sequences according to clusters and store it in files/groups/
# #grep -r ">" $FASTA_PATH/$f_org_name/ | parallel --max-procs $n_threads --colsep ":>" --recend "\n" echo {2} | awk -F"$seqID_delimiter" '{print $5}' |  sort | uniq > files/genes/$f_org_name/ALL_CLUSTERS
# grep -r ">" $FASTA_PATH/$f_org_name/ | awk -F ":>" '{print $2}' | awk -F"$seqID_delimiter" '{print $5}' |  sort | uniq > files/genes/$f_org_name/ALL_CLUSTERS

# #rm $TEMP_PATH/$f_org_name.*.*
# #rm files/genes/$5/*.gtf_slice

# #rm $3*
# #rm files/bed/$5*

# #exit 1
