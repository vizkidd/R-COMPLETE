#!/bin/bash
#$ -cwd
#$ -V
#$ -l data

# # INTERNAL SCRIPT FOR finding orthologous genes
# $1 - reference organism
# $2 - gene list (to look for orthologous genes)
# $3 - output file for gene list from OrthoDB
# $4 - ANNO_FILE

##FUNCTIONS
function remove_duplicate_genes() {
  #$1 - gene_list
  #$2 - TEMP_PATH
  cat $1 | sort | uniq > $1.tmp
  cat $1.tmp > $1
  rm -f $1.tmp
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


f_org_name=$1
#readarray gene_list < $2
gene_list=($(cat $2 | grep -v -w -i "gene"))
ANNO_FILE=$4

#export -f remove_duplicate_genes
#export -f check_param

source fasta_functions.sh

TEMP_PATH=$(grep -i -w "temp_path" parameters.txt | check_param) 
ORTHODB_PATH_PREFIX=$(grep -i -w "orthodb_path_prefix" parameters.txt | check_param) 
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | check_param) 
GENE_SEARCH_MODE=$(grep -i -w "gene_search_mode" parameters.txt | check_param)
n_threads=$(grep -i -w "max_concurrent_jobs" parameters.txt | check_param)
ODB_ID_DELIM="||"

MODE=""
if [[ $GENE_SEARCH_MODE=="HARD" ]]; then
  MODE="-w"
fi

mkdir -p files/genes/$f_org_name/odb/

#Check if files exist
if [[ ( ! -s "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz || ! -s "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz ) && ! -s "$ORTHODB_PATH_PREFIX"_species.tab.gz ]] ; then
  >&2 color_FG_Bold $Red "2. OrthoDB files missing/corrupt"
  >&2 color_FG_Bold $Red "2. Require: "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz, "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz, "$ORTHODB_PATH_PREFIX"_species.tab.gz"
  exit -1
fi

if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz ; then
  if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz ; then
    exit -1
  else
    ODB_FILE="$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz
    >&1 echo $(color_FG $Yellow "2. Selected Fixed ODB (User) file : ")$(color_FG_BG_Bold $White $BG_Purple "$ODB_FILE")
  fi
  >&2 color_FG_Bold $Red "2. Error in ODB files! Rerun merge_OG2genes_OrthoDB.sh"
  exit -1
else
  ODB_FILE="$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
  >&1 echo $(color_FG $Yellow "2. Selected Fixed ODB file : ")$(color_FG_BG_Bold $White $BG_Purple "$ODB_FILE")
fi

readarray refs < $REF_ORGS

s_names=($(awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo "${gene_list[@]}")))
genes_strip=($(awk -F'_' '{print $(echo '$1')}' <(echo "${s_names[@]}")))
lookup_genes=($(echo "${genes_strip[@]}" "${gene_list[@]}" | sort | uniq | awk 'NF' ))
org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')

org_ID=$(zgrep -i -P "\t\b($org_name)\b\t" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
if [[ -z $org_ID ]] ; then
  org_ID=$(zgrep -i -P "($org_name)" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
fi

if [[ ! -z $org_ID ]] ; then

  for ref in "${refs[@]}"; do 
    ref_org_name=$(echo $ref | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
    ref_org_ID=""
    ref_org_ID=$(zgrep -i -P "\t\b($ref_org_name)\b\t" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
    if [[ -z $ref_org_ID ]] ; then
      ref_org_ID=$(zgrep -i -P "($ref_org_name)" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
    fi
    if [[ ! -z $ref_org_ID && $ref_org_ID != $f_org_name ]] ; then
      ref_org_IDs+=($ref_org_ID)
    fi
  done

zgrep -f <(printf -- '%s\n' "${ref_org_IDs[@]}") $ODB_FILE | gzip -c > $TEMP_PATH/$f_org_name/odb.clusters.gz
zcat -f $TEMP_PATH/$f_org_name/odb.clusters.gz | awk '{print $2}' | sed 's/,/\n/g' | grep -w "$org_ID" - | grep -i $MODE -f <(printf -- '%s\n' "${lookup_genes[@]}") - | awk '{split($0,a,/\|\|/); print a[2]}' | sort -u > $3
>&1 echo $(color_FG $Green "2. DONE: Found Orthologous genes : ")$(color_FG_BG_Bold $White $BG_Purple "$3")

zgrep "$org_ID" $TEMP_PATH/$f_org_name/odb.clusters.gz | awk '{split($3,a,","); for(key in a){print $1"\t"a[key]}}' | grep -i $MODE -f <(cat $2 $3 | grep -v -w -i "gene") | sort -u | awk '{if($1 in a){a[$1]=a[$1]","$2}else{a[$1]=$2;} } END{for(key in a){print key"\t"a[key]}}' > files/genes/$f_org_name/odb.final_map

>&1 echo $(color_FG $Green "2. DONE: Mapped gene names to clusters : ")$(color_FG_BG_Bold $White $BG_Purple "files/genes/$f_org_name/odb.final_map")
rm $TEMP_PATH/$f_org_name/odb.clusters.gz

else
  >&2 color_FG_Bold $Red "2. $org_name not found in OrthoDB files"
fi
#>&1 color_FG $Green "2. Checking OrthoDB....Done! "
exit 0