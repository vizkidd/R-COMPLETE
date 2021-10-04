#!/bin/bash
##INTERNAL SCRIPT FOR RUNNING IN PARALLEL
# $1 - gene name
# $2 - gene name
# $3 - OrthoDB organism ID
# $4 - f_org_name
# $5 - s_name
# $6 - ANNO_FILE

f_org_name=$4
s_name=$5
ref_org=$1  # =$(</dev/stdin)
gene_name=$2
org_ID=$3
ANNO_FILE=$6
#printf $1 $2 $3 $4 $5 $6

TEMP_PATH=$(grep -i -w "temp_path" parameters.txt | awk -F'=' '{print $2}') 
ORTHODB_PATH=$(grep -i -w "orthodb_files_path" parameters.txt | awk -F'=' '{print $2}')
ORTHODB_PREFIX=$(grep -i -w "orthodb_prefix" parameters.txt | awk -F'=' '{print $2}') 
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'=' '{print $2}') 

ref_org_name=$(echo $ref_org | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
ref_org_ID=$(grep -i -P "\t\b($ref_org_name)\b\t" $ORTHODB_PATH/"$ORTHODB_PREFIX"_species.tab | awk '{print $2}')
if [[ $ref_org_ID ==  "" ]] ; then
  ref_org_ID=$(grep -i -P "($ref_org_name)" "$ORTHODB_PATH/$ORTHODB_PREFIX"_species.tab | awk '{print $2}')
fi
grep -i -w $ref_org_ID $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab | grep -i -P "\b($gene_name)\b\t" | awk '{print $1}' > "$TEMP_PATH/$ref_org.$f_org_name.$s_name.ref_ID"
#printf $ref_org_ID $ref_org_name $ORTHODB_PATH
cat "$TEMP_PATH/$ref_org.$f_org_name.$s_name.ref_ID"
if [[ $(wc -l "$TEMP_PATH/$ref_org.$f_org_name.$s_name.ref_ID" | awk '{print $1}') !=  0 ]] ; then
  #printf "$ref_org_name\t$ref_org_ID\n"
  grep -w -f "$TEMP_PATH/$ref_org.$f_org_name.$s_name.ref_ID" $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | awk '{print $1}' | sort | uniq > "$TEMP_PATH/$ref_org.$f_org_name.$s_name.clusters"
  #cat "$TEMP_PATH/$ref_org.$f_org_name.$s_name.clusters"
  if [[ $(wc -l "$TEMP_PATH/$ref_org.$f_org_name.$s_name.clusters" | awk '{print $1}') !=  0 ]] ; then
  grep -w -f "$TEMP_PATH/$ref_org.$f_org_name.$s_name.clusters" $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | grep -w $org_ID | awk '{print $2}' | sort | uniq > "$TEMP_PATH/$ref_org.$f_org_name.$s_name.genes"
  if [[ $(wc -l "$TEMP_PATH/$ref_org.$f_org_name.$s_name.genes" | awk '{print $1}') !=  0 ]] ; then
  grep -w -f "$TEMP_PATH/$ref_org.$f_org_name.$s_name.genes" $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab | awk '{print $4}' > "$TEMP_PATH/$ref_org.$f_org_name.$s_name.names"
  if [[ $(wc -l "$TEMP_PATH/$ref_org.$f_org_name.$s_name.names" | awk '{print $1}') !=  0 ]] ; then
  #cat "$TEMP_PATH/$ref_org.$f_org_name.$s_name.names"
  grep -i -f "$TEMP_PATH/$ref_org.$f_org_name.$s_name.names" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > "$TEMP_PATH/$ref_org.$f_org_name.$s_name" # files/genes/$f_org_name/$file_out.rna_list
fi
fi
fi
fi
#rm $TEMP_PATH/$ref_org.$f_org_name.$s_name.*