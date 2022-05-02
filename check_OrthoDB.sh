#!/bin/bash
#$ -cwd
#$ -V
#$ -l data

# # INTERNAL SCRIPT FOR RUNNING IN PARALLEL
# $1 - reference organism
# $2 - gene list
# $3 - OrthoDB organism ID
# $4 - f_org_name
# $5 - ANNO_FILE

f_org_name=$4
ref_org=$1  # =$(</dev/stdin)
readarray gene_list < $2
org_ID=$3
ANNO_FILE=$5
#printf $1 $2 $3 $4 $5 $6

if [[ $ref_org == $f_org_name || -z $ref_org ]]; then
  #SAME ORGANISM or empty variable, exit
  exit 0
fi

TEMP_PATH=$(grep -i -w "temp_path" parameters.txt | awk -F'=' '{print $2}') 
ORTHODB_PATH=$(grep -i -w "orthodb_files_path" parameters.txt | awk -F'=' '{print $2}')
ORTHODB_PREFIX=$(grep -i -w "orthodb_prefix" parameters.txt | awk -F'=' '{print $2}') 
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'=' '{print $2}') 
GENES_META=$(grep -i -w "genes_meta" parameters.txt | awk -F'=' '{print $2}')
n_threads=$(nproc --all)

#echo $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab $ORTHODB_PATH/"$ORTHODB_PREFIX"_species.tab
#echo $s_name $gene_name

#Check if files exist
if [[ ! -s $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab || ! -s $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab || ! -s $ORTHODB_PATH/"$ORTHODB_PREFIX"_species.tab ]] ; then
  echo "OrthoDB files missing"
  echo "Require: _OG2genes.tab, _genes.tab, _species.tab"
  exit -1
fi

name_strip=($(echo ${gene_full[@]} | sed 's/.$//'))
s_names=$(parallel --max-procs $n_threads "echo {} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' " ::: ${gene_list[@]})
genes_strip=$(parallel --max-procs $n_threads "echo {} | awk -F'_' '{print $(echo '$1')}' " ::: ${s_names[@]})
lookup_genes=($(echo "${s_names[@]}" "${genes_strip[@]}" "${name_strip[@]}" | sort | uniq ))

##Saving the perl exp to files so we can do a file grep
#parallel "printf %s\\\n \\\b\({}\)\\\b\\\t " ::: ${lookup_genes[@]} > files/genes/$f_org_name/genes.pexp 
#echo ${lookup_genes[@]} | awk '{print "^"$0 }' > files/genes/$f_org_name/genes.rexp 
#parallel "printf %s\\\n \\\b{}" ::: ${lookup_genes[@]} > files/genes/$f_org_name/genes.rexp 

ref_org_name=$(echo $ref_org | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')

#echo $ref_org_name

#Get the OrgID for the organism from OrthoDB
ref_org_ID=$(grep -i -P "\t\b($ref_org_name)\b\t" $ORTHODB_PATH/"$ORTHODB_PREFIX"_species.tab | awk '{print $2}')
if [[ $ref_org_ID ==  "" ]] ; then
  ref_org_ID=$(grep -i -P "($ref_org_name)" "$ORTHODB_PATH/$ORTHODB_PREFIX"_species.tab | awk '{print $2}')
fi

#echo $ref_org_ID
#echo "here1"
##FIND THE GeneIDs from OrthoDB for the SUBJECT reference organism
#awk "IGNORECASE = 1;/^$ref_org_ID/" $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab | awk '{print $1"\t"$4}' > "$TEMP_PATH/$ref_org.$f_org_name.tmp"
if [[ ! -s "$TEMP_PATH/$ref_org.$f_org_name.tmp" ]] ; then
grep -i -w $ref_org_ID $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab | awk '{print $1"\t"$4}' > "$TEMP_PATH/$ref_org.$f_org_name.tmp"
fi
#grep -i -w $ref_org_ID $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab | awk '{print $1"\t"$4}' | grep -i -P "\b($gene_name)\b\t" | awk '{print $1}' > "$TEMP_PATH/$ref_org.$f_org_name.$s_name.ref_ID"
#parallel --max-procs $n_threads "printf %s\\\n \\\b{}" ::: ${lookup_genes[@]} | grep -E -f - <(awk '{print $1"\t"$4}' <(grep -i -w $ref_org_ID $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab)) > "$TEMP_PATH/$ref_org.$f_org_name.ref_ID"
if [[ ! -s "files/genes/$f_org_name/refs/$ref_org.$f_org_name.ref_ID" ]] ; then
parallel --max-procs $n_threads -j0 "printf %s\\\n \\\b{}" ::: ${lookup_genes[@]} | grep -E -f - "$TEMP_PATH/$ref_org.$f_org_name.tmp" >  "files/genes/$f_org_name/refs/$ref_org.$f_org_name.ref_ID" #"$TEMP_PATH/$ref_org.$f_org_name.ref_ID"
fi
##IF the gene exists in OrthoDB,
if [[ $(awk 'END{print NR;}' "files/genes/$f_org_name/refs/$ref_org.$f_org_name.ref_ID") !=  0 && ! -s "files/genes/$f_org_name/refs/$ref_org.$f_org_name.clusters" ]] ; then
  #printf "$ref_org_name\t$ref_org_ID\n"
  
  ##Find the orthologous clusters in which the (SUBJECT) genes participate and extract them
  #time awk '{print $1}' "$TEMP_PATH/$ref_org.$f_org_name.ref_ID" | grep -w -f - $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | sort | uniq > "$TEMP_PATH/$ref_org.$f_org_name.clusters"
  parallel --pipe-part --block 10M --max-procs $n_threads --recend "\n" -j $n_threads -a $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab "grep -w -f <(awk '"'{print $1}'"' files/genes/$f_org_name/refs/$ref_org.$f_org_name.ref_ID"  | awk '{print $1}' | sort | uniq > "files/genes/$f_org_name/refs/$ref_org.$f_org_name.clusters"
fi
  #awk '{print $1}' "$TEMP_PATH/$ref_org.$f_org_name.ref_ID" | parallel --max-procs $n_threads "grep -w {} $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | sort | uniq" >> "$TEMP_PATH/$ref_org.$f_org_name.clusters"
  #awk '{print $1}' "$TEMP_PATH/$ref_org.$f_org_name.ref_ID" | parallel --progress --max-procs $n_threads "awk '/^{}/' $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | sort | uniq" >> "$TEMP_PATH/$ref_org.$f_org_name.clusters"
  #cat "$TEMP_PATH/$ref_org.$f_org_name.ref_ID" | parallel --colsep '\t' --progress --max-procs $n_threads "awk '/^{1}/' $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | sort | uniq" >> "$TEMP_PATH/$ref_org.$f_org_name.clusters"
  
  #cat "$TEMP_PATH/$ref_org.$f_org_name.$s_name.clusters"
#echo "here2"
  ##In respective clusters, Get the GeneIDs of (QUERY) genes in the clusters where (SUBJECT) gene is also participating in
  if [[ $(awk 'END{print NR;}' "files/genes/$f_org_name/refs/$ref_org.$f_org_name.clusters") !=  0 &&  ! -s "files/genes/$f_org_name/refs/$ref_org.$f_org_name.genes" ]] ; then
  #time awk '{print $1}' "$TEMP_PATH/$ref_org.$f_org_name.clusters" | grep -w -f - $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | grep -w $org_ID | sort | uniq > "$TEMP_PATH/$ref_org.$f_org_name.genes"
  parallel --pipe-part --block 10M --max-procs $n_threads --recend "\n" -j $n_threads -a $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab "grep -w $org_ID | grep -w -f files/genes/$f_org_name/refs/$ref_org.$f_org_name.clusters " | sort | uniq > "files/genes/$f_org_name/refs/$ref_org.$f_org_name.genes"
fi
  #awk '{print $1}' "$TEMP_PATH/$ref_org.$f_org_name.clusters" | parallel --max-procs $n_threads  "grep -w {} $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | grep -w $org_ID | sort | uniq" >> "$TEMP_PATH/$ref_org.$f_org_name.genes"
  #awk '{print $1}' "$TEMP_PATH/$ref_org.$f_org_name.clusters" | parallel --progress --max-procs $n_threads  "awk '/^{}/' $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | grep -w $org_ID | sort | uniq" >> "$TEMP_PATH/$ref_org.$f_org_name.genes"
  #cat "$TEMP_PATH/$ref_org.$f_org_name.clusters" | parallel --colsep '\t' --progress --max-procs $n_threads  "awk '/^{1}/' $ORTHODB_PATH/"$ORTHODB_PREFIX"_OG2genes.tab | grep -w $org_ID | sort | uniq" >> "$TEMP_PATH/$ref_org.$f_org_name.genes"
  
# echo "here3" 
  #Convert GeneIDs to gene name
  if [[ $(awk 'END{print NR;}' "files/genes/$f_org_name/refs/$ref_org.$f_org_name.genes") !=  0 && ! -s "files/genes/$f_org_name/refs/$ref_org.$f_org_name.names" ]] ; then
  parallel --pipe-part --block 10M --max-procs $n_threads --colsep '\t' --recend "\n" -j $n_threads -a $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab "grep -w -f <(cut -f 2 files/genes/$f_org_name/refs/$ref_org.$f_org_name.genes) " | awk '{a[$1]=$4} END{for ( key in a ) print key"\t"a[key]}' > "files/genes/$f_org_name/refs/$ref_org.$f_org_name.names"
 fi
  #awk '{print $2}' "$TEMP_PATH/$ref_org.$f_org_name.genes" | parallel --max-procs $n_threads "grep -w {} $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab" >> "$TEMP_PATH/$ref_org.$f_org_name.names"
  #awk '{print $2}' "$TEMP_PATH/$ref_org.$f_org_name.genes" | parallel --progress --max-procs $n_threads "awk '/^{}/' $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab" >> "$TEMP_PATH/$ref_org.$f_org_name.names"
  #cat "$TEMP_PATH/$ref_org.$f_org_name.genes" | parallel --colsep '\t' --progress --max-procs $n_threads "awk '/^{2}/' $ORTHODB_PATH/"$ORTHODB_PREFIX"_genes.tab" >> "$TEMP_PATH/$ref_org.$f_org_name.names"
  
# echo "here4" 
  #Check if these genes are present in the GTF/GFF(annotation) file
  if [[ $(awk 'END{print NR;}' "files/genes/$f_org_name/refs/$ref_org.$f_org_name.names" | awk '{print $1}') !=  0 ]] ; then
  #cat "$TEMP_PATH/$ref_org.$f_org_name.names"
  #perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq
  #awk '{print $4}' "$TEMP_PATH/$ref_org.$f_org_name.names" | grep -i -f - $ANNO_FILE > "$TEMP_PATH/$ref_org.$f_org_name.gtf_slice" # files/genes/$f_org_name/$file_out.rna_list
  # && ! -s $GENES_META/$f_org_name/{2}.$ref_org.gtf_slice
  awk '{ a=$2; gsub(/[[:punct:]]/,"_",$2); print a"\t"$2;}' "files/genes/$f_org_name/refs/$ref_org.$f_org_name.names" | sort | uniq | parallel --max-procs $n_threads --colsep '\t' "if [[ ! -z {1} && ! -z {2} ]] ; then grep -i {1} $ANNO_FILE > $GENES_META/$f_org_name/{2}.$ref_org.gtf_slice  
fi; "
# && ! -s $GENES_META/$f_org_name/{2}.$ref_org.rna_list
  awk '{ a=$2; gsub(/[[:punct:]]/,"_",$2); print a"\t"$2;}' "files/genes/$f_org_name/refs/$ref_org.$f_org_name.names" | sort | uniq | parallel --max-procs $n_threads --colsep '\t' "if [[ ! -z {1} && ! -z {2} ]] ; then grep -w {1} $GENES_META/$f_org_name/{2}.$ref_org.gtf_slice | perl -lne 'print @m if @m=(/((?:transcript_id)\s+\S+)/g);' | cut -f 2 -d ' ' | cut -f 2 -d '\"' | sort | uniq > $GENES_META/$f_org_name/{2}.$ref_org.rna_list 
   fi; " #$TEMP_PATH/$ref_org.$f_org_name.{2};
  # | perl -lne 'print @m if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sort | uniq > $TEMP_PATH/$ref_org.$f_org_name.{2}
#echo "here5"
fi

#rm $TEMP_PATH/$ref_org.$f_org_name.$s_name.*
#cat "$TEMP_PATH/$ref_org.$f_org_name.$s_name"
#awk '{print $1"\t"$2"\t"$3}' "$TEMP_PATH/$ref_org.$f_org_name.$s_name.names" "$TEMP_PATH/$ref_org.$f_org_name.$s_name"
exit 0