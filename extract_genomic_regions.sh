#!/bin/bash
#$ -cwd
#$ -V
#$ -l data
#$ -pe smp 128
#$ -l m_mem_free 8G

# $1 = genome fasta file
# $2 = formatted gtf/gff3(UNTESTED) file
# $3 = basename for bedfiles
# $4 = input gene list
# $5 = org name(eg x_laevis)

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
remove_duplicate_genes() {
  #$1 - gene_list
  #$2 - TEMP_PATH
  cat $1 | sort | uniq > $1.tmp
  cat $1.tmp > $1
  rm -f $1.tmp
}

store_RNA_list() {
  #$1 - Full Gene name(Unstripped)
  #$2 - Safe Gene name(Stripped - without punctuation) - for file storage
  #$3 - {(Genes metadata path (For storing *.gtf_slice and *.rna_list))/(Name of the organism)}
  zgrep -i -w "$1" $3/$2.gtf_slice | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > $3/$2.rna_list
}

store_GTF_block() {
  #$1 - File to store the gene name in
  #$2 - Genes metadata path (For storing *.gtf_slice and *.rna_list)
  #$3 - Number of threads
  #$4 - TEMP PATH

   #echo $1 $2
   #printf '%s\n' "$@"
   # tmp_fifo=$(mktemp -p $4)
   local tmp_file=$(mktemp -u $4/XXXXXXXX.tmp)
   #mkfifo $tmp_fifo
   cat > $tmp_file 
   # cat $tmp_fifo
   #IFS= read -a gtf_block <<< $(cat)
   #gtf_block=($(cat))
   #echo "${gtf_block[@]}"
   #printf '%s' "${gtf_block[@]}"
   #echo $1 $2
   #exit
   local g_name=($(cat $tmp_file | zgrep -Po 'gene_name "\K[^"]+' | sed 's/ //g' | sort | uniq))
   #echo ${g_name[@]}
   if [[ "${#g_name[@]}" > 1 ]] ; then
      zcat -f $tmp_file | store_GTF_blocks $1 $2 $3 $4 -nt  #parallel --pipe --recstart '' --recend '>' --compress --compress-program gzip -N 1 -L 1 --tmpdir $4 --max-procs $3 "store_GTF_block $1 $2 $3 $4"
      rm -f $tmp_file
      return 0;
   elif [[ "${#g_name[@]}" == 1 ]] ; then
      ##s_name=$(parallel --keep-order --max-procs $5 "echo {} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' " ::: ${g_name[@]})
      local s_name=$(echo $g_name | awk -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1') #($(awk -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo ${g_name[@]})))
      cat $tmp_file | sed 's/>//g' | gzip -c > $2/$3/$s_name.gtf_slice
      #cat | gzip -c > $2/$3/$s_name.gtf_slice
      #mv -f $tmp_file $2/$s_name.gtf_slice
      rm -f $tmp_file
      if [[ -s $2/$s_name.gtf_slice ]]; then 
           echo $g_name
           echo $g_name >> $1; 
           store_RNA_list $g_name $s_name $2 ;
         else
           rm -f $2/$s_name.gtf_slice $tmp_file ;
      fi
   fi
  #rm -f $tmp_file
}

store_GTF_blocks() {
  #$1 - File to store the gene name in
  #$2 - Genes metadata path (For storing *.gtf_slice and *.rna_list)
  #$3 - Number of threads
  #$4 - TEMP PATH
  #$5 - -nt/-ct (no temp file/cat to temp file)
  #$6 - IF -ct is given then, $6 is the temp file

  # g_name=($(zgrep -Po 'gene_name "\K[^"]+' $1 | sort | uniq))
  # ##s_name=$(parallel --keep-order --max-procs $5 "echo {} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' " ::: ${g_name[@]})
  # #s_name=($(awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo ${g_name[@]})))
  # ##g_name=($(cat $1 | grep -Po 'gene_name "\K[^"]+' | sort | uniq | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;'))
  # ###g_name=$(cat $1 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | sort | uniq | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
  # #echo "${g_name[@]}" 
  # #echo "----------"
  # if [[ "${#g_name[@]}" > 1 ]] ; then
  #   #echo "parallel" 
  #   #parallel --max-procs $5 --link "zgrep -w {1} $1 | zcat -  $3/$4/{2}.gtf_slice 1> >(gzip -c > $3/$4/{2}.gtf_slice) 2> /dev/null ; 
  #   #touch $3/$4/{2}.gtf_slice
     #cat | parallel --pipe --recstart '\n' --recend '>' --compress --compress-program gzip -N 1 --tmpdir $6 --max-procs $5 "store_GTF_block $2 $3 $4 $6" #--compress --compress-program gzip --tmpdir $6  #--recstart '\n'
   if [[ $5 == "-nt" ]]; then
     cat | parallel --pipe --recstart '' --recend '>' --compress --compress-program gzip --tmpdir $4 -N 1 --max-procs 1 "store_GTF_block $1 $2 $3 $4" #--recstart '\n' -1 -N 1 -L 1
     return 0
   elif [[ $5 == "-ct" ]]; then
      #local tmp_file=$(mktemp -u $4/XXXXXXXX.tmp)
      ##mkfifo $tmp_fifo
      #cat > $tmp_file
      parallel --pipe-part --recstart '' --recend '>' --compress --compress-program gzip -a $6 --tmpdir $4 --max-procs 1 "store_GTF_blocks $1 $2 $3 $4 -nt " #--block "-$3"

      #rm -f $tmp_file
      return 0
   fi
     
  # else
  # #echo "one" 
  # #cat $1 >> $3/$4/${s_name[@]}.gtf_slice
  #   cat $1 | store_GTF_block $3/$4/ $2
  # fi

}

# store_GTF_blocks() {
#   #$1 - Part(Slice) of GTF as a filename
#   #$2 - File to store the gene name in
#   #$3 - Genes metadata path (For storing *.gtf_slice and *.rna_list)
#   #$4 - Name of the organism
#   #$5 - Number of threads
#   #$6 - TEMP PATH

#   #cat $1
#   g_name=($(zgrep -Po 'gene_name "\K[^"]+' $1 | sort | uniq))
#   ##s_name=$(parallel --keep-order --max-procs $5 "echo {} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' " ::: ${g_name[@]})
#   s_name=($(awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo ${g_name[@]})))
#   ##g_name=($(cat $1 | grep -Po 'gene_name "\K[^"]+' | sort | uniq | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;'))
#   ###g_name=$(cat $1 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | sort | uniq | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
#   #echo "${g_name[@]}" 
#   #echo "----------"
#   if [[ "${#g_name[@]}" > 1 ]] ; then
#     #echo "parallel" 
#     #parallel --max-procs $5 --link "zgrep -w {1} $1 | zcat -  $3/$4/{2}.gtf_slice 1> >(gzip -c > $3/$4/{2}.gtf_slice) 2> /dev/null ; 
#     #touch $3/$4/{2}.gtf_slice
#     parallel --max-procs $5 --link "zgrep -w {1} $1 | zcat -f -  $3/$4/{2}.gtf_slice 1>  >(gzip -f -c > $3/$4/{2}.gtf_slice) 2> /dev/null ; 
#      if [[ -s $3/$4/{2}.gtf_slice ]]; then 
#        echo {1} >> $2; 
#        store_RNA_list {1} {2} $3/$4/ ;
#      else
#        rm -f $3/$4/{2}.gtf_slice;
#      fi" ::: ${g_name[@]} ::: ${s_name[@]} #--compress --compress-program gzip --tmpdir $6 
#   else
#   #echo "one" 
#   #cat $1 >> $3/$4/${s_name[@]}.gtf_slice
#     cat $1 | zcat -f -  $3/$4/{2}.gtf_slice 1> >(gzip -c > $3/$4/${s_name[@]}.gtf_slice) 2> /dev/null ;
#     if [[ -s $3/$4/${s_name[@]}.gtf_slice ]]; then 
#       echo ${g_name[@]} >> $2
#       store_RNA_list ${g_name[@]} ${s_name[@]} $3/$4/
#     else
#       rm -f $3/$4/${s_name[@]}.gtf_slice
#     fi
#   fi

# }

# # store_GTF_blocks() {
# #   #$1 - Part(Slice) of GTF as a filename
# #   #$2 - File to store the gene name in
# #   #$3 - Genes metadata path (For storing *.gtf_slice and *.rna_list)
# #   #$4 - Name of the organism
# #   #$5 - Number of threads
# #   #$6 - TEMP PATH
# # #GTF_block=(`< $1`)
# # ##IFS=$'\n' read -r -a GTF_block < "$1"
# # #echo "${GTF_block[@]}"

# # #cat $1
# # #echo "------"
# # #return 0
# # group_sep_count=$(zgrep -c ">" $1)
# # #echo "here"
# # #echo "$group_sep_count"
# # #echo $1 $2 $3 $4 $5 $6
# # if [[ $group_sep_count > 1 && ! -z $group_sep_count ]]; then
# #   #echo "count not 1"
# #   parallel --max-procs $5 -a $1 --pipe-part --tmpdir $6 --recstart '\n' --recend '>' --block -1 --compress --compress-program gzip --cat "store_GTF_blocks {} $2 $3 $4 $5 $6"
# #   return 0
# # else
# #   #cat $1
# #   #g_name=($(zgrep -Po 'gene_name "\K[^"]+' $1 | sort | uniq))
# #   ##s_name=$(parallel --keep-order --max-procs $5 "echo {} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' " ::: ${g_name[@]})
# #   #s_name=($(awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo ${g_name[@]})))
# #   ##g_name=($(cat $1 | grep -Po 'gene_name "\K[^"]+' | sort | uniq | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;'))
# #   ###g_name=$(cat $1 | perl -lne 'print "@m" if @m=(/((?:gene_name)\s+\S+)/g);' | sort | uniq | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
# #   ##echo "${g_name[@]}"
# #   #if [[ "${#g_name[@]}" > 1 ]] ; then
# #   #  echo "parallel" 
# #   #  parallel --max-procs $5 --link --tmpdir $6  "zgrep -w {1} $1 | gzip -c  > $3/$4/{2}.gtf_slice; 
# #   #   if [[ -s $3/$4/{2}.gtf_slice ]]; then 
# #   #     echo {1} >> $2; 
# #   #     store_RNA_list {1} {2} $3/$4/ ;
# #   #   else
# #   #     rm -f $3/$4/{2}.gtf_slice;
# #   #   fi" ::: ${g_name[@]} ::: ${s_name[@]} #--compress --compress-program gzip
# #   #else
 
# #   g_name=($(zgrep -Po 'gene_name "\K[^"]+' $1 | uniq ))
# #   s_name=($(echo ${g_name[@]} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1'))
# #   #echo "one" 
# #   #cat $1 >> $3/$4/${s_name[@]}.gtf_slice
# #   cat $1 | gzip -c >> $3/$4/${s_name[@]}.gtf_slice
# #   if [[ -s $3/$4/${s_name[@]}.gtf_slice ]]; then 
# #     echo ${g_name[@]} >> $2
# #     store_RNA_list ${g_name[@]} ${s_name[@]} $3/$4/
# #   else
# #     rm -f $3/$4/${s_name[@]}.gtf_slice
# #   fi
# #   #fi

# # fi
# # return 0
# # }

map_names_to_ODB_clusters() {
  #1 - Organism name
  #2 - Metadata output path
  
  local f_org_name=$1
  local GENES_META=$2

  echo "Additional : Mapping gene names to clusters (files/genes/$f_org_name/odb/odb.full_map)"

  if [[ ! -s  files/genes/$f_org_name/odb/odb.genes_map || ! -s files/genes/$f_org_name/odb/odb.clusters_map ]]; then
    echo "Checking OrthoDB failed, odb.genes_map or odb.cluster_map NOT FOUND in files/genes/$f_org_name/odb/"
    return 0
  fi

  #cat files/genes/$f_org_name/odb/*.genes | sort | uniq > files/genes/$f_org_name/odb/odb.clusters_map
  #cat files/genes/$f_org_name/odb/*.names | sort | uniq > files/genes/$f_org_name/odb/odb.genes_map

  #awk '{print $1"\t"$4}' "$ORTHODB_PATH/$ORTHODB_PREFIX"_genes.tab | grep "$org_ID" | grep -w -f $5 | sort -k1 | uniq >> files/genes/$f_org_name/odb/odb.genes_map
    
    ##find $GENES_META/$f_org_name/*.gtf_slice --empty --delete
    #find $GENES_META/$f_org_name/ -type f -empty -print -delete

  local gene_list=($(find $GENES_META/$f_org_name/*.gtf_slice -exec basename {} \; | cut -d '.' -f 1 | sort | uniq ))

   ##MAP GENE NAMES to clusters
   awk 'fname != FILENAME { fname = FILENAME; idx++ } FNR > 1 && idx == 1 { gene_name[$1] = $2 } FNR > 1 && idx == 2 { clusters[$2]=$1 }END{for(key in gene_name) print key"\t"clusters[key]"\t"gene_name[key];}' files/genes/$f_org_name/odb/odb.genes_map files/genes/$f_org_name/odb/odb.clusters_map > files/genes/$f_org_name/odb/odb.group_map
  ##COLLAPSE GENE name rows (in case a field has multiple records)  
  #parallel --max-procs $n_threads "printf '\(%s\)\{1\}\n' {}" ::: ${gene_list[@]} | grep -i -f - -o -n -G files/genes/$f_org_name/odb/odb.group_map | awk -F':' '{if($1 in a)a[$1]=a[$1]","$2;else a[$1]=$2 ;}END{for(key in a)print key"\t"a[key];}'  | sort -k1 -n > files/genes/$f_org_name/odb/genes.match_map
  printf '\(%s\)\{1\}\n' "${gene_list[@]}" | grep -i -f - -o -n -G files/genes/$f_org_name/odb/odb.group_map | awk -F':' '{if($1 in a)a[$1]=a[$1]","$2;else a[$1]=$2 ;}END{for(key in a)print key"\t"a[key];}'  | sort -k1 -n > files/genes/$f_org_name/odb/genes.match_map

  ##FINAL map, GENE list -> CLUSTER list -> SEARCH_GROUP list (second awk is to colllapse rows based on clusters)
  awk 'START{j=0;k=0;} fname != FILENAME { fname = FILENAME; idx++; }  FNR > 0 && idx == 1 { j=FNR; genes[j]=$3; search_group[j]=$3; clusters[j]=$2;  }  FNR > 0 && idx == 2 { k=FNR; if(match(k,$1)==0) search_group[$1]=$2; }END{if(j<k) l=k;else l=j; for(i=1; i<l;i++) print clusters[i]"\t"genes[i]"\t"search_group[i];} ' files/genes/$f_org_name/odb/odb.group_map files/genes/$f_org_name/odb/genes.match_map | awk '{if($1 in a)a[$1]=a[$1]","$2;else a[$1]=$2 ;}END{for(key in a)print key"\t"a[key];}' > files/genes/$f_org_name/odb/odb.full_map

  echo "DONE: Mapped gene names to clusters (files/genes/$f_org_name/odb/odb.full_map)
"
}

# check_OrthoDB() { #DEPRECATED
#   #1 - genelist to loopkup in OrthoDB
#   #2 - Organism name
#   #3 - FULL GTF File(not slice)
#   #4 - Metadata output path
#   #5 - Full gene list (union of lists 1, 2.1, 2.2)
#  # proc_list=()
 
#  #org_name=$(echo $2 | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
#  ANNO_FILE=$3
#  GENES_META=$4
# #ORTHODB_PATH_PREFIX=$(grep -i -w "orthodb_path_prefix" parameters.txt | awk -F'==' '{print $2}') 
# REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'==' '{print $2}') 

# #org_name=$(echo $3 | awk '{ gsub(/[[:punct:]]/, " ", $0); print $1" "$2; } ;')
# #if [[ -s files/genes/$f_org_name/$file_out.rna_list ]] ; then
# #  rm files/genes/$f_org_name/$file_out.rna_list
# #fi

# # if [[ ! -s "$ORTHODB_PATH_PREFIX"_OG2genes.tab.gz || ! -s "$ORTHODB_PATH_PREFIX"_genes.tab.gz || ! -s "$ORTHODB_PATH_PREFIX"_species.tab.gz ]] ; then
# #   echo "OrthoDB files missing/corrupt"
# #   echo "Require: $ORTHODB_PATH_PREFIX_OG2genes.tab.gz, $ORTHODB_PATH_PREFIX_genes.tab.gz, $ORTHODB_PATH_PREFIX_species.tab.gz"
# #   return -1
# # fi

# # if [[ ! -s "$ORTHODB_PATH_PREFIX"_genes.tab ]] ; then
# #   zcat "$ORTHODB_PATH_PREFIX"_genes.tab.gz > "$ORTHODB_PATH_PREFIX"_genes.tab
# # fi
# # if [[ ! -s "$ORTHODB_PATH_PREFIX"_OG2genes.tab ]] ; then
# #   zcat "$ORTHODB_PATH_PREFIX"_OG2genes.tab.gz > "$ORTHODB_PATH_PREFIX"_OG2genes.tab
# # fi

# # readarray refs < $REF_ORGS
# #readarray gene_list < $1
# org_ID=$(zgrep -i -P "\t\b($org_name)\b\t" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
# #gene_ID=$(zcat $ORTHODB_PATH/odb10v1_genes.tab.gz | grep $org_ID | grep -i -P "\b($1)\b\t" | awk '{print $1}')
# if [[ -z $org_ID ]] ; then
#   org_ID=$(zgrep -i -P "($org_name)" "$ORTHODB_PATH_PREFIX"_species.tab.gz | awk '{print $2}')
# fi

# if [[ ! -z $org_ID ]] ; then
#   #echo $1 $org_ID $org_name $s_name $ANNO_FILE  
#   #parallel -j1 --compress --pipepart -a "$REF_ORGS" --recstart '\n' --block -1 "./check_OrthoDB.sh {} $1 $org_ID $f_org_name $s_name"
#   #echo "ext_para_check"

#   parallel --max-procs 1 "if [[ $2 != {} ]] ; then
#     ./jobhold.sh check_ortho ./check_OrthoDB.sh {} $1 $org_ID $2 $ANNO_FILE
#   fi" ::: ${refs[@]}  #parallel -N 251 --pipe #-j ${#refs[@]} #$(echo $n_threads/${#refs[@]} + 1 | bc) ##echo $2 {} ;

# #mkdir -p files/genes/$f_org_name/refs


#   #cat $TEMP_PATH/*.$f_org_name.$s_name > files/genes/$f_org_name/$file_out.rna_list
#   #rm $TEMP_PATH/*.$f_org_name.$s_name
# else
#   echo "$org_name not found in OrthoDB files"
# fi
# #for id in ${proc_list[@]}
# #do
# #  echo $id
# #  wait $id
# #done
# }

# check_gene() { #DEPRECATED, just here for reference
#   gene_full=$1
#   f_org_name=$2

#   ANNO_FILE=$3
#   GENES_META=$4


#   if [[ $gene_full !=  "" ]] ; then
#         gene=$(echo $gene_full | sed 's/.$//')
#     file_out=$(echo $gene_full | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')

#   #   if [[ -s files/genes/$f_org_name/1.list && -s files/genes/$f_org_name/2.1.list && -s files/genes/$f_org_name/2.2.list && -s files/genes/$f_org_name/2.list ]] ; then
#   #   if ! grep -q -w "$gene" files/genes/$f_org_name/1.list files/genes/$f_org_name/2.1.list files/genes/$f_org_name/2.2.list files/genes/$f_org_name/2.list ; then
#   #     return 2
#   #   fi
#   # fi

#     #echo $gene_full #$gene #$f_org_name $file_out $bed_prefix
#   #Saving gene info list in a file
#     zgrep -i "$gene" $ANNO_FILE | gzip -c > $GENES_META/$f_org_name/$file_out.gtf_slice
#     if [ -s $GENES_META/$f_org_name/$file_out.gtf_slice ]; 
#     then
#       #grep -i "$gene" $2 | awk -F'\t' '{print $9}' | awk -F';' '{print $3}' | awk '{print $2}' | sed 's/"//g' | sort | uniq > files/genes/$5/$gene.rna_list
#       if grep -q -w -i "$gene_full" $GENES_META/$f_org_name/$file_out.gtf_slice ; then
#         echo "1 - Gene name has perfect match ($gene_full)"
#         echo $gene_full >> files/genes/$f_org_name/1.list
#         grep -i -w "$gene_full" $GENES_META/$f_org_name/$file_out.gtf_slice | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > $GENES_META/$f_org_name/$file_out.rna_list
#       elif grep -q -i "$gene_full" $GENES_META/$f_org_name/$file_out.gtf_slice ; then
#         #grep -i "$gene_full" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
#           echo "2.1 - Gene name(full) matches a supergroup or subgroup, will check OrthoDB ($gene_full)"
#           echo $gene_full >> files/genes/$f_org_name/2.1.list
#           echo $gene_full >> files/genes/$f_org_name/2.list
#       elif grep -q -i "$gene" $GENES_META/$f_org_name/$file_out.gtf_slice ; then
#         #grep -i "$gene_full" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
#           echo "2.2 - Gene name(stripped) matches a supergroup or subgroup, will check OrthoDB ($gene)"
#           echo $gene_full $gene >> files/genes/$f_org_name/2.2.list
#           echo $gene_full >> files/genes/$f_org_name/2.list
#       fi
#     fi
#   fi
# }


# #      else
# #        echo "3 - Gene name has a partial/no match, saving partial results"
# #        grep -i "$gene" files/genes/$f_org_name/$file_out.gtf_slice | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$f_org_name/$file_out.rna_list
# #      fi
# #    
# #    if [ ! -s files/genes/$f_org_name/$file_out.rna_list ]; then
# #      grep -i "$gene_full" files/genes/$f_org_name/$file_out.gtf_slice | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$f_org_name/$file_out.rna_list
# #    fi
# #    
# #    cat files/genes/$f_org_name/$file_out.rna_list
# #    
# #    if [ -s files/genes/$f_org_name/$file_out.rna_list ] ; then
# #      echo $gene_full
# #      ## GET GTF stats - UTR lengths, CDS_count, exon_count (both fully and partially annotated info (eg, UTR info is hidden in exon and CDS info and is extracted using exon boundaries vs CDS boundaries))
# #      Rscript extract_gtf_info.R files/genes/$f_org_name/$file_out.gtf_slice $f_org_name $file_out
# #      #time parallel -j ${#transcript_regions[@]} "get_fasta $gene_full {} $bed_prefix $GENOME_FILE $f_org_name $FASTA_PATH $file_out 'files/genes/$f_org_name/$file_out.gtf_slice' $TEMP_PATH " ::: ${transcript_regions[@]} #$8
# #      parallel "get_fasta $gene_full {} $bed_prefix $GENOME_FILE $f_org_name $FASTA_PATH $file_out $ANNO_FILE $TEMP_PATH " ::: ${transcript_regions[@]} #parallel -N 251 --pipe parallel  #time #-j ${#transcript_regions[@]} 
# #fi

get_fasta() {
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
local GENES_META=$9

#echo $(grep -i -w -f files/genes/$org_name/$s_name.rna_list $base_bed_$reg.bed)
#echo $gene_full $reg $base_bed $genome_fa $org_name $FASTA_PATH $s_name $GTF_PATH $TEMP_PATH 
#flank_len=${10}
#echo $flank_len

if [ ! -s $FASTA_PATH/$s_name.$reg ]; then #$FASTA_PATH/$file_out.cds
>&2 echo $s_name $reg
# if grep -q -w -i "$s_name" files/genes/$org_name/gtf_stats.csv ; then
#   return 2
# fi

grep -i -w -f $GENES_META/$org_name/$s_name.rna_list "$base_bed"_"$reg.bed" > $TEMP_PATH/"$s_name"_"$reg.bed"

##TO get flanks 
#grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3"

##TO find MEDIAN values
#grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | sort -n | awk '{arr[NR]=$1} END {if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}'

##TO find MEAN
#grep -i -f files/genes/some_org/cat1.rna_list files/genes/some_org/cat1.gtf_slice | grep -i "utr" | grep -i "three\|3" | awk '{print $5-$4}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'

##EXTEND regions using bedtools slop
if [[ "$reg" == "3utr" ]]; then
  ##SHIFT by 3bp to omit stop codons in 3UTRs
  #bedtools shift -s 3 -i $TEMP_PATH/$5_"$1"_$2.bed -g $TEMP_PATH/$5_genomeFile.txt ##WONT WORK, and dont need to because stop codons are part of CDS
  ##GET AVERAGE UTR LENGTH for getting FLANK LEN
  #grep -w -i -f files/genes/$5/$7.rna_list files/genes/$5/$7.gtf_slice | grep -i "utr" | grep -i "three" | awk '{print $5-$4}' > files/genes/$5/$7.3_flank_lens
  ##CHECK if *.*_flank_lens files are empty otherwise we have to extraction the UTR info manually as it wasnt provided in the annotation
  #if [[ ! -s files/genes/$5/$7.3_flank_lens ]];
  #then  
  #  ##IF empty then get flank length from calculated UTR lengths
  #  flank_len=$(grep -w -i "$7" files/genes/$5/predicted_utr_lens.csv | awk -F',' '{ sum += $9; n++ } END { if (n > 0) print sum / n; }')
  #else
  #  flank_len=$(cat files/genes/$5/$7.3_flank_lens | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }') ##USING MEAN for now
  #fi
  local flank_len=$(grep -w -i "$s_name" files/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($9 > 3) sum += $9; n++ } END { if (n > 0) print sum / n; }') ##3' UTR length must be > 3
  if [[ $flank_len == "" ]]; then
    return 2
  fi

  if [ -s $TEMP_PATH/"$s_name"_"$reg.bed" ]; then 
  	grep -i "$reg" $TEMP_PATH/"$s_name"_"$reg.bed" > $TEMP_PATH/"$s_name"_"$reg"_FLANK.bed
  	##REMOVE the line because we stored flanking regions in a seperate bed file
  	#sed -i'' "/$(echo $reg)_FLANK/d" "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed"
  	#bedtools slop -s -r $flank_len -l 0 -i "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.tmp -g "$TEMP_PATH"/"$org_name"_genomeFile.txt > "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  else
    ## if bed file is empty then we will have to take flanks from CDS
  	bedtools flank -s -r $flank_len -l 0 -i $TEMP_PATH/"$s_name"_cds.bed -g files/genes/$org_name/genomeFile.txt > $TEMP_PATH/"$s_name"_"$reg"_FLANK.bed
  	sed -i'' "s/cds/3utr_FLANK/g" $TEMP_PATH/"$s_name"_"$reg"_FLANK.bed
  fi
  bedtools getfasta -s -split -fi $genome_fa -bed $TEMP_PATH/"$s_name"_"$reg"_FLANK.bed -nameOnly -fullHeader > "$FASTA_PATH/$s_name.$reg.tmp" #NOTUSING name+ because BLAST doesn't like long IDs (and we dont want chromosome, genomic position info)
elif [[ "$reg" == "5utr" ]]; then
  ##GET AVERAGE UTR LENGTH for getting FLANK LEN
  #grep -w -i -f files/genes/$5/$7.rna_list files/genes/$5/$7.gtf_slice | grep -i "utr" | grep -i "five" | awk '{print $5-$4}' > files/genes/$5/$7.5_flank_lens
  ##CHECK if *.*_flank_lens files are empty otherwise we have to extraction the UTR info manually as it wasnt provided in the annotation
  #  if [[ ! -s files/genes/$5/$7.3_flank_lens ]];
  #then  
  #  ##IF empty then get flank length from calculated UTR lengths
  #  flank_len=$(grep -w -i "$7" files/genes/$5/predicted_utr_lens.csv | awk -F',' '{ sum += $8; n++ } END { if (n > 0) print sum / n; }')
  #else
  #  flank_len=$(cat files/genes/$5/$7.5_flank_lens | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }') ##USING MEAN for now
  #fi
  local flank_len=$(grep -w -i "$s_name" files/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($8 > 0) sum += $8; n++ } END { if (n > 0) print sum / n; }') ##5' UTR length must be > 0
  if [[ $flank_len == "" ]]; then
    return 2
  fi
  if [ -s $TEMP_PATH/"$s_name"_"$reg.bed" ]; then 
  	grep -i "$reg" $TEMP_PATH/"$s_name"_"$reg.bed" > $TEMP_PATH/"$s_name"_"$reg"_FLANK.bed
  	##REMOVE the line because we stored flanking regions in a seperate bed file
  	#sed -i'' "/$(echo $reg)_FLANK/d" "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed"
  	#bedtools slop -s -l $flank_len -r 0 -i "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.tmp -g "$TEMP_PATH"/"$org_name"_genomeFile.txt > "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  else
    ## if bed file is empty then we will have to take flanks from CDS
  	bedtools flank -s -l $flank_len -r 0 -i $TEMP_PATH/"$s_name"_cds.bed -g files/genes/"$org_name"/genomeFile.txt > $TEMP_PATH/"$s_name"_"$reg"_FLANK.bed
  	sed -i'' "s/cds/5utr_FLANK/g" $TEMP_PATH/"$s_name"_"$reg"_FLANK.bed
  fi
  bedtools getfasta -s -split -fi $genome_fa -bed $TEMP_PATH/"$s_name"_"$reg"_FLANK.bed -nameOnly -fullHeader > "$FASTA_PATH/$s_name.$reg.tmp" #NOTUSING name+ because BLAST doesn't like long IDs
fi

##CDS
if [[ -s $TEMP_PATH/"$s_name"_"$reg.bed" ]]; then
    bedtools getfasta -s -split -fi $genome_fa -bed $TEMP_PATH/"$s_name"_"$reg.bed" -nameOnly -fullHeader > "$FASTA_PATH/$s_name.$reg.tmp" #NOTUSING name+ because BLAST doesn't like long IDs
  #else
  #  echo $TEMP_PATH/"$s_name"_"$reg.bed is empty..."
fi

#python label_sequenceIDs.py $org_name $FASTA_PATH/$s_name.$reg.tmp $gene_full $GTF_PATH $(grep -i "filter" parameters.txt | awk -F'=' '{print $reg}') ##$org $file #$filename

#while IFS='>' read -r rec;
#do
#echo $FASTA_PATH/$s_name.$reg

#echo "sh label_sequenceIDs.sh $org_name $gene_full $GTF_PATH $FASTA_PATH/$s_name.$reg.tmp $FASTA_PATH/$s_name.$reg"
#parallel -j1 --compress --pipepart -a "$FASTA_PATH/$s_name.$reg.tmp" --recstart '>' --block -1 "sh label_sequenceIDs.sh $org_name $gene_full $GTF_PATH $FASTA_PATH/$s_name.$reg.tmp $FASTA_PATH/$s_name.$reg" #-word_size 5 -evalue 1e-25

#read -r line seq <<< $(echo "$rec")
#printf "%s\n%s" $line $seq | parallel -j1 "sh label_sequenceIDs.sh $org_name $gene_full $GTF_PATH $FASTA_PATH/$s_name.$reg.tmp $FASTA_PATH/$s_name.$reg" #-word_size 5 -evalue 1e-25
#done< "$FASTA_PATH/$s_name.$reg.tmp"
##sh label_sequenceIDs.sh $org_name $gene_full $GTF_PATH $FASTA_PATH/$s_name.$reg.tmp $FASTA_PATH/$s_name.$reg #-word_size 5 -evalue 1e-25
#echo $flank_len

#rm $FASTA_PATH/$s_name.$reg.tmp

fi
}
####

##ENTRYPOINT

#Check if scripts exist
if [[ ! -s extract_transcript_regions.py || ! -s extract_gtf_info.R || ! -s check_OrthoDB.sh || ! -s label_sequenceIDs.sh || ! -s parameters.txt ]] ; then
  echo "Missing scripts &/or parameters.txt!"
  echo "Need : extract_transcript_regions.py, extract_gtf_info.R, check_OrthoDB.sh, label_sequenceIDs.sh, parameters.txt"
  exit -1
fi

GENOMES_PATH=$(grep -i -w "genomes_path" parameters.txt | awk -F'==' '{print $2}') 
ANNOS_PATH=$(grep -i -w "annos_path" parameters.txt | awk -F'==' '{print $2}') 
BED_PATH=$(grep -i -w "bed_path" parameters.txt | awk -F'==' '{print $2}') 
FASTA_PATH=$(grep -i -w "fasta_path" parameters.txt | awk -F'==' '{print $2}') 
TEMP_PATH=$(grep -i -w "temp_path" parameters.txt | awk -F'==' '{print $2}') 
CLEAN_EXTRACT=$(grep -i -w "clean_extract" parameters.txt | awk -F'==' '{print $2}') 
GENE_SEARCH_MODE=$(grep -i -w "gene_search_mode" parameters.txt | awk -F'==' '{print $2}')
ORTHODB_PATH_PREFIX=$(grep -i -w "orthodb_path_prefix" parameters.txt | awk -F'==' '{print $2}') 
#flank_len=2000 #$(grep -i -w "flank_len" parameters.txt | awk -F'==' '{print $2}')   #2000
REMOVE_DOWNLOADS=$(grep -i -w "remove_downloads" parameters.txt | awk -F'==' '{print $2}') 
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'==' '{print $2}') 
GENES_META=$(grep -i -w "genes_meta" parameters.txt | awk -F'==' '{print $2}') 
seqID_delimiter=$(grep -i -w "seqID_delimiter" parameters.txt | awk -F'==' '{print $2}') 
LABEL_SEQS=$(grep -i -w "label_sequence_IDs" parameters.txt | awk -F'==' '{print $2}') 
#PY2_PATH=$(grep -i -w "python2_path" parameters.txt | awk -F'==' '{print $2}')
PY3_PATH=$(grep -i -w "python3_path" parameters.txt | awk -F'==' '{print $2}')
PY3_PATH="${PY3_PATH/#\~/$HOME}"
#n_threads=$(nproc --all)
n_threads=$(grep -i -w "max_concurrent_jobs" parameters.txt | awk -F'==' '{print $2}')
#transcript_regions=("cds" "exons" "noncodingexons" "3utr" "5utr")

#if [[ -z $n_threads || $n_threads > 30 ]]; then
#if [[ -z "$n_threads" || "$n_threads" -gt 31 ]]; then
#  n_threads=10
#fi

export -f get_fasta
#export -f check_gene
#export -f check_OrthoDB
export -f store_GTF_blocks
export -f store_RNA_list
export -f map_names_to_ODB_clusters
export -f remove_duplicate_genes
export -f store_GTF_block

GENOME_FILE=$1
ANNO_FILE=$2
f_org_name=$5
bed_prefix=$3

 transcript_regions=("cds" "3utr" "5utr") #("cds" "exon" "3utr" "5utr" "transcript") # "codingexons" "codingintrons" "noncodingexons" "noncodingintrons" "introns") #("cds" "exons" "3utr" "5utr" "codingexons" "codingintrons" "noncodingexons" "noncodingintrons" "introns")
#org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0); print $1" "$2; } ;')
orthodb_lookup_genes=()

#echo $org_name
echo "Command : $0 $@"

#NOT cheking if file exists because we do preliminary processing in get_genomes_ensembl.R
# if [[ ! -s $GENOME_FILE ]] ; then
#   if [[ -s $GENOME_FILE.gz ]]; then
#     GENOME_FILE=$(echo $GENOME_FILE.gz)
#   elif [[ -s $GENOME_FILE.fa ]]; then
#     GENOME_FILE=$(echo $GENOME_FILE.fa)
#   else
#     echo "$GENOME_FILE doesn't exist"
#     exit -1;
#   fi
# fi

# if [[ ! -s $ANNO_FILE ]] ; then
#   if [[ -s $ANNO_FILE.gz ]]; then
#     ANNO_FILE=$(echo $ANNO_FILE.gz)
#   elif [[ -s $ANNO_FILE.fa ]]; then
#     ANNO_FILE=$(echo $ANNO_FILE.fa)
#   else
#     echo "$ANNO_FILE doesn't exist"
#     exit -1;
#   fi
# fi

if [[ $CLEAN_EXTRACT ==  "TRUE" ]] ; then
  rm -rf $FASTA_PATH/$5
  rm -rf $GENES_META/$f_org_name
  rm -rf files/genes/$f_org_name
  rm -rf $TEMP_PATH/$f_org_name
  rm -f $GENOME_FILE.fai
  rm -f files/genes/$f_org_name/gtf_stats.csv
fi

mkdir -p $FASTA_PATH/$f_org_name
#mkdir files/genes
#mkdir files/genes/$f_org_name
mkdir -p files/genes/$f_org_name/odb/
mkdir -p $TEMP_PATH/$f_org_name
#mkdir $BED_PATH
mkdir -p $GENES_META/$f_org_name/
#rm logs/job_status.o
#rm $bed_prefix*

##Unzipping genome because samtools wont accept gzip compression(only bgzip)
#if [[ ${GENOME_FILE##*.} ==  "gz" ]] ; then
#gfile_name=${GENOME_FILE%.*}
#gunzip -c -f $GENOME_FILE > $gfile_name
##zcat $GENOME_FILE > $gfile_name
#GENOME_FILE=$gfile_name
#fi

genome_pipe="$TEMP_PATH/$f_org_name/genome_pipe"
rm -f $genome_pipe
#if [[ ! -p $genome_pipe ]]; then
mkfifo $genome_pipe
#fi

if [[ ${GENOME_FILE##*.} ==  "gz" ]] ; then
  gfile_name=${GENOME_FILE%.*}
  if [[ ! -s $gfile_name || $CLEAN_EXTRACT == "TRUE" ]] ; then
    zcat -f $GENOME_FILE > $genome_pipe &
    nohup cat $genome_pipe > $gfile_name &> /dev/null &
    genome_proc_id=$(echo $!)
  fi
  GENOME_FILE=$gfile_name
fi

##NOT unpacking GTF gz files because its faster in downstream programs, just use the ones which support gzipped files
# if [[ ${ANNO_FILE##*.} ==  "gz" ]] ; then
# afile_name=${ANNO_FILE%.*}
# gunzip -c -f $ANNO_FILE > $afile_name
# #zcat $ANNO_FILE > $afile_name
# ANNO_FILE=$afile_name
# fi

#if [[ ${ANNO_FILE##*.} !=  "gtf" ]] ; then
if [[ $(echo $ANNO_FILE | grep -q -i "gtf") ]] ; then
  file_name=${ANNO_FILE%.*}
  #gunzip $ANNO_FILE
  #gffread $ANNO_FILE -T -O -E -o $ANNOS_PATH/"$5".gtf
  nohup zcat $ANNO_FILE | gffread - -T -O -E | gzip -c - > $ANNOS_PATH/"$5".gtf.gz &> /dev/null &
  anno_proc_id=$(echo $!)
  ANNO_FILE=$ANNOS_PATH/"$5".gtf.gz
  #echo $ANNO_FILE
fi

#if [[ ! -s $GENOME_FILE.bz2 ]]; then
#  zcat $GENOME_FILE | pbzip2 -r -c  > $GENOME_FILE.bz2 #| bgzip -c --threads $n_threads --index --index-name $GENOME_FILE.bgz.gzi > $GENOME_FILE.bz2
#fi
##zcat $GENOME_FILE | pbzip2 -r -c  > $GENOME_FILE.bz2 #| bgzip -c --threads $n_threads --index --index-name $GENOME_FILE.bgz.gzi > $GENOME_FILE.bz2
#GENOME_FILE=$(echo "$GENOME_FILE.bz2")

#echo $gfile_name
#echo $afile_name
#echo $GENOME_FILE
#echo $ANNO_FILE

##SKIPPING verifying GENOME index file because it takes the same amount of time to generate it
#if [[ -s $GENOME_FILE.fai ]]; then
# if [[ $(awk 'END{print NR;}' $GENOME_FILE.fai ) != $(zgrep -c -o ">" $GENOME_FILE ) ]] ; then
#   echo "Updating Genomic Index..."
#   rm $GENOME_FILE.fai
#   #samtools faidx --fai-idx $GENOME_FILE.fai <(zcat $GENOME_FILE) 
#  samtools faidx --mark-strand sign -@ $n_threads --fai-idx $GENOME_FILE.fai $GENOME_FILE
# fi
#else
#    #samtools faidx --fai-idx $GENOME_FILE.fai <(zcat $GENOME_FILE) 
#    samtools faidx --mark-strand sign -@ $n_threads --fai-idx $GENOME_FILE.fai $GENOME_FILE
#fi

#samtools faidx --fai-idx $GENOME_FILE.fai <(pbzip2 -c -k -d $GENOME_FILE) #$GENOME_FILE #<(bgzip -c --threads $n_threads $GENOME_FILE) 

#UNCOMMENT
echo "0.1 Building Reference Genome Index (Samtools faidx)..."
if [[ ! -s files/genes/$f_org_name/genomeFile.txt ]] ; then
  if [[ ! -s $GENOME_FILE.fai ]] ; then
    # if [[ ${GENOME_FILE##*.} ==  "gz" ]] ; then
    #   time samtools faidx --fai-idx $GENOME_FILE.fai <(zcat $GENOME_FILE)
    # else
    #   time samtools faidx --fai-idx $GENOME_FILE.fai <(cat $GENOME_FILE)
    # fi
    #time samtools faidx --fai-idx $GENOME_FILE.fai < $genome_pipe
    time samtools faidx --fai-idx $GENOME_FILE.fai $genome_pipe
  fi
  awk -v OFS='\t' {'print $1,$2'} $GENOME_FILE.fai > files/genes/$f_org_name/genomeFile.txt
fi

#rm $3*
#rm $bed_prefix*
#lsof -R -p $$ -u $USER | wc -l
#printf %s"\t"%s "File handles count.." $(lsof -R -p $$ -u $USER | wc -l)
#printf %s"\t"%s "File handles count(proc).." $(ls -la /proc/$$/fd | wc -l)

while [[ $(lsof -R -p $$ | wc -l) -gt $(($(ulimit -Sn)-200)) ]]; #limit file descriptors for child & parent process to 1024-200 #$(lsof -R -p $$ -u $USER | wc -l) > $(($(ulimit -Sn)-200))
do
  printf %s"\t"%s"\n" "Waiting for file handles to close.." $(lsof -R -p $$ | wc -l)
  sleep 5
done

##DONT need to generate bed files, they are now generated in get_genomes_ensembl.R
# if [[ -s $ANNO_FILE ]] ; then
#  if ls $bed_prefix* 1> /dev/null 2>&1; then ##Check if BED files exist
#   if [[ $(ls -1 $bed_prefix* | parallel --max-procs 10 "if [[ ! -s {} ]] ; then echo 0; fi" | sort | uniq) == 0 ]] ; then
#    ##Check if BED files are NOT empty
#     rm -f $bed_prefix*
#     #$PY3_PATH extract_transcript_regions.py -i $ANNO_FILE --gtf  -o $bed_prefix
#   fi; 
#  else
#   #$PY3_PATH extract_transcript_regions.py -i $ANNO_FILE --gtf  -o $bed_prefix
#  fi 
# else
#   exit 2
# fi


#readarray gene_list < $4
# if [[ -z "$(ls -A $GENES_META/$f_org_name/)" || $CLEAN_EXTRACT == "TRUE" ]]; then
#    gene_list=($(cat $4 | sort | uniq | grep -v -w -i "gene"))
# else
#    temp_list=($(cat $4 | sort | uniq | grep -v -w -i "gene"))
#    existing_list=($(find $GENES_META/$f_org_name/*.gtf_slice -exec basename {} \; | cut -d '.' -f 1 | sort | uniq))
#    gene_list=($(echo ${temp_list[@]/${existing_list[@]}}))
# fi

gene_list=($(cat $4 | sort | uniq | grep -v -w -i "gene"))

#s_names=$(parallel --max-procs $n_threads "echo {} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' " ::: ${gene_list[@]})
s_names=($(awk -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo ${gene_list[@]})))

echo "1. Checking Gene Names & Splitting GTFs for parallel processing..."

rm -f files/genes/$f_org_name/1.list
rm -f files/genes/$f_org_name/2.list
rm -f files/genes/$f_org_name/full.list
rm -f files/genes/$f_org_name/final.list
rm -f $TEMP_PATH/$f_org_name/*

#Updated implementation of check_gene, the function is now deprecated
#time parallel --max-procs $n_threads "check_gene {} $f_org_name $ANNO_FILE $GENES_META " ::: ${gene_list[@]} #parallel -N 251 --pipe  #-j ${#gene_list[@]} #$(($n_threads * 2)) #get_fasta is called from check_gene function after preliminary checks 
#./jobhold.sh check_gene

#for i in "${gene_list[@]}"; do 
#  eexp_gene="$eexp_gene -e $i"; #creating expression string -e <gene> to pass to grep
#done
#echo "${eexp_gene[@]}"

#eexp_gene=$(printf -- '-e %s ' "${gene_list[@]}")
eexp_gene=$(printf -- '%s\n' "${gene_list[@]}")

#echo ${#gene_list[@]}/$n_threads | bc

#time echo "${eexp_gene[@]}" | zgrep -i -C 0 --group-separator='>' -f - $ANNO_FILE | parallel --max-procs $n_threads --pipe --line-buffer --tmpdir $TEMP_PATH/$f_org_name/ --recstart '\n' --recend '>' -N $(echo ${#gene_list[@]}/$n_threads | bc) --compress --compress-program gzip --cat "store_GTF_blocks {} files/genes/$f_org_name/1.list $GENES_META $f_org_name $n_threads $TEMP_PATH/$f_org_name/" ## -N $(echo ${#gene_list[@]}/$n_threads | bc) 
#time zgrep -i -C 0 --group-separator='>' -e sync -e sybu -e dazl -e anln -e vg1 -e acot7 $ANNO_FILE | parallel --max-procs $n_threads --pipe --tmpdir $TEMP_PATH/$f_org_name/ --recstart '\n' --recend '>' -N $(echo ${#gene_list[@]}/$n_threads | bc) --compress --compress-program gzip --cat "store_GTF_blocks {} files/genes/$f_org_name/1.list $GENES_META $f_org_name $n_threads $TEMP_PATH/$f_org_name/" ## -N 1

#PLAIN - PASSED < 10mins
# Code : time echo "${eexp_gene[@]}" | zgrep -i -C 0 --group-separator='>' -f - $ANNO_FILE | parallel --max-procs $n_threads --pipe --tmpdir $TEMP_PATH/$f_org_name/ --recstart '\n' --recend '>' -N 1 --compress --compress-program gzip --cat "store_GTF_blocks {} files/genes/$f_org_name/1.list $GENES_META $f_org_name $n_threads $TEMP_PATH/$f_org_name/" ## -N 1 # -N $(echo ${#gene_list[@]}/$n_threads | bc)
# no N - 8m25s
# -N $(echo ${#gene_list[@]}/$n_threads | bc) - 7m51s
# -N 1 -9m16s
#RECURSION _FAILED > 30mins
# Code : time echo "${eexp_gene[@]}" | zgrep -i -C 0 --group-separator='>' -f - $ANNO_FILE | parallel --max-procs $n_threads --pipe --tmpdir $TEMP_PATH/$f_org_name/ --recstart '\n' --recend '>' --compress --compress-program gzip -N $(echo ${#gene_list[@]}/$n_threads | bc) --cat "store_GTF_blocks {} files/genes/$f_org_name/1.list $GENES_META $f_org_name $n_threads $TEMP_PATH/$f_org_name/" ## -N 1 #--pipe -N $(echo ${#gene_list[@]}/$n_threads | bc)
# no N - & block =$(echo $group_sep_count/$5 | bc) --> FAILED
# -N $(echo ${#gene_list[@]}/$n_threads | bc) & block =$(echo $group_sep_count/$5 | bc)
# -N $(echo ${#gene_list[@]}/$n_threads | bc) & block =-1 -->

MODE=""
if [[ $GENE_SEARCH_MODE=="HARD" ]]; then
  MODE="-w"
fi
#zgrep -i -A 0 --group-separator='>' -e "adrm1" -e "sync" ../../mrna_loc/files/annos/danio_rerio.gtf.gz | csplit --suppress-matched - '/>/' '{*}'
wait $anno_proc_id
time echo "${eexp_gene[@]}" | zgrep -i $MODE -A 0 --group-separator='>' -f - $ANNO_FILE | parallel --max-procs $n_threads --pipe --tmpdir $TEMP_PATH/$f_org_name/ --recstart '' --recend '>' --compress --compress-program gzip --cat "store_GTF_blocks files/genes/$f_org_name/1.list $GENES_META/$f_org_name/ $n_threads $TEMP_PATH/$f_org_name/ -ct {}" ## -N 1 # -N $(echo ${#gene_list[@]}/$n_threads | bc) --recstart '\n' --cat {} -N $(echo ${#gene_list[@]}/$n_threads | bc)

#time parallel --max-procs $n_threads --tmpdir $TEMP_PATH/$f_org_name/ --recstart "-e" --recend "\n" -X "echo {} -----------" ::: "${eexp_gene[@]}"

#time parallel --max-procs $n_threads --tmpdir $TEMP_PATH/$f_org_name/ --compress --compress-program gzip --recstart "-e" --recend "\n" -X "grep -i -C 0 --group-separator='>' {} <(zcat $ANNO_FILE)" ::: "${eexp_gene[@]}"
#| parallel --max-procs $n_threads --tmpdir $TEMP_PATH/$f_org_name/ --recstart '\n' --recend '>' --compress --compress-program gzip --cat "store_GTF_blocks {} files/genes/$f_org_name/1.list $GENES_META $f_org_name $n_threads $TEMP_PATH/$f_org_name/" ##-N 1

#Remove duplicate genes
#remove_duplicate_genes files/genes/$f_org_name/1.list
echo ${gene_list[@]/($(cat files/genes/$f_org_name/1.list)))} > files/genes/$f_org_name/2.list

##Looking up genes which were not present in GTF annotations (because of differing gene names between organisms)
#awk '{gsub(" ","\n");print;}' files/genes/$f_org_name/*.list files/genes/$f_org_name/2.list | sort | uniq > files/genes/$f_org_name/full.list
cat files/genes/$f_org_name/1.list files/genes/$f_org_name/2.list > files/genes/$f_org_name/full.list

#######################################################################################################

echo "2. Checking OrthoDB for missing & orthologous genes..." #(log:$TEMP_PATH/$f_org_name/check_ortho.[o/e])

##if [[ -s files/genes/$f_org_name/2.list ]]; then
#  time check_OrthoDB files/genes/$f_org_name/full.list $f_org_name $ANNO_FILE $GENES_META $4 ;
##fi

if [[ ! -s files/genes/$f_org_name/odb.list || ! -s files/genes/$f_org_name/final.list ]] ; then
  #time ./jobhold.sh check_ortho ./check_OrthoDB.sh $f_org_name files/genes/$f_org_name/full.list $ANNO_FILE ;
  time ./check_OrthoDB.sh $f_org_name files/genes/$f_org_name/full.list $ANNO_FILE #1> $TEMP_PATH/$f_org_name/check_ortho.o 2> $TEMP_PATH/$f_org_name/check_ortho.e 
fi
nohup bash -c "map_names_to_ODB_clusters $f_org_name $GENES_META" & #&> $TEMP_PATH/$f_org_name/check_ortho.o &
odb_proc_id=$(echo $!)

##REFRESH geene list based on file names
cat files/genes/$f_org_name/1.list files/genes/$f_org_name/2.list files/genes/$f_org_name/odb.list > files/genes/$f_org_name/full.list
#gene_list=($(find $GENES_META/$f_org_name/*.gtf_slice | parallel "basename {}" | cut -d '.' -f 1 | sort | uniq ))
gene_list=($(cat files/genes/$f_org_name/full.list | sort | uniq | grep -v -w -i "gene")) #| awk  -v s_var='_' '{gsub(/[[:punct:]]/,s_var);}1'))

#s_names=$(parallel --keep-order --max-procs $n_threads "echo {} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' " ::: ${gene_list[@]})
s_names=($(awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo ${gene_list[@]})))

## GET GTF stats - UTR lengths, CDS_count, exon_count (both fully and partially annotated info (eg, UTR info is hidden in exon and CDS info and is extracted using exon boundaries vs CDS boundaries)) 

#echo ${s_names[@]} ${gene_list[@]}
exit
#######################################################################################################

echo "3. Extracting Transcript Stats from GTF_Slices...(log:$TEMP_PATH/$f_org_name/get_GTF_info.[o/e])"

#echo "With parallel"
#time parallel --link --tmpdir $TEMP_PATH --compress --max-procs $(echo $n_threads/$(nproc) + 1 | bc) "if [[ -s $GENES_META/$f_org_name/{1}.rna_list && ! -z {1} && ! -z {2} && -s $GENES_META/$f_org_name/{1}.gtf_slice && ! -s $GENES_META/$f_org_name/{1}.gtf_info ]] ; then
#    Rscript extract_gtf_info_parallel.R $GENES_META/$f_org_name/{1}.gtf_slice $f_org_name {2} $GENES_META/$f_org_name/{1}.gtf_info ;
#  fi" ::: ${s_names[@]} ::: ${gene_list[@]} 1>> $TEMP_PATH/get_GTF_info.o 2>> $TEMP_PATH/get_GTF_info.e #$ANNO_FILE #./jobhold.sh get_GTF_info 
##rm $GENES_META/$f_org_name/*.gtf_info
if [[ ! -s files/genes/$f_org_name/gtf_stats.csv || ! -s files/genes/$f_org_name/final.list ]] ; then
  #echo "With R"
  #time ./jobhold.sh get_GTF_info Rscript --vanilla --verbose extract_gtf_info.R $GENES_META/$f_org_name/ $f_org_name files/genes/$f_org_name/gtf_stats.csv ; #1> $TEMP_PATH/get_GTF_info.o 2> $TEMP_PATH/get_GTF_info.e
  time Rscript --vanilla --verbose extract_gtf_info.R $GENES_META/$f_org_name/ $f_org_name files/genes/$f_org_name/gtf_stats.csv 1> $TEMP_PATH/$f_org_name/get_GTF_info.o 2> $TEMP_PATH/$f_org_name/get_GTF_info.e
fi
#find $GENES_META/$f_org_name/*.gtf_info -exec sed 1d {} \; > files/genes/$f_org_name/gtf_stats.csv
# parallel --max-procs $n_threads "if [[ -s $GENES_META/$f_org_name/{}.gtf_info ]] ; then
#   sed 1d $GENES_META/$f_org_name/{}.gtf_info ;
# fi" ::: ${gene_list[@]} > files/genes/$f_org_name/gtf_stats.csv
if [[ ! -s files/genes/$f_org_name/gtf_stats.csv || ! -s files/genes/$f_org_name/final.list ]] ; then
  >2& echo "ERROR: Extraction of transcript stats failed..."
  >2& echo "Check $TEMP_PATH/$f_org_name/get_GTF_info.[o/e]"
  >2& echo "Remove $TEMP_PATH/$f_org_name/ & $GENES_META/$f_org_name/ and re-run the pipeline"
fi

sed 1d files/genes/$f_org_name/gtf_stats.csv | awk -F',' '{print $1"\n"}' | sort | uniq > files/genes/$f_org_name/final.list

gene_list=($(cat files/genes/$f_org_name/final.list))
#s_names=$(parallel --keep-order --max-procs $n_threads "echo {} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' " ::: ${gene_list[@]})
s_names=($(awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' <(echo ${gene_list[@]})))
deletion_list=($(find $GENES_META/$f_org_name/*.gtf_slice -exec basename {} \; | cut -d '.' -f 1)/${s_names[@]}})

echo ${deletion_list[@]}
exit

#DELETE GENEs which are not in the final list
#echo "Number of Genes in FINAL LIST:"
parallel --max-procs $n_threads --tmpdir $TEMP_PATH/$f_org_name/ --compress --compress-program gzip -X "(cd $GENES_META/$f_org_name/ && rm -f {}.gtf_slice {}.rna_list )" ::: "${deletion_list[@]}"

#######################################################################################################

echo "4. Fetching sequences from Genome..." #(log:$TEMP_PATH/$f_org_name/get_FASTA.[o/e])
echo "Transcript Regions : ${transcript_regions[@]}"

wait $genome_proc_id

time parallel --max-procs $n_threads "if [[ ! -z {1} && ! -z {2} && -s $GENES_META/$f_org_name/{1}.rna_list ]] ; then
 get_fasta {1} {2} $bed_prefix $GENOME_FILE $f_org_name $FASTA_PATH/$f_org_name $ANNO_FILE $TEMP_PATH/$f_org_name $GENES_META ; 
fi" ::: ${gene_list[@]} ::: ${transcript_regions[@]} #1> $TEMP_PATH/$f_org_name/get_FASTA.o 2> $TEMP_PATH/$f_org_name/get_FASTA.e

if [[ $LABEL_SEQS ==  "TRUE" ]] ; then
  echo "4.1 Labelling sequences..."
  wait $odb_proc_id

  #Get CLUSTERS for FINAL LIST OF SELECTED GENES
  #printf -- '%s\n' "${gene_list[@]}" | grep -i -w -f - files/genes/$f_org_name/odb/odb.full_map | sort -k 1 | uniq > files/genes/$f_org_name/odb/odb.final_map 

  time parallel --link --max-procs $n_threads "printf -- %s\t%s {1} {2}" ::: ${gene_list[@]} ::: ${s_names[@]} | parallel --colsep "\t" --max-procs $n_threads "printf -- %s\t%s\t%s {1} {2} {3}" ::: ${transcript_regions[@]} #| "./label_sequenceIDs.sh $f_org_name {1} $ANNO_FILE $FASTA_PATH/$f_org_name/{2}.$t_r $FASTA_PATH/$f_org_name/{2}.$t_r.tmp files/genes/$f_org_name/odb/odb.full_map" 
else
  time parallel --link --max-procs $n_threads "mv -f $FASTA_PATH/$f_org_name/{1}.{2}.tmp $FASTA_PATH/$f_org_name/{1}.{2}" ::: ${gene_list[@]} ::: ${s_names[@]}
fi

#######################################################################################################

echo "5. Generating Metadata and Cleaning up..."
#find files/genes/$5 -iname "*.rna_list" -empty | awk -F'/' '{print $NF}' | sed 's/.rna_list//g' > $FASTA_PATH/$5/MISSING_GENES
#grep -v -i -f $FASTA_PATH/$5/MISSING_GENES $4 | sort > $FASTA_PATH/$5/AVAILABLE_GENES
sed 1d files/genes/$f_org_name/gtf_stats.csv | awk -F',' '{print $2}' | sort | uniq >  files/genes/$f_org_name/AVAILABLE_GENES
grep -v -i -f files/genes/$f_org_name/AVAILABLE_GENES $4 | sort | uniq > files/genes/$f_org_name/MISSING_GENES

find $FASTA_PATH/$f_org_name/ -type f -name "*.fai" -exec rm -f {} +

#rm -f $TEMP_PATH/$f_org_name/*

wait $odb_proc_id

# ##Populate clusters(files/genes/<org>/ALL_CLUSTERS) and group FASTA sequences according to clusters and store it in files/groups/
# #grep -r ">" $FASTA_PATH/$f_org_name/ | parallel --max-procs $n_threads --colsep ":>" --recend "\n" echo {2} | awk -F"$seqID_delimiter" '{print $5}' |  sort | uniq > files/genes/$f_org_name/ALL_CLUSTERS
# grep -r ">" $FASTA_PATH/$f_org_name/ | awk -F ":>" '{print $2}' | awk -F"$seqID_delimiter" '{print $5}' |  sort | uniq > files/genes/$f_org_name/ALL_CLUSTERS


# if [[ $REMOVE_DOWNLOADS ==  "FALSE" ]] ; then
  
#   if [[ ! -s $GENOMES_PATH/$5.fa.gz ]]; then
#     mv $GENOME_FILE $GENOMES_PATH/$5.fa
#     gzip -f -c --fast $GENOMES_PATH/$5.fa # --rsyncable
#   fi
#   #bgzip -i $GENOMES_PATH/$5.fa
#   #bgzip -i $ANNOS_PATH/$5.gtf
#   if [[ ! -s $ANNOS_PATH/$5.gtf.gz ]]; then
#     mv $ANNO_FILE $ANNOS_PATH/$5.gtf
#     gzip -f -c --fast $ANNOS_PATH/$5.gtf # --rsyncable
#   fi
# else
#   rm $GENOME_FILE
#   rm $ANNO_FILE
# fi

# #rm $TEMP_PATH/$f_org_name.*.*
# #rm files/genes/$5/*.gtf_slice

# #rm $3*
# #rm files/bed/$5*

# #exit 1
