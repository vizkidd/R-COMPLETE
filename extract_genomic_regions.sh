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
check_OrthoDB() {
  #1 - genelist to loopkup in OrthoDB
  #2 - Organism name
  #3 - FULL GTF File(not slice)
  #4 - Metadata output path
  #5 - Full gene list (union of lists 1, 2.1, 2.2)
 # proc_list=()
 org_name=$(echo $2 | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
 ANNO_FILE=$3
 GENES_META=$4
ORTHODB_PATH=$(grep -i -w "orthodb_files_path" parameters.txt | awk -F'=' '{print $2}') 
ORTHODB_PREFIX=$(grep -i -w "orthodb_prefix" parameters.txt | awk -F'=' '{print $2}') 
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'=' '{print $2}') 

#org_name=$(echo $3 | awk '{ gsub(/[[:punct:]]/, " ", $0); print $1" "$2; } ;')
echo "Checking OrthoDB..."
#if [[ -s files/genes/$f_org_name/$file_out.rna_list ]] ; then
#  rm files/genes/$f_org_name/$file_out.rna_list
#fi
readarray refs < $REF_ORGS
#readarray gene_list < $1
org_ID=$(grep -i -P "\t\b($org_name)\b\t" "$ORTHODB_PATH/$ORTHODB_PREFIX"_species.tab | awk '{print $2}')
#gene_ID=$(zcat $ORTHODB_PATH/odb10v1_genes.tab.gz | grep $org_ID | grep -i -P "\b($1)\b\t" | awk '{print $1}')
if [[ -z $org_ID ]] ; then
  org_ID=$(grep -i -P "($org_name)" "$ORTHODB_PATH/$ORTHODB_PREFIX"_species.tab | awk '{print $2}')
fi
if [[ ! -z $org_ID ]] ; then
  #echo $1 $org_ID $org_name $s_name $ANNO_FILE  
  #parallel -j1 --compress --pipepart -a "$REF_ORGS" --recstart '\n' --block -1 "./check_OrthoDB.sh {} $1 $org_ID $f_org_name $s_name"
  time parallel --max-procs 1 "./jobhold.sh check_ortho ./check_OrthoDB.sh {} $1 $org_ID $2 $ANNO_FILE" ::: ${refs[@]}  #parallel -N 251 --pipe #-j ${#refs[@]}
mkdir files/genes/$f_org_name/refs
cat files/genes/$f_org_name/refs/*.genes | sort | uniq > files/genes/$f_org_name/odb.clusters_map
cat files/genes/$f_org_name/refs/*.names | sort | uniq > files/genes/$f_org_name/odb.genes_map

#awk '{print $1"\t"$4}' "$ORTHODB_PATH/$ORTHODB_PREFIX"_genes.tab | grep "$org_ID" | grep -w -f $5 | sort -k1 | uniq >> files/genes/$f_org_name/odb.genes_map

 ##MAP GENE NAMES to clusters
 awk 'fname != FILENAME { fname = FILENAME; idx++ } FNR > 1 && idx == 1 { gene_name[$1] = $2 } FNR > 1 && idx == 2 { clusters[$2]=$1 }END{for(key in gene_name) print key"\t"clusters[key]"\t"gene_name[key];}' files/genes/$f_org_name/odb.genes_map files/genes/$f_org_name/odb.clusters_map > files/genes/$f_org_name/odb.group_map
##COLLAPSE GENE name rows (in case a field has multiple records)
parallel --max-procs $n_threads "printf '\(%s\)\{1\}\n' {}" ::: ${gene_list[@]} | grep -i -f - -o -n -G files/genes/$f_org_name/odb.group_map | awk -F':' '{if($1 in a)a[$1]=a[$1]","$2;else a[$1]=$2 ;}END{for(key in a)print key"\t"a[key];}'  | sort -k1 -n > files/genes/$f_org_name/genes.match_map

##FINAL map, GENE list -> CLUSTER list -> SEARCH_GROUP list (second awk is to colllapse rows based on clusters)
awk 'START{j=0;k=0;} fname != FILENAME { fname = FILENAME; idx++; }  FNR > 0 && idx == 1 { j=FNR; genes[j]=$3; search_group[j]=$3; clusters[j]=$2;  }  FNR > 0 && idx == 2 { k=FNR; if(match(k,$1)==0) search_group[$1]=$2; }END{if(j<k) l=k;else l=j; for(i=1; i<l;i++) print clusters[i]"\t"genes[i]"\t"search_group[i];} ' files/genes/$f_org_name/odb.group_map files/genes/$f_org_name/genes.match_map | awk '{if($1 in a)a[$1]=a[$1]","$2;else a[$1]=$2 ;}END{for(key in a)print key"\t"a[key];}' > files/genes/$f_org_name/odb.final_map
 

  find $GENES_META/$f_org_name/*.*.gtf_slice | parallel --max-procs $n_threads " basename {}" | cut -d '.' -f 1 | sort | uniq | parallel --max-procs $n_threads " find $GENES_META/$f_org_name/{}*.gtf_slice | xargs cat | sort -k2,3 -n -r > $GENES_META/$f_org_name/{}.gtf_slice" 
  rm $GENES_META/$f_org_name/*.*.gtf_slice; 
  find $GENES_META/$f_org_name/*.*.rna_list | parallel --max-procs $n_threads " basename {}" | cut -d '.' -f 1 | sort | uniq | parallel --max-procs $n_threads " find $GENES_META/$f_org_name/{}*.rna_list | xargs cat | sort | uniq > $GENES_META/$f_org_name/{}.rna_list" 
  rm $GENES_META/$f_org_name/*.*.rna_list; 
  
  #cat $TEMP_PATH/*.$f_org_name.$s_name > files/genes/$f_org_name/$file_out.rna_list
  #rm $TEMP_PATH/*.$f_org_name.$s_name
else
  echo "$org_name not found in OrthoDB files"
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

  ANNO_FILE=$3
  GENES_META=$4


  if [[ $gene_full !=  "" ]] ; then
        gene=$(echo $gene_full | sed 's/.$//')
    file_out=$(echo $gene_full | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')

  #   if [[ -s files/genes/$f_org_name/1.list && -s files/genes/$f_org_name/2.1.list && -s files/genes/$f_org_name/2.2.list && -s files/genes/$f_org_name/genes.odb ]] ; then
  #   if ! grep -q -w "$gene" files/genes/$f_org_name/1.list files/genes/$f_org_name/2.1.list files/genes/$f_org_name/2.2.list files/genes/$f_org_name/genes.odb ; then
  #     return 2
  #   fi
  # fi

    #echo $gene_full #$gene #$f_org_name $file_out $bed_prefix
  #Saving gene info list in a file
    grep -i "$gene" $ANNO_FILE > $GENES_META/$f_org_name/$file_out.gtf_slice
    if [ -s $GENES_META/$f_org_name/$file_out.gtf_slice ]; 
    then
      #grep -i "$gene" $2 | awk -F'\t' '{print $9}' | awk -F';' '{print $3}' | awk '{print $2}' | sed 's/"//g' | sort | uniq > files/genes/$5/$gene.rna_list
      if grep -q -w -i "$gene_full" $GENES_META/$f_org_name/$file_out.gtf_slice ; then
        echo "1 - Gene name has perfect match ($gene_full)"
        echo $gene_full >> files/genes/$f_org_name/1.list
        grep -i -w "$gene_full" $GENES_META/$f_org_name/$file_out.gtf_slice | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > $GENES_META/$f_org_name/$file_out.rna_list
      elif grep -q -i "$gene_full" $GENES_META/$f_org_name/$file_out.gtf_slice ; then
        #grep -i "$gene_full" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
          echo "2.1 - Gene name(full) matches a supergroup or subgroup, will check OrthoDB ($gene_full)"
          echo $gene_full >> files/genes/$f_org_name/2.1.list
          echo $gene_full >> files/genes/$f_org_name/genes.odb
      elif grep -q -i "$gene" $GENES_META/$f_org_name/$file_out.gtf_slice ; then
        #grep -i "$gene_full" $ANNO_FILE | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$5/$file_out.rna_list
          echo "2.2 - Gene name(stripped) matches a supergroup or subgroup, will check OrthoDB ($gene)"
          echo $gene_full $gene >> files/genes/$f_org_name/2.2.list
          echo $gene_full >> files/genes/$f_org_name/genes.odb
      fi
    fi
  fi
}


#      else
#        echo "3 - Gene name has a partial/no match, saving partial results"
#        grep -i "$gene" files/genes/$f_org_name/$file_out.gtf_slice | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$f_org_name/$file_out.rna_list
#      fi
#    
#    if [ ! -s files/genes/$f_org_name/$file_out.rna_list ]; then
#      grep -i "$gene_full" files/genes/$f_org_name/$file_out.gtf_slice | perl -lne 'print "@m" if @m=(/((?:transcript_id)\s+\S+)/g);' | awk '{print $2}' | sed 's/[";]//g' | sort | uniq > files/genes/$f_org_name/$file_out.rna_list
#    fi
#    
#    cat files/genes/$f_org_name/$file_out.rna_list
#    
#    if [ -s files/genes/$f_org_name/$file_out.rna_list ] ; then
#      echo $gene_full
#      ## GET GTF stats - UTR lengths, CDS_count, exon_count (both fully and partially annotated info (eg, UTR info is hidden in exon and CDS info and is extracted using exon boundaries vs CDS boundaries))
#      Rscript extract_gtf_info.R files/genes/$f_org_name/$file_out.gtf_slice $f_org_name $file_out
#      #time parallel -j ${#transcript_regions[@]} "get_fasta $gene_full {} $bed_prefix $GENOME_FILE $f_org_name $FASTA_PATH $file_out 'files/genes/$f_org_name/$file_out.gtf_slice' $TEMP_PATH " ::: ${transcript_regions[@]} #$8
#      parallel "get_fasta $gene_full {} $bed_prefix $GENOME_FILE $f_org_name $FASTA_PATH $file_out $ANNO_FILE $TEMP_PATH " ::: ${transcript_regions[@]} #parallel -N 251 --pipe parallel  #time #-j ${#transcript_regions[@]} 
#fi

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

gene_full=$1
s_name=$(echo $gene_full | awk '{ gsub(/[[:punct:]]/, "_", $0) } 1;')
reg=$2
base_bed=$3
genome_fa=$4
org_name=$5
FASTA_PATH=$6
GTF_PATH=$7
TEMP_PATH=$8
GENES_META=$9

#echo $(grep -i -w -f files/genes/$org_name/$s_name.rna_list $base_bed_$reg.bed)
#echo $gene_full $reg $base_bed $genome_fa $org_name $FASTA_PATH $s_name $GTF_PATH $TEMP_PATH 
#flank_len=${10}
#echo $flank_len

if [ ! -s $FASTA_PATH/$org_name/$s_name.$reg ]; then #$FASTA_PATH/$org_name/$file_out.cds

# if grep -q -w -i "$s_name" files/genes/$org_name/gtf_stats.csv ; then
#   return 2
# fi

grep -i -w -f $GENES_META/$org_name/$s_name.rna_list "$base_bed"_"$reg.bed" > "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed"

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
  flank_len=$(grep -w -i "$s_name" files/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($9 > 3) sum += $9; n++ } END { if (n > 0) print sum / n; }') ##3' UTR length must be > 3
  if [[ $flank_len == "" ]]; then
    return 2
  fi

  if [ -s "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed" ]; then 
  	grep -i "$reg" "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed" > "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  	##REMOVE the line because we stored flanking regions in a seperate bed file
  	#sed -i'' "/$(echo $reg)_FLANK/d" "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed"
  	#bedtools slop -s -r $flank_len -l 0 -i "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.tmp -g "$TEMP_PATH"/"$org_name"_genomeFile.txt > "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  else
    ## if bed file is empty then we will have to take flanks from CDS
  	bedtools flank -s -r $flank_len -l 0 -i "$TEMP_PATH"/"$org_name"_"$s_name"_cds.bed -g "$TEMP_PATH"/"$org_name"_genomeFile.txt > "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  	sed -i'' "s/cds/3utr_FLANK/g" "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  fi
  bedtools getfasta -s -split -fi $genome_fa -bed "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed -nameOnly -fullHeader >> "$FASTA_PATH/$org_name/$s_name.$reg.tmp" #NOTUSING name+ because BLAST doesn't like long IDs (and we dont want chromosome, genomic position info)
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
  flank_len=$(grep -w -i "$s_name" files/genes/$org_name/gtf_stats.csv | awk -F',' '{ if ($8 > 0) sum += $8; n++ } END { if (n > 0) print sum / n; }') ##5' UTR length must be > 0
  if [[ $flank_len == "" ]]; then
    return 2
  fi
  if [ -s "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed" ]; then 
  	grep -i "$reg" "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed" > "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  	##REMOVE the line because we stored flanking regions in a seperate bed file
  	#sed -i'' "/$(echo $reg)_FLANK/d" "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed"
  	#bedtools slop -s -l $flank_len -r 0 -i "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.tmp -g "$TEMP_PATH"/"$org_name"_genomeFile.txt > "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  else
    ## if bed file is empty then we will have to take flanks from CDS
  	bedtools flank -s -l $flank_len -r 0 -i "$TEMP_PATH"/"$org_name"_"$s_name"_cds.bed -g "$TEMP_PATH"/"$org_name"_genomeFile.txt > "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  	sed -i'' "s/cds/5utr_FLANK/g" "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed
  fi
  bedtools getfasta -s -split -fi $genome_fa -bed "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg"_FLANK.bed -nameOnly -fullHeader >> "$FASTA_PATH/$org_name/$s_name.$reg.tmp" #NOTUSING name+ because BLAST doesn't like long IDs
fi

if [ -s "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed" ]; then
    bedtools getfasta -s -split -fi $genome_fa -bed "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed" -nameOnly -fullHeader >> "$FASTA_PATH/$org_name/$s_name.$reg.tmp" #NOTUSING name+ because BLAST doesn't like long IDs
  else
    echo "$TEMP_PATH"/"$org_name"_"$s_name"_"$reg.bed is empty..."
fi

#python label_sequenceIDs.py $org_name $FASTA_PATH/$org_name/$s_name.$reg.tmp $gene_full $GTF_PATH $(grep -i "filter" parameters.txt | awk -F'=' '{print $reg}') ##$org $file #$filename

#while IFS='>' read -r rec;
#do
#echo $FASTA_PATH/$org_name/$s_name.$reg

#echo "sh label_sequenceIDs.sh $org_name $gene_full $GTF_PATH $FASTA_PATH/$org_name/$s_name.$reg.tmp $FASTA_PATH/$org_name/$s_name.$reg"
#parallel -j1 --compress --pipepart -a "$FASTA_PATH/$org_name/$s_name.$reg.tmp" --recstart '>' --block -1 "sh label_sequenceIDs.sh $org_name $gene_full $GTF_PATH $FASTA_PATH/$org_name/$s_name.$reg.tmp $FASTA_PATH/$org_name/$s_name.$reg" #-word_size 5 -evalue 1e-25

#read -r line seq <<< $(echo "$rec")
#printf "%s\n%s" $line $seq | parallel -j1 "sh label_sequenceIDs.sh $org_name $gene_full $GTF_PATH $FASTA_PATH/$org_name/$s_name.$reg.tmp $FASTA_PATH/$org_name/$s_name.$reg" #-word_size 5 -evalue 1e-25
#done< "$FASTA_PATH/$org_name/$s_name.$reg.tmp"
##sh label_sequenceIDs.sh $org_name $gene_full $GTF_PATH $FASTA_PATH/$org_name/$s_name.$reg.tmp $FASTA_PATH/$org_name/$s_name.$reg #-word_size 5 -evalue 1e-25
#echo $flank_len

#rm $FASTA_PATH/$org_name/$s_name.$reg.tmp

fi
}
####

##ENTRYPOINT

#Check if scripts exist
if [[ ! -s extract_transcript_regions.py || ! -s extract_gtf_info.R || ! -s check_OrthoDB.sh || ! -s label_sequenceIDs.sh ]] ; then
  echo "Missing scripts!"
  echo "Need : extract_transcript_regions.py, extract_gtf_info.R, check_OrthoDB.sh, label_sequenceIDs.sh"
  exit -1
fi

GENOMES_PATH=$(grep -i -w "genomes_path" parameters.txt | awk -F'=' '{print $2}') 
ANNOS_PATH=$(grep -i -w "annos_path" parameters.txt | awk -F'=' '{print $2}') 
FASTA_PATH=$(grep -i -w "fasta_path" parameters.txt | awk -F'=' '{print $2}') 
TEMP_PATH=$(grep -i -w "temp_path" parameters.txt | awk -F'=' '{print $2}') 
CLEAN_EXTRACT=$(grep -i -w "clean_extract" parameters.txt | awk -F'=' '{print $2}') 
ORTHODB_PATH=$(grep -i -w "orthodb_files_path" parameters.txt | awk -F'=' '{print $2}') 
ORTHODB_PREFIX=$(grep -i -w "orthodb_prefix" parameters.txt | awk -F'=' '{print $2}') 
#flank_len=2000 #$(grep -i -w "flank_len" parameters.txt | awk -F'=' '{print $2}')   #2000
REMOVE_DOWNLOADS=$(grep -i -w "remove_downloads" parameters.txt | awk -F'=' '{print $2}') 
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'=' '{print $2}') 
GENES_META=$(grep -i -w "genes_meta" parameters.txt | awk -F'=' '{print $2}') 
seqID_delimiter=$(grep -i -w "seqID_delimiter" parameters.txt | awk -F'=' '{print $2}') 
LABEL_SEQS=$(grep -i -w "label_sequence_IDs" parameters.txt | awk -F'=' '{print $2}') 
#PY2_PATH=$(grep -i -w "python2_path" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH=$(grep -i -w "python3_path" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH="${PY3_PATH/#\~/$HOME}"
n_threads=$(nproc --all)
#transcript_regions=("cds" "exons" "noncodingexons" "3utr" "5utr")

export -f get_fasta
export -f check_gene
export -f check_OrthoDB
#export -f parallel_checkODB

GENOME_FILE=$1
ANNO_FILE=$2
f_org_name=$5
bed_prefix=$3

 transcript_regions=("cds" "exons" "3utr" "5utr") # "codingexons" "codingintrons" "noncodingexons" "noncodingintrons" "introns") #("cds" "exons" "3utr" "5utr" "codingexons" "codingintrons" "noncodingexons" "noncodingintrons" "introns")
#org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0) } 1;')
org_name=$(echo $f_org_name | awk '{ gsub(/[[:punct:]]/, " ", $0); print $1" "$2; } ;')
orthodb_lookup_genes=()

#echo $org_name
echo "Command : $0 $@"

if [[ ! -s $GENOME_FILE ]] ; then
  if [[ -s $GENOME_FILE.gz ]]; then
    GENOME_FILE=$(echo $GENOME_FILE.gz)
  elif [[ -s $GENOME_FILE.fa ]]; then
    GENOME_FILE=$(echo $GENOME_FILE.fa)
  else
    echo "$GENOME_FILE doesn't exist"
    exit -1;
  fi
fi

if [[ ! -s $ANNO_FILE ]] ; then
  if [[ -s $ANNO_FILE.gz ]]; then
    ANNO_FILE=$(echo $ANNO_FILE.gz)
  elif [[ -s $ANNO_FILE.fa ]]; then
    ANNO_FILE=$(echo $ANNO_FILE.fa)
  else
    echo "$ANNO_FILE doesn't exist"
    exit -1;
  fi
fi

if [[ $CLEAN_EXTRACT ==  "TRUE" ]] ; then
  rm $GENOME_FILE.fai
  rm -rf $FASTA_PATH/$5
  rm -rf $GENES_META/$f_org_name
  rm -rf files/genes/$f_org_name
fi

mkdir $FASTA_PATH/$f_org_name
mkdir files/genes
mkdir files/genes/$f_org_name
mkdir files/genes/$f_org_name/refs/
#mkdir files/bed
mkdir $GENES_META/$f_org_name/
#rm logs/job_status.o
#rm $bed_prefix*

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

if [[ -s $GENOME_FILE.fai ]]; then
  if [[ $(awk 'END{print NR;}' $GENOME_FILE.fai ) != $(grep -c -o ">" $GENOME_FILE ) ]] ; then
    echo "Updating Genomic Index..."
    rm $GENOME_FILE.fai
    samtools faidx $GENOME_FILE
  fi
else
  samtools faidx $GENOME_FILE
fi

awk -v OFS='\t' {'print $1,$2'} $GENOME_FILE.fai > $TEMP_PATH/$5_genomeFile.txt
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
if [[ -s $ANNO_FILE ]] ; then
  if [[ $(ls -1 $bed_prefix* | parallel --max-procs 10 "if [[ ! -s {} ]] ; then echo 0; fi" | sort | uniq) == 0 ]] ; 
  then 
    #echo "Some BED files of trancsript regions are EMPTY"; 
    rm -f $bed_prefix*
    $PY3_PATH extract_transcript_regions.py -i $ANNO_FILE --gtf  -o $bed_prefix
  fi
else
  exit 2
fi


#if [[ "$7" == "TRUE" ]]; then

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

#readarray gene_list < $4
gene_list=$(cat $4 | sort | uniq | grep -v -w -i "gene")
s_names=$(parallel --max-procs $n_threads "echo {} | awk  -v s_var='_' '{ gsub(/[[:punct:]]/,s_var);}1' " ::: ${gene_list[@]})

echo "Checking Gene Names..."

rm files/genes/$f_org_name/1.list
rm files/genes/$f_org_name/2.1.list
rm files/genes/$f_org_name/2.2.list
rm files/genes/$f_org_name/genes.odb
rm files/genes/$f_org_name/full.list

time parallel --max-procs $n_threads "check_gene {} $f_org_name $ANNO_FILE $GENES_META " ::: ${gene_list[@]} #parallel -N 251 --pipe  #-j ${#gene_list[@]} #$(($n_threads * 2)) #get_fasta is called from check_gene function after preliminary checks 
#./jobhold.sh check_gene
##Looking up genes which were not present in GTF annotations (because of differing gene names between organisms)
awk '{gsub(" ","\n");print;}' files/genes/$f_org_name/*.list files/genes/$f_org_name/genes.odb | sort | uniq > files/genes/$f_org_name/full.list
#if [[ -s files/genes/$f_org_name/genes.odb ]]; then
  check_OrthoDB files/genes/$f_org_name/full.list $f_org_name $ANNO_FILE $GENES_META $4 ;
#fi

##REFRESH geene list based on file names
gene_list=($(find $GENES_META/$f_org_name/*.gtf_slice | parallel "basename {}" | cut -d '.' -f 1 | sort | uniq ))

## GET GTF stats - UTR lengths, CDS_count, exon_count (both fully and partially annotated info (eg, UTR info is hidden in exon and CDS info and is extracted using exon boundaries vs CDS boundaries)) 
rm files/genes/$f_org_name/gtf_stats.csv
time parallel --max-procs 2 "if [[ -s $GENES_META/$f_org_name/{}.rna_list && ! -z {} && -s $GENES_META/$f_org_name/{}.gtf_slice && ! -s $GENES_META/$f_org_name/{}.gtf_info ]] ; then
    ./jobhold.sh get_GTF_info Rscript extract_gtf_info.R $GENES_META/$f_org_name/{}.gtf_slice $f_org_name {} $GENES_META/$f_org_name/{}.gtf_info ;
fi" ::: ${gene_list[@]} 1>> $TEMP_PATH/get_GTF_info.o 2>> $TEMP_PATH/get_GTF_info.e #$ANNO_FILE
# find files/genes/xenopus_tropicalis/*.gtf_info -exec sed 1d {} \;
find $GENES_META/$f_org_name/*.gtf_info -exec sed 1d {} \; > files/genes/$f_org_name/gtf_stats.csv

gene_list=($(awk -F',' '{print $2"\n"}' files/genes/$f_org_name/gtf_stats.csv | sort | uniq
))

parallel --max-procs $n_threads "if [[ ! -z {1} && ! -z {2} && -s $GENES_META/$f_org_name/{1}.rna_list ]] ; then
 get_fasta {1} {2} $bed_prefix $GENOME_FILE $f_org_name $FASTA_PATH $ANNO_FILE $TEMP_PATH $GENES_META ; 
fi" ::: ${gene_list[@]} ::: ${transcript_regions[@]}

if [[ $LABEL_SEQS ==  "TRUE" ]] ; then
parallel --max-procs 10 "if [[ ! -z {1} && ! -z {2} && -s $FASTA_PATH/$f_org_name/{1}.{2}.tmp ]] ; then 
  ./label_sequenceIDs.sh $f_org_name {1} $ANNO_FILE $FASTA_PATH/$f_org_name/{1}.{2} $FASTA_PATH/$f_org_name/{1}.{2}.tmp files/genes/$f_org_name/odb.final_map ;
fi" ::: ${gene_list[@]} ::: ${transcript_regions[@]}
fi

#find files/genes/$5 -iname "*.rna_list" -empty | awk -F'/' '{print $NF}' | sed 's/.rna_list//g' > $FASTA_PATH/$5/MISSING_GENES
#grep -v -i -f $FASTA_PATH/$5/MISSING_GENES $4 | sort > $FASTA_PATH/$5/AVAILABLE_GENES
sed 1d files/genes/$5/gtf_stats.csv | awk -F',' '{print $2}' | sort | uniq >  files/genes/$5/AVAILABLE_GENES
grep -v -i -f files/genes/$5/AVAILABLE_GENES $4 | sort | uniq > files/genes/$5/MISSING_GENES

find $FASTA_PATH/$f_org_name/ -empty -delete
find $GENES_META/$f_org_name/ -empty -delete
find files/genes/$f_org_name/ -empty -delete
find $FASTA_PATH/$f_org_name/ -type f -name "*.fai" -exec rm -f {} +

##Populate clusters(files/genes/<org>/ALL_CLUSTERS) and group FASTA sequences according to clusters and store it in files/groups/
#grep -r ">" $FASTA_PATH/$f_org_name/ | parallel --max-procs $n_threads --colsep ":>" --recend "\n" echo {2} | awk -F"$seqID_delimiter" '{print $5}' |  sort | uniq > files/genes/$f_org_name/ALL_CLUSTERS
grep -r ">" $FASTA_PATH/$f_org_name/ | awk -F ":>" '{print $2}' | awk -F"$seqID_delimiter" '{print $5}' |  sort | uniq > files/genes/$f_org_name/ALL_CLUSTERS
#parallel --max-procs $n_threads "printf '%s\t%s\n' {1} {2}" :::: <(grep -H -f files/genes/$f_org_name/ALL_CLUSTERS -r $FASTA_PATH/$f_org_name/ | awk -F':>' -v s_delim="$seqID_delimiter" '{split($2,a,s_delim); print $1"\t"$2"\t"a[5]}')  ::: $(echo ${transcript_regions[@]} | awk '{print;}') | parallel  --max-procs 1 --colsep '\t' --recend '\n'  "if [[ -s {1} && ! -z {2} && ! -z {1} && ! -z {3} ]] ; then samtools faidx {1}  {2} >> $GROUPS_PATH/{3}.{4} ; fi"
#parallel --max-procs $n_threads " printf '%s\t%s\n' {1} {2}" :::: <(grep -H -f files/genes/$f_org_name/ALL_CLUSTERS -r $FASTA_PATH/$f_org_name/ | awk -F':>' -v s_delim="$seqID_delimiter" '{split($2,a,s_delim); n=split($1,b,"."); print $1"\t"$2"\t"a[5]"\t"b[n]'}) | parallel  --max-procs 1 --colsep '\t' --recend '\n'  "if [[ -s {1} && ! -z {2} && ! -z {1} && ! -z {3} ]] ; then samtools faidx {1}  {2} >> $GROUPS_PATH/{3}.{4} ; fi"

if [[ $REMOVE_DOWNLOADS ==  "FALSE" ]] ; then
  mv $GENOME_FILE $GENOMES_PATH/$5.fa
  mv $ANNO_FILE $ANNOS_PATH/$5.gtf
  #bgzip -i $GENOMES_PATH/$5.fa
  #bgzip -i $ANNOS_PATH/$5.gtf
  gzip -f --fast $GENOMES_PATH/$5.fa # --rsyncable
  gzip -f --fast $ANNOS_PATH/$5.gtf # --rsyncable
else
  rm $GENOME_FILE
  rm $ANNO_FILE
fi

#rm $TEMP_PATH/$f_org_name.*.*
#rm files/genes/$5/*.gtf_slice

#rm $3*
#rm files/bed/$5*

#exit 1
