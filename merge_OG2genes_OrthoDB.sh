#!/bin/bash

##Optional : Give genelist to reduce the ODB file size to user specified genes

##NOTE: I should give link for downloading Fixed_ODB file

###ENTRYPOINT

source functions.sh

ORTHODB_PATH_PREFIX=$(grep -i -w "orthodb_path_prefix" parameters.txt | check_param) 
CLEAN_EXTRACT=$(grep -i -w "clean_extract" parameters.txt | check_param) 
n_threads=$(grep -i -w "max_concurrent_jobs" parameters.txt | check_param)

if [[ $CLEAN_EXTRACT == "TRUE" ]]; then
#	rm -f "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz
	rm -f "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
fi

if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz &> /dev/null ; then
	time join -t $'\t' -1 2 -2 1 <(zcat "$ORTHODB_PATH_PREFIX"_OG2genes.tab.gz | sort --parallel=$n_threads -k 2) <(zcat "$ORTHODB_PATH_PREFIX"_genes.tab.gz | sort --parallel=$n_threads -k 1) | awk -F "\t" '{if($2 in a)a[$2]=a[$2]","$1"||"$5;else a[$2]=$1"||"$5 ;}END{for(key in a)print key"\t"a[key];}' |  awk '{split($2,a,","); delete c; for(key in a){split(a[key],b,/\|\|/); if(b[2] in c==0){c[b[2]]=0;} } e=""; for(gene in c){ if(length(e)==0){e=gene;}else{e=gene","e;} } print $1"\t"$2"\t"e }' | gzip -c > "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz 
else
	>&1 color_FG $Green "File exists!"	
fi
>&1 echo $(color_FG $Green "Fixed ODB file stored in : ")$(color_FG_BG_Bold $White $BG_Purple "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz)

if [ $# -eq 1 ]
  then
    if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz &> /dev/null ; then
		time zgrep -f <(cat $1 | sort | uniq | grep -v -w -i "gene") "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz | gzip -c > "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
	else
		>&1 color_FG $Green "File exists!"	
	fi
    >&1 echo $(color_FG $Green "(User gene list) Fixed ODB file stored in : ")$(color_FG_BG_Bold $White $BG_Purple "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz)
fi

