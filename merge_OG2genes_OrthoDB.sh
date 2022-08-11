#!/bin/bash

##Optional : Give genelist to reduce the ODB file size to user specified genes

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

###ENTRYPOINT

export -f check_param

ORTHODB_PATH_PREFIX=$(grep -i -w "orthodb_path_prefix" parameters.txt | check_param) 
n_threads=$(grep -i -w "max_concurrent_jobs" parameters.txt | check_param)

if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz ; then
	join -t $'\t' -1 2 -2 1 <(zcat "$ORTHODB_PATH_PREFIX"_OG2genes.tab.gz | sort --parallel=$n_threads -k 2) <(zcat "$ORTHODB_PATH_PREFIX"_genes.tab.gz | sort --parallel=$n_threads -k 1) | awk -F "\t" '{if($2 in a)a[$2]=a[$2]","$1"||"$5;else a[$2]=$1"||"$5 ;}END{for(key in a)print key"\t"a[key];}'  | gzip -c > "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz
fi
echo "Fixed ODB file stored in : $ORTHODB_PATH_PREFIX_OGgenes_fixed.tab.gz"

if [ $# -eq 1 ]
  then
    if ! gunzip -t "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz ; then
		zgrep -f <(cat $1 | sort | uniq | grep -v -w -i "gene") "$ORTHODB_PATH_PREFIX"_OGgenes_fixed.tab.gz | gzip -c > "$ORTHODB_PATH_PREFIX"_OGgenes_fixed_user.tab.gz
	fi
    echo "(User gene list) Fixed ODB file stored in : $ORTHODB_PATH_PREFIX_OGgenes_fixed_user.tab.gz"
fi

