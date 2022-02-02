#!/bin/bash

##TAKES the genelist file as input
##REQUIRES /data/meyer/software/maketree/makeTree.sh
if [ $# -eq 0 ]
  then
    echo "Give genelist as input(exiting)."
    exit
fi

genes=($(cat $1 | awk '{print $1}' | sort | uniq))

ALN_PATH=$(grep -i -w "alignments_path" parameters.txt | awk -F'=' '{print $2}') 

parallel -j ${#genes[@]} "/data/meyer/software/maketree/makeTree.sh $ALN_PATH/{}.aln" ::: ${genes[@]}