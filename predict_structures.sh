#!/bin/bash

gene_list=$1 #"files/genelist.txt"
reference_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'=' '{print $2}') 
ALN_PATH=$(grep -i -w "alignments_path" parameters.txt | awk -F'=' '{print $2}') 
FT_PATH=$(grep -i -w "fasttree_path" parameters.txt | awk -F'=' '{print $2}') 
FASTA_PATH=$(grep -i -w "fasta_path" parameters.txt | awk -F'=' '{print $2}')
RNADEC_PATH=$(grep -i -w "rnadecoder_path" parameters.txt | awk -F'=' '{print $2}')
PY2_PATH=$(grep -i -w "python2_path" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH=$(grep -i -w "python3_path" parameters.txt | awk -F'=' '{print $2}')

readarray genes < $gene_list

#parallel -j ${#genes[@]} "$FT_PATH -nt -gtr -quote -gamma -bionj  -slow  $ALN_PATH/{}.aln > $ALN_PATH/{}.tree" ::: ${genes[@]}

#parallel -j ${#genes[@]} "$TRANSAT_PATH -v -d --mantissa --indep_pvalues -a $ALN_PATH/{}.aln -t $ALN_PATH/{}.tree -c $ALN_PATH/{}.ann -out $ALN_PATH/{}.out   " ::: ${genes[@]}

#parallel -j ${#genes[@]} "python $TOOLS/rnadecoder/bin/fasta_to_col.py -f $ALN_PATH/{}.aln -a $ALN_PATH/{}.ann -o $ALN_PATH/{}.col -cs $TOOLS/rnadecoder/bin/fasta_ann2rnadecoder_col.pl" ::: ${genes[@]}

#parallel -j ${#genes[@]} "python $TOOLS/rnadecoder/bin/create_rnadec_xml.py -c $ALN_PATH/{}.col -t $ALN_PATH/{}.tree -o $ALN_PATH/{}.fold.xml -m fold" ::: ${genes[@]}
#parallel -j ${#genes[@]} "python $TOOLS/rnadecoder/bin/create_rnadec_xml.py -c $ALN_PATH/{}.col -t $ALN_PATH/{}.tree -o $ALN_PATH/{}.scan.xml -m scan" ::: ${genes[@]}

##CREATING directories and moving files for proper file strucutre
#rm -rf $ALN_PATH/rnadecoder_output/ 
#mkdir $ALN_PATH/rnadecoder_output/scan/
#mkdir $ALN_PATH/rnadecoder_output/fold/
#parallel -j ${#genes[@]} "mkdir $ALN_PATH/rnadecoder_output/scan/{}" ::: ${genes[@]}
#parallel -j ${#genes[@]} "mkdir $ALN_PATH/rnadecoder_output/fold/{}" ::: ${genes[@]}
#parallel -j ${#genes[@]} "mv $ALN_PATH/{}.fold.xml $ALN_PATH/rnadecoder_output/fold/" ::: ${genes[@]}
#parallel -j ${#genes[@]} "mv $ALN_PATH/{}.scan.xml $ALN_PATH/rnadecoder_output/scan/" ::: ${genes[@]}

#parallel -j ${#genes[@]} "$TOOLS/rnadecoder/bin/RNA-decoder $ALN_PATH/rnadecoder_output/fold/{}.fold.xml " ::: ${genes[@]}
#parallel -j ${#genes[@]} "$TOOLS/rnadecoder/bin/RNA-decoder $ALN_PATH/rnadecoder_output/scan/{}.scan.xml " ::: ${genes[@]}

##CREATING directories and moving files for proper file strucutre
rm -rf files/rnadecoder_output/ 
rm -rf files/transat_output/ 
mkdir files/rnadecoder_output/
mkdir files/rnadecoder_output/fold/
mkdir files/rnadecoder_output/scan/
mkdir files/transat_output/
parallel -j ${#genes[@]} "mkdir files/rnadecoder_output/scan/{}" ::: ${genes[@]}
parallel -j ${#genes[@]} "mkdir files/rnadecoder_output/fold/{}" ::: ${genes[@]}

proc_list=()

#parallel -j ${#genes[@]} "$FT_PATH -nt -gtr -quote -gamma -bionj  -slow  $ALN_PATH/{}.aln > $ALN_PATH/{}.tree" ::: ${genes[@]}
##CREATING ANNOTATION
parallel -j ${#genes[@]} "Rscript annotate_bases.R $ALN_PATH/{}.aln" ::: ${genes[@]}
##CREATING ANNOTATION for RNADECODER (TRANSAT accepts 0 for noncoding regions but RNADECODER only accepts 3)
parallel -j ${#genes[@]} "sed 's/0/3/g' $ALN_PATH/{}.ann > $ALN_PATH/{}.rd.ann " ::: ${genes[@]}
parallel -j ${#genes[@]} "$PY3_PATH $RNADEC_PATH/fasta_to_col.py -f $ALN_PATH/{}.aln -a $ALN_PATH/{}.rd.ann -o $ALN_PATH/{}.col -cs $TOOLS/rnadecoder/bin/fasta_ann2rnadecoder_col.pl" ::: ${genes[@]}

parallel -j ${#genes[@]} "$PY3_PATH $RNADEC_PATH/create_rnadec_xml.py -c $ALN_PATH/{}.col -t $ALN_PATH/{}.tree -o files/rnadecoder_output/fold/{}/{}.fold.xml -m fold" ::: ${genes[@]}
parallel -j ${#genes[@]} "$PY3_PATH $RNADEC_PATH/create_rnadec_xml.py -c $ALN_PATH/{}.col -t $ALN_PATH/{}.tree -o files/rnadecoder_output/scan/{}/{}.scan.xml -m scan" ::: ${genes[@]}

#parallel -j ${#genes[@]} "mv files/{}.fold.xml $ALN_PATH/rnadecoder_output/fold/" ::: ${genes[@]}
#parallel -j ${#genes[@]} "mv files/{}.scan.xml $ALN_PATH/rnadecoder_output/scan/" ::: ${genes[@]}
for gene_name in "${genes[@]}"
do
gene=$(echo "$gene_name" | xargs)
j_name=$(echo "$gene"_fold)
nohup ./jobhold.sh "$j_name do_rnadec.sh files/rnadecoder_output/fold/$gene/$gene.fold.xml" &>> logs/job_status.o &
proc_list+=("$!")
j_name=$(echo $(echo "$gene" | xargs)_scan)
nohup ./jobhold.sh "$j_name do_rnadec.sh files/rnadecoder_output/scan/$gene/$gene.scan.xml" &>> logs/job_status.o &
proc_list+=("$!")
j_name=$(echo $(echo "$gene" | xargs)_tsat)
nohup ./jobhold.sh "$j_name do_transat.sh $ALN_PATH/$gene.aln $ALN_PATH/$gene.tree $ALN_PATH/$gene.ann files/transat_output/$gene "  &>> logs/job_status.o &
proc_list+=("$!")
done

for proc_id in "${proc_list[@]}"
do
	echo $proc_id
	wait $proc_id #|| true
done

exit 0