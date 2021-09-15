#!/bin/bash

##ALIGNMENT script to align and stitch different transcript regions

make_db() {
echo "$1"
if [ -s "$1" ]; then
	makeblastdb -in "$1" -dbtype nucl  -hash_index || true # -parse_seqids -out $2
else
	echo "$1 empty..."
fi
}

extract_sequences() {
	#rm files/gene.tmp
	gene=$1
	region=$2
	#if [ ! -s $BLASTDB_PATH/$gene.$region ]; then
	if [[ ! -s $(blastdb_path -db $BLASTDB_PATH/$gene.$region) ]]; then
		cat $(grep -r -i -l "$gene" $FASTA_PATH | grep ".$region") > $BLASTDB_PATH/$gene.$region
		make_db $BLASTDB_PATH/$gene.$region
	fi
	
	grep -w "$gene" files/HSP.csv | awk -F',' '{print $1"\n"$2}' | grep -f $reference_ORGS | sort | uniq | sed "s/cds/$region/g" > $ALN_PATH/$gene.ref
	grep -w "$gene" files/HSP.csv | awk -F',' '{print $1"\n"$2}' | grep -v -f $reference_ORGS | sort | uniq | sed "s/cds/$region/g" > $ALN_PATH/$gene.org
	seqkit grep -n -f $ALN_PATH/$gene.ref $BLASTDB_PATH/$gene.$region > $ALN_PATH/$gene.ref_orgs.$region
	seqkit grep -n -f $ALN_PATH/$gene.org $BLASTDB_PATH/$gene.$region > $ALN_PATH/$gene.orgs.$region
	cat $ALN_PATH/$gene.ref_orgs.$region $ALN_PATH/$gene.orgs.$region > $ALN_PATH/$gene.$region
	#cat $(find files/fasta -iname "*.cds" | grep -f $reference_ORGS | grep "$gene") > $ALN_PATH/$gene.ref_orgs.cds
	#cat $(find files/fasta -iname "*.cds" | grep -v -f $reference_ORGS | grep "$gene") > $ALN_PATH/$gene.orgs.cds
}

align_cds() {
	if [[ $(wc -l $reference_ORGS) > 2 ]]; then
	##IF wc -l ref_orgs.txt > 1; create_profile alignment
		#java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/$gene.ref_orgs.cds -out_NT $ALN_PATH/$gene.ref_orgs_NT.cds -out_AA $ALN_PATH/$gene.ref_orgs_AA.cds
		#java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/$gene.orgs.cds -out_NT $ALN_PATH/$gene.orgs_NT.cds -out_AA $ALN_PATH/$gene.orgs_AA.cds
		#java -jar $MACSE_PATH -prog alignTwoProfiles -p1 $ALN_PATH/$gene.ref_orgs_NT.cds -p2 $ALN_PATH/$gene.orgs_NT.cds -out_NT $ALN_PATH/"$gene"_NT.cds.aln -out_AA $ALN_PATH/"$gene"_AA.cds.aln
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/{}.ref_orgs.cds -out_NT $ALN_PATH/{}.ref_orgs_NT.cds -out_AA $ALN_PATH/{}.ref_orgs_AA.cds" ::: ${genes[@]}
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/{}.orgs.cds -out_NT $ALN_PATH/{}.orgs_NT.cds -out_AA $ALN_PATH/{}.orgs_AA.cds" ::: ${genes[@]}
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog alignTwoProfiles -p1 $ALN_PATH/{}.ref_orgs_NT.cds -p2 $ALN_PATH/{}.orgs_NT.cds -out_NT $ALN_PATH/{}_NT.cds.aln -out_AA $ALN_PATH/{}_AA.cds.aln" ::: ${genes[@]}
	else 	
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/{}.cds -out_NT $ALN_PATH/{}_NT.cds.aln -out_AA $ALN_PATH/{}_AA.cds.aln" ::: ${genes[@]}
		#java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/$gene.cds -out_NT $ALN_PATH/"$gene"_NT.cds.aln -out_AA $ALN_PATH/"$gene"_AA.cds.aln
	fi

	##GAP REMOVAL
	#if [[ -s $ALN_PATH/"$gene"_NT.cds.aln ]]; then
	#	Rscript remove_Gaps.R $ALN_PATH/"$gene"_NT.cds.aln
	#fi
	parallel -j ${#genes[@]} "Rscript remove_Gaps.R $ALN_PATH/{}_NT.cds.aln" ::: ${genes[@]}
}

align_3utr() {

	parallel -j ${#genes[@]} "mafft-qinsi --inputorder --treeout --maxiterate 1000  --localpair --leavegappyregion --thread -1 $ALN_PATH/{}.3utr > $ALN_PATH/{}_NT.3utr.aln" ::: ${genes[@]}

	##GAP REMOVAL
	#if [[ -s $ALN_PATH/"$gene"_NT.cds.aln ]]; then
	#	Rscript remove_Gaps.R $ALN_PATH/"$gene"_NT.cds.aln
	#fi
	parallel -j ${#genes[@]} "Rscript remove_Gaps.R $ALN_PATH/{}_NT.3utr.aln" ::: ${genes[@]}
}

align_5utr() {

	parallel -j ${#genes[@]} "mafft-qinsi --inputorder --treeout --maxiterate 1000  --localpair --leavegappyregion --thread -1 $ALN_PATH/{}.5utr > $ALN_PATH/{}_NT.5utr.aln" ::: ${genes[@]}

	##GAP REMOVAL
	#if [[ -s $ALN_PATH/"$gene"_NT.cds.aln ]]; then
	#	Rscript remove_Gaps.R $ALN_PATH/"$gene"_NT.cds.aln
	#fi
	parallel -j ${#genes[@]} "Rscript remove_Gaps.R $ALN_PATH/{}_NT.5utr.aln" ::: ${genes[@]}
}

###ENTRYPOINT

gene_list="files/genelist.txt"
reference_ORGS=$(grep -i "ref_orgs" parameters.txt | awk -F'=' '{print $2}') 
ALN_PATH=$(grep -i "alignments_path" parameters.txt | awk -F'=' '{print $2}') 
MACSE_PATH=$(grep -i "macse_path" parameters.txt | awk -F'=' '{print $2}') 
MAFFT_PATH=$(grep -i "mafft_path" parameters.txt | awk -F'=' '{print $2}') 
BLASTDB_PATH=$(grep -i "blastdb_path" parameters.txt | awk -F'=' '{print $2}')
FASTA_PATH=$(grep -i "fasta_path" parameters.txt | awk -F'=' '{print $2}')

mkdir $ALN_PATH

if [ $MAFFT_PATH == "" ]; then
	MAFFT_PATH=$(which mafft)
fi
if [ $MACSE_PATH == "" ]; then
	MACSE_PATH=$(which macse)
fi

if [ $MAFFT_PATH == "" ] || [ $MACSE_PATH == "" ]; then
	echo "PLEASE install MACSE and MAFFT and give the path in parameters.txt"
fi

#rm $ALN_PATH/*

while IFS= read -r gene
do
	extract_sequences $gene cds
	extract_sequences $gene 3utr
	extract_sequences $gene 5utr
done < "$gene_list"

readarray genes < $gene_list

##ALIGN 
#align_cds
align_3utr
align_5utr
