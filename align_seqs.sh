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
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog trimNonHomologousFragments -out_NT $ALN_PATH/{}_NT.ref_orgs.cds.trimmed -seq $ALN_PATH/{}.ref_orgs.cds -out_AA $ALN_PATH/{}_AA.ref_orgs.cds.trimmed -out_mask_detail $ALN_PATH/{}_NT.ref_orgs.cds.mask -out_trace $ALN_PATH/{}_NT.ref_orgs.cds.trace -out_trim_info $ALN_PATH/{}_NT.ref_orgs.cds.trim_info" ::: ${genes[@]}
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog trimNonHomologousFragments -out_NT $ALN_PATH/{}_NT.orgs.cds.trimmed -seq $ALN_PATH/{}.orgs.cds -out_AA $ALN_PATH/{}_AA.orgs.cds.trimmed -out_mask_detail $ALN_PATH/{}_NT.orgs.cds.mask -out_trace $ALN_PATH/{}_NT.orgs.cds.trace -out_trim_info $ALN_PATH/{}_NT.orgs.cds.trim_info" ::: ${genes[@]}
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/{}_NT.ref_orgs.cds.trimmed -out_NT $ALN_PATH/{}.ref_orgs_NT.cds -out_AA $ALN_PATH/{}.ref_orgs_AA.cds" ::: ${genes[@]}
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/{}_NT.orgs.cds.trimmed -out_NT $ALN_PATH/{}.orgs_NT.cds -out_AA $ALN_PATH/{}.orgs_AA.cds" ::: ${genes[@]}
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog alignTwoProfiles -p1 $ALN_PATH/{}.ref_orgs_NT.cds -p2 $ALN_PATH/{}.orgs_NT.cds -out_NT $ALN_PATH/{}_NT.cds.aln.fm -out_AA $ALN_PATH/{}_AA.cds.aln.fm" ::: ${genes[@]}
	else 	
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog trimNonHomologousFragments -out_NT $ALN_PATH/{}_NT.cds.trimmed -seq $ALN_PATH/{}.cds -out_AA $ALN_PATH/{}_AA.cds.trimmed -out_mask_detail $ALN_PATH/{}_NT.cds.mask -out_trace $ALN_PATH/{}_NT.cds.trace -out_trim_info $ALN_PATH/{}_NT.cds.trim_info" ::: ${genes[@]}
		parallel -j ${#genes[@]} "java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/{}_NT.cds.trimmed -out_NT $ALN_PATH/{}_NT.cds.aln.fm -out_AA $ALN_PATH/{}_AA.cds.aln.fm" ::: ${genes[@]}
		#java -jar $MACSE_PATH -prog alignSequences -seq $ALN_PATH/$gene.cds -out_NT $ALN_PATH/"$gene"_NT.cds.aln -out_AA $ALN_PATH/"$gene"_AA.cds.aln
	fi

	##GAP REMOVAL
	#if [[ -s $ALN_PATH/"$gene"_NT.cds.aln ]]; then
	#	Rscript remove_Gaps.R $ALN_PATH/"$gene"_NT.cds.aln
	#fi
	
	##CONVERTING frameshift mutations(!)(introduced by MACSE) to gaps(-)
	parallel -j ${#genes[@]} "sed 's/!/-/g' $ALN_PATH/{}_NT.cds.aln.fm > $ALN_PATH/{}_NT.cds.aln" ::: ${genes[@]}
	parallel -j ${#genes[@]} "Rscript remove_Gaps.R $ALN_PATH/{}_NT.cds.aln" ::: ${genes[@]}
}

align_3utr() {

	parallel -j ${#genes[@]} "$MAFFT_PATH  --inputorder --maxiterate 1000  --localpair --leavegappyregion --thread -1 $ALN_PATH/{}.3utr > $ALN_PATH/{}_NT.3utr.aln" ::: ${genes[@]} #--treeout #--treein $ALN_PATH/{}.cds.mtree

	##GAP REMOVAL
	#if [[ -s $ALN_PATH/"$gene"_NT.cds.aln ]]; then
	#	Rscript remove_Gaps.R $ALN_PATH/"$gene"_NT.cds.aln
	#fi
	parallel -j ${#genes[@]} "Rscript remove_Gaps.R $ALN_PATH/{}_NT.3utr.aln" ::: ${genes[@]}
}

align_5utr() {

	parallel -j ${#genes[@]} "$MAFFT_PATH  --inputorder --maxiterate 1000  --localpair --leavegappyregion --thread -1 $ALN_PATH/{}.5utr > $ALN_PATH/{}_NT.5utr.aln" ::: ${genes[@]} #--treeout #--treein $ALN_PATH/{}.cds.mtree

	##GAP REMOVAL
	#if [[ -s $ALN_PATH/"$gene"_NT.cds.aln ]]; then
	#	Rscript remove_Gaps.R $ALN_PATH/"$gene"_NT.cds.aln
	#fi
	parallel -j ${#genes[@]} "Rscript remove_Gaps.R $ALN_PATH/{}_NT.5utr.aln" ::: ${genes[@]}
}

stitch_alns() {
	for gn in ${genes[@]}
	do
		echo $gn
		min_seqs=$(printf 3utr,$(grep ">" $ALN_PATH/"$gn"_NT.3utr.aln | wc -l)"\n"cds,$(grep ">" $ALN_PATH/"$gn"_NT.cds.aln | wc -l)"\n"5utr,$(grep ">" $ALN_PATH/"$gn"_NT.5utr.aln | wc -l) | sort -k2  | head -n 1)
		reg=$(echo $min_seqs | awk -F',' '{print $1}')
		tr_list=($(grep ">" $ALN_PATH/"$gn"_NT."$reg".aln | sed 's/>//g' | sed "s/$reg/full/g"))
		rm $ALN_PATH/"$gn".aln
		rm $ALN_PATH/"$gn".lens
		#printf "five_start,five_end,cds_start,cds_end,three_start,three_end\n" > $ALN_PATH/"$gn".lens
		#printf "five_len,cds_len,three_len\n" > $ALN_PATH/"$gn".lens
		for tr in ${tr_list[@]}
		do
			five_seq=$(seqkit grep -n -p $(echo $tr | sed "s/full/5utr/g") $ALN_PATH/"$gn"_NT.5utr.aln | grep ">" -v | sed "s/ //g" )
			cds_seq=$(seqkit grep -n -p $(echo $tr | sed "s/full/cds/g") $ALN_PATH/"$gn"_NT.cds.aln | grep ">" -v | sed "s/ //g" )
			three_seq=$(seqkit grep -n -p $(echo $tr | sed "s/full/3utr/g") $ALN_PATH/"$gn"_NT.3utr.aln | grep ">" -v | sed "s/ //g" )
			five_len=$(echo $five_seq | wc -m )
			cds_len=$(echo $cds_seq | wc -m )
			three_len=$(echo $three_seq | wc -m )
			if [ "$five_seq" != "" ] && [ "$cds_seq" != "" ] && [ "$three_seq" != "" ]; then
				##printf -- 1,"$five_len",$(($five_len + 1)),$(($five_len + 1 + $cds_len )),$(( $five_len + 1 + $cds_len +1 )),$(( $five_len + 1 + $cds_len + $three_len ))"""\n"
				printf -- ">$tr\n$five_seq$mrna_regions_delimiter$cds_seq$mrna_regions_delimiter$three_seq\n" | sed "s/ //g" >> $ALN_PATH/"$gn".aln #| sed 's/!/-/g' >> $ALN_PATH/"$gn".aln
				##printf -- "$(echo $five_seq | wc -m )","$(echo $cds_seq | wc -m )","$(echo $three_seq | wc -m )""\n" >> $ALN_PATH/"$gn".lens
				##printf -- 1,"$five_len",$(($five_len + 1)),$(($five_len + 1 + $cds_len )),$(( $five_len + 1 + $cds_len +1 )),$(( $five_len + 1 + $cds_len + $three_len ))"""\n" >> $ALN_PATH/"$gn".lens
				#printf -- "$five_len,$cds_len,$three_len\n" >> $ALN_PATH/"$gn".lens
			fi
		done

		#awk '/^>/ {print (NR>1?"\n":"")$0;; next} {printf "%s",$0;} END{print "";}' $ALN_PATH/"$gn".tmp > $ALN_PATH/"$gn".aln
		#rm $ALN_PATH/"$gn".tmp
		#printf -- "$five_len,$cds_len,$three_len\n"
		tar -czvf $ALN_PATH/$gn.tar $ALN_PATH/"$gn".3utr $ALN_PATH/"$gn"_AA*.trimmed $ALN_PATH/"$gn"_AA.cds.aln.fm $ALN_PATH/"$gn"_NT*.trimmed $ALN_PATH/"$gn"_NT.*.cds.trim_info $ALN_PATH/"$gn"_NT.*.cds.mask $ALN_PATH/"$gn"_NT.*.cds.trace $ALN_PATH/"$gn".cds $ALN_PATH/"$gn".5utr $ALN_PATH/"$gn"_NT.*.aln $ALN_PATH/"$gn"_NT.cds.* $ALN_PATH/"$gn".*.tree $ALN_PATH/"$gn".*.mtree $ALN_PATH/"$gn"_NT.cds.aln.fm $ALN_PATH/"$gn".orgs.* $ALN_PATH/"$gn".ref_orgs.* $ALN_PATH/"$gn".ref_orgs_NT.* $ALN_PATH/"$gn".orgs_NT.* $ALN_PATH/"$gn".ref_orgs_AA.* $ALN_PATH/"$gn".orgs_AA.* #$ALN_PATH/"$gn"_AA.*.aln
		rm $ALN_PATH/"$gn".3utr $ALN_PATH/"$gn"_AA*.trimmed $ALN_PATH/"$gn"_NT*.trimmed $ALN_PATH/"$gn"_AA.cds.aln.fm $ALN_PATH/"$gn"_NT.*.cds.trim_info $ALN_PATH/"$gn"_NT.*.cds.mask $ALN_PATH/"$gn"_NT.*.cds.trace $ALN_PATH/"$gn"_NT.cds.* $ALN_PATH/"$gn".cds $ALN_PATH/"$gn".5utr $ALN_PATH/"$gn"_NT.*.aln $ALN_PATH/"$gn".*.tree $ALN_PATH/"$gn".*.mtree $ALN_PATH/"$gn"_NT.cds.aln.fm $ALN_PATH/"$gn".orgs.* $ALN_PATH/"$gn".ref_orgs.* $ALN_PATH/"$gn".ref_orgs_NT.* $ALN_PATH/"$gn".orgs_NT.* $ALN_PATH/"$gn".ref_orgs_AA.* $ALN_PATH/"$gn".orgs_AA.* #$ALN_PATH/"$gn"_AA.*.aln
	done

	##DOING this because its not a lot of MSA programs support special chars like frameshift mutations(!) or stitched alignments(5utr|cds|3utr) (which was something I was trying to do)
	#parallel -j ${#genes[@]} "sed 's/ //g' $ALN_PATH/{}.aln > $ALN_PATH/{}.aln" ::: ${genes[@]}
	#parallel -j ${#genes[@]} "Rscript remove_Gaps.R $ALN_PATH/{}.aln" ::: ${genes[@]}

}

###ENTRYPOINT

gene_list=$1 #"files/genelist.txt"
reference_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'=' '{print $2}') 
ALN_PATH=$(grep -i -w "alignments_path" parameters.txt | awk -F'=' '{print $2}') 
MACSE_PATH=$(grep -i -w "macse_path" parameters.txt | awk -F'=' '{print $2}') 
MAFFT_PATH=$(grep -i -w "mafft_path" parameters.txt | awk -F'=' '{print $2}') 
BLASTDB_PATH=$(grep -i -w "blastdb_path" parameters.txt | awk -F'=' '{print $2}')
FASTA_PATH=$(grep -i -w "fasta_path" parameters.txt | awk -F'=' '{print $2}')
FT_PATH=$(grep -i -w "fasttree_path" parameters.txt | awk -F'=' '{print $2}') 
mrna_regions_delimiter=$(grep -i -w "mrna_regions_delimiter" parameters.txt | awk -F'=' '{print $2}')

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
align_cds
parallel -j ${#genes[@]} "$FT_PATH -nt -gtr -quote -gamma -bionj  -slow  $ALN_PATH/{}_NT.cds.aln > $ALN_PATH/{}.cds.tree" ::: ${genes[@]}
#parallel -j ${#genes[@]} "ruby newick2mafft.rb 1 $ALN_PATH/{}.cds.tree > $ALN_PATH/{}.cds.mtree" ::: ${genes[@]}
align_3utr
align_5utr

stitch_alns

