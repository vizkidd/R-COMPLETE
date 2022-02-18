#!/bin/bash
#$ -N "existing_loc"
#$ -cwd
#$ -V
#$ -l data
#$ -l longrun
#$ -o "logs/run_pipeline_ensembl.o"
#$ -e "logs/run_pipeline_ensembl.e"
#$ -m ea

#- STARTS the pipeline
#- the script takes only one parameter - file containing the gene list 

##ENTRYPOINT
# $1 is genelist.txt
if [ $# -eq 0 ]
  then
    echo "Give genelist as input(exiting)."
    exit
fi

echo $1

#source /home/vsuresh/.guix-profile/etc/profile

GENOMES_PATH=$(grep -i -w "genomes_path" parameters.txt | awk -F'=' '{print $2}')
ANNOS_PATH=$(grep -i -w "annos_path" parameters.txt | awk -F'=' '{print $2}')
FASTA_PATH=$(grep -i -w "fasta_path" parameters.txt | awk -F'=' '{print $2}')
TEMP_PATH=$(grep -i -w "temp_path" parameters.txt | awk -F'=' '{print $2}')
BLASTDB_PATH=$(grep -i -w "blastdb_path" parameters.txt | awk -F'=' '{print $2}')
REF_ORGS=$(grep -i -w "ref_orgs" parameters.txt | awk -F'=' '{print $2}')
#PY2_PATH=$(grep -i -w "python2_path" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH=$(grep -i -w "python3_path" parameters.txt | awk -F'=' '{print $2}')
PY3_PATH="${PY3_PATH/#\~/$HOME}"

source fasta_functions.sh

mkdir files
#mkdir files/genes
mkdir files/all2all
mkdir files/oneway
mkdir files/all2all_final
mkdir $ANNOS_PATH
mkdir files/genes
mkdir files/bed
mkdir $GENOMES_PATH
mkdir $FASTA_PATH
mkdir $TEMP_PATH

##OMIT for now
#rm $1
#awk -F',' 'BEGIN{OFS = "\n"} NR>1 {print $3}'  files/collab_Junker/count_tables/danio_vegetal_pole_ENSDART.csv | sed 's/.$//' > $1 ##Removing the last character in all genes
#awk -F',' 'BEGIN{OFS = "\n"} NR>1 {print $3}'  files/collab_Junker/count_tables/danio_vegetal_pole_ENSDART.csv  > files/genelist.txt

#dos2unix genelist.txt

##/gnu/store/71wrbjz4xx53mhd40qfzm5czfzsc89sx-profile/bin/Rscript fetch_fasta_biomaRt.R $1 ##very OLD CODE

#RIGHT NOW, we manually run the code for non-ensembl orgs
#(cd $GENOMES_PATH && curl  -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_10.0_genome.fasta.gz) #&
#(cd $ANNOS_PATH && curl -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_9.1_Xenbase.v10-lift.gff3) #&
#gffread $ANNOS_PATH/XENTR_9.1_Xenbase.v10-lift.gff3 -T -O -E -o $ANNOS_PATH/XENTR_9.1_Xenbase.v10-lift.gtf
#gunzip -d $GENOMES_PATH/XENTR_10.0_genome.fasta.gz
#samtools faidx $GENOMES_PATH/XENTR_10.0_genome.fasta
#mv $ANNOS_PATH/XENTR_9.1_Xenbase.v10-lift.gtf $ANNOS_PATH/xenopus_tropicalis.gtf
#./extract_genomic_regions.sh "$GENOMES_PATH/XENTR_10.0_genome.fasta" "$ANNOS_PATH/xenopus_tropicalis.gtf" "files/bed/XT_genes" $1 xenopus_tropicalis

#(cd $GENOMES_PATH && curl -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_genome.fa.gz) #&
#(cd $ANNOS_PATH && curl -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_Xenbase.gff3) #&
#sed "s/gene=;/gene=vds.S;/g" $ANNOS_PATH/XENLA_9.2_Xenbase.gff3 > $ANNOS_PATH/XENLA_9.2_Xenbase_fixed.gff3
#gffread $ANNOS_PATH/XENLA_9.2_Xenbase_fixed.gff3 -T -O -E -o $ANNOS_PATH/XENLA_9.2_Xenbase.gtf
#gunzip -d $GENOMES_PATH/XENLA_9.2_genome.fa.gz
#samtools faidx $GENOMES_PATH/XENLA_9.2_genome.fa
#mv $ANNOS_PATH/XENLA_9.2_Xenbase.gtf $ANNOS_PATH/xenopus_laevis.gtf
#./extract_genomic_regions.sh "$GENOMES_PATH/XENLA_9.2_genome.fa" "$ANNOS_PATH/xenopus_laevis.gtf" "files/bed/XL_genes" $1 xenopus_laevis

##/gnu/store/71wrbjz4xx53mhd40qfzm5czfzsc89sx-profile/bin/Rscript get_genomes_ensembl.R files/genomes files/annos $1  $(echo "$(pwd)/files/fasta")
time Rscript get_genomes_ensembl.R $1 || false

delete_empty_orgs $FASTA_PATH $1 files/UNAVAILABLE_ORGS

find $FASTA_PATH -name MISSING_GENES| awk -F'/' '{print $(NF-1)","$0}' > files/MISSING_GENES_FINAL
find $FASTA_PATH -name AVAILABLE_GENES| awk -F'/' '{print $(NF-1)","$0}' > files/AVAILABLE_GENES_FINAL
find $FASTA_PATH/* -type d |  sort | uniq | awk -F'/' '{print $NF}' > files/available_orgs.txt

count_seqs4genes $FASTA_PATH cds files/cds_BS.txt
count_seqs4genes $FASTA_PATH 3utr files/3utr_BS.txt
count_genes4orgs $FASTA_PATH $1 files/gene_counts_BS.txt

#rm files/temp/*.bed
#rm files/temp/*.tmp
#rm valid_ENSDART.txt
#rm selected_ENSDART.txt
#rm all_ENSDART.txt

find $FASTA_PATH -type f -name ".*" -execdir rm -f {} + #Deleting residues 
find $FASTA_PATH -type f -name "*.fai" -execdir rm -f {} +

#awk -F',' 'NR>1 {print $1}' files/collab_Junker/danio_vegetal_pole_ENSDART.csv  > valid_ENSDART.txt ## these transcripts were found in the zebrafish cell
#select_transcripts valid_ENSDART.txt selected_ENSDART.txt all_ENSDART.txt $FASTA_PATH/danio_rerio/
#awk -F'\t' 'NR>1 {print $1}' files/collab_Junker/Xtropicalis_vegetal_isoforms.txt  > valid_ENSDART.txt ## these transcripts were found in the xenopus tropicalis cell
#select_transcripts valid_ENSDART.txt selected_ENSDART.txt all_ENSDART.txt $FASTA_PATH/xenopus_tropicalis/ 
#awk -F'\t' 'NR>1 {print $1}' files/collab_Junker/Xlaevis_vegetal_pole_isoforms.txt  > valid_ENSDART.txt ## these transcripts were found in the xenopus laevis cell
#select_transcripts valid_ENSDART.txt selected_ENSDART.txt all_ENSDART.txt $FASTA_PATH/xenopus_laevis/

#printf "danio_rerio\nxenopus_tropicalis\nxenopus_laevis\n" > files/reference_ORGS.txt

count_seqs4genes $FASTA_PATH/ cds files/cds_AS.txt
count_seqs4genes $FASTA_PATH/ 3utr files/3utr_AS.txt
count_genes4orgs $FASTA_PATH $1 files/gene_counts_AS.txt

##select organisms to search orthologs for
ls $FASTA_PATH > files/selected_ORGS.txt
find files/genes -iname "gtf_stats.csv" -exec sed 1d {} \; > files/gtf_stats.csv

##ID Alignments - GENERATE numeric ids for FASTA IDS (because they are long and downstream analysis have difficulty taking long names) 
#x=0; grep ">" -H -r $FASTA_PATH | awk -F'>' -v x=$x '{ x = ++x; print substr($1, 1, length($1)-1) "\t" $2 "\t" x; }' > files/rna_ids.txt ##(full_filename  fasta_id  numeric_id)
index_fastaIDs files/rna_ids.txt $FASTA_PATH

time ./find_orthologs.sh files/selected_ORGS.txt $1 #100 ##This also selects the transcripts

##CHANGING FASTA IDs to numeric IDs because some programs dont work well with long FASTA IDs
#parallel -j0 -N 252 --pipe parallel -a files/rna_ids.txt -j 0 "fastaID_to_number" 
#parallel -j0 -a files/rna_ids.txt "fastaID_to_number" 

time ./align_seqs.sh $1

time ./predict_structures.sh $1

rm $TEMP_PATH/*
rm files/bed/*

exit



##MASK Stop codons to NNN. Check the bum-tail analogy as to why I am doing this
#for file in files/fasta/*/*.3utr;
#do
#	#echo $file
#	mask_stops_3utr $file
#done

#for file in files/fasta/xenopus_*/*.cds;
#do
#	#echo $file
#	mask_stops_cds $file
#done

count_seqs4genes files/fasta/ cds files/cds_OS.txt
count_seqs4genes files/fasta/ 3utr files/3utr_OS.txt
count_genes4orgs files/fasta $1 files/gene_counts_OS.txt

##./make_blast_db.sh files/fasta $1 files/ #Collate sequences from all organisms for each gene


##MASK Stop codons to NNN. Check the bum-tail analogy as to why I am doing this
#for file in files/fasta/*/*.3utr;
#do
#	#echo $file
#	mask_stops $file
#done

find files/fasta | grep -I -i -f $1 | awk -F'/' '{print $(NF-1),$(NF)}' | awk -F'.' '{print $1}' | sort | uniq > files/AVAILABLE_GENES_FINAL

rm files/stats.txt
/gnu/store/71wrbjz4xx53mhd40qfzm5czfzsc89sx-profile/bin/Rscript gene_stats.R >> files/stats.txt
