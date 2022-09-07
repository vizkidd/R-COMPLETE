#!/bin/bash
#$ -N "s_vegetal"
#$ -cwd
#$ -V
#$ -l data
# $ -l longrun
#$ -o "logs/s_vegetal.o" 
#$ -e "logs/s_vegetal.e" 

source functions.sh

GENOMES_PATH=$(grep -i -w "genomes_path" parameters.txt | check_param)
ANNOS_PATH=$(grep -i -w "annos_path" parameters.txt | check_param)

mkdir logs

# ###USER GENOMES AND ANNOTATIONS ARE DOWNLOADED HERE###########
# if [[ ! -s $GENOMES_PATH/xenopus_tropicalis.fa.gz ]]; then
# 	(cd $GENOMES_PATH && curl  -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_10.0_genome.fasta.gz) #&
# 	mv -f $GENOMES_PATH/XENTR_10.0_genome.fasta.gz $GENOMES_PATH/xenopus_tropicalis.fa.gz
# fi
# if [[ ! -s $ANNOS_PATH/xenopus_tropicalis.gtf ]]; then
# 	(cd $ANNOS_PATH && curl -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_9.1_Xenbase.v10-lift.gff3) #&
# 	gffread $ANNOS_PATH/XENTR_9.1_Xenbase.v10-lift.gff3 -T -O -E -o $ANNOS_PATH/XENTR_9.1_Xenbase.v10-lift.gtf
# 	mv $ANNOS_PATH/XENTR_9.1_Xenbase.v10-lift.gtf $ANNOS_PATH/xenopus_tropicalis.gtf
# fi

# if [[ ! -s $GENOMES_PATH/xenopus_laevis.fa.gz ]]; then
# 	(cd $GENOMES_PATH && curl -OJL https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz) #&
# 	mv -f $GENOMES_PATH/XENLA_10.1_genome.fa.gz $GENOMES_PATH/xenopus_laevis.fa.gz
# fi
# if [[ ! -s $ANNOS_PATH/xenopus_laevis.gtf ]]; then
# 	(cd $ANNOS_PATH && curl -OJL https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_GCF.gff3) #&
# 	#sed "s/gene=;/gene=vds.S;/g" $ANNOS_PATH/XENLA_9.2_Xenbase.gff3 > $ANNOS_PATH/XENLA_9.2_Xenbase_fixed.gff3
# 	#gffread $ANNOS_PATH/XENLA_9.2_Xenbase_fixed.gff3 -T -O -E -o $ANNOS_PATH/XENLA_9.2_Xenbase.gtf
# 	gffread $ANNOS_PATH/XENLA_10.1_GCF.gff3 -T -O -E -o $ANNOS_PATH/XENLA_10.1_GCF.gtf
# 	mv $ANNOS_PATH/XENLA_10.1_GCF.gtf $ANNOS_PATH/xenopus_laevis.gtf
# fi
# ###############################################################

# ./merge_OG2genes_OrthoDB.sh files/genelist.txt || false
Rscript get_genomes_ensembl.R files/genelist.txt 
