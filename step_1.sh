#!/bin/bash
#$ -N "get_genomes"
#$ -cwd
#$ -V
#$ -l data
#$ -l longrun
#$ -o "logs/get_genomes.o"
#$ -e "logs/get_genomes.e

#(cd files/genomes && curl  -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_10.0_genome.fasta.gz) #&
#(cd files/annos && curl -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xentr10.0/XENTR_9.1_Xenbase.v10-lift.gff3) #&
#gffread files/annos/XENTR_9.1_Xenbase.v10-lift.gff3 -T -O -E -o files/annos/XENTR_9.1_Xenbase.v10-lift.gtf
#gunzip -d files/genomes/XENTR_10.0_genome.fasta.gz
#samtools faidx files/genomes/XENTR_10.0_genome.fasta
#mv files/annos/XENTR_9.1_Xenbase.v10-lift.gtf files/annos/xenopus_tropicalis.gtf
#./extract_genomic_regions.sh "files/genomes/XENTR_10.0_genome.fasta" "files/annos/xenopus_tropicalis.gtf" "files/bed/XT_genes" files/genelist.txt xenopus_tropicalis files/fasta

#(cd files/genomes && curl -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_genome.fa.gz) #&
#(cd files/annos && curl -OJL http://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.2/XENLA_9.2_Xenbase.gff3) #&
#sed "s/gene=;/gene=vds.S;/g" files/annos/XENLA_9.2_Xenbase.gff3 > files/annos/XENLA_9.2_Xenbase_fixed.gff3
#gffread files/annos/XENLA_9.2_Xenbase_fixed.gff3 -T -O -E -o files/annos/XENLA_9.2_Xenbase.gtf
#gunzip -d files/genomes/XENLA_9.2_genome.fa.gz
#samtools faidx files/genomes/XENLA_9.2_genome.fa
#mv files/annos/XENLA_9.2_Xenbase.gtf files/annos/xenopus_laevis.gtf
#./extract_genomic_regions.sh "files/genomes/XENLA_9.2_genome.fa" "files/annos/xenopus_laevis.gtf" "files/bed/XL_genes" files/genelist.txt xenopus_laevis files/fasta

##/gnu/store/71wrbjz4xx53mhd40qfzm5czfzsc89sx-profile/bin/Rscript get_genomes_ensembl.R files/genomes files/annos $1  $(echo "$(pwd)/files/fasta")
Rscript get_genomes_ensembl.R files/genelist.txt || false #$(echo "$(pwd)/files/fasta")
