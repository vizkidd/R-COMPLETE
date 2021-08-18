#!/bin/bash
#$ -N "find_orthologs"
#$ -cwd
#$ -V
#$ -l data
#$ -l longrun
#$ -o "logs/find_orthologs.o"
#$ -e "logs/find_orthologs.e


delete_empty_orgs() {
	#For each org
	#files/fasta - $1
	#files/genelist.txt - $2
	#files/UNAVAILABLE_ORGS -$3
	rm $3
	for org in $1/*
	do
		if [[ $(grep -c -v -i -f $org/MISSING_GENES $2) == 0 ]]; 
		then
		 	echo $org "is empty!...Removing"
			echo $org >> $3
		 	rm -rf $org
		fi 
	done
}
select_transcripts() {
#1 - transcript list file
#2 - intermediate output (transcripts from CURRENT organism)
#3 - final output (transcripts from ALL organisms)
#4 - organism fasta path
OFS=$'\n' grep -I -i -r -f $1 $4 > $2 ##Get the sequence IDs and file paths of all the valid transcripts
cat $2 >> $3
./select_transcripts.sh $2 $4
}

echo $1

source /home/vsuresh/.guix-profile/etc/profile

GENOMES_PATH=$(grep -i "genomes_path" parameters.txt | awk -F'=' '{print $2}')
ANNOS_PATH=$(grep -i "annos_path" parameters.txt | awk -F'=' '{print $2}')
FASTA_PATH=$(grep -i "fasta_path" parameters.txt | awk -F'=' '{print $2}')
TEMP_PATH=$(grep -i "temp_path" parameters.txt | awk -F'=' '{print $2}')
BLASTDB_PATH=$(grep -i "blastdb_path" parameters.txt | awk -F'=' '{print $2}')
REF_ORGS=$(grep -i "ref_orgs" parameters.txt | awk -F'=' '{print $2}')

delete_empty_orgs $FASTA_PATH $1 files/UNAVAILABLE_ORGS

find $FASTA_PATH -name MISSING_GENES| awk -F'/' '{print $(NF-1)","$0}' > files/MISSING_GENES_FINAL
find $FASTA_PATH -name AVAILABLE_GENES| awk -F'/' '{print $(NF-1)","$0}' > files/AVAILABLE_GENES_FINAL
find $FASTA_PATH/* -type d |  sort | uniq | awk -F'/' '{print $NF}' > files/available_orgs.txt

#rm files/temp/*.bed
#rm files/temp/*.tmp
#rm valid_ENSDART.txt
#rm selected_ENSDART.txt
#rm all_ENSDART.txt

find $FASTA_PATH -type f -name ".*" -execdir rm -f {} + #Deleting residues 
find $FASTA_PATH -type f -name "*.fai" -execdir rm -f {} +

awk -F',' 'NR>1 {print $1}' files/collab_Junker/danio_vegetal_pole_ENSDART.csv  > valid_ENSDART.txt ## these transcripts were found in the zebrafish cell
select_transcripts valid_ENSDART.txt selected_ENSDART.txt all_ENSDART.txt $FASTA_PATH/danio_rerio/
awk -F'\t' 'NR>1 {print $1}' files/collab_Junker/Xtropicalis_vegetal_isoforms.txt  > valid_ENSDART.txt ## these transcripts were found in the xenopus tropicalis cell
select_transcripts valid_ENSDART.txt selected_ENSDART.txt all_ENSDART.txt $FASTA_PATH/xenopus_tropicalis/ 
awk -F'\t' 'NR>1 {print $1}' files/collab_Junker/Xlaevis_vegetal_pole_isoforms.txt  > valid_ENSDART.txt ## these transcripts were found in the xenopus laevis cell
select_transcripts valid_ENSDART.txt selected_ENSDART.txt all_ENSDART.txt $FASTA_PATH/xenopus_laevis/

#printf "danio_rerio\nxenopus_tropicalis\nxenopus_laevis\n" > files/reference_ORGS.txt

##select organisms to search orthologs for
ls $FASTA_PATH > files/selected_ORGS.txt

./find_orthologs.sh files/selected_ORGS.txt $1 #100 ##This also selects the transcripts
