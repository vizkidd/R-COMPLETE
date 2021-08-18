##IF wc -l ref_orgs.txt > 1; create_profile alignment
##ELSE continue;
java -jar macse.jar -prog alignSequences -seq 

if [[ wc -l ref_orgs.txt > 2 ]]; then
	find files/fasta 
else continue;

##ALIGNMENT script to align and stitch different transcript regions
