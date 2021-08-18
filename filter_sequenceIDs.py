from Bio import SeqIO
from StringIO import StringIO
import sys
from os.path import exists
#import subprocess
#import multiprocessing
import re, mmap
from difflib import SequenceMatcher

def log_results(result):
    if(result is not None):
    	return(result)

def log_e(e):
  print("ERROR: "+e)

def contains_digits(d):
	_digits = re.compile(r'\d')
	return bool(_digits.search(d))

def filter_gene_name(rec,gene,gene_name,filter_val):
	#print(gene, gene_name)
	#print("anno:"+gene_name)
	#print("SG:"+gene)		
	if(gene_name[-2:] != "--" and filter_val.lower() != "none"):
		if(contains_digits(gene)):
			gene_group = gene[:re.search(r"\d", gene).start()]
		else:
			match = SequenceMatcher(None, gene, gene_name).find_longest_match(0, len(gene), 0, len(gene_name))
			gene_group=gene[match.a: match.a + match.size]
		if(contains_digits(gene_name)):
			anno_group = gene_name[:re.search(r"\d", gene_name).start()]
		else:
			match = SequenceMatcher(None, gene_name, gene).find_longest_match(0, len(gene_name), 0, len(gene))
			anno_group=gene_name[match.a: match.a + match.size]
		#print(gene_group)
		#print(anno_group)
		if(filter_val.lower() == "strict"):
			filter_strict = (gene_name.lower() == gene.lower() or gene_name.lower() == gene[:-1].lower() or (gene_group.lower() == anno_group.lower() and (len(gene_group)>1 and len(anno_group)>1)))
			if(filter_strict):
				#print(rec)
				return(rec)
		if(filter_val.lower() == "moderate"):
			#print("here1")
			filter_moderate = (gene.lower() in gene_name.lower() or gene[:-1].lower() in gene_name.lower() or (gene_group.lower() in gene_name.lower() and (len(gene_group)>1 and len(anno_group)>1)))
			if(filter_moderate):
				#print("here2")
				#print(rec)
				return(rec)
		#if(filter_val.lower() == "none"):
		#	return(rec)	
	else:
		return(rec)	
	#print(rec)
	return(None)



org=sys.argv[1]
out_file=sys.argv[2]
gene=sys.argv[3] ##Gene name from search group (files/genelist.txt), this is usedd to search GTFs for gene names
gene_name=sys.argv[4] ##Gene name from GTF
transcript_id=sys.argv[5]
filter_val="none" ##filter=strict/moderate/none
with open("parameters.txt",'r') as fh:
	value=[ keys for keys in fh.readlines() if ("genename_filter" in keys) ]
	filter_val=str(value[0].split("=")[1]).rstrip('\n')

seqID_delimiter="::"
with open("parameters.txt",'r') as fh:
	value=[ keys for keys in fh.readlines() if ("seqID_delimiter" in keys) ]
	seqID_delimiter=str(value[0].split("=")[1]).rstrip('\n')


#print(org, transcript_id, gene, gene_name, filter_val, out_file)
new_rec=[]
#existing_rec=[]
#print(org, file)
#print(sys.stdin)
rec=SeqIO.parse(sys.stdin,"fasta")
#print(rec)

for index, record in enumerate(rec):
	#print(record)
	record.id = record.id + seqID_delimiter + gene + seqID_delimiter + org + seqID_delimiter + gene_name
	record.description = record.description + seqID_delimiter + gene + seqID_delimiter + org + seqID_delimiter + gene_name
	#new_rec.append(log_results(filter_gene_name(record,gene,gene_name,filter_val)))	
	new_rec.append(filter_gene_name(record,gene,gene_name,filter_val))


#if exists(out_file):
#	with open(out_file, 'r+') as f:
#		for record in SeqIO.parse(f, 'fasta'):
#			#print(record)
#			existing_rec.append(record)
#	f.close()
#print(new_rec)
new_rec=[x for x in new_rec if x is not None]
if new_rec is not None:
#	existing_rec.extend(new_rec)
#	print(sys.argv)
#	print(new_rec)
#	print(existing_rec)
	with open(out_file, 'a+') as f:
		SeqIO.write(new_rec,f,"fasta")
	f.close()


#quit()

# def _old_get_gene_name(rec,gtf,filter_val):
# 	#print("NAME:" + multiprocessing.current_process().name)
# 	with open(gtf, 'r+') as f:
# 		gtf_data = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
# 		try:
# 			rec_clean = re.compile(r"\b"+rec.id.split("||")[0] + r"\b.*", re.M|re.I) #re.compile(rec.id.split("||")[0] + ".*", re.M|re.I|)
# 			#print(rec_clean)
# 			matched_lines = re.findall(rec_clean, gtf_data)
# 			#print(matched_lines)
# 			gene_name=re.sub("""[\\s";]""","",re.findall(r'gene_name [^;]+', "".join(matched_lines))[0].replace("gene_name",""))
# 			#print(re.findall(r'gene_name [^;]+', "".join(matched_lines)))
# 			#if (gene.lower() not in gene_name.lower()):
# 			#	gene_name = gene_name+"--"	
# 		except:
# 			#IF gene name not found then add trailing '--' to gene-group name
# 			#print("NOT found")
# 			gene_name=gene+"--"
# 		finally:
# 			rec.id = rec.id + "::" + gene + "::" + org + "::" + gene_name
# 			rec.description = rec.description + "::" + gene + "::" + org + "::" + gene_name
# 			#rec.id = rec.id + "::" + org
# 			#rec.description = rec.description + "::" + org
# 			#print("anno:"+gene_name)
# 			#print("SG:"+gene)
			
# 			if(gene_name[-2:] != "--"):
# 				if(contains_digits(gene)):
# 					gene_group = gene[:re.search(r"\d", gene).start()]
# 				else:
# 					match = SequenceMatcher(None, gene, gene_name).find_longest_match(0, len(gene), 0, len(gene_name))
# 					gene_group=gene[match.a: match.a + match.size]
# 				if(contains_digits(gene_name)):
# 					anno_group = gene_name[:re.search(r"\d", gene_name).start()]
# 				else:
# 					match = SequenceMatcher(None, gene_name, gene).find_longest_match(0, len(gene_name), 0, len(gene))
# 					anno_group=gene_name[match.a: match.a + match.size]
# 				#print(gene_group)
# 				#print(anno_group)
# 				if(filter_val.lower() == "strict"):
# 					filter_strict = (gene_name.lower() == gene.lower() or gene_name.lower() == gene[:-1].lower() or (gene_group.lower() == anno_group.lower() and (len(gene_group)>1 and len(anno_group)>1)))
# 					if(filter_strict):
# 						#print(rec)
# 						return(rec)
# 				if(filter_val.lower() == "moderate"):
# 					filter_moderate = (gene.lower() in gene_name.lower() or gene[:-1].lower() in gene_name.lower() or (gene_group.lower() in gene_name.lower() and (len(gene_group)>1 and len(anno_group)>1)))
# 					if(filter_moderate):
# 						#print(rec)
# 						return(rec)		
# 			#print(tmp_rec)
# 			return(None)

###OLD CODE WAS HERE
###FILTERS
##filter_strict = lambda gene_name, gene : True if gene_name.lower() == gene.lower() else False
##filter_moderate = lambda gene_name, gene : True if gene.lower() in gene_name.lower() else False
#
#pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
####OLDER CODE WAS HERE
#for rec in SeqIO.parse(file, "fasta"):
#	pool.apply_async(get_gene_name,(rec,gene,anno_name,filter_val.lower()), callback=log_results)
#	#print(get_gene_name(rec,gtf,filter_val.lower()))
#pool.close()
#pool.join()
##print(new_rec)
##print(new_rec[0])
#SeqIO.write(new_rec,file,"fasta")

#####OLD CODE


####OLDER_CODE
#with open(gtf, 'r+') as f:
#	#gtf_data = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) #f.readlines() #f.read()#
#	#print(gtf_data)
#	#for rec in SeqIO.parse(file, "fasta"):
#			##TRYING to get the name of the gene IN the annotation
#			##COULDNT EXECUTE PERL CMD IN BASH FROM PYTHON. SRY about the mess xD
#			#," | ", "perl -lne ", "'print ",'"@m"'," if @m=(/((?:gene_name)\s+\S+)/g);' | awk '{print $2}' | ","sed ",""'s/[";]//g'""," | sort | uniq "])
#			rec_clean = rec.id.split("_")[0]
#			cmd_txt = "".join(["grep -i -w ", '"',rec_clean,'"'," ",gtf])
#			#print(cmd_txt)
#			#print(sys.path[0])
#			bash_cmd = subprocess.Popen(cmd_txt, stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=sys.path[0], shell=True)
#			output, errors = bash_cmd.communicate()
#			bash_cmd.wait()
#			print(output)
#			pattern = re.compile("gene_name .*", re.MULTILINE)
#			names = pattern.findall(output)
#			print(names)
#			##CONVERTING FROM LIST TO SET AND THEN BACK AGAIN TO LIST to keep only unique gene name. Crude way but works
#			gene_name = re.sub("""[";]""","","".join(list(set(names))).split(" ")[1])
#			#print(gene_name)
#	#gtf_data.close()
#f.close()
