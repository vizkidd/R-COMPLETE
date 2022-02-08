##Script source and credits: hongqin
##https://github.com/hongqin/Simple-reciprocal-best-blast-hit-pairs

import sys, re
import pickle
import csv
from Graph import Graph
from Graph import Vertex
#import multiprocessing

# ===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(argsCount):
	if len(sys.argv) < argsCount:
		print("Orthologous Gene Finder")
		print("By Lee Bergstrand (modified by Viz)\n")
		print("Please refer to source code for documentation\n")
		print("Usage: " + sys.argv[0] + " BLASTOUTPUT1 BLASTOUTPUT2 outfile \n")
		print("Examples:" + " python RBH-v1.py files/all2all/danio_rerio-xenopus_tropicalis.out files/all2all/xenopus_tropicalis-danio_rerio.out files/all2all/dr-xt.out")
		print("Note: BLAST output format is 6")
		print('-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames (qcovhsp) sstrand qlen slen qseq sseq nident positive"')
		print("This script is an improvement and modification of two scripts: (Thanks and citations)")
		print("*https://github.com/hongqin/Simple-reciprocal-best-blast-hit-pairs  :: This script only supports one hit queries (but amazing work), I modified it to support the graph based BackBLAST algorithm")
		print("*https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST  :: This script only supports proteomes (but amazing work), I use the HSP based graph propogation from this script")
		sys.exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

# -------------------------------------------------------------------------------------------------

##ENTRYPOINT

argsCheck(4)
debug = 9
minIdent=25 #default #taking from parameters file
with open("parameters.txt",'r') as fh:
	value=[ keys for keys in fh.readlines() if ("minIdent" in keys) ]
	minIdent=float(value[0].split("=")[1])

seqID_delimiter="::"
with open("parameters.txt",'r') as fh:
	value=[ keys for keys in fh.readlines() if ("seqID_delimiter" in keys) ]
	seqID_delimiter=str(value[0].split("=")[1]).rstrip('\n')

#Commenting this comment, we're cheking if frames are 1/1 or 3/3 #frames = "1/1" CAN'T DO FRAMES BECAUSE IT WONT ACCOUNT FOR THE PRESCENCE/ABSCENCE OF START/STOP CODONS IN CDS SEQS
infl1 = sys.argv[1]
infl2 = sys.argv[2]
outfile = sys.argv[3]

BLASTGraph = Graph()
forward_list = list()
final_list = list()
check_list = list()
#parse first BLAST results
FL1 = open(infl1, 'r')
#print(FL1)
D1 = {} #dictionary for BLAST file ONE
for Line in FL1:
	if ( Line[0] != '#' ):
		Line.strip()
		Elements = re.split('\t', Line)
		if(float(Elements[2]) >= minIdent): # and Elements[14] == frames):
			queryId = Elements[0]
			subjectId = Elements[1]
			query_gene = str(queryId).split(seqID_delimiter)[1]
			subject_gene = str(subjectId).split(seqID_delimiter)[1]
			if ( query_gene == subject_gene ):
				forward_list.append(Elements)
				#print(queryId, subjectId)
				BLASTGraph.addVertex(queryId)
				BLASTGraph.addVertex(subjectId)
				BLASTGraph.addEdge(queryId, subjectId, float(Elements[15]))
				if ( not ( queryId in D1.keys() ) ):
					D1[queryId] = list()
				if ( not ( subjectId in D1[queryId] ) ):
					D1[queryId].append(subjectId)
					#BLASTGraph.addEdge(queryId, subjectId, Elements[15]) #HSPs


##IMPROVE GRAPH
##
## >>> dir(Graph)
#['__contains__', '__doc__', '__init__', '__iter__', '__module__', 'addEdge', 'addVertex', 'getVertex', 'getVertices']
#>>> dir(Vertex)
#['__doc__', '__init__', '__module__', '__str__', 'addNeighbor', 'getConnections', 'getId', 'getWeight']
##

#print(BLASTGraph.getVertices())

#for vert in BLASTGraph.getVertices():
#	if ( BLASTGraph.getVertex(vert).getConnections() ):
#		print(BLASTGraph.getVertex(vert).getConnections()[0])

#exit()

if (debug): D1.keys() 

#parse second BLAST results
FL2 = open(infl2, 'r')
#print(FL2)
D2 = {}
for Line in FL2:
	if ( Line[0] != '#' ):
		Line.strip()
		Elements = re.split('\t', Line)
		if(float(Elements[2]) >= minIdent): # and Elements[14] == frames):
			queryId = Elements[0]
			subjectId = Elements[1]
			query_gene = str(queryId).split(seqID_delimiter)[1]
			subject_gene = str(subjectId).split(seqID_delimiter)[1]
			if ( query_gene == subject_gene ):
				final_list.append(Elements)
				#print(queryId, subjectId)
				BLASTGraph.addEdge(queryId, subjectId, float(Elements[15]))     #HSPs
				if ( not ( queryId in D2.keys() ) ):
					D2[queryId] = list()
				if ( not ( subjectId in D2[queryId] ) ):
					D2[queryId].append(subjectId)

if (debug): D2.keys() 

BackBlastOutput = list(forward_list)
final_list.append(forward_list)
#print(BLASTGraph)

backup_list = list()

for hit in forward_list:
	queryId = BLASTGraph.getVertex(hit[0])
	subjectId = BLASTGraph.getVertex(hit[1])
	query_gene = str(queryId).split(seqID_delimiter)[1]
	subject_gene = str(subjectId).split(seqID_delimiter)[1]
	#print(hit)
	#print(query_gene, subject_gene)
	if ( query_gene == subject_gene ): ##Remove hits between different genes
		#print("-------------------------------------")
		topBackHitScore = -1
		# Find the top score of the best reciprocal BLAST hit.
		for backHit in subjectId.getConnections():
			#print(backHit)
			backHitScore = subjectId.getWeight(backHit)  # The edge weight between the subject and its reciprocal BLAST hit is the BLAST score.
			#print(backHitScore)
			if backHitScore >= topBackHitScore:
				topBackHitScore = backHitScore

		#print(topBackHitScore)
		# Check if the query is the best reciprocal BLAST hit for the subject.
		deleteHit = False
		if queryId in subjectId.getConnections():
			BackHitToQueryScore = float(subjectId.getWeight(queryId))  # The edge weight between the subject and the query is the reciprocal BLAST score.
			#print(query_gene, subject_gene)
			#print(topBackHitScore)
			#print(BackHitToQueryScore)
			if BackHitToQueryScore < float(topBackHitScore):
			# If the query is not the best reciprocal BLAST hit simply delete it from the BackBlastOutput.
				deleteHit = True
				#print("Delete hit 1 (low score)")
		else: 
			deleteHit = True  # If the query is not a reciprocal BLAST hit simply delete it from the BackBlastOutput.
		#	#print(subjectId.getConnections())
		#	#print(hit in backup_list)
		#	#print("Delete hit 2 (no connections)")
	else:
		deleteHit = True
		

	if deleteHit:
		#print("Delete hit")
		#print(BackBlastOutput[BackBlastOutput.index(hit)])
		#print(BackHitToQueryScore, topBackHitScore)
		#print(BackBlastOutput[BackBlastOutput.index(hit)])
		backup_list.append(BackBlastOutput.pop(BackBlastOutput.index(hit)))  # Delete the forward BLAST hit from BackBlastOutput.
	#else:
	#	print(BackBlastOutput[BackBlastOutput.index(hit)])

print(len(forward_list))
print(len(BackBlastOutput))

#    else:
#		print(BackBlastOutput[BackBlastOutput.index(hit)])
#print(BackBlastOutput)

#Now, pick the share pairs

#SharedPairs={}
#for id1 in D1.keys():
#	value1 = D1[id1]
#	for each_hit in value1:
#		if ( each_hit in D2.keys() ):
#			if ( id1 in D2[each_hit] ) : #a shared best reciprocal pair
#				if ( not ( each_hit in SharedPairs.keys() ) ):
#					SharedPairs[each_hit] = list()
#				if ( not ( id1 in SharedPairs[each_hit] ) ):
#					SharedPairs[each_hit].append(id1)

#if (debug): SharedPairs 

#outfl = open( outfile, 'w')

#for k1 in SharedPairs.keys():
#	line = k1 + '\t' + ",".join(SharedPairs[k1]) + '\n'
#	outfl.write(line)
#	
#outfl.close()


#with open(outfile+'.pickle', 'wb') as handle:
#    pickle.dump(SharedPairs, handle, protocol=pickle.HIGHEST_PROTOCOL)

for row in BackBlastOutput:
	check_list.append(row[0])
	check_list.append(row[1])

try:
	writeFile = open(outfile, "w")
	print(">> Output file created.")
	print(">> Writing Data...")
	for row1 in final_list:
		if(row1[0] in check_list or row1[1] in check_list):
			writeFile.write("\t".join(row1))
	writeFile.close()
except IOError:
	print(">> Failed to create " + outfile)
	sys.exit(1)

FL1.close()
FL2.close()

print("Done. RBH from", sys.argv[1], "and", sys.argv[2], "are in", sys.argv[3])
