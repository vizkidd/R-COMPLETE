from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys

##CODE can be used to modify sequence 'motifs' in certain positions. Here I am using it to ambiguate stop codons to NNN (3utr or cds sequence)

##ENTRYPOINT

parser = argparse.ArgumentParser()
parser.add_argument('-f','--file', help='input FASTA file ', required=True,action='store',type=str)
parser.add_argument('-rf','--frame_start', type=int, help='Reading frame [1-based] [(1)-3 or position where you want to start reading]', required=False,action='store')
parser.add_argument('-s','--size', help='input frame size (codons=3) ', required=True,action='store',type=int)
parser.add_argument('-p','--pos', type=str, help='Comma seperated positions of frame to mask (NOTE : [1-based index] [position of frame, not the nt base] [giving 0 for index would take the last postion] [accepts a single positions & negative values for indices])', required=False,action='store')
parser.add_argument('-m','--mask', type=str, help='Mask character [give single character]', required=False,action='store')
parser.add_argument('-a','--add', type=str, help='Add masking character(s) at the position [True/(False)]', required=False,action='store')
parser.add_argument('-r','--rep', type=str, help='Repeat mask character to fit frame length? [(True)/False]', required=False,action='store')
parser.add_argument('-mf','--mask_func', type=str, help='Mask function? [lower/upper]', required=False,action='store')
parser.add_argument('-cm','--comp_mask', type=str, help='Comma seperated mask(s) to compare with [mask length should be same as frame length]', required=False,action='store')
parser.add_argument('-cmm','--comp_mask_method', type=str, help='Compare mask method [(positive or +)/negative or -]', required=False,action='store')
parser.add_argument('-o','--out_file', type=str, help='Output file (if not given, printed to STD OUT)', required=False,action='store')
args = parser.parse_args()

#print(vars(args))

if (args.pos is None and args.comp_mask is None):
	parser.error('Need --pos and/or --comp')

if (args.mask is None and args.rep is not None):
	parser.error('--rep requires --mask')

if (args.mask is None and args.mask_func is None):
	parser.error('Give either --mask or --mask_func')

if (args.mask is not None and args.mask_func is not None):
	parser.error('Give either --mask or --mask_func. Not both')

if (args.mask_func is not None and (args.mask_func.lower() != "lower" and args.mask_func.lower() != "upper")):
	print(args.mask_func)
	parser.error('--mask_func only takes lower/upper as values')

if (args.add is not None and args.pos is None):
	parser.error('--add only works with --pos ')

if (args.add is not None and args.mask_func is not None):
	parser.error('Cannot use both --add and --mask_func (--add only works with --pos)')

if (args.comp_mask_method is not None and args.comp_mask is None):
	parser.error('Mask method needs masks to compare (--comp_mask_method needs --comp_mask)')

file=args.file
if(args.frame_start is None):
	frame_start=1
else:
	frame_start=args.frame_start
frame_size=args.size
if(args.pos is not None):
	pos=args.pos.split(',') ##can give multiple positions (comma-seperated values) NOTE : [1-based index] [position of frame, not the nt base] [giving 0 for index would take the last postion]
else:
	pos = None
mask_method=True
if(args.comp_mask_method is not None):
	if(args.comp_mask_method == "negative" or args.comp_mask_method == "-"):
		mask_method=False
mask_char=args.mask
mask_func = args.mask_func
rep_mask_char=args.rep
add_char=False
if(args.add is not None):
	if(args.add.lower() == "true" or args.add.lower() == "t"):
		add_char=True
#compare_mask=args.comp
if(args.comp_mask is not None):
	masks=args.comp_mask.split(',')
else:
	masks = None
out_file=args.out_file

if(mask_func is None and rep_mask_char is not None):
	if(rep_mask_char.lower() == "true" or rep_mask_char.lower() == "t"):
		mask_char = mask_char * frame_size
		#print(mask_char)

#print(rep_mask_char)
#print(mask_char)
#print(masks)

new_rec=[]

for rec in SeqIO.parse(file, "fasta"):
	split_seq = [rec.seq[i:i+frame_size] for i in range(frame_start-1, len(rec.seq), frame_size)]
	#split_seq = ([[split_seq[int(j)]=mask_char] for i,j in enumerate(pos)])

	#print(split_seq[len(split_seq)-1])

	if(pos is not None):
		for i,j in enumerate(pos):
			#if(int(j)==0):
			#	j=len(split_seq)
			#print(split_seq[int(j)-1])
			if(masks is not None):
				if(mask_method): ##POSITIVE COMPARE MASK	
					if(str(split_seq[int(j)-1]) in masks):
						if(mask_char is not None):
							if(add_char):
								split_seq.insert(int(j)-1, Seq(mask_char))
							else:						
								split_seq[int(j)-1] = mask_char; ## replace
						else: ##Apply mask function
							if(mask_func.lower() == "lower"):
								split_seq[int(j)-1] = str(split_seq[int(j)-1]).lower()
							else:
								split_seq[int(j)-1] = str(split_seq[int(j)-1]).upper()
				else: 		 ##NEGATIVE COMPARE MASK	
					if(str(split_seq[int(j)-1]) not in masks):
						if(mask_char is not None):
							if(add_char):
								split_seq.insert(int(j)-1, Seq(mask_char))
							else:						
								split_seq[int(j)-1] = mask_char; ## replace
						else: ##Apply mask function
							if(mask_func.lower() == "lower"):
								split_seq[int(j)-1] = str(split_seq[int(j)-1]).lower()
							else:
								split_seq[int(j)-1] = str(split_seq[int(j)-1]).upper()
			else:
				if(mask_char is not None):
						if(add_char):
							split_seq.insert(int(j)-1, Seq(mask_char))
						else:
							split_seq[int(j)-1] = mask_char; ## replace or add
				else: ##Apply mask function
					if(mask_func.lower() == "lower"):
						split_seq[int(j)-1] = str(split_seq[int(j)-1]).lower()
					else:
						split_seq[int(j)-1] = str(split_seq[int(j)-1]).upper() # convert to zero-based index
	elif(masks is not None):
		for i,j in enumerate(split_seq):
			#print(split_seq[i])
			if(mask_method):	 ##POSITIVE COMPARE MASK			
				if(j in masks):
					if(mask_char is not None):
						split_seq[int(i)] = mask_char;
					else: ##Apply mask function
						if(mask_func.lower() == "lower"):
							#print(j)
							split_seq[int(i)] = str(split_seq[int(i)]).lower()
						else:
							split_seq[int(i)] = str(split_seq[int(i)]).upper()      
			else:			##NEGATIVE COMPARE MASK			
				if(j not in masks):
					if(mask_char is not None):
						split_seq[int(i)] = mask_char;
					else: ##Apply mask function
						if(mask_func.lower() == "lower"):
							#print(j)
							split_seq[int(i)] = str(split_seq[int(i)]).lower()
						else:
							split_seq[int(i)] = str(split_seq[int(i)]).upper()      

	#print(split_seq[-1],split_seq[0],rec.seq[0:frame_start-1] )
	rec.seq = rec.seq[0:frame_start-1] + Seq(''.join(str(codon) for codon in split_seq))
	if(out_file is None):
		print(">" + 	rec.id)
		print(rec.seq)
	new_rec.append(rec)
if(out_file is not None):
	with open(out_file, "w") as handle:
		SeqIO.write(new_rec,handle,"fasta")
