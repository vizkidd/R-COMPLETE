##CREATES CODON ANNOTATION REQUIRED bY RNADECODER and TRANSAT

library(seqinr)
library(stringi)
library(fs)

set.seed(123)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) Path to Alignment(MSA) file(FASTA only annotated with 'mrna_regions_delimiter' (in parameters.txt))", call.=FALSE)
}

aln_file <- args[1]
#lens_file <- args[2]

#aln_file <- "files/alns/plin.aln"
if(file.info(aln_file)$size == 0 || !file.exists(aln_file)){
  stop(paste(aln_file,"doesn't exist...!"), call.=FALSE)
}
param_file <- "parameters.txt"
param_table <- read.table(param_file,sep="=")
alignments_path <- as.character(param_table[which(param_table=="alignments_path"),c(2)])
mrna_regions_delimiter <- as.character(param_table[which(param_table=="mrna_regions_delimiter"),c(2)])

#lens_values <- unique(read.table(lens_file,header = T,sep = ","))

#if(nrow(lens_values) > 1){
#  stop(paste(aln_file,"not aligned properly or the lengths in the .lens file are errornous...!"), call.=FALSE)
#}
aln <- read.fasta(aln_file,as.string = T)
indices <- c()
indices <- lapply(aln, function(x){grepRaw(mrna_regions_delimiter,x,fixed = T,value = F,all=T)})
aln_length <-unique(unlist(lapply(aln, function(x){stri_length(x)}))) ##SUBTRACT by 2 to account for the two '|'s(mrna_regions_delimiter) that were introduced to break the mRNA into CDS and UTR regions
if(length(aln_length)==1){ ##IF all alignments are of same length
if(!is.null(indices)){
  if(all(unlist(lapply(indices,function(x){length(x)==2})))){ ##ONLY 2 indices should be present (one which breaks 5utr and cds and other which seperates cds and 3utr)
    if(length(unique(unlist(indices)))==2){ ##ALL indices shoud be the same otherwise there is a problem with the alignment
      mrna_sites <- unique(unlist(indices))
      
      codon_annotation <- stri_flatten(c(rep(0,mrna_sites[1]-1),mrna_regions_delimiter,rep(1:3,(mrna_sites[2]-mrna_sites[1]-1)/3),mrna_regions_delimiter,rep(0,aln_length-mrna_sites[2])))
      ##NOW remove the '|'s(mrna_regions_delimiter) from alignments and codon annotation
      codon_annotation <- gsub(mrna_regions_delimiter,"",codon_annotation,fixed = T)
      transcript_names <- names(aln)
      transcript_stats <- sapply(aln,function(x){print(summary(x))})
      #hist(unlist(transcript_stats["GC",]))
      transcript_seqs <- lapply(aln,function(x){gsub(mrna_regions_delimiter,"",getSequence(x,as.string = T),fixed = T)})
      write.fasta(sequences = transcript_seqs,names = transcript_names,file.out = aln_file,open = "w")
      writeLines(text=c(">codon",codon_annotation),sep = "\n",con = file(path(alignments_path,paste(path_ext_remove(basename(aln_file)),".ann",sep = "")),open="wt"))
    }
  }
}
}

##Check this piece of code to get varants from an MSA:https://www.biostars.org/p/9975/
