library(ape)
library(Biostrings)

set.seed(123)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) Path to Alignment file(FASTA only) ", call.=FALSE)
}

aln_file <- args[1]
#aln_file <- "files/alns/plin.aln"
if(file.info(aln_file)$size == 0 || !file.exists(aln_file)){
  stop(paste(aln_file,"doesn't exist...!"), call.=FALSE)
}
param_file <- "parameters.txt"
param_table <- read.table(param_file,sep="=")
gap_threshold <- as.numeric(param_table[which(param_table=="gap_threshold"),c(2)])
plot_out_path <- as.character(param_table[which(param_table=="plot_path"),c(2)])


pdf(file =paste(plot_out_path,paste(basename(aln_file),".pdf",sep=""),sep="/"),title = basename(aln_file))
#aln <- readDNAMultipleAlignment("files/alns/dazl_NT.aln",format = "fasta")
aln <- read.FASTA(aln_file)
#print(aln_file)
#print(aln)
plot.new()
print(text(.5, .5, "BEFORE GAP REMOVAL"))
print(checkAlignment(as.matrix(aln)))
#print(del.colgapsonly(as.matrix(aln),threshold = gap_threshold))
aln <- as.DNAbin(as.alignment(del.colgapsonly(as.matrix(aln),threshold = gap_threshold)))
aln <- as.DNAbin(as.alignment(del.rowgapsonly(as.matrix(aln),threshold = 1)))
#print(aln)
plot.new()
print(text(.5, .5, "AFTER GAP REMOVAL"))
print(checkAlignment(as.matrix(aln)))
dev.off()
write.dna(aln,aln_file, colsep = "",nbcol=-1, format = "fasta")