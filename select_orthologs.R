#.libPaths(file.path("/mnt","d","Work","r_libs_ubuntu"))
#.libPaths()
suppressMessages(require(data.table))
suppressMessages(require(readr))
##Find orthologous sequences

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) gene name (2) path to BLAST (tablar) output (3) output directory path (4) threshold (0-100) (0 - for everything, 100 - only exactly matching sequences", call.=FALSE)
}
set.seed(123)

gene <- args[1]
table_path <- args[2]
out_path <- args[3]
thres <- as.numeric(args[4])
#print(thres)
#gene <- "anln"
#tables_path <- "files/tables"
#out_path <- "files/orths"

write(paste(gene, ":(select_orthologs.R)"), stdout())
write(paste(gene, ":(select_orthologs.R)"), stderr())

gene_table <- read.table(table_path, header = F)

gene_table <- gene_table[order(gene_table$V3,decreasing = T),]

gene_table <- gene_table[!gene_table$V3<thres,]

#gene_orths <- c(gene_table$V1, gene_table$V2)
gene_orths <- c(levels(factor(gene_table$V1)), levels(factor(gene_table$V2)))
#gene_orths <- levels(factor(gene_orths))

#print(gene_table)
#print(gene_orths)

write(gene_orths, file = file.path(out_path, paste(gene,".orths",sep = "")))
