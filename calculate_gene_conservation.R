suppressMessages(require(tidyverse))
#library(tidyselect)
require(dplyr)
suppressMessages(require(stringr))
suppressMessages(require(stringi))
suppressMessages(require(data.table))
suppressMessages(require(gplots))
suppressMessages(require(showtext))
#library(ggplot2)
#library(plotly)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) blast results path (2) gene (3) identity tables output directory (4) [calculated] threshold (0-100) (0 - for everything, 100 - only exactly matching sequences (5) plots output directory (6) available orgs list (7) reference organisms list (8) organisms info directory", call.=FALSE)
}

set.seed(123)

font_add_google("Gochi Hand", "gochi")
#font_add_google("Schoolbell", "bell")

## Automatically use showtext to render text
showtext_auto()

param_file <- "parameters.txt"
param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^") #Convert multibyte seperator to one byte sep #read.table(param_file,sep="==")
delimiter <- as.character(param_table[which(param_table=="seqID_delimiter"),c(2)])

blast_results_path <- args[1]
gene <- args[2]
identities_out_path <- args[3]
thres_out <- args[4]
plots_out_path <- as.character(param_table[which(param_table=="plot_path"),c(2)])
aorgs_path <- args[5]
sorgs_path <- as.character(param_table[which(param_table=="ref_orgs"),c(2)])
oinfo_out_path <- args[6]
num_transcripts <- args[7]

write(paste(gene, ":(calculate_gene_conservation.R)"), stdout()) ##DBG messages
write(paste(gene, ":(calculate_gene_conservation.R)"), stderr()) ##DBG messages
#gene <- "rnf4"
#blast_results_path <- "rnf4"
#oinfo_out_path <- "files/orgs"
#identities_out_path <- "files/idents"
#thres_out <- "files/gene_thresholds.txt"
#plots_out_path <- "files/plots"
#aorgs_path <- file.path("files","available_orgs.txt")
#sorgs_path <- file.path("files","selected_ORGS.txt")

dir.create(identities_out_path)
dir.create(oinfo_out_path)
dir.create(plots_out_path)
print(blast_results_path)
if(file.exists(blast_results_path)){
blast_results <- read.table(file.path(blast_results_path))
}else{

    quit()
}

blast_results <- blast_results %>% mutate(from_org=unlist(lapply(stri_split_fixed(blast_results$V1,delimiter,n=4,tokens_only = T),function(x){return(x[3])}))) #(from_org=unlist(stri_split_fixed(blast_results$V1,delimiter,n=4))[3])
blast_results <- blast_results %>% mutate(to_org=unlist(lapply(stri_split_fixed(blast_results$V2,delimiter,n=4,tokens_only = T),function(x){return(x[3])}))) #(to_org=unlist(stri_split_fixed(blast_results$V2,delimiter,n=4))[3])

# done automatically => blast_results <- transform(blast_results, V3 = as.numeric(V3), V4 = as.numeric(V4)) ##Convert seq identity and blast score columns to numeric type
#head(blast_results)

orgs.available <- factor(scan(aorgs_path, character()))
orgs.selected <- factor(scan(sorgs_path, character()))

blast_results$V1 <- factor(blast_results$V1)
blast_results$V2 <- factor(blast_results$V2)
blast_results$from_org <- factor(blast_results$from_org)
blast_results$to_org <- factor(blast_results$to_org)

min_seq_identity <- data.frame(matrix(numeric(), nrow = length(orgs.selected), ncol = length(orgs.available)))
max_seq_identity <- data.frame(matrix(numeric(), nrow = length(orgs.selected), ncol = length(orgs.available)))
rownames(min_seq_identity) <- orgs.selected
rownames(max_seq_identity) <- orgs.selected
colnames(min_seq_identity) <- orgs.available
colnames(max_seq_identity) <- orgs.available

min_seq_identity[is.na(min_seq_identity)] <- 101
max_seq_identity[is.na(max_seq_identity)] <- -1


for (s_org in orgs.selected) {
  for (a_org in orgs.available) {
    if(s_org!=a_org){
    min_seq_identity[s_org,a_org] <- min(blast_results[which(blast_results$"from_org" == s_org & blast_results$"to_org" == a_org),"V3"] )
    max_seq_identity[s_org,a_org] <- max(blast_results[which(blast_results$"from_org" == s_org & blast_results$"to_org" == a_org),"V3"] )
    }
  }
}

min_seq_identity <- apply(min_seq_identity, MARGIN = 2 , function(x) replace(x, is.infinite(x), 0))
max_seq_identity <- apply(max_seq_identity, MARGIN = 2 , function(x) replace(x, is.infinite(x), 0))

if(identical(class(max_seq_identity),"numeric") && identical(class(min_seq_identity),"numeric")){  #If only one organism is selected then we have to adjust datatypes
  max_seq_identity <- t(as.data.frame(max_seq_identity))
  rownames(max_seq_identity) <- orgs.selected
  min_seq_identity <- t(as.data.frame(min_seq_identity))
  rownames(min_seq_identity) <- orgs.selected
}

#heatmap(as.matrix(max_seq_identity), col = heat.colors(256))
#heatmap(as.matrix(min_seq_identity), col = heat.colors(256))
#heatmap(as.matrix(max_seq_identity - min_seq_identity) , col = heat.colors(256))

#data <- expand.grid(X=rownames(max_seq_identity), Y=colnames(max_seq_identity))
#data$Z <- apply(data, MARGIN = 1, function(x){
#  max_seq_identity[x[1],x[2]]
#})

#p <- ggplot(data, aes(X, Y, fill= Z)) +  geom_tile()

#ggplotly(p, tooltip="text")
#head(min_seq_identity)
#head(max_seq_identity)
if(nrow(max_seq_identity)>1){
png(filename =paste(file.path(plots_out_path,gene),"max.png",sep = "."),width = 15, height = 15 , units = "in", res = 100)
heatmap.2(max_seq_identity, trace = "row", tracecol = "black", margins = c(12,12), scale = "none",cexRow=1.2,cexCol=1.2)
dev.off()}
if(nrow(min_seq_identity)>1){
png(filename =paste(file.path(plots_out_path,gene),"min.png",sep = "."),width = 15, height = 15 , units = "in", res = 100)
heatmap.2(min_seq_identity, trace = "row", tracecol = "black", margins = c(12,12),scale = "none",cexRow=1.2,cexCol=1.2)
dev.off()}

discardable_orgs <- orgs.available[which(is.na(match(orgs.available, levels(blast_results$to_org))))] # We have no sequences of the gene for these orgs
saveable_orgs <- orgs.available[which(!is.na(match(orgs.available, levels(blast_results$to_org))))]

nref_orgs <- c()
#nref_orgs <- orgs.selected[which(is.na(match(orgs.selected, levels(blast_results$to_org))))]

for (s_org in orgs.selected) {
  if(max_seq_identity[which(rownames(max_seq_identity)==s_org),which(colnames(max_seq_identity)==s_org)]==0){
    print(paste(gene,": is not present in : ", s_org))
    nref_orgs <- c(nref_orgs, s_org)
  }
}

if(any(is.na(match(orgs.selected, levels(blast_results$from_org))))){ ## Are all reference species blasted 
  print(paste(gene,": Either gene missing/Please rerun BLAST for ", gene," (",blast_results_path,")"))
  #quit(status = 2)
}

if(any(is.na(match(orgs.available, levels(blast_results$to_org))))){ ## Do blast results contain all the available species? 
    print(paste(gene, ": Species not blasted against: ")) # because gene not available
    print(discardable_orgs)
}

if(any(is.na(match(orgs.selected, levels(blast_results$to_org))))){
  print(paste("WARNING: ",gene, " not present in all reference species:"))
  print(nref_orgs)
} ##Between reference species

max_seq_identity[,discardable_orgs] <- NA
max_seq_identity[nref_orgs,] <- NA
min_seq_identity[,discardable_orgs] <- NA
min_seq_identity[nref_orgs,] <- NA
max_seq_identity <- as.data.frame(max_seq_identity)
min_seq_identity <- as.data.frame(min_seq_identity)

semi_ortho_orgs <- unique(colnames(max_seq_identity)[which(max_seq_identity == 0, arr.ind = T)[,"col"]])

semi_orths <- max_seq_identity[, c(semi_ortho_orgs)]

for (a_org in names(semi_orths)) {
  if(all(is.na(match(semi_orths[, (names(semi_orths) %in% a_org)],0)))){
    discardable_orgs <- c(discardable_orgs, a_org)
    semi_orths <- semi_orths[, !(names(semi_orths) %in% a_org)]
    semi_ortho_orgs <- semi_ortho_orgs[!semi_ortho_orgs %in% a_org]
  }
}

print("These organisms are semi-orthologous (ie, they are not orthologous to all reference species):")
print(semi_ortho_orgs)
print(semi_orths)

write.table(semi_orths, file = file.path(identities_out_path,paste(gene,".semiorgs",sep = "")),quote = F, row.names = T, col.names = T)

##Correct isoform representation, for all reference species
##As in, we do not care about the similarity of the reference sequence and it's isoforms. So we make it 0
##Having these values also gives rise to ridiculous scores for poorly conserved genes (camk2g scores - before correction = 100, after correction = 48) 
#min_seq_identity[which(min_seq_identity == 100, arr.ind = T)] <- 0
#min_seq_identity[,(names(min_seq_identity) %in% orgs.selected)] <- 0
for (s_org in orgs.selected) {
  min_seq_identity[s_org,s_org] <- 0
}

if(length(semi_ortho_orgs) > 0){
  max_seq_identity <- max_seq_identity[, !(names(max_seq_identity) %in% semi_ortho_orgs)]
  min_seq_identity <- min_seq_identity[,!(names(min_seq_identity) %in% semi_ortho_orgs)]
}

if(length(discardable_orgs) > 0){
  max_seq_identity <- max_seq_identity[,!(names(max_seq_identity) %in% discardable_orgs)]
  min_seq_identity <- min_seq_identity[,!(names(min_seq_identity) %in% discardable_orgs)]
}

##Omit NAs
min_gene_conservation <- as.integer(max(min(max_seq_identity,na.rm = T),max(min_seq_identity,na.rm = T)))
print(min_gene_conservation)
print(max_seq_identity)
low_orthology_orgs <- unique(names(max_seq_identity)[which(apply(max_seq_identity,MARGIN=c(1,2),function(x){ x < min_gene_conservation }), arr.ind = T)[,"col"]])

if(!identical(low_orthology_orgs, character(0))){
  print(paste("These organisms did not pass the minimum threshold(",min_gene_conservation,"):"))
  print(all_of(low_orthology_orgs))
}

saveable_orgs <- anti_join(as.data.frame(list(saveable_orgs), col.names=c("orgs")),as.data.frame(list(semi_ortho_orgs), col.names=c("orgs")))
saveable_orgs <- as.vector(t(saveable_orgs))

#values stored to files

if(length(discardable_orgs) > 0){
fwrite(list(discardable_orgs), file = file.path(oinfo_out_path,paste(gene,".dorgs",sep = "")))
}
if(length(saveable_orgs) > 0){
fwrite(list(saveable_orgs), file = file.path(oinfo_out_path,paste(gene,".sorgs",sep = "")))
}
if(length(low_orthology_orgs) > 0){
fwrite(list(low_orthology_orgs), file = file.path(oinfo_out_path,paste(gene,".lorgs",sep = "")))
}
if(length(nref_orgs) > 0){
fwrite(list(nref_orgs), file = file.path(oinfo_out_path,paste(gene,".nref",sep = "")))
}
write.table(min_seq_identity, file = file.path(identities_out_path,paste(gene,".min",sep = "")), quote = F, row.names = T, col.names = T)
write.table(max_seq_identity, file = file.path(identities_out_path,paste(gene,".max",sep = "")), quote = F, row.names = T, col.names = T)

fileConn<-file(thres_out,open = "at")
writeLines(paste(gene, min_gene_conservation,length(nref_orgs)==0,num_transcripts,sep = ","), fileConn)
close(fileConn)

