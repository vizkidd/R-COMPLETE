require(dplyr)
require(tidyverse)
require(tools)

#args = commandArgs(trailingOnly=TRUE)

#if (length(args)==0) {
#  stop("", call.=FALSE)
#}
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  #stop("Give the (1) folder path for storing genomes (2) folder path for storing annotations (3) the filename with gene symbols(without header)", call.=FALSE)
  stop("Give the (1) folder path for storing genomes (2) folder path for storing annotations (3) the filename with gene symbols(without header) (4) FASTA OUTPUT PATH(absolute path)", call.=FALSE)
}

set.seed(123)

thresholds <- read.csv("files/gene_thresholds.txt",sep = ",", header = F, col.names = c("gene","score","all_refs","num_transcripts"))
thresholds <- thresholds[order(thresholds$all_refs,thresholds$score,decreasing = T),]

#available.orgs <- factor(scan("files/available_orgs.txt", character()))
#org.counts <- data.frame(counts=rep(0, length(available.orgs)))
#rownames(org.counts) <- available.orgs

param_file <- "parameters.txt"
param_table <- read.table(param_file,sep="=")
fasta_path <- as.character(param_table[which(param_table=="fasta_path"),c(2)])
gene_drop_thresh <- as.numeric(param_table[which(param_table=="gene_drop_thresh"),c(2)])
ref_orgs <- as.character(param_table[which(param_table=="ref_orgs"),c(2)])

ref_orgs_list <- factor(scan(ref_orgs, character()))
orgs.selected <- dir(fasta_path) #factor(scan("files/reference_ORGS.txt", character()))

gene_list <- gsub('[[:punct:] ]+','_', factor(scan(args[1], character())))
genes.selected <- thresholds[which(!is.na(match(thresholds$all_refs,TRUE)) & thresholds$score > gene_drop_thresh),] #gsub('[[:punct:] ]+','_', factor(scan(args[2], character())))
genes.discarded <- thresholds[which(!is.na(match(thresholds$all_refs,FALSE)) | thresholds$score < gene_drop_thresh),] #gene_list[-match(genes.selected,gene_list)] 

ref <- data.frame(matrix(0, nrow=nrow(genes.selected), ncol=length(orgs.selected),
                       dimnames=list(genes.selected$gene, orgs.selected)),
                stringsAsFactors=F)
#rownames(ref) <- genes.selected


org.set <- c()
full_list <- c()
for (s_file in list.files(path="files/orgs/",pattern = "\\.sorgs$", full.names = T)) {
  if(!is.na(match(file_path_sans_ext(basename(s_file)),genes.selected$gene))){
      print(file_path_sans_ext(basename(s_file)))
      org_list <- c()
      org_list <- data.frame(orgs=scan(s_file, character()))
      full_list <- c(full_list, list(org_list))
      if(nrow(org_list) > 0){
        print(s_file)
        print(orgs.selected[which(orgs.selected %in% org_list$orgs)])
        ref[file_path_sans_ext(basename(s_file)),c(orgs.selected[which(orgs.selected %in% org_list$orgs)])] <- 1
        print(nrow(org_list))
        if(length(org.set)==0){
          org.set <- org_list 
        } else{
          org.set %>% inner_join(org_list)
        }
      }
    }
}

ref$sums <- rowSums(ref)

print(ref)

#print(paste("WARNING: ", " Best genes between organisms"))
#print(genes.selected)
#print(ref[which(ref$"sums"==3),])

org.counts <- data.frame(table(unlist(full_list)))
org.ranks <- org.counts[order(-org.counts$Freq),]
names(org.ranks) <- c("orgs", "Freq")
print(org.ranks)
best_orgs <- org.set
org.set <- org.ranks[sort(match(org.set$orgs, org.ranks$orgs)),]
print(org.set)
print(org.ranks[c(match(orgs.selected,org.ranks$orgs)),])
#print(paste("WARNING: ",min(org.ranks[c(match(orgs.selected,org.ranks$orgs)),]$Freq), " genes can be analyzed effectively"))


#if(!all(is.na(match(orgs.selected, best_orgs$orgs)))){ 
#  print("WARNING: Not all organisms are present in best orgs......adding them to the list anyways")
#  org.ranks[which(org.ranks$Var1 %in% orgs.selected[which(is.na(match(orgs.selected, best_orgs$orgs)))]),]$Var1
#}

#org.set2 <- c()
#for (s_file in list.files(path="files/orgs/",pattern = "\\.sorgs$", full.names = T)) {
#  if(!is.na(match(file_path_sans_ext(basename(s_file)),gsub('[[:punct:] ]+','_',genes.selected)))){
#    print(file_path_sans_ext(basename(s_file)))
#    org_list <- c()
#    org_list <- data.frame(orgs=scan(s_file, character()))
#    if(nrow(org_list) > 0){
#      print(s_file)
#      print(org_list)
#      print(nrow(org_list))
#      if(length(org.set2)==0){
#        org.set2 <- org_list 
#      } else{
#        #org.set3 %>% full_join(org_list)
#        org.set2 <- merge(org.set2, org_list, all=T)
#      }
#    }
#  }
#}

low_orthology_orgs <- c()
for (l_file in list.files(path="files/orgs/",pattern = "\\.lorgs$", full.names = T)) {
  if(!is.na(match(file_path_sans_ext(basename(l_file)),gsub('[[:punct:] ]+','_',genes.selected$gene)))){
    org_list <- c()
    org_list <- data.frame(orgs=scan(l_file, character()))
    if(nrow(org_list) > 0){
      print(l_file)
      #print(org_list)
      print(nrow(org_list))
      if(length(low_orthology_orgs)==0){
        low_orthology_orgs <- org_list 
      } else{
        low_orthology_orgs %>% inner_join(org_list)
        #low_orthology_orgs <- merge(low_orthology_orgs, org_list, all=T)
      }
    }
  }
}

discardable_orgs <- c()
for (l_file in list.files(path="files/orgs/",pattern = "\\.dorgs$", full.names = T)) {
  if(!is.na(match(file_path_sans_ext(basename(l_file)),gsub('[[:punct:] ]+','_',genes.selected$gene)))){
    org_list <- c()
    org_list <- data.frame(orgs=scan(l_file, character()))
    if(nrow(org_list) > 0){
      print(l_file)
      #print(org_list)
      print(nrow(org_list))
      if(length(discardable_orgs)==0){
        discardable_orgs <- org_list 
      } else{
        discardable_orgs %>% inner_join(org_list)
        #discardable_orgs <- merge(discardable_orgs, org_list, all=T)
      }
    }
  }
}


org.set1 <- anti_join(org.set,low_orthology_orgs)
org.set1 <- anti_join(org.set1, discardable_orgs)

write.table(org.set1, "files/oneway/oneway.orglist",na = "", col.names = F, row.names = F,quote = F)
write.table(org.set1$orgs, "files/oneway/SET",na = "", col.names = F, row.names = F,quote = F)
if( nrow(org.set1) > 30){
  org.sets <- split(org.set1, rep(1:3, length.out = nrow(org.set1), each = ceiling(nrow(org.set1)/3)))
  write.table(org.sets$`1`$orgs, "files/oneway/SET1",na = "", col.names = F, row.names = F,quote = F)
  write.table(org.sets$`2`$orgs, "files/oneway/SET2",na = "", col.names = F, row.names = F,quote = F)
  write.table(org.sets$`3`$orgs, "files/oneway/SET3",na = "", col.names = F, row.names = F,quote = F)
} 

#write.table(genes.selected, "files/GENES_SELECTED",na = "", col.names = F, row.names = F,quote = F)


