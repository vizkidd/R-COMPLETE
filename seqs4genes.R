.libPaths("D:/Work/r_libs")
.libPaths("/mnt/d/Work/r_libs_ubuntu")
.libPaths("/data/meyer/viz/r_libs")
library(tidyverse)
library(reshape2)
library(ggplot2)
library(dplyr)

counts_BS <- read.csv("files/cds_OS.txt", header = F,sep = " ")
#rnames <- counts_BS$V1
#counts_BS <- as.data.frame(counts_BS[,2])
#rownames(counts_BS) <- rnames
#colnames(counts_BS) <- "V1"
counts_BS <- cbind(counts_BS, V3=rep("OS", nrow(counts_BS)))
#counts_BS$V2<- scale(counts_BS$V2)

counts_AS <- read.csv("files/cds_AS.txt", header = F, sep = " ")
#rnames <- counts_AS$V1
#counts_AS <- as.data.frame(counts_AS[,2])
#rownames(counts_AS) <- rnames
#colnames(counts_AS) <- "V1"
counts_AS <- cbind(counts_AS, V3=rep("AS", nrow(counts_AS)))
#counts_AS$V2<- scale(counts_AS$V2)

counts <- merge(counts_AS, counts_BS, all = T, by="V1")

counts <- counts[order(counts$V2.y, decreasing = T),]

counts <- cbind(counts,counts$V2.y-counts$V2.x)

colnames(counts) <- c("gene","AS","OS","diff")

print(counts)

if(all(counts$diff==0)){
  print("Per gene, No transcripts were lost after transcript selection (AS) and ortholog selection (OS)!")
}else{
  print(counts[which(counts$diff!=0), 1:3])
  print("Per gene, These numeber of transcripts were lost after transcript selection (AS) and ortholog selection (OS)!")
}

