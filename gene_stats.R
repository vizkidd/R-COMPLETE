#.libPaths("D:/Work/r_libs")
#.libPaths("/mnt/d/Work/r_libs_ubuntu")
.libPaths("/data/meyer/viz/r_libs")
library(tidyverse)
library(reshape2)
library(ggplot2)

set.seed(123)

gene_orgs <- read.csv("files/AVAILABLE_GENES_FINAL", header = F, sep=" ")
gene_orgs <- split(gene_orgs,f = factor(gene_orgs$V1))
print(gene_orgs %>% reduce(inner_join, by="V2") %>% summarize(V2) %>% arrange(V2))

print("These genes are available in all organisms")

file_paths <- read.csv("files/MISSING_GENES_FINAL", header = F)
file_paths <- split(file_paths$V2, file_paths$V1)
org_list <- c()

org_list <- lapply(file_paths, FUN = function(x){
  factor(as.vector(t(read.table(as.character(x),header = F))))
})

org_list <- melt(org_list)

stats <- apply(org_list,2, table)

print(stats)

print("$value - Number of organisms containing the gene")
print("$L1 - Number of genes in each organism")

number_of_bar <- nrow(stats$value)

stats$value <- as.data.frame(cbind(stats$value, seq(1,number_of_bar)))

label_data <- stats$value

angle <-  90 - 360 * (label_data$V2-0.5) /number_of_bar

label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

p <- ggplot(stats$value, aes(x=as.factor(V2), y=V1)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("skyblue", 0.7)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-80,50) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=V2, y=V1+10, label=rownames(label_data), hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

p


##GENES For each ORGANISM
counts_BS <- read.csv("files/gene_counts_OS.txt", header = F)
#rnames <- counts_BS$V1
#counts_BS <- as.data.frame(counts_BS[,2])
#rownames(counts_BS) <- rnames
#colnames(counts_BS) <- "V1"
#counts_BS <- cbind(counts_BS, V3=rep("OS", nrow(counts_BS)))


counts_AS <- read.csv("files/gene_counts_AS.txt", header = F)
#rnames <- counts_AS$V1
#counts_AS <- as.data.frame(counts_AS[,2])
#rownames(counts_AS) <- rnames
#colnames(counts_AS) <- "V1"
#counts_AS <- cbind(counts_AS, V3=rep("AS", nrow(counts_AS)))

#counts <- full_join(counts_AS,counts_BS)

#ggplot(counts, aes(fill=V3, y=V2, x=V1)) + 
#  geom_bar(position="dodge", stat="identity")

counts <- merge(counts_AS, counts_BS, all = T, by="V1")

#counts <- counts[order(counts$V2.y, decreasing = T),]

counts <- cbind(counts,counts$V2.y-counts$V2.x)

colnames(counts) <- c("org","AS","OS","diff")

print(counts)

if(all(counts$diff==0)){
  print("Per organism, No transcripts were lost during transcript selection (AS) and ortholog selection (OS)!")
}else{
  print(counts[which(counts$diff!=0), 1:3])
  print("Per organism, These number of transcripts were lost during transcript selection (AS) and ortholog selection (OS)!")
}

#SEQUENCES for each GENE
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

