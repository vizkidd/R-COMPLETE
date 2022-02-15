suppressMessages(require(rtracklayer))
suppressMessages(require(GenomicRanges))
suppressMessages(require(Rgb))
suppressMessages(require(stringi))

###CODE BASED ON IM's INPUT
feature_length <- function(transcript_id,gtf_data,feature="CDS"){
  ###gtf_data SHOULD BE READ USING Rgb::read.gtf
  ###With the followiing options, 
  ###attr = c("split"),features = c("CDS")
  #gtf_sel <- gtf_data[which(gtf_data$transcript_id == transcript_id & gtf_data$feature == feature),]
  #print(transcript_id)
  gtf_sel <- gtf_data[which(stri_detect(fixed = transcript_id,str = gtf_data$attributes)),]
  gtf_sel <- gtf_sel[which(gtf_sel$feature == feature & gtf_sel$frame==0),]
  #print(gtf_sel[order(gtf_sel$start),])
  
  return(sum(gtf_sel[order(gtf_sel$start),]$end - gtf_sel[order(gtf_sel$start),]$start))
}

feature_count <- function(transcript_id,gtf_data,feature="CDS"){
  ###gtf_data SHOULD BE READ USING Rgb::read.gtf
  ###With the followiing options, 
  ###attr = c("split"),features = c("CDS")
  
  #gtf_sel <- gtf_data[which(gtf_data$transcript_id == transcript_id & gtf_data$feature == feature),]
  gtf_sel <- gtf_data[which(stri_detect(fixed = transcript_id,str = gtf_data$attribute)),]
  gtf_sel <- gtf_sel[which(gtf_sel$feature == feature),]
  #print(transcript_id)
  #print(gtf_sel)
  return(length(gtf_sel))
  #return(max(gtf_sel$exon_number))
}

get_gtf_attributes <- function(gtf_data,attribute){
  #attribute_mat <- c()
  #matching_rows<-grep(x=gtf_data$attributes,pattern=attribute)
  return(unlist(lapply(strsplit(gtf$attributes, split="; "), function(x){
    attr_present <- grepl(x=x,pattern = attribute,ignore.case = T)
    if(any(attr_present)){
      return(unlist(strsplit(x = x[which(attr_present)], split=" "))[2])
    }else{
      return("")
    }
  })))
  #unlist(lapply(strsplit(unlist(lapply(strsplit(gtf$attributes, split="; "), function(x){x[grep(x=x,pattern=attribute)]})),split=" "), function(x){return(x[2])}))
  
  #return(gtf_data[which(stri_detect(fixed = attribute,str = gtf_data$attributes)),])
}

subset_gtf_attribute <- function(gtf_data,attribute,value){
  return(gtf_data[grep(x=strsplit(unlist(lapply(strsplit(gtf$attributes, split="; "), function(x){x[grep(x=x,pattern=attribute)]})),split=" "),pattern = value,fixed = T),])
}

set.seed(123)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) GTF (*.gtf_slice) (2) Org_name (3) Search group of gene(gene name used to grep)", call.=FALSE)
}

slice_file <- args[1] #"files/genes/xenopus_tropicalis/gdf.gtf_slice"
org_name <- args[2]
search_group <- args[3]
if(file.info(slice_file)$size==0){
  stop("GTF Slice is empty, most likely gene not found")
}

gtf <- read.gtf(file = slice_file ,attr = c("intact")) #import.gff(slice_file,format="gtf") #import.gff("../mrna_loc/files/annos/xenopus_tropicalis.gtf",format="gtf")
gtf$ranges <- IRanges(gtf$start, gtf$end)
gtf$gene_id <- get_gtf_attributes(gtf,"gene_id")
gtf$transcript_id <- get_gtf_attributes(gtf,"transcript_id")
gtf$gene_name <- get_gtf_attributes(gtf,"gene_name")

gtf_gr <- GRanges(
  seqnames=Rle(gtf$seqname),
  ranges=IRanges(gtf$start, gtf$end),
  strand=Rle(droplevels(gtf$strand)),
  gene_name=gtf$gene_name,
  gene_id=gtf$gene_id,
  transcript_id=gtf$transcript_id,
  feature=gtf$feature,
  score=gtf$score,
  attribute=gtf$attributes
)

#tx <- subset(gtf, type == "mRNA")
exon <- gtf_gr[gtf_gr$feature=="exon",] #subset(gtf, type == "exon")
exon <- exon[strand(exon)=="+",] #subset(exon, strand == "+")

cds <- gtf_gr[gtf_gr$feature=="CDS",] #subset(gtf, type == "CDS")
cds <- cds[strand(cds)=="+",] #subset(cds, strand == "+")
#cds <- range(split(cds,cds$transcript_id))  #range(multisplit(cds, cds$Parent))

if(length(cds)==0 || length(exon)==0 || is.na(cds) || is.na(exon)){ ##Probably the mrna is not in the + strand so we can safely discard it
  stop(paste("No + strand info (or) region of gtf missing in the gene : ", search_group))
}

gene_cds <- multisplit(cds,as.list(cds$gene_id))
rna_cds <- multisplit(cds,as.list(cds$transcript_id))

gene_exon <- multisplit(exon,as.list(exon$gene_id))
rna_exon <- multisplit(exon,as.list(exon$transcript_id))

gene <- c()
rna <- c()

gene <- rna_exon
rna <- rna_cds

##GET only rna info, omiited because we use gtf slices instead of full files
#gene <- gene[which(grepl(pattern="gene",names(gene),ignore.case = T)==T)]
#rna <- rna[which(grepl(pattern="rna",names(rna),ignore.case = T)==T)]

gene_names <- names(gene)
rna_names <- names(rna)

gene_list <- range(gene)
rna_list <- range(rna)
rna_boundry <- c()
gene_boundry <- c()

#for (i in 1:length(rna)) {
#  seqlevels(rna[[i]]) <- sort(seqlevels(rna[[i]]))
#  rna[[i]] <- unique(sort(rna[[i]], by=~start+end))
#  for(j in 2:length(rna[[i]])){
#    prev_range <- c()
#    curr_range <- c()
#    prev_range <- attr(rna[[i]][j-1],"range")
#    curr_range <- attr(rna[[i]][j],"range")
#    print(names(rna)[[i]])
#    print(start(curr_range)-end(prev_range))
#  }
#}

for (i in 1:length(rna_list)) {
  rna_boundry<- append(rna_boundry,(end(rna_list[i])-start(rna_list[i])))
}
rna_boundry <- unlist(rna_boundry)

for (i in 1:length(gene_list)) {
  gene_boundry<- append(gene_boundry,(end(gene_list[i])-start(gene_list[i])))
}
gene_boundry <- unlist(gene_boundry)

boundary_coverage <- c()
utr_len <- c() #exon-cds to get length of UTR(3' + 5' UTR)
for (i in 1:length(gene)){
  r_id <- unique(gene[[i]]$transcript_id)
  if(!is.na(rna_boundry[r_id]) && !is.na(gene_boundry[r_id])){
    #print(paste(unique(gene[[i]]$gene_id),unique(sort(gene[[i]]$transcript_id)),sep = "->"))
    #print(sum(gene_boundry[r_id])/sum(rna_boundry[r_id]))
    gtf_rna <- subset(gtf_gr, gtf_gr$transcript_id==r_id)
    
    exon_range <- attr(range(gene[[r_id]]),"ranges")
    rna_range <- attr(range(rna[[r_id]]),"ranges")
    five_len <- start(rna_range)-start(exon_range)
    three_len <- end(exon_range)-end(rna_range)
    exon_cnt <- feature_count(gtf_data = gtf_rna,transcript_id = r_id, feature = "exon")
    cds_cnt <- feature_count(gtf_data = gtf_rna,transcript_id = r_id, feature = "CDS")
    
    if(five_len==0){ five_len <- sum(width(gtf_rna[grep(gtf_rna$feature, pattern = "five")])) }
    if(three_len <= 3){ three_len <- sum(width(gtf_rna[grep(gtf_rna$feature, pattern = "three")])) }
    #  if(sum(gene_boundry[r_id])/sum(rna_boundry[r_id]) == 1){ ## if  exon/sum(CDS) == 1, then UTR annotation isnt available in exons, we search in the GTF for UTR info
#    #print(r_id)
#    
#    #print(paste(exon_range,rna_range,five_len,three_len))
#  #  if(five_len!=0 && three_len > 3){ ##3 because stop codons can be annotated along with CDS
#    if(five_len==0){ five_len <- sum(width(gtf_rna[grep(gtf_rna$feature, pattern = "five")])) }
#    if(three_len <= 3){ three_len <- sum(width(gtf_rna[grep(gtf_rna$feature, pattern = "three")])) }
#  #  }
#  }
    utr_len <- rbind(utr_len, data.frame(org=org_name,search_group=search_group,gene_name=unique(gene[[i]]$gene_name),gene_id=unique(gene[[i]]$gene_id),rna_id=unique(gene[[i]]$transcript_id),total_exon_len=sum(gene_boundry[unique(gene[[i]]$transcript_id)]),total_cds_len=sum(rna_boundry[unique(sort(gene[[i]]$transcript_id))]), five_len=five_len,three_len=three_len, exon_count=exon_cnt,cds_count=cds_cnt))
  }
  #print(gene_boundry[unique(gene[[i]]$gene_id)]/sum(rna_boundry[unique(sort(gene[[i]]$transcript_id))]))
  #boundary_coverage<-append(boundary_coverage,list(data.frame(gene_coverage=gene_boundry[unique(gene[[i]]$gene_id)]/sum(rna_boundry[unique(sort(gene[[i]]$transcript_id))]),transcript_count=length(unique(sort(gene[[i]]$transcript_id))))))
}

#OUTLIER removal, removing stuff from 3rd quartiles
pdf(file = NULL)
utr_len<-utr_len[which(is.na(match(utr_len$five_len,boxplot(utr_len$five_len,outline=F,plot=F)$out))),]
utr_len<-utr_len[which(is.na(match(utr_len$three_len,boxplot(utr_len$three_len,outline=F,plot=F)$out))),]
dev.off()

write.table(utr_len, paste("files/genes/",org_name,"/gtf_stats.csv",sep = ""), sep = ",", quote = F,row.names = F, col.names = !file.exists(paste("files/genes/",org_name,"/gtf_stats.csv",sep = "")), append = T)
#write.csv(utr_len,file = paste("files/genes/",org_name,"/predicted_utr_lens.csv",sep = ""),quote = F,row.names = F,append = T)

#seqlevels(gene_list) <- sort(seqlevels(gene_list))
#seqlevels(rna_list) <- sort(seqlevels(rna_list))
#gene_list <- unique(sort(gene_list, by=~start+end))
#rna_list <- unique(sort(rna_list, by=~start+end))

#inter_gene <- c()
#inter_rna <- c()

#for (i in 2:length(gene_list)) {
#  prev_gene <- attr(gene_list[i-1],"range")
#  curr_gene <- attr(gene_list[i],"range")
#  #print(i)
#  print(names(gene_list)[i])
#  print(start(curr_gene)-end(prev_gene))
#}

#for (i in 2:length(rna_list)) {
#  prev_rna <- attr(rna_list[i-1],"range")
#  curr_rna <- attr(rna_list[i],"range")
#  #print(i)
#  print(names(rna_list)[i])
#  print(start(curr_rna)-end(prev_rna))
#}

#utrs <- psetdiff(tx, cds[tx$ID])
