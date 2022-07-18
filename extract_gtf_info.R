suppressMessages(require(rtracklayer))
suppressMessages(require(GenomicRanges))
suppressMessages(require(Rgb))
suppressMessages(require(stringi))
suppressMessages(require(parallel))
suppressMessages(require(dataPreparation))
suppressMessages(require(dplyr))
suppressMessages(require(bedr))
#suppressMessages(require(future))

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
  return(unlist(lapply(strsplit(gtf_data$attributes, split="; "), function(x){
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

parallel_extract_info <- function(slice_file, strandedness="+") {
  if(!file.exists(slice_file) || file.info(slice_file)$size==0){
    print(paste("GTF Slice is empty, most likely gene not found",slice_file))
    return()
  }
  #print(slice_file)
  #gtf_con <- gzfile(slice_file,"r")
  gtf <- read.gtf(file = slice_file ,attr = c("intact")) #import.gff(slice_file,format="gtf") #import.gff("../mrna_loc/files/annos/xenopus_tropicalis.gtf",format="gtf")
  #gtf <- na.omit(gtf)
  #close(gtf_con)
  
  gtf$ranges <- IRanges(gtf$start, gtf$end)
  gtf$gene_id <- get_gtf_attributes(gtf,"gene_id")
  gtf$transcript_id <- get_gtf_attributes(gtf,"transcript_id")
  gtf$gene_name <- get_gtf_attributes(gtf,"gene_name")
  g_name <- unique(gtf$gene_name)
  
  #if(file.exists(output_file) && file.info(output_file)$size>0){
  #  print(paste("Output file exists(",output_file,")"))
  #  return()
  #}
  
  #print(paste(g_name, slice_file, output_file))
  
  gtf_gr <- GRanges(
    seqnames=Rle(gtf$seqname),
    ranges=IRanges(gtf$start, gtf$end),
    strand=Rle(droplevels(gtf$strand)),
    gene_name=gtf$gene_name,
    gene_id=gtf$gene_id,
    transcript_id=gtf$transcript_id,
    feature=gtf$feature,
    score=gtf$score
    #attribute=gtf$attributes
  )
  
  ##ONLY POSITIVE STRAND GENES
  #return(gtf_gr)
  
  gtf_gr <- gtf_gr[strand(gtf_gr)==strandedness,]
  
  if(is.null(gtf_gr) || length(gtf_gr)==0 || is.na(gtf_gr)){ ##Probably the mrna is not in the + strand so we can safely discard it
    print(paste("No + strand info (or) region of gtf missing for the gene : ", g_name))
    return()
  }
  
  ##tx <- subset(gtf, type == "mRNA")
  exon <- gtf_gr[gtf_gr$feature=="exon",] #subset(gtf, type == "exon")
  #exon <- exon[strand(exon)==strandedness,] #subset(exon, strand == "+")
  
  cds <- gtf_gr[gtf_gr$feature=="CDS",] #subset(gtf, type == "CDS")
  #cds <- cds[strand(cds)==strandedness,] #subset(cds, strand == "+")
  ##cds <- range(split(cds,cds$transcript_id))  #range(multisplit(cds, cds$Parent))
  
  if(length(cds)==0 || length(exon)==0 || is.na(cds) || is.na(exon)){ ##Probably the mrna is not in the + strand so we can safely discard it
    print(paste("Missing CDS/exon features for the gene : ", g_name))
    return()
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
  rna_boundry <- sapply(rna_list,end)-sapply(rna_list,start)
  names(rna_boundry) <- names(rna_list)
  ##rna_boundry <- sapply((sapply(seq_along(rna_list), function(i){
  ##  return(end(rna_list[i])-start(rna_list[i]))
  ##})), drop)
  #for (i in 1:length(rna_list)) {
  #  rna_boundry<- append(rna_boundry,(end(rna_list[i])-start(rna_list[i])))
  #}
  #rna_boundry <- unlist(rna_boundry)
  gene_boundry <- sapply(gene_list,end)-sapply(gene_list,start)
  names(gene_boundry) <- names(gene_list)
  ##gene_boundry <- sapply(sapply(seq_along(gene_list), function(i){
  ##  return(end(gene_list[i])-start(gene_list[i]))
  ##}), drop)
  #for (i in 1:length(gene_list)) {
  #  gene_boundry <- append(gene_boundry,(end(gene_list[i])-start(gene_list[i])))
  #}
  #gene_boundry <- unlist(gene_boundry)
  
  boundary_coverage <- c()
  utr_len <- c() #exon-cds to get length of UTR(3' + 5' UTR)
  utr_len <- bind_rows(lapply(seq_along(gene),function(i){
  #for (i in 1:length(gene)){
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
      #utr_len <- rbind(utr_len, data.frame(gene_name=unique(gene[[i]]$gene_name),gene_id=unique(gene[[i]]$gene_id),rna_id=unique(gene[[i]]$transcript_id),total_exon_len=sum(gene_boundry[unique(gene[[i]]$transcript_id)]),total_cds_len=sum(rna_boundry[unique(sort(gene[[i]]$transcript_id))]), five_len=five_len,three_len=three_len, exon_count=exon_cnt,cds_count=cds_cnt)) #search_group=search_group
      return(data.frame(gene_name=unique(gene[[i]]$gene_name),gene_id=unique(gene[[i]]$gene_id),transcript_id=unique(gene[[i]]$transcript_id),total_exon_len=sum(gene_boundry[unique(gene[[i]]$transcript_id)]),total_cds_len=sum(rna_boundry[unique(sort(gene[[i]]$transcript_id))]), five_len=five_len,three_len=three_len, exon_count=exon_cnt,cds_count=cds_cnt)) #search_group=search_group
      }
    #print(gene_boundry[unique(gene[[i]]$gene_id)]/sum(rna_boundry[unique(sort(gene[[i]]$transcript_id))]))
    #boundary_coverage<-append(boundary_coverage,list(data.frame(gene_coverage=gene_boundry[unique(gene[[i]]$gene_id)]/sum(rna_boundry[unique(sort(gene[[i]]$transcript_id))]),transcript_count=length(unique(sort(gene[[i]]$transcript_id))))))
  }))
  
  #OUTLIER removal, removing stuff from 3rd quartiles
  #pdf(file = NULL)
  #utr_len<-utr_len[which(is.na(match(utr_len$five_len,boxplot(utr_len$five_len,outline=F,plot=F)$out))),]
  #utr_len<-utr_len[which(is.na(match(utr_len$three_len,boxplot(utr_len$three_len,outline=F,plot=F)$out))),]
  #dev.off()
  #print("here")
  #print(utr_len)
  if(nrow(utr_len) > 1){
    utr_len <- as.data.frame(remove_sd_outlier(utr_len,cols=c("five_len","three_len"), verbose = F))
  }
  #print(head(gtf_gr))
  #bed_data <- unique(data.frame(seqnames=as.vector(seqnames(gtf_gr)),start=start(gtf_gr),end=end(gtf_gr)))#,feature=gtf_gr$feature,strand=gtf_gr$strand))
  ##CLEANUP
  #rm("gtf","gtf_gr","gene_boundry","rna_boundry","exon,cds","gene_list","rna_list","gene_exon","rna_exon","gene_cds","rna_cds","utr_len")
  #suppressMessages(gc(verbose = F,full = T))
  #names(utr_len) <- g_name
  #return(list(utr_len=utr_len,bed_data=bed_data))
  return(list(utr_len=utr_len,bed_data=as.data.frame(gtf_gr)))
  #write.table(utr_len, gzfile(output_file, open = "w"), sep = ",", quote = F,row.names = F, col.names = T)
  ##return(0)
  #write.table(utr_len, paste("files/genes/",org_name,"/",output_file,sep = ""), sep = ",", quote = F,row.names = F, col.names = !file.exists(paste("files/genes/",org_name,"/gtf_stats.csv",sep = "")), append = T)
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
}

##ENTRYPOINT

set.seed(123)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) GTF_Slice files path (2) Org_name (3) Output File (is in csv,NOT Gzipped)", call.=FALSE)
}

param_file <- "parameters.txt"

if(!file.exists(param_file) || file.info(param_file)$size < 0){
  stop("ERROR: parameters.txt is missing and is required")
}

#read.table(textConnection(gsub("==", "@", readLines(param_file))),sep="@")

param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^") #Convert multibyte seperator to one byte sep #read.table(param_file,sep="==")

max_concurrent_jobs <- as.numeric(param_table[which(param_table=="max_concurrent_jobs"),c(2)])
BED_PATH <- param_table[which(param_table=="bed_path"),c(2)]
TRANSCRIPT_ID_DELIM <- param_table[which(param_table=="transcript_delimiter"),c(2)]
STRAND <- param_table[which(param_table=="strand"),c(2)]
n_cores <- detectCores(all.tests = TRUE, logical = TRUE)
if(is.na(n_cores) || n_cores > max_concurrent_jobs){
  n_cores=max_concurrent_jobs
}

slice_file_path <- args[1] #"files/genes/xenopus_tropicalis/gdf.gtf_slice"
org_name <- args[2]
#gene_list_file <- args[3]
output_file <- args[3]

if (stri_isempty(slice_file_path) || stri_isempty(org_name) || stri_isempty(output_file)) { #stri_isempty(gene_list_file) ||
  stop("One/Many parameters are empty")
}

if(file.exists(output_file) && file.info(output_file)$size>0){
  print(paste("Output file exists(",output_file,")"))
  stop("Output file exists")
}

#mcmapply has some issues and R doesnt not work with GNU parallel (Rscript)
# https://stat.ethz.ch/pipermail/r-help/2021-March/470501.html
##CAREFUL with mcmapply return value (do not return(NULL) or return(0) etc, unless you are returning  data (return(data)) just do a plain return())

utr_len_list <- mcmapply(FUN=parallel_extract_info, list.files(path=slice_file_path, pattern = "*.gtf_slice",full.names = T), MoreArgs=list(strandedness=STRAND), USE.NAMES = F, mc.cores = n_cores-1, SIMPLIFY = F, mc.cleanup=F) # MoreArgs=list( output_path=output_path, org_name=org_name) MoreArgs=list(strandedness=STRAND, transcript_delim=TRANSCRIPT_ID_DELIM) mc.preschedule = T

utr_len_list[sapply(utr_len_list, is.null)] <- NULL

save("utr_len_list", file="utr_len_list.RData")

utr_len_df <- c()
#utr_len_df <- bind_rows(utr_len_list)
utr_len_df <- bind_rows(lapply(utr_len_list, function(x){
  return(x[["utr_len"]]) #return(x[[1]])
}))

if(nrow(utr_len_df) > 0 || length(utr_len_df) > 0 || !is.null(utr_len_df)){
  utr_len_df <- utr_len_df %>% mutate(org=org_name)
  write.table(utr_len_df, file(output_file, open = "w"), sep = ",", quote = F,row.names = F, col.names = T) #gzfile(output_file, open = "w")

  bed_df <- c()
  bed_df <- bind_rows(lapply(utr_len_list, function(x){
    #tmp_lens <- x[["utr_len"]]
    tmp_bed <- x[["bed_data"]] #return(x[[2]])
    #utr_len_df$transcript_id_ext <- paste(utr_len_df$transcript_id, utr_len_df$feature,sep = transcript_delim)
    #print(head(tmp_lens))
    #print(head(tmp_bed))
    tmp_bed$transcript_id[tolower(tmp_bed$feature)=="gene"] <- unique(tmp_bed$gene_name)
    #print(tmp_bed[,c("seqnames","start","end")])
    gtf_bed <- c()
    gtf_bed <- convert2bed(tmp_bed[,c("seqnames","start","end")], check.sort = F, check.valid=F, check.merge=F, check.chr = any(grepl("chr",tmp_bed[,c("seqnames")])))
    #gtf_bed <- convert2bed(org_gtf[,c("seqname","start","end")], check.sort = F, check.chr = any(grepl("chr",org_gtf[,c("seqname")])))
    #gtf_bed <- bedr.sort.region(gtf_bed, check.chr = any(grepl("chr",org_gtf[,c("seqname")])))
    
    gtf_bed$name <- paste(tmp_bed$transcript_id, tmp_bed$feature,sep = TRANSCRIPT_ID_DELIM) #utr_len_df$transcript_id_ext
    gtf_bed$feature <- tmp_bed$feature
    gtf_bed$strand <- tmp_bed$strand
    #gtf_bed <- gtf_bed %>% mutate(feature=tmp_lens$feature) %>% mutate(strand=tmp_lens$strand)
    return(gtf_bed)
  }))
  
  #gtf_bed <- bedr.sort.region(gtf_bed, check.chr = any(grepl("chr",org_gtf[,c("seqname")])))
  if(nrow(bed_df) > 0 || !is.null(bed_df)){
    bed_split <- split(bed_df,as.factor(bed_df$feature))
    BED_PATH_PREFIX <- file.path(BED_PATH,org_name)
    write.table(bed_split$CDS, file = paste(BED_PATH_PREFIX,"cds.bed",sep="_"),quote = F,sep = "\t", row.names = F, col.names = F)
    write.table(bed_split$five_prime_utr, file = paste(BED_PATH_PREFIX,"5utr.bed",sep="_"),quote = F,sep = "\t", row.names = F, col.names = F)
    write.table(bed_split$three_prime_utr, file = paste(BED_PATH_PREFIX,"3utr.bed",sep="_"),quote = F,sep = "\t", row.names = F, col.names = F)
    write.table(bed_split$exon, file = paste(BED_PATH_PREFIX,"exon.bed",sep="_"),quote = F,sep = "\t", row.names = F, col.names = F)
    write.table(bed_split$transcript, file = paste(BED_PATH_PREFIX,"transcript.bed",sep="_"),quote = F,sep = "\t", row.names = F, col.names = F)
  }
}