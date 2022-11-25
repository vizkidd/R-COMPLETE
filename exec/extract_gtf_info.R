suppressMessages(require(rtracklayer))
suppressMessages(require(GenomicRanges))
suppressMessages(require(Rgb))
suppressMessages(require(stringi))
suppressMessages(require(parallel))
#suppressMessages(require(dataPreparation))
suppressMessages(require(dplyr))
suppressMessages(require(bedr))
suppressMessages(require(purrr))


###FUNCTIONS

#Get Lengths of features like CDS/exons/utrs
feature_length <- function(transcript_id,gtf_data,feature="CDS"){
  gtf_sel <- gtf_data[which(gtf_data$transcript_id == transcript_id),]
  gtf_sel <- unique(gtf_sel[which(gtf_sel$feature == feature),])
  return(sum(gtf_sel[order(gtf_sel$start),]$end - gtf_sel[order(gtf_sel$start),]$start))
}

#Get counts of features, like number of CDS/exon Blocks per transcript
feature_count <- function(transcript_id,gtf_data,feature="CDS"){
  gtf_sel <- gtf_data[which(gtf_data$transcript_id == transcript_id),]
  gtf_sel <- unique(gtf_sel[which(gtf_sel$feature == feature),])
  return(length(gtf_sel))
}

#Split GTF attributes into induvidual columns
get_gtf_attributes <- function(gtf_data,attribute){
  return(unlist(lapply(strsplit(gtf_data$attributes, split="; "), function(x){
    attr_present <- grepl(x=x,pattern = attribute,ignore.case = T)
    if(any(attr_present)){
      #return(unlist(strsplit(x = x[which(attr_present)], split=" "))[2])
      return(sub(x=unlist(strsplit(x = x[which(attr_present)], split=" "))[2],replacement = "",pattern = "[[:punct:]]$"))
    }else{
      return("")
    }
  })))
}

#subset_gtf_attribute <- function(gtf_data,attribute,value){
#  return(gtf_data[grep(x=strsplit(unlist(lapply(strsplit(gtf$attributes, split="; "), function(x){x[grep(x=x,pattern=attribute)]})),split=" "),pattern = value,fixed = T),])
#}

extract_info <- function(gtf, strandedness="+") {

  gtf_split <- base::split(gtf, as.factor(gtf$gene_name))

  #save(file = "gtf_split.RData",list="gtf_split")

  return(mclapply(gtf_split, function(gtf_x){
    gtf_gr <- GRanges(
      seqnames=Rle(gtf_x$seqname),
      ranges=IRanges(gtf_x$start, gtf_x$end),
      strand=Rle(droplevels(gtf_x$strand)),
      gene_name=gtf_x$gene_name,
      gene_id=gtf_x$gene_id,
      transcript_id=gtf_x$transcript_id,
      feature=gtf_x$feature,
      score=gtf_x$score
      #attribute=gtf$attributes
    )

    g_name <- unique(gtf_x$gene_name)
    #print(g_name)

    if (strandedness=="+" || strandedness=="-") {
      gtf_gr <- gtf_gr[strand(gtf_gr)==strandedness,]
      if(is.null(gtf_gr) || length(gtf_gr)==0 || all(is.na(gtf_gr))){ ##Probably the mrna is not in the + strand so we can safely discard it
        message(paste("No",strandedness,"strand info (or) region of gtf missing for the gene : ", g_name))
        return()
      }
    }

    exon <- gtf_gr[gtf_gr$feature=="exon",]
    cds <- gtf_gr[gtf_gr$feature=="CDS",]

    if(length(cds)==0 || length(exon)==0 || is.na(cds) || is.na(exon)){ ##Probably the mrna is not in the strand we want so we can safely discard it
      message(paste("Missing CDS/exon features for the gene : ", g_name))
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

    gene_names <- names(gene)
    rna_names <- names(rna)

    gene_list <- range(gene)
    rna_list <- range(rna)
    rna_boundry <- c()
    gene_boundry <- c()

    rna_boundry <- sapply(rna_list,end)-sapply(rna_list,start)
    names(rna_boundry) <- names(rna_list)
    gene_boundry <- sapply(gene_list,end)-sapply(gene_list,start)
    names(gene_boundry) <- names(gene_list)

    boundary_coverage <- c()
    utr_len <- c() #exon-cds to get length of UTR(3' + 5' UTR)
    utr_len <- bind_rows(lapply(seq_along(gene),function(i){
      r_id <- unique(gene[[i]]$transcript_id)
      if(!is.na(rna_boundry[r_id]) && !is.na(gene_boundry[r_id])){
        #print(paste(unique(gene[[i]]$gene_id),unique(sort(gene[[i]]$transcript_id)),sep = "->"))
        #print(sum(gene_boundry[r_id])/sum(rna_boundry[r_id]))
        gtf_rna <- subset(gtf_gr, gtf_gr$transcript_id==r_id)

        exon_range <- attr(range(gene[[r_id]]),"ranges")
        rna_range <- attr(range(rna[[r_id]]),"ranges")
        five_len <- start(unique(rna_range))-unique(start(exon_range))
        three_len <- end(unique(exon_range))-unique(end(rna_range))
        #cds_len <- feature_length(gtf_data = gtf_rna,transcript_id = r_id, feature = "CDS")
        exon_cnt <- feature_count(gtf_data = gtf_rna,transcript_id = r_id, feature = "exon")
        cds_cnt <- feature_count(gtf_data = gtf_rna,transcript_id = r_id, feature = "CDS")

        if(five_len==0){ five_len <- sum(width(unique(gtf_rna[grep(gtf_rna$feature, pattern = "five")]))) }
        if(three_len <= 3){ three_len <- sum(width(unique(gtf_rna[grep(gtf_rna$feature, pattern = "three")]))) }

        return(bind_rows(list(gene_name=unique(gene[[i]]$gene_name),gene_id=unique(gene[[i]]$gene_id),transcript_id=unique(gene[[i]]$transcript_id),total_cds_len= sum(rna_boundry[unique(sort(gene[[i]]$transcript_id))]) , five_len=five_len, three_len=three_len, exon_count=exon_cnt,cds_count=cds_cnt )) ) #,g.exon_start=start(unique(exon_range)),g.exon_end=end(unique(exon_range)),total_exon_len=sum(gene_boundry[unique(gene[[i]]$transcript_id)])
      }
    }))

    #OUTLIER removal, removing stuff from 3rd quartiles
    #pdf(file = NULL)
    #utr_len<-utr_len[which(is.na(match(utr_len$five_len,boxplot(utr_len$five_len,outline=F,plot=F)$out))),]
    #utr_len<-utr_len[which(is.na(match(utr_len$three_len,boxplot(utr_len$three_len,outline=F,plot=F)$out))),]
    #dev.off()

    #if(nrow(utr_len) > 1){
    #  utr_len <- as.data.frame(remove_sd_outlier(utr_len,cols=c("five_len","three_len"), verbose = F))
    #}

    return(list(utr_len=utr_len,bed_data=as.data.frame(gtf_gr),gene_name=g_name))
  }, mc.cores = n_cores))
}

get_bed <- function(bed_df, BED_PATH_PREFIX, TRANSCRIPT_REGIONS) {
  if(nrow(bed_df) > 0 || !is.null(bed_df)){
    bed_split <- base::split(bed_df,as.factor(bed_df$feature))
    lapply(match(TRANSCRIPT_REGIONS,names(bed_split)), function(x){
      bed_region <- names(bed_split)[x]
      write.table(bed_split[[x]], file = paste(BED_PATH_PREFIX,paste(bed_region,"bed",sep="."),sep="_"),quote = F,sep = "\t", row.names = F, col.names = F)
    })
    write.table(bed_df, file = paste(BED_PATH_PREFIX,"bed",sep="."),quote = F,sep = "\t", row.names = F, col.names = T)
  }
}

get_stats_parallel <- function(slice_file_path, STRAND, n_cores, org_name, output_file, TRANSCRIPT_ID_DELIM, BED_PATH_PREFIX, TRANSCRIPT_REGIONS) {
  #mcmapply has some issues and R doesnt not work with GNU parallel (Rscript)
  # https://stat.ethz.ch/pipermail/r-help/2021-March/470501.html
  ##CAREFUL with mcmapply return value (do not return(NULL) or return(0) etc, unless you are returning  data (return(data)) just do a plain return())
  gtf_data <- mcmapply(FUN=function(slice_file){
    if(!file.exists(slice_file) || file.info(slice_file)$size==0){
      message(paste("GTF Slice is empty, most likely gene not found",slice_file))
      return()
    }

    gtf <- invisible(read.gtf(file = slice_file ,attr = c("intact"),quiet=T))

    #gtf$ranges <- IRanges(gtf$start, gtf$end)
    gtf$gene_id <- get_gtf_attributes(gtf,"gene_id")
    gtf$transcript_id <- get_gtf_attributes(gtf,"transcript_id")
    gtf$gene_name <- tolower(get_gtf_attributes(gtf,"gene_name"))
    #g_name <- unique(gtf$gene_name)
    return(gtf)
  }, list.files(path=slice_file_path, pattern = "*.gtf_slice",full.names = T), USE.NAMES = F, mc.cores = n_cores, SIMPLIFY = F)
  gtf_data <- bind_rows(gtf_data)

  if(nrow(gtf_data) == 0){
    stop(paste("GTF data is empty : ",org_name, " - Possibly no genes were found. (or gene_name attribute is missing from GTF)"))
  }

  utr_len_list <- extract_info(gtf_data,strandedness = STRAND)

  if(length(utr_len_list) == 0){
    return()
  }
  #save(file = "utr_len_list.RData",list="utr_len_list")

  #utr_len_list  <- lapply(utr_len_list, function(x){
  #  return(unlist(x,use.names = F,recursive = F))
  #})
  utr_len_list[sapply(utr_len_list, is.null)] <- NULL

  utr_len_df  <- bind_rows(mclapply(utr_len_list, function(x){
    #for each gene
    gene_stats <- as.data.frame(x[[1]])

    gene_stats[is.na(gene_stats)] <- 0

    ##get flank lengths
    gene_stats <- gene_stats %>% mutate(five_flank=ceiling(mean(unique(gene_stats$five_len))))
    gene_stats <- gene_stats %>% mutate(three_flank=ceiling(mean(unique(gene_stats$three_len))))

    #do variance correction
    gene_stats$five_flank <- gene_stats$five_flank + ceiling(mean(unique(abs(gene_stats$five_len-gene_stats$five_flank))))
    gene_stats$three_flank <- gene_stats$three_flank + ceiling(mean(unique(abs(gene_stats$three_len-gene_stats$three_flank))))

    gene_stats$transcript_length.estimated <- gene_stats$five_flank +  gene_stats$total_cds_len + gene_stats$three_flank
    gene_stats$transcript_length.annotated <- gene_stats$five_len + gene_stats$total_cds_len + gene_stats$three_len

    gene_stats$five_flank[gene_stats$five_flank==0] <- gene_stats$transcript_length.annotated
    gene_stats$three_flank[gene_stats$three_flank==0] <- gene_stats$transcript_length.annotated


    return(gene_stats)
  },  mc.cores = n_cores))

  if(nrow(utr_len_df) > 0 || length(utr_len_df) > 0 || !is.null(utr_len_df)){
    utr_len_df <- utr_len_df %>% mutate(org=org_name)
    utr_len_df$gene_name <- tolower(utr_len_df$gene_name)
    write.table(utr_len_df, file(output_file, open = "w"), sep = ",", quote = F,row.names = F, col.names = T,na = "-")

    bed_df <- bind_rows(mclapply(seq_along(utr_len_list), function(x){
      #print(x)
      tmp_bed <- utr_len_list[[x]][[2]]
      tmp_bed <- base::split(tmp_bed,as.factor(tmp_bed$gene_name))
      gtf_beds <- bind_rows(lapply(tmp_bed, function(x){
        if(any(!is.na(match(levels(factor(x$feature)),"gene")))){
          x$transcript_id[tolower(x$feature)=="gene"] <- unique(x$gene_name)
        }
        x$feature[tolower(x$feature)=="three_prime_utr"] <- "3utr"
        x$feature[tolower(x$feature)=="five_prime_utr"] <- "5utr"
        x$feature[tolower(x$feature)=="cds"] <- "cds"
        gtf_bed <- c()
        gtf_bed <- invisible(convert2bed(x[,c("seqnames","start","end")], check.sort = F, check.valid=F, check.merge=F, check.chr = any(grepl("chr",x[,c("seqnames")])), verbose = F))

        gtf_bed$name <- paste(x$transcript_id, x$feature,sep = TRANSCRIPT_ID_DELIM)
        gtf_bed$feature <- x$feature
        gtf_bed$strand <- x$strand
        gtf_bed <- gtf_bed %>% mutate(gene_name=unique(x$gene_name))
        #print(unique(x$gene_name))
        return(gtf_bed)
      }))

      return(gtf_beds)
    }, mc.cores = n_cores))

    get_bed(bed_df=bed_df, BED_PATH_PREFIX=BED_PATH_PREFIX, TRANSCRIPT_REGIONS=TRANSCRIPT_REGIONS)
  }
}

check_param <- function(param_table,param_id,optional=F,CAST_FUN=as.character,create_dir=F){
  if(is.character(param_table)){
    param_table <- read.table(textConnection(gsub("==", "^", readLines(param_table))),sep="^", header = T) #Convert multibyte seperator to one byte sep
  }
  param_index <- which(param_table==param_id)
  if(is.null(param_index) || is.na(param_index)){
    stop(paste("Parameter :",param_id,"is missing!"))
  }
  param_value <- param_table[param_index,c(2)]
  #print(CAST_FUN(param_value))
  if(create_dir){
    dir.create(CAST_FUN(param_value),showWarnings = F,recursive = T)
  }
  if(!stringi::stri_isempty(param_value) || optional){
    return(CAST_FUN(param_value))
  }else{
    stop(paste("Parameter :",param_id,"is empty and is not optional!"))
  }
}
##ENTRYPOINT

set.seed(123)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) GTF_Slice files path (2) Org_name (3) Output File (is in csv,NOT Gzipped) (4) Parameter file", call.=FALSE)
}

slice_file_path <- args[1]
org_name <- args[2]
output_file <- args[3]
param_file <- args[4] #"parameters.txt"

if(!file.exists(param_file) || file.info(param_file)$size <= 4){
  stop("ERROR: parameters.txt is missing and is required")
}

param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^", header = T) #Convert multibyte seperator to one byte sep #read.table(param_file,sep="==")

max_concurrent_jobs <<- check_param(param_table,"max_concurrent_jobs",optional=T,CAST_FUN=as.numeric)
if(is.na(max_concurrent_jobs) || length(max_concurrent_jobs) == 0){
  max_concurrent_jobs <- parallel::detectCores(all.tests = T, logical = T)
}
BED_PATH <<- tools::file_path_as_absolute(check_param(param_table,"bed_path",optional=F,CAST_FUN=as.character,create_dir=T))
TRANSCRIPT_ID_DELIM <<-  check_param(param_table,"transcript_delimiter",optional=F,CAST_FUN=as.character)
TRANSCRIPT_REGIONS <<- tolower(gsub("[[:space:]]","",x = unlist(stringi::stri_split( check_param(param_table,"transcript_regions",optional=F,CAST_FUN=as.character) ,fixed = ","))))
CLEAN_EXTRACT <<- check_param(param_table,"clean_extract",optional=F,CAST_FUN=as.logical)
STRAND <<- check_param(param_table,"strand",optional=F,CAST_FUN=as.character)
n_cores <<- detectCores(all.tests = TRUE, logical = TRUE)
if(is.na(n_cores) || n_cores > max_concurrent_jobs){
  n_cores <<- max_concurrent_jobs
}

dir.create(path = file.path(BED_PATH,org_name),recursive = T,showWarnings = F)
BED_PATH_PREFIX <- file.path(BED_PATH,org_name,org_name)

if (stri_isempty(slice_file_path) || stri_isempty(org_name) || stri_isempty(output_file)) {
  stop("One/Many parameters are empty")
}

if (CLEAN_EXTRACT==T) {
  unlink(x = output_file,force = T,expand = T)
  unlink(x = paste(BED_PATH_PREFIX,"bed",sep="."),force = T,expand = T)
}

utr_len_df <- c()
bed_df <- c()

#if(!file.exists(paste(BED_PATH_PREFIX,"bed",sep=".")) || file.info(paste(BED_PATH_PREFIX,"bed",sep="."))$size<4){ #file.exists(output_file) || file.info(output_file)$size<4
#  get_stats_parallel(slice_file_path=slice_file_path, STRAND=STRAND, n_cores=n_cores, org_name=org_name, output_file=output_file, TRANSCRIPT_ID_DELIM=TRANSCRIPT_ID_DELIM, BED_PATH_PREFIX=BED_PATH_PREFIX, TRANSCRIPT_REGIONS=TRANSCRIPT_REGIONS)
#}else{
tryCatch({
 #print(paste("Output file exists(",output_file,")...Checking BED files!"))

  utr_len_df <- read.table(file = file(output_file, open = "r"), sep = ",",header = T)

  bed_df <- read.table(file = file(paste(BED_PATH_PREFIX,"bed",sep="."), open = "r"), sep = "\t",header = T)
  bed_genes <- tolower(bed_df$gene_name) #sapply(stri_split(bed_df[bed_df$feature=="gene",c("name")], fixed = TRANSCRIPT_ID_DELIM), function(x){return(x[[1]])})
  if(!all(!is.na(match(unique(bed_genes),unique(utr_len_df$gene_name)))) || !all(!is.na(match(unique(utr_len_df$gene_name),unique(bed_genes))))){
    print("Some genes were missing, re-extracting")
    get_stats_parallel(slice_file_path=slice_file_path, STRAND=STRAND, n_cores=n_cores, org_name=org_name, output_file=output_file, TRANSCRIPT_ID_DELIM=TRANSCRIPT_ID_DELIM, BED_PATH_PREFIX=BED_PATH_PREFIX, TRANSCRIPT_REGIONS=TRANSCRIPT_REGIONS)
  }else{
    get_bed(bed_df=bed_df, BED_PATH_PREFIX=BED_PATH_PREFIX, TRANSCRIPT_REGIONS=TRANSCRIPT_REGIONS)
  }

  print("Files Checked!")
},error=function(cond){
  get_stats_parallel(slice_file_path=slice_file_path, STRAND=STRAND, n_cores=n_cores, org_name=org_name, output_file=output_file, TRANSCRIPT_ID_DELIM=TRANSCRIPT_ID_DELIM, BED_PATH_PREFIX=BED_PATH_PREFIX, TRANSCRIPT_REGIONS=TRANSCRIPT_REGIONS)
})

quit()



