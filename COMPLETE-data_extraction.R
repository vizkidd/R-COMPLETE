suppressMessages(require(biomartr))
suppressMessages(require(dplyr))
suppressMessages(require(biomaRt))
suppressMessages(require(parallel))
suppressMessages(require(tools)) ##FOR removing file extensions (gz)
suppressMessages(require(RCurl))
suppressMessages(require(curl))
suppressMessages(require(purrr))
suppressMessages(require(fs))
suppressMessages(require(processx))
suppressMessages(require(ps))
suppressMessages(require(stringi))
suppressMessages(require(data.table))
suppressMessages(require(stringr))
suppressMessages(require(seqinr))
suppressMessages(require(GenomicRanges))
suppressMessages(require(tictoc))

##Credits to https://nbisweden.github.io/workshop-RNAseq/2011/lab_download.html for GTF download through biomart with biomaRt

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) Gene List", call.=FALSE)
}

###FUNCTIONS
print_toc <- function(clk){
  return(toc.outmsg(clk$tic,clk$toc,clk$msg))
}

fetch_stat <- function(x,col_name){
  tmp_stat <- x[col_name][!is.na(x[col_name])]
  if (length(tmp_stat)==0) {
    tmp_stat <- 0
  }
  return(tmp_stat)
}

check_mart_dataset <- function(org){
  #mart_err <- paste("Dataset not found for :",org)
  #class(mart_err) <- c("try-error", class(mart_err))
  #split_org <- stri_split_fixed(org,"_",2,simplify = T)
  split_org <- stri_split(org,fixed="_",simplify = T)
  if(stri_isempty(split_org[2])){
    stop(paste(org,"name not in proper format, eg danio_rerio\n"))
  }
  mart.dataset <- grep(x = org.meta.list$dataset, pattern=regex(split_org[2],ignore_case = T),fixed=F, value = T)
  if(length(mart.dataset)==0){
    stop(paste("\nDataset not found in BIOMART for :",org,"\nYou can provide the organism in user data\n"))
  }else{
    return(mart.dataset)
  }
}

mart_connect <- function(MART_FUN=NULL,args=c(),verbose=T){
  if (!is.null(MART_FUN)) {
    #print(args)
    time_out <- 0
    while (time_out < 600) { #max time out is 10mins, code will fibonacci to it stepping up 5 seconds
      time_out <- time_out+5
      catch_value <- tryCatch(do.call(MART_FUN,args),
                              error=function(cond){
                                if(verbose){
                                  message(cond)
                                  #print(class(cond))
                                  message(paste("\nWaiting for",time_out,"s & trying again...\n"))
                                }
                                Sys.sleep(time_out)
                              })
      if(!any(grepl(pattern = "error|exception|try|fail|timeout",ignore.case = T,x = class(catch_value))) && !is.null(catch_value)){
        #print("Done!")
        return(catch_value)
      }
    }
    stop("\nEnsembl is unresponsive, please check connection\n")
  }else{
    return(0)
  }
}

extract_stats <- function(gtf_data) {
  
  #print(g_name)
  #print(paste(g_name, slice_file, output_file))
  strandedness=STRAND
  gtf_split <- base::split(gtf_data, as.factor(gtf_data$external_gene_name))
  return(mclapply(gtf_split, function(gtf_x){
    g_name <- unique(gtf_x$external_gene_name)
    print(g_name)
    gtf_gr <- unique(GRanges(
      seqnames=Rle(gtf_x$ensembl_transcript_id),
      ranges=IRanges(gtf_x$transcript_start, gtf_x$transcript_end),
      strand=Rle(droplevels(factor(gtf_x$strand))),
      gene_name=gtf_x$external_gene_name,
      gene_id=gtf_x$ensembl_gene_id,
      #transcript_id=gtf_x$ensembl_transcript_id,
      feature="transcript" #gtf_x$feature,
      #score=gtf_x$score
      #attribute=gtf$attributes
    ))
    
    if (strandedness=="+" || strandedness=="-") {
      gtf_gr <- gtf_gr[strand(gtf_gr)==strandedness,]  
      if(is.null(gtf_gr) || length(gtf_gr)==0 || all(is.na(gtf_gr))){ ##Probably the mrna is not in the + strand so we can safely discard it
        message(paste("No",strandedness,"strand info (or) region of gtf missing for the gene : ", g_name))
        return()
      }
    }
    
    transcript_stats <- lapply(base::split(gtf_x, as.factor(gtf_x$ensembl_transcript_id)), function(x){
      five_start <- fetch_stat(x,"5_utr_start")
      five_end <- fetch_stat(x,"5_utr_end")
      five_len <- sum( unique(five_end)-unique(five_start) ) #+ 1
      three_start <- fetch_stat(x,"3_utr_start")
      three_end <- fetch_stat(x,"3_utr_end")
      three_len <- sum( unique(three_end)-unique(three_start) ) #+ 1
      cds_start <- min(x["cds_start"])
      cds_end <- max(x["cds_end"])
      cds_len <- sum(x["cds_end"]-x["cds_start"]) # + 1
      exon_len <- sum(x["exon_chrom_end"]-x["exon_chrom_start"])
      
      
      stats_df <- bind_rows(list(transcript_id=as.character(unique(x["ensembl_transcript_id"])),cds_count=nrow(x),exon_count=length(unlist(unique(x["exon_chrom_end"]-x["exon_chrom_start"]))),transcript_length=as.numeric(unique(x["transcript_length"])),five_len=five_len,three_len=three_len, cds_start=cds_start,cds_end=cds_end,total_cds_len=cds_len,total_exon_len=exon_len,g.exon_start=unique(x["exon_chrom_start"]),g.exon_end=unique(x["exon_chrom_end"]),g.transcript_start=as.numeric(unique(x["transcript_start"])),g.transcript_end=as.numeric(unique(x["transcript_end"])),chromosome_name=as.character(unique(x["chromosome_name"])) ) ) #,g.three_start=three_start,g.three_end=three_end, g.five_start=five_start,g.five_end=five_end
      ##Remove rows with len==3, for the same transcript it is the coordinates for stop codons & start codons
      #stats_df <- stats_df[stats_df$three_len>3,]
      #stats_df <- stats_df[stats_df$five_len>3,] 
      #print(stats_df)
      return(stats_df)
    })
    transcript_stats <- bind_rows(transcript_stats)
    #transcript_stats <- transcript_stats[match(seqnames(gtf_gr)@values,transcript_stats$transcript_id),]
    gtf_gr <- invisible(full_join(as.data.frame(gtf_gr),transcript_stats, by = c("seqnames"="transcript_id") ))
    #gtf_gr$five_len[gtf_gr$five_len<=4] <- ceiling(mean(gtf_gr$five_len)) #correcting for Start codons
    #gtf_gr$three_len[gtf_gr$three_len<=4] <- ceiling(mean(gtf_gr$three_len)) #correcting for stop codons
    gtf_gr[is.na(gtf_gr)] <- 0
    ##gtf_gr <- remove_sd_outlier(as.data.frame(gtf_gr),cols=c("five_len","three_len"), verbose = F)
    gtf_gr <- gtf_gr %>% mutate(five_flank=ceiling(mean(unique(gtf_gr$five_len)))) #correcting for Start codons
    gtf_gr <- gtf_gr %>% mutate(three_flank=ceiling(mean(unique(gtf_gr$three_len)))) #correcting for stop codons
    #gtf_gr <- remove_sd_outlier(as.data.frame(gtf_gr),cols=c("five_len","three_len"), verbose = F)
    
    #do variance correction
    gtf_gr$five_flank <- gtf_gr$five_flank + ceiling(mean(unique(abs(gtf_gr$five_len-gtf_gr$five_flank))))
    gtf_gr$three_flank <- gtf_gr$three_flank + ceiling(mean(unique(abs(gtf_gr$three_len-gtf_gr$three_flank))))
    #gtf_gr$five_flank <- gtf_gr$five_flank + ceiling(mean(unique(gtf_gr$five_flank-gtf_gr$five_len)))
    #gtf_gr$three_flank <- gtf_gr$three_flank + ceiling(mean(unique(gtf_gr$three_flank-gtf_gr$three_len)))
    
    
    return(gtf_gr[,c("gene_name","gene_id","seqnames","total_exon_len","total_cds_len", "five_len","three_len", "exon_count","cds_count","g.exon_start","g.exon_end","transcript_length","five_flank","three_flank","g.transcript_start","g.transcript_end","chromosome_name","strand")])
  },mc.silent = T,mc.cores = numWorkers)) #,mc.cleanup = T
}

add_to_process <- function(p_cmd,p_args=list(),verbose=T){
  process_list[sapply(process_list, function(x){
    return(!x$is_alive())
  })] <<- NULL
  tryCatch(expr = function(){
    process_list[sapply(process_list, is.null)] <<- NULL
  }, error = function(){
    if(is.null(process_list)){
      process_list <<- c()
    }
  })
  if (verbose) {
    print(paste("Adding process to list...(",length(process_list),")",sep=""))  
  }
  if(length(process_list)>=numWorkers || ps::ps_num_fds() >= 250){
    #for (p_id in seq_along(process_list)) {
    #  if(process_list[[p_id]]$is_alive()){
    #    print(paste("Process Q Full...Waiting for a process to end(",length(process_list),")",sep=""))
    #    save(process_list, files="proces_list.RData")
    #    process_list[[p_id]]$wait(timeout=-1)
    #    process_list[[p_idx]] <<- NULL
    #    break;
    #  }
    if (verbose) {
      print(paste("Process Q Full...Waiting for a process to end(",length(process_list),")",sep=""))
    }
    #save(process_list, file="proces_list.RData")
    dead_procs <- c()
    for (x in seq_along(process_list)) {
      if(process_list[[x]]$is_alive()){
        process_list[[x]]$wait(timeout=-1)
        break;
      }else{
        dead_procs <- c(dead_procs,x)
      }
    }
    
    #lapply(seq_along(process_list), function(x){
    #  if(process_list[[x]]$is_alive()){
    #    process_list[[x]]$wait(timeout=-1)
    #    break;
    #  }else{
    #    dead_procs <- c(dead_procs,x)
    #  }
    #})
    
    if (length(dead_procs)>0) {
      dead_procs <- unique(dead_procs)
      process_list[[dead_procs]] <<- NULL 
    }
  }
  
  proc <- process$new(command=p_cmd,args = p_args, supervise = TRUE,stdout = "",stderr = "2>&1" ) #stderr = T, stdout =  T
  #proc$wait(timeout=-1)
  process_list <<- append(process_list,proc)
  
  return(proc)
}

check_files <-function(fasta_path,org,genes){
  # if(CLEAN_EXTRACT){
  #   print(paste(fasta_path,": Check FAILED! (CLEAN_EXTRACT : TRUE)"))
  #   return(FALSE)
  # }
  if(dir.exists(fasta_path) && length(dir(fasta_path,all.files = F)) > 0 && length(grep(list.files(path=fasta_path), pattern="tmp", invert=T, value=TRUE)) > 0){
    missing_genes <- c()
    available_genes <- c()
    if(file.exists(paste("files/genes/",org,"/","MISSING_GENES",sep=""))){
      missing_genes <- gsub('[[:punct:] ]+','_', factor(scan(paste("files/genes/",org,"/","MISSING_GENES",sep=""), character())))
    }
    if(file.exists(paste("files/genes/",org,"/","AVAILABLE_GENES",sep=""))){
      available_genes <- gsub('[[:punct:] ]+','_', factor(scan(paste("files/genes/",org,"/","AVAILABLE_GENES",sep=""), character())))
    }
    files_in_dir <- unique(sapply(list.files(fasta_path,no.. = T,recursive = F), FUN=function(x){stri_split(str = x, fixed='.',simplify = T)[1]})) 
    missing_genes <- missing_genes[is.na(match(missing_genes, available_genes))] 
    print(paste("Genes in Dir:",length(files_in_dir),", Missing:",length(missing_genes),", Available:",length(available_genes),", User Genes:",length(genes)))
    if((all(!is.na(match(files_in_dir,available_genes))) && all(!is.na(match(available_genes,files_in_dir)))) && all(is.na(match(missing_genes, files_in_dir)))){
      print(paste(fasta_path,": Check PASSED!"))
      return(TRUE)
    }
  }
  print(paste(fasta_path,": Check FAILED!"))
  return(FALSE)
}

get_gtf_mart <- function(org, genes){
  ##Credits to https://nbisweden.github.io/workshop-RNAseq/2011/lab_download.html for GTF download through biomart with biomaRt
  gtf_attributes <- c("chromosome_name",
                      "strand",
                      "ensembl_gene_id",
                      "ensembl_transcript_id",
                      "external_gene_name",
                      "start_position",
                      "end_position",
                      "exon_chrom_start",
                      "exon_chrom_end",
                      "transcript_version",
                      "transcript_length",
                      "5_utr_start",
                      "5_utr_end",
                      "transcript_start",
                      "transcription_start_site",
                      "transcript_end",
                      "cds_start",
                      "cds_end",
                      "3_utr_start",
                      "3_utr_end")
  mart.dataset <-tryCatch(check_mart_dataset(org),error=function(cond){
    stop(cond)
  })
  using.mart.data <- mart_connect(useMart,args = list(ENSEMBL_MART, mart.dataset))
  
  ##Splitting genes to the recomended number of queries for biomart through biomaRt
  #bm_gtf <- bind_rows(lapply(split(genes, ceiling(seq_along(genes)/500)), function(split_genes){
  #  return(mart_connect(getBM,args=list(mart=using.mart.data,attributes=gtf_attributes,uniqueRows=T, useCache=F, filters = c("external_gene_name"), values = split_genes, curl=curl_handle)))
  #}))
  bm_gtf <- mart_connect(getBM,args=list(mart=using.mart.data,attributes=gtf_attributes,uniqueRows=T, useCache=F, filters = c("external_gene_name"), values = genes, curl=curl_handle))
  bm_gtf <- dplyr::arrange(bm_gtf,chromosome_name,start_position)
  
  bm_gtf$strand[bm_gtf$strand==1] <- "+"
  #bm_gtf$strand[bm_gtf$strand== -1] <- "-"
  bm_gtf$strand[is.na(match(bm_gtf$strand,"+"))] <- "-"
  #bm_gtf$source_name <- "Ensembl"
  return(bm_gtf)
}

get_FASTA_mart <- function(org,gtf_stats, fasta_path){
  
  dir.create(path=fasta_path,showWarnings=F,recursive=T)
  seq_attributes <- c("5utr","coding","3utr") #"external_gene_name"
  mart.dataset <-tryCatch(check_mart_dataset(org),error=function(cond){
    stop(cond)
  })
  using.mart.data <- mart_connect(useMart,args = list(ENSEMBL_MART, mart.dataset))
  
  check_rows_idx <- unlist(sapply(seq_along(1:nrow(gtf_stats)), FUN=function(x){
    #print(gtf_stats[x,]<=3)
    if(any(na.omit(gtf_stats[x,c("five_len","three_len","five_flank","three_flank")]<=3))){ #any(is.na(gtf_stats[x,])) , not checking for CDS==NA because I was able to obtain coding sequences even when cds_len==NA
      return(x)
    }
  }))
  invalid_stats <- gtf_stats[check_rows_idx,]
  valid_transcripts <- unique(gtf_stats[-check_rows_idx,]$transcript_id)
  
  if(length(valid_transcripts) > 0){
    bm_seq <- mclapply(seq_attributes, function(x){return(mart_connect(getBM,args=list(mart=using.mart.data,attributes=c("ensembl_transcript_id",x),uniqueRows=T, useCache=F, filters = c("ensembl_transcript_id"), values = valid_transcripts, curl=curl_handle)))},mc.cores = numWorkers)
    
    bm_df <- bm_seq %>% purrr::reduce(full_join, by = "ensembl_transcript_id")
    unavailable_transcripts <- unique( unlist(apply(bm_df,MARGIN = 1, FUN = function(x){
      row_check <- grepl("unavailable",x=x,ignore.case = T)
      if(all(row_check) || row_check[2]==T){ #Removing transcripts which do not have any regions or CDS
        return(as.character(x["ensembl_transcript_id"]))
      }
    })) )
    if(length(unavailable_transcripts) > 0){
      bm_df <- bm_df[is.na(match(bm_df$ensembl_transcript_id,unavailable_transcripts)),]
      invalid_stats <- invalid_stats[is.na(match(invalid_stats$transcript_id,unavailable_transcripts)),]
    }
    
    bm_df <- bm_df %>% mutate(flanking_5utr=F)
    bm_df <- bm_df %>% mutate(flanking_3utr=F)
    
    #invalid_stats <- base::split(invalid_stats, as.factor(invalid_stats$gene_name))
    
    bm_cds <- mart_connect(getBM,args=list(mart=using.mart.data,attributes= c("ensembl_transcript_id","coding"),uniqueRows=T, useCache=F, filters =c("ensembl_transcript_id"), values = invalid_stats$transcript_id, curl=curl_handle))
    
    unavailable_transcripts <- unique( bm_cds$ensembl_transcript_id[grepl("unavailable",x=bm_cds$coding,ignore.case = T)] ) #Removing transcripts which do not have any regions or CDS
    if(length(unavailable_transcripts) > 0){
      invalid_stats <- invalid_stats[is.na(match(invalid_stats$transcript_id,unavailable_transcripts)),]
      bm_cds <- bm_cds[is.na(match(bm_cds$ensembl_transcript_id,unavailable_transcripts)),]
    }
    
    transcripts_with_5utr <- unique( na.omit(invalid_stats$transcript_id[invalid_stats[,grep("five",colnames(invalid_stats),ignore.case = T,value = F)] > 3]) )
    transcripts_with_3utr <-unique( na.omit(invalid_stats$transcript_id[invalid_stats[,grep("three",colnames(invalid_stats),ignore.case = T,value = F)] > 3]) )
    bm_5utr <- mart_connect(getBM,args=list(mart=using.mart.data,attributes= c("ensembl_transcript_id","5utr"),uniqueRows=T, useCache=F, filters =c("ensembl_transcript_id"), values = transcripts_with_5utr, curl=curl_handle))
    bm_3utr <- mart_connect(getBM,args=list(mart=using.mart.data,attributes= c("ensembl_transcript_id","3utr"),uniqueRows=T, useCache=F, filters =c("ensembl_transcript_id"), values = transcripts_with_3utr , curl=curl_handle))
    
    bm_5utr <- bm_5utr %>% mutate(flanking_5utr=F)
    bm_3utr <- bm_3utr %>% mutate(flanking_3utr=F)
    
    transcripts_without_utr <- unique( c( bm_5utr$ensembl_transcript_id[grepl("unavailable",x=bm_5utr[,c("5utr")],ignore.case = T)], bm_3utr$ensembl_transcript_id[grepl("unavailable",x=bm_3utr[,c("3utr")],ignore.case = T)] ) ) #We need flanks for unavilable utrs
    if(length(transcripts_without_utr) > 0){
      no_5utr <- intersect(transcripts_without_utr,transcripts_with_5utr)
      no_3utr <- intersect(transcripts_without_utr,transcripts_with_3utr)
      invalid_stats <- invalid_stats[!is.na(match(invalid_stats$transcript_id,unique( c(transcripts_without_utr, no_5utr, no_3utr) ) ) ),]
      #bm_5utr <- bm_5utr[is.na(match(bm_5utr$ensembl_transcript_id,transcripts_without_utr)),]
      #bm_3utr <- bm_3utr[is.na(match(bm_3utr$ensembl_transcript_id,transcripts_without_utr)),]
      bm_5utr <- bm_5utr[is.na(match(bm_5utr$ensembl_transcript_id,no_5utr)),]
      bm_3utr <- bm_3utr[is.na(match(bm_3utr$ensembl_transcript_id,no_3utr)),]
    }#else{
    #  invalid_stats <- invalid_stats[is.na(match(invalid_stats$transcript_id,unique( c(transcripts_with_5utr, transcripts_with_3utr) ) ) ),]
    #}
  }
  missing_any_flank_info <- unique(c(invalid_stats$transcript_id[invalid_stats$five_flank==0],invalid_stats$transcript_id[invalid_stats$three_flank==0])) #unique( invalid_stats$transcript_id[apply(X = (invalid_stats[,unique(c(grep("five",colnames(invalid_stats),ignore.case = T,value = F),grep("three",colnames(invalid_stats),ignore.case = T,value = F)))] ==0) , MARGIN = 1, FUN = any)] ) #missing_any_flank_info <- unique( invalid_stats$transcript_id[apply(X = (invalid_stats[,unique(c(grep("five",colnames(invalid_stats),ignore.case = T,value = F),grep("three",colnames(invalid_stats),ignore.case = T,value = F)))] ==0) , MARGIN = 1, FUN = all)] )
  
  flank_stats <- invalid_stats #invalid_stats[is.na(match(invalid_stats$transcript_id,missing_any_flank_info)),]  
  if (nrow(flank_stats) > 0) {
    if (length(missing_any_flank_info) > 0) { ##All these transcripts have neither UTR Lengths nor Flank lengths info (possibly because the gene does not have isoforms or info is missing), so I set the utr flanks to transcript lengths and correct for variance
      arbitrary_flanks <- flank_stats[!is.na(match(flank_stats$transcript_id,missing_any_flank_info)),]  
      arbitrary_flanks <- arbitrary_flanks[order(arbitrary_flanks$gene_name),]
      arb_flank_values <- bind_rows(lapply(base::split(arbitrary_flanks,as.factor(arbitrary_flanks$gene_name)), function(x){
        x$five_flank[ which(x$five_flank==0) ] <- mean(unique(arbitrary_flanks$transcript_length))
        x$three_flank[ which(x$three_flank==0) ] <- mean(unique(arbitrary_flanks$transcript_length))
        #Do variance correction
        return( data.frame(gene_name=x$gene_name,transcript_id=x$transcript_id,five_flank= x$five_flank + ceiling(mean(unique(abs(x$five_len-x$five_flank)))), three_flank=x$three_flank + ceiling(mean(unique(abs(x$three_len-x$three_flank)))) ) )
      }))
      arbitrary_flanks <- arbitrary_flanks %>% dplyr::select(-c("five_flank","three_flank")) %>% full_join(arb_flank_values, by = c("gene_name","transcript_id"))
      flank_stats <- flank_stats[is.na(match(flank_stats$transcript_id,missing_any_flank_info)),]
      flank_stats <- full_join(flank_stats, arbitrary_flanks)
    }
    #rm("invalid_stats") ##Cleaning up to save space, too many variables
    
    ##Calculate genomic coordinates for flanks
    flank_stats <- flank_stats %>% mutate(g.five_flank_start= abs(flank_stats$g.transcript_start-flank_stats$five_flank)-1 )
    flank_stats <- flank_stats %>% mutate(g.five_flank_end= abs(flank_stats$g.transcript_start-1) )
    flank_stats <- flank_stats %>% mutate(g.three_flank_start= abs(flank_stats$g.transcript_end+1) )
    flank_stats <- flank_stats %>% mutate(g.three_flank_end=  abs(flank_stats$g.transcript_end+flank_stats$three_flank)+1 )
    
    if (nrow(flank_stats) > 0) {
      ##Get Flanking sequences
      bm_flanks_five <- mart_connect(biomaRt::getSequence,args = list(id=unique(flank_stats$transcript_id),type="ensembl_transcript_id",seqType="coding_transcript_flank", upstream=flank_stats$five_flank , mart=using.mart.data, useCache = F) )#mart_connect(getBM,args=list(mart=using.mart.data,attributes=c("start_position","end_position","cdna"),uniqueRows=T, useCache=F, filters = c("chromosome_name","start","end","strand"), values = list(as.character(flank_stats$chromosome_name), flank_stats$g.five_flank_start, flank_stats$g.five_flank_end, STRAND), curl=curl_handle))
      bm_flanks_three <- mart_connect(biomaRt::getSequence,args = list(id=unique(flank_stats$transcript_id),type="ensembl_transcript_id",seqType="coding_transcript_flank", downstream = flank_stats$three_flank , mart=using.mart.data, useCache = F))
      
      names(bm_flanks_five)[grep(pattern="coding_transcript_flank",names(bm_flanks_five))] <- "5utr"
      names(bm_flanks_three)[grep(pattern="coding_transcript_flank",names(bm_flanks_three))] <- "3utr"
      
      bm_flanks_five <- bm_flanks_five %>% mutate(flanking_5utr=T)
      bm_flanks_three <- bm_flanks_three %>% mutate(flanking_3utr=T)
      
      bm_5utr_merged <- purrr::reduce(list(bm_5utr,bm_flanks_five), dplyr::full_join)
      bm_3utr_merged <- purrr::reduce(list(bm_3utr,bm_flanks_three), dplyr::full_join)
      
      bm_seqs <- purrr::reduce(list(bm_5utr_merged,bm_cds,bm_3utr_merged), dplyr::full_join)
    }else{
      bm_seqs <- purrr::reduce(list(bm_5utr,bm_cds,bm_3utr), dplyr::full_join)  
    }
  }else{
    bm_seqs <- purrr::reduce(list(bm_5utr,bm_cds,bm_3utr), dplyr::full_join)
  }
  
  if(nrow(bm_df)>0){
    bm_df <- purrr::reduce(list(bm_seqs,bm_df), dplyr::full_join)
  }else{
   stop(paste("Error fetching FASTA for :",org))
  }
  names(bm_df)[grep(pattern="ensembl_transcript_id",names(bm_df))] <- "transcript_id"
  names(bm_df)[grep(pattern="coding",names(bm_df))] <- "cds"
  
  bm_df <- inner_join(gtf_stats[,c("gene_name","safe_gene_name","transcript_id","strand")],bm_df)
  bm_df <- unique(bm_df)
  
  mclapply(base::split(bm_df[,c("transcript_id","gene_name","safe_gene_name",TRANSCRIPT_REGIONS,"flanking_5utr","flanking_3utr","strand")],as.factor(bm_df$gene_name)), function(x){
    lapply(TRANSCRIPT_REGIONS, function(y){
      if (grepl("3utr|5utr",y,ignore.case = T)) { ##Adding "_FLANK" for transcripts with UTR Flanks
        flank_col <- grep(paste("flanking",y,sep="_"),colnames(x),ignore.case = T,value = T)
        x_unique <- unique(x[,c("transcript_id","safe_gene_name",y,flank_col,"strand")])
        seq_names <- paste(as.character(x_unique[,"transcript_id"]), y,sep = TRANSCRIPT_ID_DELIM) 
        seq_names[x_unique[,flank_col]==T] <- paste(seq_names[x_unique[,flank_col]==T], "FLANK",sep = "_")
      }else{
        x_unique <- unique(x[,c(y,"transcript_id","safe_gene_name","strand")])
        seq_names <- paste(as.character(x_unique[,"transcript_id"]), y,sep = TRANSCRIPT_ID_DELIM)   
      }
      seq_names <- paste(seq_names,"(",x_unique[,"strand"],")",sep = "")   
      write.fasta(as.list(x_unique[,y]), seq_names,file.out=paste(fasta_path,"/",unique(x_unique[,"safe_gene_name"]),".",y,".tmp",sep="") )
    })
  },mc.cores = numWorkers, mc.silent = T)
  
  
  
  return(gtf_stats[!is.na(match(gtf_stats$transcript_id,bm_df$transcript_id)),])
}

fetch_fasta <- function(org_row) {
  
  org <- org_row["name"]
  tic(msg=paste("Processed:",org))
  org_fasta_path <- file.path(OUT_PATH ,org) 
  
  if(!CLEAN_EXTRACT && check_files(org_fasta_path,org,genes)){
    print(print_toc(toc(quiet = T)))
    return(org_row)
  }
  
  odb_list <- paste("files/genes/",org,"/odb.list",sep = "")
  odb_gene_map <- paste("files/genes/",org,"/odb.final_map",sep = "")
  gtf_stats_file <- paste("files/genes/",org,"/gtf_stats.csv",sep = "")
  
  if(CLEAN_EXTRACT){
    unlink(org_fasta_path,recursive = T,force=T,expand = T)
    unlink(odb_list,force=T,expand = T)
    unlink(odb_gene_map,force=T,expand = T)
    unlink(gtf_stats_file,force=T,expand = T)
  }
  
  gtf_data <- get_gtf_mart(org, genes)
  
  if(!file.exists(odb_list) || file.info(odb_list)$size == 0){
    proc <- do.call(add_to_process,list(p_cmd = c("./check_OrthoDB.sh"),p_args = c(org,gene_list, odb_list, odb_gene_map))) #time
    proc$wait(timeout=-1)
  }
  if(file.exists(odb_list) && file.info(odb_list)$size > 0){
    odb_list_genes <- factor(scan(odb_list, character())) 
    odb_list_genes <- odb_list_genes[grep("gene",tolower(odb_list_genes), invert = T, fixed = T)]
    odb_list_data <- get_gtf_mart(org, odb_list_genes)
    gtf_data <- unique(merge(odb_list_data,gtf_data))
  }else{
    message(paste("ODB gene list could not be found for : ",org))
  }
  
  gtf_stats <- extract_stats(gtf_data)
  gtf_stats <- bind_rows(gtf_stats)
  names(gtf_stats)[grep(pattern="seqnames",names(gtf_stats))] <- "transcript_id"
  gtf_stats$safe_gene_name <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(gtf_stats$gene_name))
  gtf_stats <- gtf_stats[gtf_stats$gene_name!=0 | gtf_stats$gene_id!=0,] # cleaning up
  
  
  if(nrow(gtf_stats)==0){
    message(paste("No genes passed filters for : ",org))
    print(print_toc(toc(quiet = T)))
    return()
  }
  
  gtf_stats <- suppressMessages(get_FASTA_mart(org,gtf_stats,org_fasta_path))
  
  if(nrow(gtf_stats)==0){
    message(paste("Error fetching FASTA for : ",org))
    print(print_toc(toc(quiet = T)))
    return()
  }
  
  local_procs <- c()
  local_procs <- unlist(mclapply(base::split(gtf_stats,as.factor(gtf_stats$gene_name)), function(x){
    gene <- unique(x[,"gene_name"])
    safe_name <- unique(x[,"safe_gene_name"])
    return(unlist(lapply(TRANSCRIPT_REGIONS, function(y){
      output_fasta <- paste(org_fasta_path,"/",safe_name,".",y,sep="")
      input_fasta <- paste(org_fasta_path,"/",safe_name,".",y,".tmp",sep="")
      proc <- do.call(add_to_process,list(p_cmd = c("./label_sequenceIDs.sh"),p_args = c(org,gene, output_fasta,input_fasta, odb_gene_map), verbose=F))
      return(proc)
    })))
  }, mc.cores = numWorkers,mc.silent = T))
  invisible(mclapply(local_procs, function(x){
    #print(x)
    if(x$is_alive()){
      x$wait(timeout=-1)
    }
  }, mc.cores = numWorkers))
  
  names(gtf_stats)[grep(pattern="transcript_length",names(gtf_stats))] <- "transcript_length.annotated"
  gtf_stats <- gtf_stats %>% mutate(org=org)
  gtf_stats$transcript_length.estimated <- gtf_stats$five_flank + gtf_stats$total_cds_len + gtf_stats$three_flank
  gtf_stats <- gtf_stats[,c("gene_name","gene_id","transcript_id","total_cds_len","five_len","three_len","exon_count","cds_count","five_flank","three_flank","transcript_length.estimated","transcript_length.annotated","org")]
  
  write.table(x = gtf_stats,file = gtf_stats_file,quote = F,sep = ",",row.names = F,col.names = T)
  write.table(x= unique(tolower(gtf_stats$gene_name)),file=paste(file.path("files","genes",org) ,"/","AVAILABLE_GENES",sep=""),quote=F,row.names=F,col.names=F )
  write.table(x= unique(genes[is.na(match( tolower(genes),tolower(unique(gtf_stats$gene_name)) ))]) ,file=paste(file.path("files","genes",org),"/","MISSING_GENES",sep=""),quote=F,row.names=F,col.names=F )
  
  print(paste("DONE :", org))
  print(print_toc(toc(quiet = T)))
  return(org_row)
}

fetch_genome_user <- function(data){
  genome_path <- c()
  gtf_path <- c()
  #print(data)
  genome <- data["genome"]
  gtf <- data["gtf"]
  org <- data["org"]
  
  if(genome == "-" || gtf == "-"){
    return( tryCatch(fetch_fasta(org),error=function(cond){
      message(cond)
      return(NULL)
    }) ) 
    #break;
  }
  
  tic(msg = paste("Processed:",org))
  genome_path<-paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")
  gtf_path<-paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")  
  org_fasta_path <- file.path(OUT_PATH ,org)
  
  if(grepl("://|http|ftp|www",genome)){
    if(!file.exists(genome_path) || !file.info(genome_path)$size > 20 || CLEAN_EXTRACT){
      curl_fetch_disk(genome, paste(GENOMES_PATH, "/",basename(URLdecode(genome)),sep=""))
      genome_path <- paste(GENOMES_PATH, "/",basename(URLdecode(genome)),sep="")
    } # else{
    #   if(system2("gzip",args = c("-t", genome_path), wait = T,stdout = NULL, stderr = NULL) == 0){
    #     genome_path<-genome_path 
    #   }else{
    #     lapply(list.files(GENOMES_PATH, full.names = T, ignore.case = T, no.. = T, pattern = regex(stri_split_fixed(org,"_",2,simplify = T),ignore_case = T)),function(x){
    #       if(grepl(x,pattern = "fa|gz")){
    #         try(file.remove(x, showWarnings=F))
    #       }
    #     })
    #     curl_fetch_disk(genome, basename(URLdecode(genome)))
    #     genome_path <- paste(GENOMES_PATH, "/",basename(URLdecode(genome)),sep="") #basename(URLdecode(genome))
    #   }
    # }
    
  }else{
    if(file.exists(genome) && file.info(genome)$size > 0){
      genome_path <- genome
    }else{
      message(paste("User Genome not found :",org,"-",genome))
      print(print_toc(toc(quiet = T)))
      return(NULL)
    }
  }
  
  if(grepl("://|http|ftp|www",gtf)){
    if(!file.exists(gtf_path) || !file.info(gtf_path)$size > 20 || CLEAN_EXTRACT){
      curl_fetch_disk(gtf, paste(ANNOS_PATH, "/",basename(URLdecode(gtf)),sep=""))
      gtf_path <- paste(ANNOS_PATH, "/",basename(URLdecode(gtf)),sep="")
    }# else{
    #   if(system2("gzip",args = c("-t", gtf_path), wait = T,stdout = NULL, stderr = NULL) == 0){
    #     gtf_path<-gtf_path
    #   }else{
    #     lapply(list.files(ANNOS_PATH, full.names = T, ignore.case = T, no.. = T, pattern = regex(stri_split_fixed(org,"_",2,simplify = T),ignore_case = T)),function(x){
    #       if(grepl(x,pattern = "gtf|gz")){
    #         try(file.remove(x, showWarnings=F))
    #       }
    #     })
    #     curl_fetch_disk(gtf, paste(ANNOS_PATH, "/",basename(URLdecode(gtf)),sep=""))
    #     gtf_path <- paste(ANNOS_PATH, "/",basename(URLdecode(gtf)),sep="")
    #   }
    # }
    
  }else{
    if(file.exists(gtf) && file.info(gtf)$size > 0){
      gtf_path <- gtf
    }else{
      message(paste("User GTF not found :",org,"-",gtf))
      print(print_toc(toc(quiet = T)))
      return(NULL)
    }
  }
  
  # if(!is.logical(genome_path) && file_ext(genome_path) == "gz"){
  #   try(file_move(genome_path, paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")))
  #   genome_path <- paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")
  # }else if(!is.logical(genome_path)){
  #   g_ext <- file_ext(genome_path)
  #   try(file_move(genome_path, paste(GENOMES_PATH, "/",org,".",g_ext,sep = "")))
  #   genome_path <- paste(GENOMES_PATH, "/",org,".",g_ext,sep = "")
  # }
  # if(!is.logical(gtf_path) && file_ext(gtf_path) == "gz"){
  #   try(file_move(gtf_path, paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")))
  #   gtf_path <- paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")
  # }else if(!is.logical(gtf_path)){
  #   a_ext <- file_ext(gtf_path)
  #   try(file_move(gtf_path, paste(ANNOS_PATH, "/",org,".",a_ext,sep = "")))
  #   gtf_path <- paste(ANNOS_PATH, "/",org,".",a_ext,sep = "")
  # }
  
  #print(paste(genome_path,gtf_path,org))
  
  if(CLEAN_EXTRACT || !try(check_files(org_fasta_path,org,genes))){
    if (!is.logical(gtf_path) && !is.logical(genome_path)) {
      do.call(add_to_process,list(p_cmd = c("./jobhold.sh"), p_args = c(paste("extract",org,sep="_"), "./extract_genomic_regions.sh",genome_path, gtf_path, gene_list, org)))
      #return(do.call(add_to_process,list(p_cmd = c("./extract_genomic_regions.sh"), p_args = c(genome_path, gtf_path, gene_list, org))))
      print(print_toc(toc(quiet = T)))
      return(data)      
    }else{
      return(NULL)
    }
  }else{
    print(print_toc(toc(quiet = T)))
    return(data)      
  }
}

####ENTRY POINT
set.seed(123)

options(RCurlOptions = list(ssl.verifyhost=0, ssl.verifypeer=0,timeout=200,maxconnects=200,connecttimeout=200))
#connection_options <<- curlOptions(ssl.verifyhost=0, ssl.verifypeer=0,timeout=200,maxconnects=200,connecttimeout=200))

process_list <<- c()

param_file <- "parameters.txt"

if(!file.exists(param_file) || file.info(param_file)$size < 0){
  stop("ERROR: parameters.txt is missing and is required")
}

param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^", header = T) #Convert multibyte seperator to one byte sep

GENOMES_PATH <<- param_table[which(param_table=="genomes_path"),c(2)]
ANNOS_PATH <<- param_table[which(param_table=="annos_path"),c(2)]
TEMP_PATH <<- param_table[which(param_table=="temp_path"),c(2)]
GROUPS_PATH <<- param_table[which(param_table=="groups_path"),c(2)]
OUT_PATH <<- param_table[which(param_table=="fasta_path"),c(2)]
BLAST_REGION <<- tolower(as.character(param_table[which(param_table=="blast_region"),c(2)]))
USER_GENOMES <<- as.character(param_table[which(param_table=="user_genomes"),c(2)])
CLEAN_EXTRACT <<- as.logical(param_table[which(param_table=="clean_extract"),c(2)])
TRANSCRIPT_ID_DELIM <<- param_table[which(param_table=="transcript_delimiter"),c(2)]
DATA_SOURCE <<- tolower(param_table[which(param_table=="data_source"),c(2)])
TRANSCRIPT_REGIONS <<- tolower(gsub("[[:space:]]","",x = unlist(stri_split(param_table[which(param_table=="transcript_regions"),c(2)],fixed = ","))))
STRAND <<- param_table[which(param_table=="strand"),c(2)]
max_concurrent_jobs <<- as.numeric(param_table[which(param_table=="max_concurrent_jobs"),c(2)])

curl_handle <<- getCurlMultiHandle()

gene_list <<- args[1]
genes <<- factor(scan(gene_list, character())) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
genes <<- genes[grep("gene",tolower(genes), invert = T, fixed = T)]

print(paste("MAX PROCESSES:",max_concurrent_jobs))

tryCatch(numWorkers <<- detectCores(all.tests = T, logical = T), error=function(){numWorkers <<- 2})
if(is.na(numWorkers) || numWorkers > max_concurrent_jobs){
  numWorkers <<- max_concurrent_jobs
}

dir.create("files",showWarnings = F, recursive = T)
unlink(TEMP_PATH, recursive = T,force = T,expand = T)
dir.create(TEMP_PATH,showWarnings = F, recursive = T)
unlink(GROUPS_PATH, recursive = T,force = T,expand = T)
dir.create(GROUPS_PATH,showWarnings = F, recursive = T)
dir.create(OUT_PATH,showWarnings = F, recursive = T)

ENSEMBL_MART <<- "ENSEMBL_MART_ENSEMBL"
using.mart <- mart_connect(useMart,args=list(ENSEMBL_MART)) #For biomaRt
org.meta.list <<- mart_connect(listDatasets,args=list(mart=using.mart)) #For biomaRt
org.meta <- mart_connect(listGenomes,args=list(db = "ensembl", type = "all", details = T)) #For biomartr #db = tolower(GENOMES_SOURCE)

if(!stri_isempty(USER_GENOMES) && !is.null(USER_GENOMES)){
  user_data <- read.csv(USER_GENOMES,header = F)
  names(user_data) <- c("org","genome","gtf")
}

print("Transforming ODB Files...")

ODB_proc <- do.call(add_to_process,list(p_cmd = c("./merge_OG2genes_OrthoDB.sh"), p_args = c(gene_list)))
ODB_proc$wait(timeout=-1)

print("Checking and downloading transcripts...")

saved_meta <- c()

if(DATA_SOURCE=="both" || DATA_SOURCE=="user"){
  if(nrow(user_data)!=0 && !is.null(user_data)){
    saved_meta <- apply(user_data, MARGIN = 1, function(x){
      #user_proc <- fetch_genome_user(x)
      #user_proc$wait()
      tic.clear(); 
      return( tryCatch(fetch_genome_user(x),error=function(cond){
        message(cond)
        return(NULL)
      })
      ) } )
  }
}

save("saved_meta", file="saved_meta.RData")

if(DATA_SOURCE=="both" || DATA_SOURCE!="user"){
  #mclapply(org.meta$name, fetch_genome_db, mc.cores = 1)
  #lapply(org.meta$name, FUN = function(x){
  saved_meta <- list( saved_meta, apply(org.meta, MARGIN = 1, FUN = function(x){
    tic.clear(); 
    return( tryCatch(fetch_fasta(x),error=function(cond){
      message(cond)
      return(NULL)
    })
    ) } ) )
}

save("saved_meta", file="saved_meta.RData")

# print(process_list)
# for(proc in process_list){
#   if(proc$is_alive()){
#     print(proc$get_status())
#     if(ps_is_running(proc$as_ps_handle())){
#       ## if SIGCHLD is overwritten the process is lost, so we try to get pid and check if the process is running
#       proc$wait(timeout = -1) 
#     }else{
#       print(proc$print())
#       print(proc$get_status())
#       proc$interrupt()
#     }
#   }
#   #print(proc$get_exit_status())
# }

mclapply(process_list, function(x){
  if(x$is_alive()){
    x$wait(timeout=-1)
  }
}, mc.cores =  numWorkers)

saved_meta[sapply(saved_meta, is.null)] <- NULL
saved_meta <- bind_rows(saved_meta)

write.table(x = saved_meta,file = "files/org_meta.txt", quote = F,sep=",", row.names = F)

lapply(list.files(path = OUT_PATH,include.dirs=TRUE, full.names=TRUE), function(x) {
  fi <- file.info(x)
  if (fi$isdir) {
    f <- list.files(x, all.files=TRUE, recursive=TRUE, full.names=TRUE)
    sz <- sum(file.info(f)$size)
    
    #as precaution, print to make sure before using unlink(x, TRUE)
    if (sz==0L) print(x)   
  }
})

#find $FASTA_PATH -name MISSING_GENES| awk -F'/' '{print $(NF-1)","$0}' > files/MISSING_GENES_FINAL
#find $FASTA_PATH -name AVAILABLE_GENES| awk -F'/' '{print $(NF-1)","$0}' > files/AVAILABLE_GENES_FINAL
#find $FASTA_PATH/* -type d |  sort | uniq | awk -F'/' '{print $NF}' > files/available_orgs.txt
#for f_org in $FASTA_PATH/*; do 
#f_org_name=$(basename $f_org)
#parallel --max-procs $n_threads " printf '%s\t%s\n' {1} {2}" :::: <(grep -H -f files/genes/$f_org_name/ALL_CLUSTERS -r $FASTA_PATH/$f_org_name/ | awk -F'[:>]' -v s_delim="$seqID_delimiter" '{split($2,a,s_delim); n=split($1,b,"."); print $1"\t"$2"\t"a[5]"\t"b[n]'}) | parallel  --max-procs 1 --colsep '\t' --recend '\n'  "if [[ -s {1} && ! -z {2} && ! -z {1} && ! -z {3} ]] ; then samtools faidx {1}  {2} >> $GROUPS_PATH/{3}.{4} ; fi" 
#done
# ls $FASTA_PATH > files/selected_ORGS.txt
# #Coerce gtf_stats.csv of all organisms
# find files/genes -iname "gtf_stats.csv" -exec sed 1d {} \; > files/gtf_stats.csv
# ##ID Alignments - GENERATE numeric ids for FASTA IDS (because they are long and downstream analysis have difficulty taking long names) 
# #Only indexing the IDs for now because find_orthologs.sh depends on the long FASTA IDs and cannot be shorted until orthologous transcripts are obtained
# #!!!!!!!#CHANGE FASTA IDs to numeric IDs (because some programs dont work well with long FASTA IDs) ONLY BEFORE ALIGNMENT!!!!!!!!
# index_fastaIDs files/rna_ids.txt $FASTA_PATH
# time ./find_orthologs.sh files/selected_ORGS.txt $1 #100 ##This also selects the transcripts
# time ./align_seqs.sh $1
# time ./predict_structures.sh $1
# rm $TEMP_PATH/*
# Rscript gene_stats.R >> files/stats.txt

sessionInfo()