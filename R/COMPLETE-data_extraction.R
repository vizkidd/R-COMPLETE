#suppressMessages(require(biomartr))
#suppressMessages(require(dplyr))
#suppressMessages(require(biomaRt))
#suppressMessages(require(parallel))
#suppressMessages(require(tools)) ##FOR removing file extensions (gz)
#suppressMessages(require(RCurl))
#suppressMessages(require(curl))
#suppressMessages(require(purrr))
#suppressMessages(require(fs))
#suppressMessages(require(processx))
#suppressMessages(require(ps))
#suppressMessages(require(stringi))
#suppressMessages(require(data.table))
#suppressMessages(require(stringr))
#suppressMessages(require(seqinr))
#suppressMessages(require(GenomicRanges))
#suppressMessages(require(tictoc))

##Credits to https://nbisweden.github.io/workshop-RNAseq/2011/lab_download.html for GTF download through biomart with biomaRt

#args = commandArgs(trailingOnly=TRUE)

#if (length(args)==0) {
#  stop("Give the (1) Gene List\n", call.=FALSE)
#}

###FUNCTIONS

#' Pretty Print tictoc::toc() output
#'
#' This function takes the output of tictoc::toc() and
#' pretty prints it to the console
#'
#' @param clk output of function tictoc::toc()
#' @return A formatted string to print to stdout/console
print_toc <- function(clk){
  return(paste(tictoc::toc.outmsg(clk$tic,clk$toc,clk$msg),"\n\n", sep = ""))
}

# toc_log_msg <- function(tic, toc, msg, info){
#   if (is.null(msg) || is.na(msg) || length(msg) == 0)
#   {
#     outmsg <- paste0(round(toc - tic, 3), " seconds elapsed")
#   }
#   else
#   {
#     outmsg <- paste0(info, ": ", msg, ": ",
#                      round(toc - tic, 3), " seconds elapsed")
#   }
#   outmsg
# }

#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param x Data Frame with columns named as required
#' @param col_name Name of the column for which the value is to be extracted
#' @return A value which was requested from the GTF data frame
fetch_stat <- function(x,col_name){
  tmp_stat <- x[col_name][!is.na(x[col_name])]
  if (length(tmp_stat)==0) {
    tmp_stat <- 0
  }
  return(tmp_stat)
}

#' Check if a dataset exists for the organism in BIOMART (through biomaRt)
#'
#' This function accepts a name of an organism and
#' returns the name of the biomaRt dataset if available
#'
#' @param org Name of the organism
#' @return Name of the dataset in biomaRt
#' @export
check_mart_dataset <- function(org){
  # #mart_err <- paste("Dataset not found for :",org)
  # #class(mart_err) <- c("try-error", class(mart_err))
  # #split_org <- stringi::stri_split_fixed(org,"_",2,simplify = T)
  # if(is.null(COMPLETE_env$org.meta.list)){
  #   #COMPLETE_env <- new.env(parent=emptyenv())
  #   COMPLETE_env$ENSEMBL_MART <- "ENSEMBL_MART_ENSEMBL"
  #   COMPLETE_env$using.mart <- mart_connect(biomaRt::useMart,args=list(COMPLETE_env$ENSEMBL_MART)) #For biomaRt
  #   COMPLETE_env$org.meta.list <- mart_connect(biomaRt::listDatasets,args=list(mart=COMPLETE_env$using.mart))
  # }

  if(is.null(COMPLETE_env$org.meta.list)){
    stop("Reload R-COMPLETE, error in initialization")
  }

  split_org <- stringi::stri_split(gsub(pattern = "[[:punct:]]|[[:space:]]",x = org, replacement = "_"),fixed="_",simplify = T)
  if(stringi::stri_isempty(split_org[,2])){
    stop(paste(org,"name not in proper format, eg danio_rerio\n"))
  }
  mart.dataset <- grep(x = COMPLETE_env$org.meta.list$dataset, pattern=stringr::regex(split_org[2],ignore_case = T),fixed=F, value = T)
  if(length(mart.dataset)==0){
    stop(paste("Dataset not found in BIOMART for :",org,"\nYou can provide the organism in user data\n\n"))
    #return(NULL)
  }else if (length(unique(mart.dataset))==1) {
    return(unique(mart.dataset))
  }else{
    f_name <- paste(tolower(substring(split_org[1],1,1)),split_org[2],sep="")
    mart.dataset <- grep(x = mart.dataset,pattern=paste("^",f_name,sep="") ,ignore.case = T,value = T) #grep(x = mart.dataset,pattern= ,ignore.case = T,value = T)
    if (length(unique(mart.dataset))==1) {
      return(unique(mart.dataset))
    }else{
      stop(paste("\nMultiple datasets found in BIOMART for :",org,":",paste(mart.dataset,collapse = ","),"\n"))
    }
  }
}

#' Calculate Flank Values and other statistics from GTF Data
#'
#' This function calculates Flank Values (for transcripts whose UTR lengths are not known in which case Flank values are used instead of lengths to obtain the Flanking regions),
#'  CDS lengths, Number of Exon & CDS Blocks from GTF Data obtained from biomaRt. Variance correction is performed for FLANKS in case the variance is huge (or 0) for the UTR lengths .This function also checks for strandedness and only keeps the genes which are present in a strand (if requested by the user).
#' The GTF Data should have these columns c("chromosome_name", "strand", "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name",
#' "start_position","end_position","exon_chrom_start","exon_chrom_end","transcript_version","transcript_length","5_utr_start","5_utr_end",
#' "transcript_start","transcription_start_site","transcript_end","cds_start","cds_end","3_utr_start","3_utr_end")
#'
#' @note For calculating stats of actual GTF file, refer to exec/extract_gtf_info.R in the github repo. Variance correction is performed per-gene (seperately for 3' and 5' UTRs) which is, UTR FLANK LENGTH + abs(variance(transcript UTR lengths))
#'
#' @examples
#'     gtf_data <- get_gtf_mart(org, genes)
#'     gtf_stats <- calculate_stats(gtf_data)
#'     gtf_stats <- dplyr::bind_rows(gtf_stats)
#'     names(gtf_stats)[grep(pattern="seqnames",names(gtf_stats))] <- "transcript_id"
#'
#' @param gtf_data GTF data obtained from biomaRt for an organism
#' @param allow_strand Only allow the specified strand. ("+","-",Default - "" (or) " " (or) "*")
#' @param n_threads Number of Threads
#' @return Transcript Statistics from GTF data
#' @export
calculate_stats <- function(gtf_data, allow_strand="", n_threads=tryCatch(parallel::detectCores(all.tests = T, logical = T), error=function(cond){return(2)})) {

  #print(g_name)
  #print(paste(g_name, slice_file, output_file))
  strandedness=allow_strand

  req_columns <- c("external_gene_name",
                   "ensembl_gene_id",
                   "ensembl_transcript_id",
                   "strand",
                   "transcript_start",
                   "transcript_end")

  if (any(is.na(match(req_columns, names(gtf_data))))) {
    stop(paste("Missing columns, Require :",paste(req_columns,collapse = ",")))
  }

  gtf_split <- base::split(gtf_data, as.factor(gtf_data$external_gene_name))
  return(parallel::mclapply(gtf_split, function(gtf_x){
    g_name <- unique(gtf_x$external_gene_name)
    #print(g_name)
    gtf_gr <- unique(GenomicRanges::GRanges(
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
      if(is.null(gtf_gr) || length(gtf_gr)==0 || all(is.na(gtf_gr))){ ##Probably the mrna is not in the required strand so we can safely discard it
        #message(paste("No",strandedness,"strand info (or) region of gtf missing for the gene : ", g_name))
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


      stats_df <- dplyr::bind_rows(list(transcript_id=as.character(unique(x["ensembl_transcript_id"])),cds_count=nrow(x),exon_count=length(unlist(unique(x["exon_chrom_end"]-x["exon_chrom_start"]))),transcript_length=as.numeric(unique(x["transcript_length"])),five_len=five_len,three_len=three_len, cds_start=cds_start,cds_end=cds_end,total_cds_len=cds_len,total_exon_len=exon_len,g.exon_start=unique(x["exon_chrom_start"]),g.exon_end=unique(x["exon_chrom_end"]),g.transcript_start=as.numeric(unique(x["transcript_start"])),g.transcript_end=as.numeric(unique(x["transcript_end"])),chromosome_name=as.character(unique(x["chromosome_name"])) ) ) #,g.three_start=three_start,g.three_end=three_end, g.five_start=five_start,g.five_end=five_end
      ##Remove rows with len==3, for the same transcript it is the coordinates for stop codons & start codons
      #stats_df <- stats_df[stats_df$three_len>3,]
      #stats_df <- stats_df[stats_df$five_len>3,]
      #print(stats_df)
      return(stats_df)
    })
    transcript_stats <- dplyr::bind_rows(transcript_stats)
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
  },mc.silent = T,mc.cores = n_threads)) #,mc.cleanup = T
}

#' Internal Function - Add process to a list of process
#'
#' This function adds process to a list and makes sure hte number of processes do not
#' exceed a user defined limit (< max_concurrent_jobs/cores/threads). Uses processx library to call
#' a child process. Process is appended to COMPLETE_env$process_list
#'
#' @param p_cmd Command to be executed
#' @param p_args Arguments to be passed to the command
#' @param verbose Print status messages?
#' @param logfile Redirect output (stdout & stderr) to this file
#' @param params_list Output of load_params()
#' @return Process ID from processx::new() which can be used for further monitoring of the process
add_to_process <- function(p_cmd,p_args=list(),verbose=F, logfile=NULL, params_list){
  if (is.null(logfile)) {
    logfile=""
  }
  if (is.null(COMPLETE_env$process_list)) {
    COMPLETE_env$process_list <- c()
  }

  if(length(COMPLETE_env$process_list) > 0){
    COMPLETE_env$process_list[sapply(COMPLETE_env$process_list, function(x){
      return(!x$is_alive())
    })] <- NULL
    tryCatch(expr = function(){
      COMPLETE_env$process_list[sapply(COMPLETE_env$process_list, is.null)] <- NULL
    }, error = function(){
      if(is.null(COMPLETE_env$process_list)){
        COMPLETE_env$process_list <- c()
      }
    })
  }
  if (verbose) {
    cat(paste("\nAdding process to list...(",length(COMPLETE_env$process_list),")\n",sep=""))
  }
  if(length(COMPLETE_env$process_list)>=params_list$numWorkers || ps::ps_num_fds() >= COMPLETE_env$max_file_handles-1){ #250
    #for (p_id in seq_along(COMPLETE_env$process_list)) {
    #  if(COMPLETE_env$process_list[[p_id]]$is_alive()){
    #    print(paste("Process Q Full...Waiting for a process to end(",length(COMPLETE_env$process_list),")",sep=""))
    #    save(COMPLETE_env$process_list, files="proces_list.RData")
    #    COMPLETE_env$process_list[[p_id]]$wait(timeout=-1)
    #    COMPLETE_env$process_list[[p_idx]] <<- NULL
    #    break;
    #  }
    if (verbose) {
      cat(paste("Process Q Full...Waiting for a process to end(",length(COMPLETE_env$process_list),")\n",sep=""))
      lapply(COMPLETE_env$process_list, function(x){ print(paste(paste(x$get_cmdline(), collapse = " "),sep="")) }) #DEBUG
    }
    #save(COMPLETE_env$process_list, file="proces_list.RData")
    dead_procs <- c()
    # for (x in seq_along(COMPLETE_env$process_list)) {
    #   if(COMPLETE_env$process_list[[x]]$is_alive()){
    #     COMPLETE_env$process_list[[x]]$wait(timeout=-1)
    #     break;
    #   }else{
    #     dead_procs <- c(dead_procs,x)
    #   }
    # }

    dead_procs <- unlist(lapply(COMPLETE_env$process_list, function(x){
      if(!x$is_alive()){
        return(x)
      }
    }))

    #lapply(seq_along(COMPLETE_env$process_list), function(x){
    #  if(COMPLETE_env$process_list[[x]]$is_alive()){
    #    COMPLETE_env$process_list[[x]]$wait(timeout=-1)
    #    break;
    #  }else{
    #    dead_procs <- c(dead_procs,x)
    #  }
    #})

    if (length(dead_procs)>0 && length(COMPLETE_env$process_list)>0){
      dead_procs <- unique(dead_procs)
      COMPLETE_env$process_list[[dead_procs]] <- NULL
    }
    if (length(COMPLETE_env$process_list)>=params_list$numWorkers || ps::ps_num_fds() >= COMPLETE_env$max_file_handles-1){
      COMPLETE_env$process_list[[1]]$wait(timeout=-1)
    }
    # if(COMPLETE_env$process_list[[1]]$is_alive()){ ## check and wait for the earliest job to complete
    #   COMPLETE_env$process_list[[1]]$wait(timeout=-1)
    # }
  }

  proc <- processx::process$new(command=p_cmd,args = p_args, cleanup = T, cleanup_tree = T,supervise = TRUE,stdout = logfile,stderr = "2>&1" ) #stderr = T, stdout =  T
  #proc$wait(timeout=-1)
  COMPLETE_env$process_list <- append(COMPLETE_env$process_list,proc)

  return(proc)
}

#' Checks FASTA files in a folder
#'
#' This function checks if all the genes (and the transcript regions (cds,3utr,5utr)) are present
#' in the given folder path. The FASTA filenames should be of the format {gene}.{transcript region} for this
#' function to work.
#'
#' @examples
#'   check_files("files/fasta/notechis_scutatus",org = "notechis_scutatus",genes = "data/genelist.txt",verbose =T, params_list = load_params("pkg_data/parameters.txt"))
#'
#' @param fasta_path Path of folder to check
#' @param org Name of the organism (format important, eg. "danio_rerio")
#' @param genes Vector of genes to check for
#' @param verbose Print check result messages?
#' @param params_list Output of load_params()
#' @return TRUE if check passed, FALSE if check failed
check_files <-function(fasta_path,org,genes, verbose=T, params_list){
  # if(params_list$CLEAN_EXTRACT){
  #   print(paste(fasta_path,": Check FAILED! (params_list$CLEAN_EXTRACT : TRUE)"))
  #   return(FALSE)
  # }

  #if (!is.vector(genes)) {
  if(length(genes) == 1 && file.exists(genes)){
    genes <- gsub('[[:punct:]]+','_', factor(scan(genes, character(), quiet = T)))#factor(scan(genes, character(),quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
    #}#else{
    #  stop(paste(genes,"does not exist!"))
    #}
  }else{
    genes <- tolower(as.vector(genes))
  }

  return(tryCatch({
    if(dir.exists(fasta_path) && length(dir(fasta_path,all.files = F)) > 0){ # && length(grep(list.files(path=fasta_path), pattern="tmp", invert=T, value=TRUE)) > 0
      missing_genes <- c()
      available_genes <- c()
      # final_genes <- c()
      # if(file.exists(paste(params_list$OUT_PATH,"/genes/",org,"/","final.list",sep=""))){
      #    final_genes <- tolower(gsub('[[:punct:]]+','_', factor(scan(paste(paste(params_list$OUT_PATH,"/genes/",org,"/","final.list",sep="")), character(), quiet = T))))
      # }
      odb_genes <- c()
      if(file.exists(paste(params_list$OUT_PATH,"/genes/",org,"/","odb.list",sep=""))){
        odb_genes <- tolower(gsub('[[:punct:]]+','_', factor(scan(paste(paste(params_list$OUT_PATH,"/genes/",org,"/","odb.list",sep="")), character(), quiet = T))))
      }
      if(file.exists(paste(params_list$OUT_PATH,"/genes/",org,"/","MISSING_GENES",sep=""))){
        missing_genes <- tolower(gsub('[[:punct:]]+','_', factor(scan(paste(params_list$OUT_PATH,"/genes/",org,"/","MISSING_GENES",sep=""), character(), quiet = T))))
      }
      if(file.exists(paste(params_list$OUT_PATH,"/genes/",org,"/","AVAILABLE_GENES",sep=""))){
        available_genes <- tolower(gsub('[[:punct:]]+','_', factor(scan(paste(params_list$OUT_PATH,"/genes/",org,"/","AVAILABLE_GENES",sep=""), character(), quiet = T))))
      }
      files_in_dir <- tolower(unique(gsub('[[:punct:]]+','_',sapply(list.files(fasta_path,no.. = T,recursive = F), FUN=function(x){stringi::stri_split(str = x, fixed='.',simplify = T)[,1]}))))
      missing_genes <- tolower(missing_genes[is.na(match(missing_genes, available_genes))])

      #final_genes <- c(genes,odb_genes)

      #if(verbose){
      # cat(paste("Org:",org,", Genes in Dir:",length(files_in_dir),", Available:",length(available_genes),", Missing:",length(missing_genes),", User Genes:",length(genes),"\n"))
      # message(paste("Org:",org,", Genes in Dir:",length(files_in_dir),", Available:",length(available_genes),", Missing:",length(missing_genes),", User Genes:",length(genes)))
      # }

      # if(length(available_genes)==0 && length(missing_genes)==0){
      #   if(verbose){
      #     message(paste("Org:",org,", Genes in Dir:",length(files_in_dir),", Available:",length(available_genes),", Missing:",length(missing_genes),", User Genes:",length(genes)," : Check FAILED!"))
      #   }
      #   return(FALSE)
      # }

      if(all(!is.na(match(intersect(available_genes,genes),files_in_dir))) && all(is.na(match(intersect(available_genes,genes),missing_genes))) ){ #all(!is.na(match(final_genes,files_in_dir))) && #&& all(is.na(match(missing_genes, files_in_dir)))  #all(which(!is.na(match(files_in_dir,available_genes)))) && all(which(!is.na(match(available_genes,files_in_dir))))) #all(!is.na(match(files_in_dir[which(!is.na(match(files_in_dir,available_genes)))],available_genes[which(!is.na(match(available_genes,files_in_dir)))]))) #all(!is.na(match(missing_genes, files_in_dir))))
        if(verbose){
          cat(paste("Org:",org,", Genes in Dir(Matching + ODB):","(",length(files_in_dir[!is.na(match(files_in_dir,genes))]),"+",length(files_in_dir[!is.na(match(files_in_dir,setdiff(odb_genes,genes)))]),")",", Genes in Dir:",length(files_in_dir),", Unavailable:",length(genes[is.na(match(genes,files_in_dir))]),", User Genes:",length(genes)," : Check PASSED!\n")) #", Unavailable:",length(files_in_dir[is.na(match(files_in_dir,c(genes,odb_genes)))])
        }

        # data.table::fwrite(x= list(unique(files_in_dir[match(files_in_dir,final_genes)])) ,file=paste(params_list$OUT_PATH,"/genes/",org ,"/AVAILABLE_GENES",sep=""),quote=F,row.names=F,col.names=F,na = "-", nThread = params_list$numWorkers )
        # data.table::fwrite(x= list(unique(files_in_dir[match(files_in_dir,final_genes)])) ,file=paste(params_list$OUT_PATH,"/genes/",org ,"/final.list",sep=""),quote=F,row.names=F,col.names=F,na = "-", nThread = params_list$numWorkers )
        # data.table::fwrite(x= list(tolower(unique(final_genes[is.na(match( tolower(final_genes),files_in_dir ))]))) ,file=paste(params_list$OUT_PATH,"/genes/",org ,"/MISSING_GENES",sep=""),quote=F,row.names=F,col.names=F ,na = "-", nThread = params_list$numWorkers )

        return(TRUE)
      }
      if(verbose){
        message(paste("Org:",org,", Genes in Dir:",length(files_in_dir),", User Genes:",length(genes)," : Check FAILED!")) #", Available:",length(available_genes),", Missing:",length(missing_genes)
        #print(org,":",setdiff(available_genes,files_in_dir)) #DEBUG
      }
      return(FALSE)
    }

    if(verbose){
      message(paste("Org:",org,", ",fasta_path," : Check FAILED!"))
    }
    return(FALSE)
  }, error=function(cond){
    message(cond)
    if(verbose){
      message(paste("Org:",org,", ",fasta_path," : Check FAILED!"))
    }
    return(FALSE)
  }))
}

#' Get GTF data from BIOMART (using biomaRt)
#'
#' This function downloads GTF data from BIOMART using biomaRt. The attributes downloaded are
#' c("chromosome_name", "strand", "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name",
#' "start_position","end_position","exon_chrom_start","exon_chrom_end","transcript_version","transcript_length",
#' "5_utr_start","5_utr_end","transcript_start","transcription_start_site","transcript_end","cds_start","cds_end","3_utr_start","3_utr_end")
#'
#' @examples
#'     gtf_data <- get_gtf_mart(org, genes)
#'
#' @param org Name of the organism (format important, eg. "danio_rerio")
#' @param gene_list Vector or File with genes to fetch the GTF data for
#' @return GTF Data
#' @export
get_gtf_mart <- function(org, gene_list){
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
  # mart.dataset <-tryCatch(check_mart_dataset(org),error=function(cond){
  #   stop(cond)
  # })
  mart.dataset <- check_mart_dataset(org)
  using.mart.data <- mart_connect(biomaRt::useMart,args = list(COMPLETE_env$ENSEMBL_MART, mart.dataset))

  mart.attributes <- biomaRt::listAttributes(using.mart.data)

  if( !all(!is.na(match(gtf_attributes,mart.attributes$name))) ){
    stop(paste("GTF Attributes were not available for :",org,"\n"))
  }

  if (length(gene_list) == 1 && file.exists(gene_list)) {
    genes <- factor(scan(gene_list, character(), quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
  }

  ##Splitting genes to the recomended number of queries for biomart through biomaRt
  #bm_gtf <- dplyr::bind_rows(lapply(split(genes, ceiling(seq_along(genes)/500)), function(split_genes){
  #  return(mart_connect(biomaRt::getBM,args=list(mart=using.mart.data,attributes=gtf_attributes,uniqueRows=T, useCache=F, filters = c("external_gene_name"), values = split_genes, curl=COMPLETE_env$curl_handle)))
  #}))
  bm_gtf <- mart_connect(biomaRt::getBM,args=list(mart=using.mart.data,attributes=gtf_attributes,uniqueRows=T, useCache=F, filters = c("external_gene_name"), values = genes, curl=COMPLETE_env$curl_handle))
  bm_gtf <- dplyr::arrange(bm_gtf,chromosome_name,start_position)

  bm_gtf$strand[bm_gtf$strand==1] <- "+"
  #bm_gtf$strand[bm_gtf$strand== -1] <- "-"
  bm_gtf$strand[is.na(match(bm_gtf$strand,"+"))] <- "-"
  #bm_gtf$source_name <- "Ensembl"
  return(bm_gtf)
}

#' Get FASTA data from BIOMART (using biomaRt)
#'
#' This function downloads FASTA data from BIOMART using biomaRt. The Data Frame from calculate_stats() is required.
#'
#' @note Columns required from biomaRt "gene_name","gene_id","transcript_id" (or) "seqnames","total_exon_len","total_cds_len", "five_len","three_len", "exon_count","cds_count","g.exon_start","g.exon_end","transcript_length","five_flank","three_flank","g.transcript_start","g.transcript_end","chromosome_name","strand"
#'
#' @examples
#'     gtf_data <- get_gtf_mart(org, genes)
#'     gtf_stats <- calculate_stats(gtf_data)
#'     gtf_stats <- dplyr::bind_rows(gtf_stats)
#'     #names(gtf_stats)[grep(pattern="seqnames",names(gtf_stats))] <- "transcript_id"
#'     #gtf_stats$safe_gene_name <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(gtf_stats$gene_name))
#'     #gtf_stats <- gtf_stats[gtf_stats$gene_name!=0 | gtf_stats$gene_id!=0,] # cleaning up
#'     params_list <- load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt"))
#'     gtf_stats <- fetch_FASTA_mart(org,gtf_stats,org_fasta_path,params_list)
#'
#' @param org Name of the organism (format important, eg. "danio_rerio")
#' @param gtf_stats Data Frame from calculate_stats()
#' @param fasta_path Output path for transcript FASTA files
#' @param params_list Output of load_params()
#' @return Data Frame with transcripts for which FASTA was obtained
#' @export
fetch_FASTA_mart <- function(org,gtf_stats, fasta_path, params_list){

  req_columns <- c("gene_name","gene_id","transcript_id","total_exon_len","total_cds_len", "five_len","three_len", "exon_count","cds_count","g.exon_start","g.exon_end","transcript_length","five_flank","three_flank","g.transcript_start","g.transcript_end","chromosome_name","strand")

  if(any(grepl(x = class(gtf_stats), pattern = "list",ignore.case = T))){
    gtf_stats <- dplyr::bind_rows(gtf_stats)
  }

  if (any(grepl(pattern="seqnames",names(gtf_stats)))) {
    names(gtf_stats)[grep(pattern="seqnames",names(gtf_stats))] <- "transcript_id"
  }

  if (any(is.na(match(req_columns, names(gtf_stats))))) {
    stop(paste("Missing columns, Require :",paste(req_columns,collapse = ",")))
  }

  if (!any(grepl(pattern="safe_gene_name",names(gtf_stats)))) {
    gtf_stats$safe_gene_name <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(gtf_stats$gene_name))
  }

  gtf_stats <- gtf_stats[gtf_stats$gene_name!=0 | gtf_stats$gene_id!=0,] # cleaning up

  dir.create(path=fasta_path,showWarnings=F,recursive=T)
  seq_attributes <- c("5utr","coding","3utr") #"external_gene_name"
  # mart.dataset <-tryCatch(check_mart_dataset(org),error=function(cond){
  #   stop(cond)
  # })
  mart.dataset <- check_mart_dataset(org)
  using.mart.data <- mart_connect(biomaRt::useMart,args = list(COMPLETE_env$ENSEMBL_MART, mart.dataset))

  if(!all(!is.na(match(seq_attributes, biomaRt::listAttributes(using.mart.data)[,c("name")])))){
    stop(paste("Atrributes :",paste(seq_attributes,collapse = ","), ": not availabe for :", org,"\n"))
  }

  check_rows_idx <- unlist(sapply(seq_along(1:nrow(gtf_stats)), FUN=function(x){
    #print(gtf_stats[x,]<=3)
    if(any(na.omit(gtf_stats[x,c("five_len","three_len","five_flank","three_flank")]<=3))){ #any(is.na(gtf_stats[x,])) , not checking for CDS==NA because I was able to obtain coding sequences even when cds_len==NA
      return(x)
    }
  }))
  invalid_stats <- gtf_stats[check_rows_idx,]
  valid_transcripts <- unique(gtf_stats[-check_rows_idx,]$transcript_id)

  bm_df <- c()

  if(length(valid_transcripts) > 0){
    bm_seq <- parallel::mclapply(seq_attributes, function(x){return(mart_connect(biomaRt::getBM,args=list(mart=using.mart.data,attributes=c("ensembl_transcript_id",x),uniqueRows=T, useCache=F, filters = c("ensembl_transcript_id"), values = valid_transcripts, curl=COMPLETE_env$curl_handle)))},mc.cores = params_list$numWorkers)

    #print(bm_seq) #DEBUG
    bm_df <- purrr::reduce(bm_seq, dplyr::full_join, by = "ensembl_transcript_id")
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
  }

  flank_stats <- invalid_stats #invalid_stats[is.na(match(invalid_stats$transcript_id,missing_any_flank_info)),]
  if (nrow(flank_stats) > 0) {

    missing_any_flank_info <- unique(c(flank_stats$transcript_id[flank_stats$five_flank==0],flank_stats$transcript_id[flank_stats$three_flank==0]))

    if (length(missing_any_flank_info) > 0) { ##All these transcripts have neither UTR Lengths nor Flank lengths info (possibly because the gene does not have isoforms or info is missing), so I set the utr flanks to transcript lengths and correct for variance

      arbitrary_flanks <- flank_stats[which(!is.na(match(flank_stats$transcript_id,missing_any_flank_info))),]
      arbitrary_flanks <- arbitrary_flanks[order(arbitrary_flanks$gene_name),]
      arb_flank_values <- dplyr::bind_rows(parallel::mclapply(base::split(arbitrary_flanks,as.factor(arbitrary_flanks$gene_name)), function(x){
        x$five_flank[ which(x$five_flank==0) ] <- ceiling(mean(unique(arbitrary_flanks$transcript_length)))
        x$three_flank[ which(x$three_flank==0) ] <- ceiling(mean(unique(arbitrary_flanks$transcript_length)))
        #Do variance correction
        return( data.frame(gene_name=x$gene_name,transcript_id=x$transcript_id,five_flank= x$five_flank + ceiling(mean(unique(abs(x$five_len-x$five_flank)))), three_flank=x$three_flank + ceiling(mean(unique(abs(x$three_len-x$three_flank)))) ) )
      }, mc.cores = params_list$numWorkers))
      arbitrary_flanks <- arbitrary_flanks %>% dplyr::select(-c("five_flank","three_flank")) %>% full_join(arb_flank_values, by = c("gene_name","transcript_id"))
      flank_stats <- flank_stats[which(is.na(match(flank_stats$transcript_id,missing_any_flank_info))),]
      flank_stats <- unique(full_join(flank_stats, arbitrary_flanks, by = c("gene_name", "gene_id", "transcript_id", "total_exon_len", "total_cds_len", "five_len", "three_len", "exon_count", "cds_count", "g.exon_start", "g.exon_end", "transcript_length", "five_flank", "three_flank", "g.transcript_start", "g.transcript_end", "chromosome_name", "strand", "safe_gene_name")))
    }

    ##Calculate genomic coordinates for flanks
    flank_stats <- flank_stats %>% mutate(g.five_flank_start= abs(flank_stats$g.transcript_start-flank_stats$five_flank)-1 )
    flank_stats <- flank_stats %>% mutate(g.five_flank_end= abs(flank_stats$g.transcript_start-1) )
    flank_stats <- flank_stats %>% mutate(g.three_flank_start= abs(flank_stats$g.transcript_end+1) )
    flank_stats <- flank_stats %>% mutate(g.three_flank_end=  abs(flank_stats$g.transcript_end+flank_stats$three_flank)+1 )

    if (nrow(flank_stats) > 0) {
      ##Get Flanking sequences
      flank_values <- unique(flank_stats[,c("transcript_id","five_flank","three_flank")])
      bm_flanks_cds <- mart_connect(biomaRt::getSequence,args = list(id=flank_values$transcript_id,type="ensembl_transcript_id",seqType="coding", mart=using.mart.data, useCache = F) )
      bm_flanks_five <- mart_connect(biomaRt::getSequence,args = list(id=flank_stats$transcript_id,type="ensembl_transcript_id",seqType="coding_transcript_flank", upstream=flank_values$five_flank , mart=using.mart.data, useCache = F) )#mart_connect(biomaRt::getBM,args=list(mart=using.mart.data,attributes=c("start_position","end_position","cdna"),uniqueRows=T, useCache=F, filters = c("chromosome_name","start","end","strand"), values = list(as.character(flank_stats$chromosome_name), flank_stats$g.five_flank_start, flank_stats$g.five_flank_end, params_list$STRAND), curl=COMPLETE_env$curl_handle))
      bm_flanks_three <- mart_connect(biomaRt::getSequence,args = list(id=flank_values$transcript_id,type="ensembl_transcript_id",seqType="coding_transcript_flank", downstream = flank_values$three_flank , mart=using.mart.data, useCache = F))

      names(bm_flanks_five)[grep(pattern="coding_transcript_flank",names(bm_flanks_five))] <- "5utr"
      names(bm_flanks_three)[grep(pattern="coding_transcript_flank",names(bm_flanks_three))] <- "3utr"

      bm_flanks_five <- bm_flanks_five %>% mutate(flanking_5utr=T)
      bm_flanks_three <- bm_flanks_three %>% mutate(flanking_3utr=T)

      transcripts_with_data <- intersect(flank_values$transcript_id,unique(c(bm_flanks_five$ensembl_transcript_id,bm_flanks_three$ensembl_transcript_id, bm_flanks_cds$ensembl_transcript_id)))

      bm_seqs <- purrr::reduce(list(bm_flanks_five,bm_flanks_cds,bm_flanks_three), full_join, by="ensembl_transcript_id")

      unavailable_transcripts <- unique( unlist(apply(bm_seqs,MARGIN = 1, FUN = function(x){
        row_check <- grepl("unavailable",x=x,ignore.case = T)
        names(row_check) <- names(bm_seqs)
        #print(row_check)
        if(all(row_check) || row_check["coding"]==T){ #Removing transcripts which do not have any regions or CDS
          return(as.character(x["ensembl_transcript_id"]))
        }
      })) )

      if(length(unavailable_transcripts)>0){
        flank_stats <- flank_stats[!is.na(match(flank_stats$transcript_id,unavailable_transcripts)),]
        #message(paste("(Some) Data missing for : ",org,":",paste(unique(unavailable_transcripts),sep=",",collapse=","),". Maybe CDS or all regions are missing for the transcripts."))
      }
      if(length(transcripts_with_data)>0){
        flank_stats <- flank_stats[!is.na(match(flank_stats$transcript_id,transcripts_with_data)),]
      }
      if (nrow(flank_stats) > 0) {
        invalid_stats <- flank_stats
      }
    }
  }
  #print("here1") #DEBUG
  #tmp1 <<- bm_df #DEBUG
  #tmp2 <<- bm_seqs #DEBUG
  #tmp3 <<- gtf_stats #DEBUG
  if(!is.null(bm_df)){
    if(nrow(bm_seqs)>0 && nrow(bm_df)>0){
      bm_df <- purrr::reduce(list(bm_seqs,bm_df), dplyr::full_join,by = c("5utr", "ensembl_transcript_id", "coding", "3utr","flanking_5utr","flanking_3utr"))
    }else if(nrow(bm_seqs)>0){
      bm_df <- bm_seqs
    }else{
      stop(paste("Error fetching FASTA for :",org,"\n"))
    }
  }else{
    if(nrow(bm_seqs)>0){
      bm_df <- bm_seqs
    }else{
      stop(paste("Error fetching FASTA for :",org,"\n"))
    }
  }
  names(bm_df)[grep(pattern="ensembl_transcript_id",names(bm_df))] <- "transcript_id"
  names(bm_df)[grep(pattern="coding",names(bm_df))] <- "cds"
  #print("here2") #DEBUG
  bm_df <- inner_join(gtf_stats[,c("gene_name","safe_gene_name","transcript_id","strand")],bm_df, by = "transcript_id")
  bm_df <- unique(bm_df)

  final_unavailable_transcripts <- unique( unlist(apply(bm_df,MARGIN = 1, FUN = function(x){
    row_check <- grepl("unavailable",x=x,ignore.case = T)
    names(row_check) <- colnames(bm_df)
    #print(row_check)
    if(all(row_check[params_list$TRANSCRIPT_REGIONS]) || row_check["cds"]==T){ #Removing transcripts which do not have any regions or CDS
      return(as.character(x["transcript_id"]))
    }
  })) )

  final_unavailable_transcripts <- unique( c(unavailable_transcripts, final_unavailable_transcripts) )
  if(length(final_unavailable_transcripts)>0){
    bm_df <- bm_df[which(is.na(match(bm_df$transcript_id,final_unavailable_transcripts))),]
    message(paste("(Some) Data missing for : ",org,": Stored in :",paste(params_list$OUT_PATH,"/genes/",org,"/non_coding.csv",sep=""),". Maybe CDS or all regions are missing for the transcripts. This could happen for non-protein coding transcripts or retained introns"))

    non_coding_data <- unique(gtf_stats[which(!is.na(match(gtf_stats$transcript_id,final_unavailable_transcripts))),c("gene_name","transcript_id")])
    data.table::fwrite(list(non_coding_data),file = paste(params_list$OUT_PATH,"/genes/",org,"/non_coding_biomart.csv",sep=""),quote = F,row.names = F,col.names = T,sep = ",",na = "-", nThread = params_list$numWorkers)
  }

  #print(head(bm_df))
  #tmp4 <<- bm_df #DEBUG
  #tmp_gtf <<- gtf_stats #DEBUG
  #tmp_missing <<- final_unavailable_transcripts #DEBUG

  parallel::mclapply(base::split(bm_df[,c("transcript_id","gene_name","safe_gene_name",params_list$TRANSCRIPT_REGIONS,"flanking_5utr","flanking_3utr","strand")],as.factor(bm_df$gene_name)), function(x){
    lapply(params_list$TRANSCRIPT_REGIONS, function(y){
      if (grepl("3utr|5utr",y,ignore.case = T)) { ##Adding "_FLANK" for transcripts with UTR Flanks
        flank_col <- grep(paste("flanking",y,sep="_"),colnames(x),ignore.case = T,value = T)
        x_unique <- unique(x[,c("transcript_id","safe_gene_name",y,flank_col,"strand")])
        seq_names <- paste(as.character(x_unique[,"transcript_id"]), y,sep = params_list$TRANSCRIPT_ID_DELIM)
        seq_names[x_unique[,flank_col]==T] <- paste(seq_names[x_unique[,flank_col]==T], "FLANK",sep = "_")
      }else{
        x_unique <- unique(x[,c(y,"transcript_id","safe_gene_name","strand")])
        seq_names <- paste(as.character(x_unique[,"transcript_id"]), y,sep = params_list$TRANSCRIPT_ID_DELIM)
      }
      seq_names <- paste(seq_names,"(",x_unique[,"strand"],")",sep = "")
      #tmp5 <<- x_unique[,y] #DEBUG
      fasta_as_c <- as.character(x_unique[,y])
      names(fasta_as_c) <- seq_names
      fasta_to_write <- Biostrings::DNAStringSet(x = fasta_as_c, use.names = T)
      #tmp6 <<- fasta_to_write #DEBUG
      #message(fasta_path)
      #seqinr::write.fasta(as.list(x_unique[,y]), seq_names,file.out=paste(fasta_path,"/",unique(x_unique[,"safe_gene_name"]),".",y,".tmp",sep="") )
      Biostrings::writeXStringSet(x = fasta_to_write,filepath = paste(fasta_path,"/",unique(x_unique[,"safe_gene_name"]),".",y,sep=""), append = F,format = "fasta")

    })
  },mc.cores = params_list$numWorkers, mc.silent = F)

  return(gtf_stats[!is.na(match(gtf_stats$transcript_id,bm_df$transcript_id)),])
}

#' Internal Function - Get FASTA data from BIOMART (using biomartr)
#'
#' This function downloads FASTA data from BIOMART using biomartr. The GENOME and GTF for the organism are downloaded and passed through the shell pipeline for extracting transcripts.
#' This function invokes external SHELL function extract_genomic_regions from fs::path_package("COMPLETE","exec","functions.sh") (just like the piepline for user data) and cannot be monitored
#'
#' @examples
#'     params_list <- load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt"))
#'     fetch_FASTA_biomartr(c(name="danio_rerio",version="106"),params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'
#' @param org_row Named vector with name of the organism (format important, eg. "danio_rerio") and other details eg, c(name="danio_rerio",version="106")
#' @param params_list Output of load_params()
#' @param gene_list Vector or File containing list of genes
#' @param verbose Print Messages?
#' @return Named Vector of the organism details
fetch_FASTA_biomartr <- function(org_row, params_list, gene_list,verbose=T){
  org <- org_row["name"]

  # if (is.null(COMPLETE_env$org.meta)) {
  #   COMPLETE_env <- new.env(parent=emptyenv())
  #   COMPLETE_env$org.meta <- mart_connect(biomartr::listGenomes,args=list(db = "ensembl", type = "all", details = T)) #For biomartr #db = tolower(GENOMES_SOURCE)
  # }

  tictoc::tic(msg=paste("Processed:",org))

  if(any(!is.na(match(COMPLETE_env$org.meta$name,org)))){
    org <- COMPLETE_env$org.meta[which(!is.na(match(COMPLETE_env$org.meta$name,org))),]
    org_name <- org$name
    genome_path<-paste(params_list$GENOMES_PATH, "/",org_name,".fa.gz",sep = "")
    gtf_path<-paste(params_list$ANNOS_PATH, "/",org_name,".gtf.gz",sep = "")
    org_fasta_path <- file.path(params_list$FASTA_OUT_PATH ,org_name)

    if( suppressMessages( biomartr::is.genome.available(db="ensembl", organism = org_name) ) ){
      if(!file.exists(gtf_path) || file.info(gtf_path)$size <= 20 || params_list$CLEAN_EXTRACT){
        gtf_ori <- biomartr::getGTF(organism = org_name, db="ensembl",path = params_list$ANNOS_PATH) #,release=org$release-1)
        if(!is.logical(gtf_ori)){
          file.rename(tools::file_path_as_absolute(gtf_ori),gtf_path)
        }
      }
      if(!file.exists(genome_path) || file.info(genome_path)$size <= 20 || params_list$CLEAN_EXTRACT){
        genome_ori <- biomartr::getGenome(organism = org_name, db="ensembl",path = params_list$GENOMES_PATH,reference = T,gunzip = F) #,release=org$release-1)
        if(!is.logical(genome_ori)){
          file.rename(tools::file_path_as_absolute(genome_ori),genome_path)
        }
      }
    }else{
      stop(paste("Organism not available :", org,"\n"))
    }

    if(length(gene_list) == 1 && file.exists(gene_list)){
      genes <- factor(scan(gene_list, character(), quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
      genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
    }else{
      genes <- tolower(as.vector(gene_list))
      tmp_gene_list <- tempfile(pattern="genelist",tmpdir = params_list$TEMP_PATH)
      data.table::fwrite(x = list(gene_list),file = tmp_gene_list ,quote = F,row.names = F,col.names = F, nThread = params_list$numWorkers)
      gene_list <- tmp_gene_list
    }

    if(params_list$CLEAN_EXTRACT || !check_files(fasta_path = org_fasta_path,org = org_name,genes = genes, verbose = verbose, params_list = params_list)){
      if ( (!is.logical(gtf_path) && !is.logical(genome_path) && file.exists(gtf_path) && file.exists(genome_path) && file.info(gtf_path)$size > 20 && file.info(genome_path)$size > 20)  && !COMPLETE_env$SKIP_USER_DATA) {
        dir.create(path = org_fasta_path,showWarnings = F,recursive = T)
        ##do.call(add_to_process,list(p_cmd = c(system.file("exec", "jobhold.sh", mustWork = T ,package = "COMPLETE")), p_args = c(param_file,paste("extract",org_name,sep="_"), system.file("exec", "extract_genomic_regions.sh", mustWork = T ,package = "COMPLETE"),genome_path, gtf_path, gene_list, org_name, param_file)))
        #do.call(add_to_process,list(p_cmd = c(system.file("exec", "jobhold.sh", mustWork = T ,package = "COMPLETE")), p_args = c(param_file,paste("extract",org_name,sep="_"), fs::path_package("COMPLETE","exec","functions.sh"),"extract_genomic_regions",genome_path, gtf_path, gene_list, org_name, param_file)))
        cat(paste("Logfile : ",params_list$TEMP_PATH,"/",org_name,".log\n",sep=""))
        do.call(add_to_process,list(p_cmd = COMPLETE_env$SHELL, p_args = c(fs::path_package("COMPLETE","exec","functions.sh"),"extract_genomic_regions",genome_path, gtf_path, gene_list, org_name, params_list$param_file), logfile=paste(params_list$TEMP_PATH,"/",org_name,".log",sep=""), params_list = params_list,verbose = verbose ))
        ##return(do.call(add_to_process,list(p_cmd = c("./extract_genomic_regions.sh"), p_args = c(genome_path, gtf_path, gene_list, org))))
        cat(print_toc(tictoc::toc(quiet = T, log = T)))
        return(c(org_row, source="r-biomartr"))
      }else{
        message(print_toc(tictoc::toc(quiet = T, log = T)))
        return(NULL)
      }
    }else{
      cat(print_toc(tictoc::toc(quiet = T, log = T)))
      return(c(org_row, source="r-biomartr"))
    }

  }else{
    message(print_toc(tictoc::toc(quiet = T, log = T)))
    stop(paste("Organism not available :", org,"\n"))
  }
}

#' Internal Function - Get FASTA data
#'
#' Main function to download FASTA data.
#'
#' @note This function wraps around other fetch_FASTA_*() function (Except for fetch_FASTA_user() which calls fetch_FASTA() if a genome or a gtf are not provided)
#'
#' @examples
#'     params_list <- load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt"))
#'     fetch_FASTA(c(name="danio_rerio",version="106"), params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'
#' @param org_row Named vector with name of the organism (format important, eg. "danio_rerio") and other details eg, c(name="danio_rerio",version="106")
#' @param params_list Output of load_params()
#' @param gene_list Vector or File containing list of genes
#' @param verbose Print Messages?
#' @return Named Vector of the organism details
#' @export
fetch_FASTA <- function(org_row, params_list, gene_list, verbose=T) {

  org <- org_row["name"]
  tictoc::tic(msg=paste("Processed:",org))
  org_fasta_path <- file.path(params_list$FASTA_OUT_PATH,org)

  #if(!is.vector(gene_list)){
  if (length(gene_list) == 1 && file.exists(gene_list)) {
    genes <- factor(scan(gene_list, character(),quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
    tmp_gene_list <- tempfile(pattern="genelist",tmpdir = params_list$TEMP_PATH)
    data.table::fwrite(x = list(gene_list),file = tmp_gene_list ,quote = F,row.names = F,col.names = F, nThread = params_list$numWorkers)
    gene_list <- tmp_gene_list
  }

  #print(check_files(fasta_path = org_fasta_path,org = org,genes = genes,params_list = params_list))

  if(!params_list$CLEAN_EXTRACT && check_files(fasta_path = org_fasta_path,org = org,genes = genes,params_list = params_list, verbose = verbose)){
    cat(print_toc(tictoc::toc(quiet = T, log = T)))
    return(c(org_row, source="r-biomaRt"))
  }

  tryCatch(check_mart_dataset(org),
           error=function(cond1){
             #message(cond)

             tryCatch({
               fetch_FASTA_biomartr(org_row = org_row, params_list = params_list, gene_list = gene_list, verbose = F)
             }, error=function(cond2){
               message(print_toc(tictoc::toc(quiet = T, log = T)))
               #message(cond1)
               stop(cond2)
             })
             message(print_toc(tictoc::toc(quiet = T, log = T)))
             #message(cond1)
             #return(org_row)
             stop(cond1)
           })

  odb_list <- paste(params_list$OUT_PATH,"/genes/",org,"/odb.list",sep = "")
  odb_gene_map <- paste(params_list$OUT_PATH,"/genes/",org,"/odb.final_map",sep = "")
  gtf_stats_file <- paste(params_list$OUT_PATH,"/genes/",org,"/gtf_stats.csv",sep = "")

  if(params_list$CLEAN_EXTRACT){
    unlink(org_fasta_path,recursive = T,force=T,expand = T)
    unlink(odb_list,force=T,expand = T)
    unlink(odb_gene_map,force=T,expand = T)
    unlink(gtf_stats_file,force=T,expand = T)
  }

  dir.create(paste(params_list$OUT_PATH,"/genes/",org, sep=""),showWarnings = F,recursive = T)

  odb_list_genes <- c()
  if(COMPLETE_env$USE_ORTHODB){
    if(params_list$CLEAN_EXTRACT || !file.exists(odb_list) || file.info(odb_list)$size == 0){
      #proc <- do.call(add_to_process,list(p_cmd = c(system.file("exec", "check_OrthoDB.sh", mustWork = T ,package = "COMPLETE")),p_args = c(org,gene_list, odb_list, odb_gene_map,param_file))) #time
      proc <- do.call(add_to_process,list(p_cmd = COMPLETE_env$SHELL, p_args = c(fs::path_package("COMPLETE","exec","functions.sh"), "check_OrthoDB",org,gene_list, odb_list, odb_gene_map,params_list$param_file, COMPLETE_env$SELECT_ALL_GENES), params_list=params_list,verbose = verbose))
      proc$wait(timeout=-1)
    }
    if(file.exists(odb_list) && file.info(odb_list)$size > 0){
      odb_list_genes <- factor(scan(odb_list, character(), quiet = T))
      odb_list_genes <- odb_list_genes[grep("gene",tolower(odb_list_genes), invert = T, fixed = T)]
      #odb_list_data <- get_gtf_mart(org, odb_list_genes)
      #gtf_data <- unique(merge(odb_list_data,gtf_data))
    }else{
      message(paste("ODB gene list could not be found for : ",org))
    }
  }

  gtf_data <- invisible( tryCatch(
    get_gtf_mart(org = org, gene_list = unique(c(genes,odb_list_genes)))
    ,error=function(cond1){
      #message(cond)
      #message(print_toc(tictoc::toc(quiet = T, log = T)))
      return(tryCatch(
        fetch_FASTA_biomartr(org_row = org_row, params_list = params_list, gene_list = genes, verbose = F)
        , error=function(cond2){
          message(print_toc(tictoc::toc(quiet = T, log = T)))
          stop(cond2)
        }))
      message(print_toc(tictoc::toc(quiet = T, log = T)))
      stop(cond1)
    }) )

  gtf_stats <- calculate_stats(gtf_data,allow_strand = params_list$STRAND, n_threads = params_list$numWorkers)
  gtf_stats <- dplyr::bind_rows(gtf_stats)
  names(gtf_stats)[grep(pattern="seqnames",names(gtf_stats))] <- "transcript_id"
  #print(head(gtf_stats))
  gtf_stats$safe_gene_name <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(gtf_stats$gene_name))
  gtf_stats <- gtf_stats[gtf_stats$gene_name!=0 | gtf_stats$gene_id!=0,] # cleaning up

  if(nrow(gtf_stats)==0){
    message(paste("No genes passed filters for : ",org))
    message(print_toc(tictoc::toc(quiet = T, log = T)))
    return()
  }

  gtf_stats <- invisible( tryCatch(fetch_FASTA_mart(org = org,gtf_stats = gtf_stats,fasta_path = org_fasta_path, params_list = params_list),error=function(cond){
    message(cond1)
    ##message(print_toc(tictoc::toc(quiet = T, log = T)))
    tryCatch(fetch_FASTA_biomartr(org_row = org_row, params_list = params_list, gene_list = genes, verbose = F), error=function(cond2){
      cat(print_toc(tictoc::toc(quiet = T, log = T)))
      stop(cond2)
    })
    stop(cond1)
  }) )

  if(nrow(gtf_stats)==0){
    message(paste("Error fetching FASTA for : ",org))
    message(print_toc(tictoc::toc(quiet = T, log = T)))
    return()
  }

  # local_procs <- c()
  # local_procs <- unlist(parallel::mclapply(base::split(gtf_stats,as.factor(gtf_stats$gene_name)), function(x){
  #   gene <- unique(x[,"gene_name"])
  #   safe_name <- unique(x[,"safe_gene_name"])
  #   return(unlist(lapply(params_list$TRANSCRIPT_REGIONS, function(y){
  #     output_fasta <- paste(org_fasta_path,"/",safe_name,".",y,sep="")
  #     input_fasta <- paste(org_fasta_path,"/",safe_name,".",y,".tmp",sep="")
  #     proc <- do.call(add_to_process,list(p_cmd = c(fs::path_package("COMPLETE","exec","functions.sh")),p_args = c("label_FASTA_files",org,gene, output_fasta,input_fasta,params_list$param_file, odb_gene_map), verbose=F, params_list=params_list))
  #     return(proc)
  #   })))
  # }, mc.cores = params_list$numWorkers,mc.silent = T))
  # invisible(parallel::mclapply(local_procs, function(x){
  #   #print(x)
  #   if(x$is_alive()){
  #     x$wait(timeout=-1)
  #   }
  #   x$kill(close_connections = TRUE)
  # }, mc.cores = params_list$numWorkers))
  #print(paste(org_fasta_path,org,genes,odb_gene_map,params_list)) #DEBUG
  
  label_FASTA_files(fasta_path = org_fasta_path,org = org,gene_list = genes,odb_gene_map = odb_gene_map,params_list = params_list, duplicates.method = "merge")
  
  names(gtf_stats)[grep(pattern="transcript_length",names(gtf_stats))] <- "transcript_length.annotated"
  gtf_stats <- gtf_stats %>% mutate(org=org)
  gtf_stats$transcript_length.estimated <- gtf_stats$five_flank + gtf_stats$total_cds_len + gtf_stats$three_flank
  gtf_stats <- gtf_stats[,c("gene_name","gene_id","transcript_id","total_cds_len","five_len","three_len","exon_count","cds_count","five_flank","three_flank","transcript_length.estimated","transcript_length.annotated","org")]
  gtf_stats$gene_name <- tolower(gtf_stats$gene_name)

  data.table::fwrite(x = list(gtf_stats),file = gtf_stats_file,quote = F,sep = ",",row.names = F,col.names = T,na = "-", nThread = params_list$numWorkers)
  #write.table(x= unique(tolower(gtf_stats$gene_name)),file=paste(params_list$OUT_PATH,"/genes/",org ,"/AVAILABLE_GENES",sep=""),quote=F,row.names=F,col.names=F,na = "-" )

  unlink( unlist(parallel::mclapply(list.dirs(path = org_fasta_path, full.names=TRUE,recursive = F), function(x) {
    if (file.info(x)$size == 0) {
      return(x) #print(x)
    }
  }, mc.cores = params_list$numWorkers, mc.silent = T, mc.preschedule = T)) ,recursive = F, force = T, expand =T)

  orgs_files <- tolower(unique(tools::file_path_sans_ext(list.files(org_fasta_path,no.. = T,recursive = F))))
  data.table::fwrite(x= list(orgs_files),file=paste(params_list$OUT_PATH,"/genes/",org ,"/AVAILABLE_GENES",sep=""),quote=F,row.names=F,col.names=F,na = "-" , nThread = params_list$numWorkers )
  data.table::fwrite(x= list(tolower(unique(genes[is.na(match( tolower(genes), orgs_files ))]))) ,file=paste(params_list$OUT_PATH,"/genes/",org ,"/MISSING_GENES",sep=""),quote=F,row.names=F,col.names=F ,na = "-", nThread = params_list$numWorkers ) #tolower(unique(gtf_stats$gene_name))

  cat(print_toc(tictoc::toc(quiet = T, log = T)))
  return(c(org_row, source="r-biomaRt"))
}

#' Get FASTA for User data
#'
#' Function to download FASTA data for user data. fetch_FASTA() is invoked if a genome and a gtf are not provided.
#'
#' @note Use "-" in Genome/GTF to check for the organism in BIOMART
#'
#' @examples
#'     params_list <- load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt"))
#'     fetch_FASTA_user(c(org="danio_rerio",genome="http://some.link",gtf="some.file",version="106"),params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'     fetch_FASTA_user(c(org="danio_rerio",genome="-",gtf="-",version="106"),params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'     fetch_FASTA_user(c(org="xenopus_laevis",genome="https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz",gtf="https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/XENLA_10.1_GCF.gff3.gz",version="106"),params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'
#' @param data Named vector with name of the organism (format important, eg. "danio_rerio"), Genome, GTF and other details eg, c(org="danio_rerio",genome="http://some.link",gtf="some.file",version="106").
#' @param params_list Output of load_params()
#' @param gene_list Filename with the list of genes
#' @param verbose Print Messages?
#' @return Named Vector of the organism details
fetch_FASTA_user <- function(data, params_list, gene_list, verbose=T){

  genome_path <- c()
  gtf_path <- c()
  #print(data)
  genome <- data["genome"]
  gtf <- data["gtf"]
  org <- data["org"]

  #if(!is.vector(gene_list)){
  if(length(gene_list) == 1 && file.exists(gene_list)){
    genes <- factor(scan(gene_list, character(),quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
    data.table::fwrite(x = list(gene_list),file = paste(params_list$TEMP_PATH,"/gene_list.txt",sep = ""),quote = F,row.names = F,col.names = F, nThread = params_list$numWorkers)
    gene_list <- paste(params_list$TEMP_PATH,"/gene_list.txt",sep = "")
  }

  if(genome == "-" || gtf == "-" || genome == " " || gtf == " " || stringi::stri_isempty(genome) || stringi::stri_isempty(gtf)){
    return( tryCatch(fetch_FASTA(org_row = c(name=as.character(org),genome="-",gtf="-"), params_list = params_list, gene_list = genes, verbose = T),error=function(cond){
      stop(cond)
      #return(NULL)
    }) )
    #break;
  }

  tictoc::tic(msg = paste("Processed:",org))
  
  if(grepl(x = basename(URLdecode(genome)), pattern="gz")){
    genome_path<-paste(params_list$GENOMES_PATH, "/",org,".",tools::file_ext(tools::file_path_sans_ext(basename(URLdecode(genome)))),sep = "")  #".gtf.gz"
  }else{
    genome_path<- paste(params_list$GENOMES_PATH, "/",org,".",tools::file_ext(basename(URLdecode(genome))),sep = "")  
  }
  
  if(grepl(x = basename(URLdecode(gtf)), pattern="gz")){
    gtf_path<-paste(params_list$ANNOS_PATH, "/",org,".",tools::file_ext(tools::file_path_sans_ext(basename(URLdecode(gtf)))),sep = "")  #".gtf.gz"
  }else{
    gtf_path<-paste(params_list$ANNOS_PATH, "/",org,".",tools::file_ext(basename(URLdecode(gtf))),sep = "")  #".gtf.gz"
  }
  org_fasta_path <- file.path(params_list$FASTA_OUT_PATH ,org)

  if(grepl("://|http|ftp|www",genome)){
    if(!file.exists(genome_path) || !file.info(genome_path)$size > 20 || params_list$CLEAN_EXTRACT){
      ret_code <- curl::curl_fetch_disk(genome, paste(params_list$GENOMES_PATH, "/",basename(URLdecode(genome)),sep=""))
      if(ret_code$status_code==404){
        message(paste("Error 404 : check genome URL :",org,"-",gtf))
        message(print_toc(tictoc::toc(quiet = T, log = T)))
        return(NULL)
      }
      #genome_path <- paste(params_list$GENOMES_PATH, "/",basename(URLdecode(genome)),sep="")
      file.rename(paste(params_list$GENOMES_PATH, "/",basename(URLdecode(genome)),sep=""),genome_path)
    } # else{
    #   if(system2("gzip",args = c("-t", genome_path), wait = T,stdout = NULL, stderr = NULL) == 0){
    #     genome_path<-genome_path
    #   }else{
    #     lapply(list.files(params_list$GENOMES_PATH, full.names = T, ignore.case = T, no.. = T, pattern = regex(stri_split_fixed(org,"_",2,simplify = T),ignore_case = T)),function(x){
    #       if(grepl(x,pattern = "fa|gz")){
    #         try(file.remove(x, showWarnings=F))
    #       }
    #     })
    #     curl::curl_fetch_disk(genome, basename(URLdecode(genome)))
    #     genome_path <- paste(params_list$GENOMES_PATH, "/",basename(URLdecode(genome)),sep="") #basename(URLdecode(genome))
    #   }
    # }

  }else{
    if(file.exists(genome) && file.info(genome)$size > 0){
      genome_path <- genome
    }else{
      message(paste("User Genome not found :",org,"-",genome))
      message(print_toc(tictoc::toc(quiet = T, log = T)))
      return(NULL)
    }
  }

  if(grepl("://|http|ftp|www",gtf)){
    if(!file.exists(gtf_path) || !file.info(gtf_path)$size > 20 || params_list$CLEAN_EXTRACT){
      ret_code <- curl::curl_fetch_disk(gtf, paste(params_list$ANNOS_PATH, "/",basename(URLdecode(gtf)),sep=""))
      #gtf_path <- paste(params_list$ANNOS_PATH, "/",basename(URLdecode(gtf)),sep="")
      if(ret_code$status_code==404){
        message(paste("Error 404 : check gtf URL :",org,"-",gtf))
        message(print_toc(tictoc::toc(quiet = T, log = T)))
        return(NULL)
      }
      #zcat -f $ANNO_FILE | gffread - -T -O -E -o - | gzip -c > $params_list$ANNOS_PATH/"$f_org_name".gtf.gz
      if(!grepl(pattern = "gtf",x = basename(URLdecode(gtf)),ignore.case = T)){

      }
      file.rename(paste(params_list$ANNOS_PATH, "/",basename(URLdecode(gtf)),sep=""),gtf_path)
    }# else{
    #   if(system2("gzip",args = c("-t", gtf_path), wait = T,stdout = NULL, stderr = NULL) == 0){
    #     gtf_path<-gtf_path
    #   }else{
    #     lapply(list.files(params_list$ANNOS_PATH, full.names = T, ignore.case = T, no.. = T, pattern = regex(stri_split_fixed(org,"_",2,simplify = T),ignore_case = T)),function(x){
    #       if(grepl(x,pattern = "gtf|gz")){
    #         try(file.remove(x, showWarnings=F))
    #       }
    #     })
    #     curl::curl_fetch_disk(gtf, paste(params_list$ANNOS_PATH, "/",basename(URLdecode(gtf)),sep=""))
    #     gtf_path <- paste(params_list$ANNOS_PATH, "/",basename(URLdecode(gtf)),sep="")
    #   }
    # }

  }else{
    if(file.exists(gtf) && file.info(gtf)$size > 0){
      gtf_path <- gtf
    }else{
      message(paste("User GTF not found :",org,"-",gtf))
      message(print_toc(tictoc::toc(quiet = T, log = T)))
      return(NULL)
    }
  }

  # if(!is.logical(genome_path) && file_ext(genome_path) == "gz"){
  #   try(file_move(genome_path, paste(params_list$GENOMES_PATH, "/",org,".fa.gz",sep = "")))
  #   genome_path <- paste(params_list$GENOMES_PATH, "/",org,".fa.gz",sep = "")
  # }else if(!is.logical(genome_path)){
  #   g_ext <- file_ext(genome_path)
  #   try(file_move(genome_path, paste(params_list$GENOMES_PATH, "/",org,".",g_ext,sep = "")))
  #   genome_path <- paste(params_list$GENOMES_PATH, "/",org,".",g_ext,sep = "")
  # }
  # if(!is.logical(gtf_path) && file_ext(gtf_path) == "gz"){
  #   try(file_move(gtf_path, paste(params_list$ANNOS_PATH, "/",org,".gtf.gz",sep = "")))
  #   gtf_path <- paste(params_list$ANNOS_PATH, "/",org,".gtf.gz",sep = "")
  # }else if(!is.logical(gtf_path)){
  #   a_ext <- file_ext(gtf_path)
  #   try(file_move(gtf_path, paste(params_list$ANNOS_PATH, "/",org,".",a_ext,sep = "")))
  #   gtf_path <- paste(params_list$ANNOS_PATH, "/",org,".",a_ext,sep = "")
  # }

  #print(paste(genome_path,gtf_path,org))

  if(params_list$CLEAN_EXTRACT || !check_files(fasta_path = org_fasta_path,org = org,genes = genes, params_list = params_list, verbose = verbose)){
    if ( (!is.logical(gtf_path) && !is.logical(genome_path)) && !COMPLETE_env$SKIP_USER_DATA) {
      dir.create(path = org_fasta_path,showWarnings = F,recursive = T)
      ##do.call(add_to_process,list(p_cmd = c(system.file("exec", "jobhold.sh", mustWork = T ,package = "COMPLETE")), p_args = c(param_file,paste("extract",org,sep="_"), system.file("exec", "extract_genomic_regions.sh", mustWork = T ,package = "COMPLETE"),genome_path, gtf_path, gene_list, org,param_file)))
      #do.call(add_to_process,list(p_cmd = c(system.file("exec", "jobhold.sh", mustWork = T ,package = "COMPLETE")), p_args = c(param_file,paste("extract",org,sep="_"), fs::path_package("COMPLETE","exec","functions.sh"),"extract_genomic_regions",genome_path, gtf_path, gene_list, org,param_file)))
      cat(paste("Logfile : ",params_list$TEMP_PATH,"/",org,".log\n",sep=""))
      do.call(add_to_process,list(p_cmd = COMPLETE_env$SHELL, p_args = c(fs::path_package("COMPLETE","exec","functions.sh"),"extract_genomic_regions",genome_path, gtf_path, gene_list, org,params_list$param_file), logfile=paste(params_list$TEMP_PATH,"/",org,".log",sep=""), params_list=params_list, verbose=verbose))
      ##return(do.call(add_to_process,list(p_cmd = c("./extract_genomic_regions.sh"), p_args = c(genome_path, gtf_path, gene_list, org))))
      cat(print_toc(tictoc::toc(quiet = T, log = T)))
      return(c(data, source="user-provided"))
    }else{
      message(print_toc(tictoc::toc(quiet = T, log = T)))
      return(NULL)
    }
  }else{
    cat(print_toc(tictoc::toc(quiet = T, log = T)))
    return(c(data, source="user-provided"))
  }
}

#' Index FASTA IDs
#'
#' Convert longer FASTA IDs into short indices. Indices are generated recursively for fasta fasta files within subfolders of the path provided.
#'
#' @param path Path with FASTA Files to index
#' @param index_out Output files with the indices of the format "file"[tab]"long_id"[tab]"index"
#' @export
index_FASTA_IDs <- function(path, index_out){
  if (!is.null(path) && !is.null(index_out)){
    processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"index_fastaIDs", index_out,path) ,spinner = T,stdout = "",stderr = "")
  }else{
    stop("Give path and output file for index")
  }
}

#' Deduplicate FASTA Sequence
#'
#' Merge, Make Unique or Delete FASTA sequences with duplicated names. IF the duplicate sequence names are CDS blocks, merge them. If the duplicate seq names are EXON blocks, make them unique. If duplicate seq names are sequence duplicates, delete them.
#'
#' @param fasta_path Path to FASTA File
#' @param duplicates.method merge/delete/make_unique. How to handle sequences with duplicate names?. For CDS/UTR blocks - merge (Concatenate sequences with duplicate names), for Exon blocks - make_unique (Sequence names/IDs are made unique), delete - for deleting all duplicate sequences (first seq is kept).
#' @param n_threads Number of Threads
#' @export
deduplicate_FASTA <- function(fasta_path, duplicates.method, n_threads=4) {
  seq_set <- Biostrings::readDNAStringSet(filepath = fasta_path,format = "fasta",use.names = T)
  seq_set <- split(seq_set,factor(names(seq_set)))
  seq_set <- dplyr::bind_rows(x = parallel::mclapply(seq_along(seq_set),function(x){
   if(stringi::stri_cmp_eq(duplicates.method,"merge")){

      merged_seq <- Biostrings::DNAString(gsub("[[:space:]]", "", paste(seq_set[[x]],collapse="")))
      merged_seq_name <- unique(names(seq_set[[x]]))
      return(data.frame(seq_name=merged_seq_name, seq=paste(merged_seq))) #return(data.frame(seq_name=merged_seq_name, seq=merged_seq)) 
    }else if(stringi::stri_cmp_eq(duplicates.method,"make_unique")){
      if(length(seq_set[[x]])==1){
        return(data.frame(seq_name=names(seq_set[[x]]), seq=paste(seq_set[[x]]))) #return(data.frame(seq_name=names(seq_set[[x]]), seq=paste(seq_set[[x]])))
      }
      uniq_seqs <- unique(paste(seq_set[[x]]))
      seq_match_val <- match(uniq_seqs,seq_set[[x]]) #match(x,uniq_seqs)

      if(any(!is.na(seq_match_val)) && length(unique(seq_match_val)) == 1){ #length(uniq_seqs) == 1
        uniq_seqs <- gsub("[[:space:]]", "", paste(uniq_seqs)) # Biostrings::DNAString(
        uniq_seq_name <- unique(names(uniq_seqs))
        names(uniq_seqs) <- uniq_seq_name
        return(data.frame(seq_name=uniq_seq_name, seq=paste(uniq_seqs))) #data.frame
      }else{
        uniq_seq_name <- names(seq_set[[x]])[seq_match_val]
        if(any(duplicated(uniq_seq_name))){
          uniq_seq_name <- paste(paste("block",1:length(uniq_seq_name),sep=""),uniq_seq_name,sep=".")
        }
        #uniq_seqs <- list(uniq_seqs)
        names(uniq_seqs) <- uniq_seq_name
        return(data.frame(seq_name=uniq_seq_name, seq=paste(uniq_seqs))) #data.frame
      }
    }else if(stringi::stri_cmp_eq(duplicates.method,"delete")){
      #x <- x[!which(duplicated(paste(x)))]
      #x <- x[!which(duplicated(names(x)))]
      x <- x[!intersect(which(duplicated(paste(seq_set[[x]]))),which(duplicated(names(seq_set[[x]]))))]
      if(nrow(seq_set[[x]]) > 0){
        return(data.frame(seq_name=names(seq_set[[x]]), seq=paste(seq_set[[x]]))) #data.frame
      }
    }
  }, mc.cores = n_threads,mc.silent = T))

  seq_set_tmp <- Biostrings::DNAStringSet(seq_set[,2], use.names = F)
  names(seq_set_tmp) <- seq_set[,1]
  seq_set <- seq_set_tmp
  return(seq_set)
}

#' Label Sequence IDs of A FASTA File (Refer ?COMPLETE_PIPELINE_DESIGN about COMPLETE.format.IDs)
#'
#' Convert short sequence IDs into longer COMPLETE format IDs
#'
#' @param fasta_path Path to FASTA File
#' @param org Name of the organism
#' @param gene Filename of gene list or Vector of gene names
#' @param odb_clusters ODB Clusters
#' @param duplicates.method merge/delete/make_unique. How to handle sequences with duplicate names?. For CDS/UTR blocks - merge (Concatenate sequences with duplicate names), for Exon blocks - make_unique (Sequence names/IDs are made unique), delete - for deleting all duplicate sequences (first seq is kept).
#' @param params_list Output of load_params()
#' @export
label_sequenceIDs <- function(fasta_path, org, gene, odb_clusters, duplicates.method, params_list) {

  if(is.null(duplicates.method) || !grepl(pattern = c("merge|delete|make_unique"), x = duplicates.method,ignore.case = T)){
    stop(paste("duplicates.method must be one of merge|delete|make_unique"))
  }

  seq_set <- deduplicate_FASTA(fasta_path=fasta_path, duplicates.method=duplicates.method, n_threads=params_list$numWorkers)

  split_seq_names <- stringi::stri_split(str = names(seq_set), fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)
  if( ncol(split_seq_names) == 1 ){ #length(COMPLETE_env$FORMAT_ID_INDEX)
    #print(paste(names(seq_set),org,x,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)) #DEBUG
    #print(names(seq_set)) #DEBUG
    names(seq_set) <- paste(names(seq_set),org,gene,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)
    Biostrings::writeXStringSet(x = seq_set,filepath = fasta_path,append = F,format = "fasta")
    #return(paste(names(seq_set),org,x,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)) #DEBUG
  }else if( ncol(split_seq_names) > length(COMPLETE_env$FORMAT_ID_INDEX) || ncol(split_seq_names) < length(COMPLETE_env$FORMAT_ID_INDEX) ){
    names(seq_set) <- paste(split_seq_names[,COMPLETE_env$FORMAT_ID_INDEX$TRANSCRIPT_ID] ,org,gene,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM) #stringi::stri_split(str = split_seq_names, fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)[,1]
    Biostrings::writeXStringSet(x = seq_set,filepath = fasta_path,append = F,format = "fasta")
  }else if(ncol(split_seq_names) != length(COMPLETE_env$FORMAT_ID_INDEX)){
    message(paste("Error : Check sequence ID format of", y))
  }
}

#' Label Sequence IDs of all FASTA Files in a folder (Refer ?COMPLETE_PIPELINE_DESIGN about COMPLETE.format.IDs)
#'
#' Convert short sequence IDs into longer COMPLETE format IDs
#'
#' @param fasta_path Path of FASTA Files
#' @param org Name of the organism
#' @param gene_list Filename of gene list or Vector of gene names
#' @param odb_gene_map Filename of OrthoDB Ortholog Clusters mapped to gene names  (Optional). If not provided, sequences cluster ID/name "ungrouped" is used
#' @param duplicates.method merge/delete/make_unique. How to handle sequences with duplicate names?. For CDS/UTR blocks - merge (Concatenate sequences with duplicate names), for Exon blocks - make_unique (Sequence names/IDs are made unique), delete - for deleting all duplicate sequences (first seq is kept).
#' @param params_list Output of load_params()
#' @export
label_FASTA_files <- function(fasta_path,org,gene_list,odb_gene_map=NULL,params_list, duplicates.method){

  if(is.null(duplicates.method) || !grepl(pattern = c("merge|delete|make_unique"), x = duplicates.method,ignore.case = T)){
    stop(paste("duplicates.method must be one of merge|delete|make_unique"))
  }

  tictoc::tic(msg = "Labelling Sequence IDs...")
  if(!is.null(odb_gene_map)){
    if(file.exists(odb_gene_map) && file.info(odb_gene_map)$size > 0){
      odb_gene_map <- read.table(file = odb_gene_map,header = F,quote = "",sep = "\t")
      #local ortho_cluster=$(grep -w $gene_name $odb_clusters | awk -F'\t' '{if (length(c) == 0){c=$1;}else{c=c","$1;}}END{print c}')
    }else{
      warning(paste(odb_gene_map,"does not exist!"))
      odb_gene_map <- NULL
    }
  }
  if (length(gene_list) == 1 && file.exists(gene_list)) {
    genes <- factor(scan(gene_list, character(),quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
  }

  furrr::future_map(genes, function(x){ #parallel::mclapply
    safe_gene <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(x))
    file_path <- paste(fasta_path,safe_gene,sep="")
    if(!is.null(odb_gene_map)){
      odb_clusters <- paste(odb_gene_map[grep(pattern = x,x = odb_gene_map[,2],ignore.case = T,value = F),1],collapse = ",")
    }else{
      odb_clusters <- "ungrouped"
    }
    furrr::future_map(list.files(path = fasta_path,pattern = safe_gene,full.names = T), function(y){ #lapply
      label_sequenceIDs(fasta_path = y, org = org, gene = x, odb_clusters = odb_clusters, params_list = params_list, duplicates.method = duplicates.method)

      # Biostrings::writeXStringSet(x = seq_set,filepath = y,append = F,format = "fasta")
    }, .options = furrr::furrr_options(seed=T, scheduling=params_list$numWorkers))

  },.options = furrr::furrr_options(seed=T, scheduling=params_list$numWorkers))#, mc.cores = params_list$numWorkers,mc.silent = T)
  cat(print_toc(tictoc::toc(quiet = T)))
}

#' Internal Function - Merge and Format OrthoDB Flat Files
#'
#' This step is essential for speeding up the extraction process. The gene information (Gene IDs and Gene Names) are merged with Ortholog Groups (Cluster IDs and Gene IDs) and the format is converted into a Tab Delimited file wrapped over a CSV of Gene IDs and Gene Names. This is performed because one Ortholog Group/Cluster can encompass multiple genes (across organisms). The final format looks like this [cluster1][tab][geneID1,geneID2..geneIDN][tab][gene_name1,gene_name2..gene_nameN].
#'
#' @note Required Flat files from OrthoDB are \*_OG2genes.tab.gz,\*_genes.tab.gz and \*_species.tab.gz. Output files are stored in paste(odb_prefix,"_OGgenes_fixed.tab.gz",sep="") and paste(odb_prefix,"_OGgenes_fixed_user.tab.gz",sep="")
#'
#' @param odb_prefix Prefix to the Flat Files from OrthoDB (Prefix of \*_OG2genes.tab.gz,\*_genes.tab.gz,\*_species.tab.gz)
#' @param quick.check Only check if files exist? (TRUE/FALSE). FALSE - When running the pipeline with a new list of genes
#' @param n_threads Number of threads
#' @param gene_list File name of gene list
#' @return TRUE if output files exist, FALSE otherwise
merge_OG2genes_OrthoDB <- function(odb_prefix,quick.check=T,n_threads=tryCatch(parallel::detectCores(all.tests = T, logical = T), error=function(cond){return(2)}),gene_list){
  tictoc::tic(msg = paste("Transforming ODB Files..."))
  if(!quick.check){
    processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"merge_OG2genes_OrthoDB",odb_prefix,!quick.check,n_threads,gene_list ) ,spinner = T,stdout = "",stderr = "")
  }
  if( (file.exists(paste(odb_prefix,"_OGgenes_fixed.tab.gz",sep="")) && file.info(paste(odb_prefix,"_OGgenes_fixed.tab.gz",sep=""))$size > 0) || (file.exists(paste(odb_prefix,"_OGgenes_fixed_user.tab.gz",sep=""))&& file.info(paste(odb_prefix,"_OGgenes_fixed_user.tab.gz",sep=""))$size > 0) ){
    cat(print_toc(tictoc::toc(quiet = T)))
    return(TRUE)
  }else{
    cat(print_toc(tictoc::toc(quiet = T)))
    return(FALSE)
  }

}

#' (1) - Extracts Sequences for Protein Coding Transcripts from Organisms
#'
#' This is the main function which calls all the other functions and performs and end-end execution of data extraction part of the pipeline. It requires a filename of a formatted parameter file and a gene list (check the github repo for an example) or fs::path_package("COMPLETE","pkg_data","parameters.txt").
#'
#' @examples
#'     COMPLETE::EXTRACT_DATA(params_list = fs::path_package("COMPLETE","pkg_data","parameters.txt"), gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"), user_data = fs::path_package("COMPLETE","pkg_data", "user_data.txt"), only.user.data = F )
#'     COMPLETE::EXTRACT_DATA(params_list = fs::path_package("COMPLETE","pkg_data","parameters.txt"), gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"), user_data = NULL, only.user.data = F )
#'     COMPLETE::EXTRACT_DATA(params_list = fs::path_package("COMPLETE","pkg_data","parameters.txt"), gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"), user_data = NULL, only.user.data = T )
#'     
#' @note If samtools/bedtools are not available in $PATH user data (genomes & GTFs) is not processed (unless "-" is used where the organism is looked-up in BIOMART using biomaRt). Files are still downloaded and saved in params_list$GENOMES_PATH and params_list$ANNOS_PATH
#'
#' @param params_list Filename of a formatted parameter file (check the github repo for an example) or Output of load_params().
#' @param gene_list Vector or File with a list of genes to extract data for(check the github repo for an example)
#' @param user_data File name or table with user-specified organisms(genomes,GTFs). File must be in CSV format and should not contain header and column names are not required for the table. Check system.file("exec", "pkg_data", "user_data.txt", mustWork = T ,package = "COMPLETE") for an example user-data file.
#' @param only.user.data Process only user data and not process all available organisms in Ensembl? (TRUE/FALSE). Default FALSE
#' @export
EXTRACT_DATA <- function(params_list, gene_list, user_data=NULL, only.user.data=F){
  set.seed(123)

  if (!curl::has_internet()) {
    stop("Check if there is an internet connection")
  }
  
  if( only.user.data && is.null(user_data) ){
    stop("only.user.data=TRUE but user_data is not provided!")
  }

  if (stringi::stri_isempty(Sys.which("samtools")) || stringi::stri_isempty(Sys.which("bedtools")) ) { # || stringi::stri_isempty(Sys.which("parallel"))
    #print("Missing samtools/bedtools in $PATH...User data will not be processed!")
    warning("Missing samtools/bedtools in $PATH...User data & biomartr genomes will not be processed! (because the shell scripts depend on these programs)")
    message("Missing samtools/bedtools in $PATH...User data & biomartr genomes will not be processed! (because the shell scripts depend on these programs)")
    COMPLETE_env$SKIP_USER_DATA <- TRUE
  }else{
    COMPLETE_env$SKIP_USER_DATA <- FALSE
  }

  if(is.character(params_list)){
    if(!file.exists(params_list) || file.info(params_list)$size < 0){
      stop("ERROR: Parameters file is missing and is required\n")
    }
    loaded_PARAMS <- load_params(params_list)
  }else{
    if(any(grepl(x = class(params_list), pattern = "COMPLETE-options"))){
      loaded_PARAMS <- params_list
    }else{
      stop("Error: params_list not valid!")
    }
  }
  print(loaded_PARAMS)

  install_parallel()

  print(paste("MAX PROCESSES:",loaded_PARAMS$numWorkers))

  if (loaded_PARAMS$CLEAN_EXTRACT) {
    unlink(x = c(paste(loaded_PARAMS$OUT_PATH,"/available_orgs.txt",sep=""),
                 paste(loaded_PARAMS$OUT_PATH,"/unavailable_orgs.txt",sep=""),
                 paste(loaded_PARAMS$OUT_PATH,"/selected_ORGS.txt",sep=""),
                 paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""),
                 paste(loaded_PARAMS$OUT_PATH,"/all_gtf_stats.csv",sep="")),force = T,expand = T)
  }

  dir.create(loaded_PARAMS$OUT_PATH,showWarnings = F, recursive = T)
  unlink(loaded_PARAMS$TEMP_PATH, recursive = T,force = T,expand = T)
  dir.create(loaded_PARAMS$TEMP_PATH,showWarnings = F, recursive = T)
  #unlink(loaded_PARAMS$GROUPS_PATH, recursive = T,force = T,expand = T)
  #dir.create(loaded_PARAMS$GROUPS_PATH,showWarnings = F, recursive = T)
  dir.create(loaded_PARAMS$FASTA_OUT_PATH,showWarnings = F, recursive = T)
  dir.create(loaded_PARAMS$GENOMES_PATH,showWarnings = F, recursive = T)
  dir.create(loaded_PARAMS$ANNOS_PATH,showWarnings = F, recursive = T)

  #COMPLETE <<- new.env(parent=emptyenv())
  #COMPLETE_env$ENSEMBL_MART <- "ENSEMBL_MART_ENSEMBL"
  #COMPLETE_env$using.mart <- mart_connect(biomaRt::useMart,args=list(COMPLETE_env$ENSEMBL_MART)) #For biomaRt
  COMPLETE_env$org.meta.list <- mart_connect(biomaRt::listDatasets,args=list(mart=COMPLETE_env$using.mart)) #For biomaRt
  COMPLETE_env$org.meta <- mart_connect(biomartr::listGenomes,args=list(db = "ensembl", type = "all", details = T)) #For biomartr #db = tolower(GENOMES_SOURCE)

  unavailable_orgs <- c()
  all_orgs <- c()
  all_orgs <-  COMPLETE_env$org.meta$name
  if(file.exists(paste(loaded_PARAMS$OUT_PATH,"/unavailable_orgs.txt",sep="")) && file.info(paste(loaded_PARAMS$OUT_PATH,"/unavailable_orgs.txt",sep=""))$size > 0){
    unavailable_orgs <- factor(scan(paste(loaded_PARAMS$OUT_PATH,"/unavailable_orgs.txt",sep=""), character(), quiet = T))
    orgs_to_fetch <- COMPLETE_env$org.meta[which(is.na(match(COMPLETE_env$org.meta$name,unavailable_orgs))),]
  }else{
    orgs_to_fetch <- COMPLETE_env$org.meta
  }

  if( (!is.null(user_data) && is.character(user_data))){ # && !COMPLETE_env$SKIP_USER_DATA
    user_data <- read.csv(user_data,header = F)
  }

  if(!is.null(user_data)){ #&& !COMPLETE_env$SKIP_USER_DATA
    names(user_data) <- c("org","genome","gtf")
    all_orgs <- c(all_orgs,user_data$org)

    if(!is.null(unavailable_orgs) && length(unavailable_orgs) > 0){
      user_data <- user_data[which(is.na(match(user_data$org,unavailable_orgs))),]
    }
    if(!is.null(user_data)  && nrow(user_data) > 0){
      orgs_to_fetch <- orgs_to_fetch[which(is.na(match(orgs_to_fetch$name,user_data$org))),]
    }

  }else{
    messsage(paste("User Data not provided or is empty!"))
  }

  if (length(gene_list) == 1 && file.exists(gene_list)) {
    genes <- factor(scan(gene_list, character(), quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
    tmp_gene_list <- tempfile(pattern="genelist",tmpdir = params_list$TEMP_PATH)
    data.table::fwrite(x = list(gene_list),file = tmp_gene_list ,quote = F,row.names = F,col.names = F, nThread = loaded_PARAMS$numWorkers)
    gene_list <- tmp_gene_list
  }

  if(COMPLETE_env$USE_ORTHODB){
    if(!merge_OG2genes_OrthoDB(odb_prefix = loaded_PARAMS$ORTHODB_PREFIX,quick.check = !loaded_PARAMS$CLEAN_EXTRACT,n_threads = loaded_PARAMS$numWorkers,gene_list = gene_list)){
      merge_OG2genes_OrthoDB(odb_prefix = loaded_PARAMS$ORTHODB_PREFIX,quick.check = F,n_threads = loaded_PARAMS$numWorkers,gene_list = gene_list)
      }
  }

  cat("Checking and downloading transcripts...\n")

  cat(paste("User Data Logs saved in : ",loaded_PARAMS$TEMP_PATH,"/*.log\n\n", sep=""))
  #message(paste("Logs saved in : ",loaded_PARAMS$TEMP_PATH,"/*.log", sep = ""))

  user_saved_meta <- c()
  saved_meta <- c()

  tictoc::tic(msg = "Total Extraction Time :")

  if(!COMPLETE_env$SKIP_USER_DATA && !is.null(user_data)){ #(loaded_PARAMS$DATA_SOURCE=="both" || loaded_PARAMS$DATA_SOURCE=="user") ){
    if(nrow(user_data)!=0 && !is.null(user_data)){
      #print(user_data) #DEBUG
      user_saved_meta <- apply(user_data, MARGIN = 1, function(x){
        #user_proc <- fetch_FASTA_user(x)
        #user_proc$wait()
        #tictoc::tic.clear();
        #print(x)
        return( tryCatch({
          fetch_FASTA_user(data = x,params_list = loaded_PARAMS, gene_list = genes);
          #cat(paste("DONE :", x["org"],"\n"));
        },error=function(cond){
          message(cond)
          return(NULL)
        })
        ) } )
    }
  }else{
    message("User Data is skipped!")
  }

  if( !only.user.data ){ #loaded_PARAMS$DATA_SOURCE=="both" || loaded_PARAMS$DATA_SOURCE!="user"){
    #parallel::mclapply(org.meta$name, fetch_genome_db, mc.cores = 1)
    #lapply(org.meta$name, FUN = function(x){
    saved_meta <- apply(orgs_to_fetch, MARGIN = 1, FUN = function(x){
      #tictoc::tic.clear();
      return( tryCatch({
        fetch_FASTA(org_row = x, params_list = loaded_PARAMS, gene_list = genes);
        #cat(paste("DONE :", x["name"],"\n"));
      },error=function(cond){
        message(cond)
        return(NULL)
      })
      ) } )
  }

  try(parallel::mclapply(COMPLETE_env$process_list, function(x){ #parallel::mclapply
    if(!is.null(x) && x$is_alive()){
      x$wait(timeout=-1)
    }
  }, mc.cores =  loaded_PARAMS$numWorkers))

  final_toc_print <- tictoc::toc(quiet = T, log = T)
  while(!is.null(final_toc_print)){
    cat(print_toc(final_toc_print))
    final_toc_print <- tictoc::toc(quiet = T, log = T)
  }

  #save(saved_meta, file ="saved_meta.RData")

  #write.table(x = bind_rows(saved_meta[[2]]),file = paste(loaded_PARAMS$OUT_PATH,"/org_meta.txt",sep=""), quote = F,sep=",", row.names = F,na = "-")
  #write.table(x = t(saved_meta[[1]]),file = paste(loaded_PARAMS$OUT_PATH,"/org_meta.txt",sep=""), quote = F,sep=",", row.names = F,col.names = F,append = T,na = "-")
  if(!COMPLETE_env$SKIP_USER_DATA && !is.null(user_saved_meta)){
    user_saved_meta[sapply(user_saved_meta, is.null)] <- NULL
    user_saved_meta <- data.frame(t(user_saved_meta))
    colnames(user_saved_meta) <- c("org","genome","gtf")
    data.table::fwrite(x = list(user_saved_meta),file = paste(loaded_PARAMS$OUT_PATH,"/org_meta.txt",sep=""), quote = F,sep=",", row.names = F,col.names = T,append = T,na = "-", nThread = loaded_PARAMS$numWorkers)
  }

  #print(saved_meta) #DEBUG
  #save(saved_meta, file ="saved_meta.RData") #DEBUG
  if(!is.null(saved_meta) && length(saved_meta) > 0){
    saved_meta[sapply(saved_meta, is.null)] <- NULL
    saved_meta <- t(saved_meta) #dplyr::bind_rows(t(saved_meta))
    data.table::fwrite(x = list(saved_meta),file = paste(loaded_PARAMS$OUT_PATH,"/org_meta.txt",sep=""), quote = F,sep=",", row.names = F,col.names = T,append = T,na = "-", nThread = loaded_PARAMS$numWorkers)
  }
  #parallel::mclapply(saved_meta, function(x){
  #  write.table(x = x,file = paste(loaded_PARAMS$OUT_PATH,"/org_meta.txt",sep=""), quote = F,sep=",", row.names = F,col.names = T,append = T,na = "-")
  #}, mc.cores = loaded_PARAMS$numWorkers, mc.silent = T, mc.preschedule = T)

  #DELETING EMPTY DIRECTORIES - THESE ORGANISMS COULD NOT BE FETCHED
  parallel::mclapply(list.dirs(path = loaded_PARAMS$FASTA_OUT_PATH, full.names=TRUE,recursive = F), function(x) { #list.files(path = loaded_PARAMS$FASTA_OUT_PATH,include.dirs=TRUE, full.names=TRUE)
    fi <- file.info(x)
    if (fi$isdir) {
      f <- list.files(x, all.files=TRUE, recursive=TRUE, full.names=TRUE)
      sz <- sum(file.info(f)$size)

      #as precaution, print to make sure before using unlink(x, TRUE)
      if (sz==0L) unlink(x,recursive = T, force = T, expand =T) #print(x)
    }
  }, mc.cores = loaded_PARAMS$numWorkers, mc.silent = T, mc.preschedule = T)

  tictoc::tic(msg= "Coercing metdata from available organisms ...")

  missing_genes_list <- parallel::mclapply(list.files(path = paste(loaded_PARAMS$OUT_PATH,"/genes/",sep=""),include.dirs=TRUE, full.names=TRUE),function(x){
    if(file.exists(paste(x,"/MISSING_GENES",sep="")) && file.info(paste(x,"/MISSING_GENES",sep=""))$size > 0 ){
      return(scan(paste(x,"/MISSING_GENES",sep=""), character(), quiet = T))
    }
  }, mc.cores =  loaded_PARAMS$numWorkers)
  missing_genes <- tolower(purrr::reduce(missing_genes_list, unique))
  available_genes_list <- parallel::mclapply(list.files(path = paste(loaded_PARAMS$OUT_PATH,"/genes/",sep=""),include.dirs=TRUE, full.names=TRUE),function(x){
    if(file.exists(paste(x,"/AVAILABLE_GENES",sep="")) && file.info(paste(x,"/AVAILABLE_GENES",sep=""))$size > 0 ){
      return(scan(paste(x,"/AVAILABLE_GENES",sep=""), character(), quiet = T))
    }
  }, mc.cores =  loaded_PARAMS$numWorkers)
  available_genes <- tolower(purrr::reduce(available_genes_list, unique))

  available_orgs <- list.dirs(path= loaded_PARAMS$FASTA_OUT_PATH, full.names = F,recursive = F) #factor(scan(paste(loaded_PARAMS$OUT_PATH,"/available_orgs.txt
  data.table::fwrite(x = list(available_orgs) ,file = paste(loaded_PARAMS$OUT_PATH,"/available_orgs.txt",sep=""), quote = F, row.names = F,col.names = F, nThread = loaded_PARAMS$numWorkers)

  if(length(available_orgs) > 0){ #file.exists(paste(loaded_PARAMS$OUT_PATH,"/available_orgs.txt",sep="")) && file.info(paste(loaded_PARAMS$OUT_PATH,"/available_orgs.txt",sep=""))$size > 0
    #available_orgs <- list.dirs(path= loaded_PARAMS$FASTA_OUT_PATH, full.names = F,recursive = F) #factor(scan(paste(loaded_PARAMS$OUT_PATH,"/available_orgs.txt",sep=""), character(), quiet = T))
    unavailable_orgs <- all_orgs[which(is.na(match(all_orgs,available_orgs)))]
    cat(paste("Unavailable Organisms :",paste(unavailable_orgs, collapse = ",")))
  }else{
    unavailable_orgs <- all_orgs
    stop("No organisms were available!. Rety with other options or a different gene list.")
  }

  data.table::fwrite(x = list(unavailable_orgs),file = paste(loaded_PARAMS$OUT_PATH,"/unavailable_orgs.txt",sep=""), quote = F, row.names = F,col.names = F, nThread = loaded_PARAMS$numWorkers)

  #if(loaded_PARAMS$CLEAN_EXTRACT || (!file.exists("files/selected_ORGS.txt") && is.na(file.info("files/selected_ORGS.txt")$size)) ){
  data.table::fwrite(x = list(list.files(path = loaded_PARAMS$FASTA_OUT_PATH,include.dirs=TRUE, full.names=F)),file = paste(loaded_PARAMS$OUT_PATH,"/selected_ORGS.txt",sep=""), quote = F, row.names = F, col.names = F, nThread = loaded_PARAMS$numWorkers)
  #}
  selected_orgs <-  factor(scan( paste(loaded_PARAMS$OUT_PATH,"/selected_ORGS.txt",sep=""), character(), quiet = T))

  all_gtf_stats <- dplyr::bind_rows(parallel::mclapply(list.files(path = paste(loaded_PARAMS$OUT_PATH,"/genes/",sep=""),include.dirs=TRUE, full.names=TRUE),function(x){
    if(file.exists(paste(x,"/gtf_stats.csv",sep="")) && file.info(paste(x,"/gtf_stats.csv",sep=""))$size > 0 ){
      tmp_tab <- read.table(file = paste(x,"/gtf_stats.csv",sep=""),header = T,sep = ",",fill = T,na.strings = "",as.is = T, colClasses = "character") #unique
      tmp_tab <- tmp_tab[stats::complete.cases(tmp_tab),]
      if(length(tmp_tab) > 0){
        return(tmp_tab)
      }else{
        return(NULL)
      }
    }
  }, mc.cores =  loaded_PARAMS$numWorkers, mc.preschedule = T))
  write.table(x = list(all_gtf_stats),file = paste(loaded_PARAMS$OUT_PATH,"/all_gtf_stats.csv",sep=""),sep = ",", quote = F, row.names = F,col.names = T,na = "-") #, nThread = loaded_PARAMS$numWorkers)

  # time ./find_orthologs.sh files/selected_ORGS.txt $1 #100 ##This also selects the transcripts
  # time ./align_seqs.sh $1
  # time ./predict_structures.sh $1
  # rm $TEMP_PATH/*
  # Rscript gene_stats.R >> files/stats.txt

  cat(print_toc(tictoc::toc(quiet = T, log = T)))

  #tictoc::tic.clearlog()
  cat(paste("Extraction Time Log:\n"))
  cat(paste(tictoc::tic.log(),collapse = "\n"))

  #sessionInfo()
}

#' Design of R-COMPLETE
#'
#' The pipeline uses R and BASH. BASH functions are invoked through R.
#' BASH functions are stored in fs::path_package("COMPLETE","exec","functions.sh")
#'
#' + REQUIRES:
#' 
#'      - Internet Connection
#'      - Referenced R packages
#'      - Linux with BASH ($SHELL must be set or /bin/bash must exist)
#'      - Parameters File (fs::path_package("COMPLETE","pkg_data","parameters.txt"))
#'      - GNU parallel (in $PATH - BASH functions)
#'      - Samtools (in $PATH - BASH functions)
#'      - Bedtools (in $PATH - BASH functions)
#'      - OrthoDB (ODB) Flat Files (>= v10.1) (Pipeline is tested with ODB v10.1) 
#'          - odb10v1_species.tab.gz - Ortho DB organism ids based on NCBI taxonomy ids (mostly species level) 
#'          - odb10v1_genes.tab.gz  -Ortho DB genes with some info 
#'          - odb10v1_OG2genes.tab.gz - OGs to genes correspondence 
#'          (OR)
#'          - odb10v1_OGgenes_fixed.tab.gz - Merged & Transformed ODB file (Done within pipeline)
#'          - odb10v1_OGgenes_fixed_user.tab.gz - Merged & Transformed ODB file BASED on user gene list (Done within pipeline)
#'
#' + PARAMETERS : 
#' 
#'           The pipeline takes a single parameter file. This design was chosen
#'             1) To expose as many options as possible to the end-user.
#'              2) The pipeline uses BASH to perform some operations (which is significantly faster than R)
#'          and the parameter file is shared between R and BASH.
#'    
#'     The file is of the format [param_id==value==comment] where param_id and value columns are CASE-SENSITIVE
#'     (because its unnecessarily hard to check and convert param types in BASH). A default/example file is in
#'     fs::path_package("COMPLETE","pkg_data","parameters.txt")
#'
#' + USER DATA : 
#' 
#'          Columns Org, genome, gtf
#'
#' + COMPLETE.format.ids : 
#' 
#'          - The Ordering of ID labels can be referred from COMPLETE_env$FORMAT_ID_INDEX
#'          - Sequences are labelled with the following long ID format of R-COMPLETE
#'          (specific to this pipeline and referred to as COMPLETE.format.ids)
#'          (seqID_delimiter & transcripID_delimiter set in parameters, "::" & "||" respectively in this context )
#'               >$transcript_id $transcripID_delimiter $transcript_region ($strand) $seqID_delimiter $seqID_delimiter $org_name $gene_name $seqID_delimiter $ortho_cluster
#'               >SOME_TRANSCRIPT||cds(+)::SOMEORG::RANDOMGENE::ORTHOLOG_CLUSTERS
#'               >ENSDART00000193157||cds(+)::danio_rerio::sulf1::18335at7898,51668at7742,360590at33208
#'
#' + FLOW :
#'
#'     1) EXTRACT_DATA() - 
#'            Extracts the transcript regions for Protein Coding Transcripts (provided in parameters, pipeline requires cds,5utr,3utr)
#'        from BIOMART and/or User provided genomes & GTFs. This functions uses biomaRt/biomartr for extracting data from BIOMART
#'        and BASH function extract_genomic_regions() for user provided data. Extraction priority/flow : User Data > biomaRt > biomartr
#'          - ODB Files are merged and transformed with BASH function merge_OG2genes_OrthoDB()
#'          - Orthologous genes are found for genes which are not present in the organism with BASH function check_OrthoDB()
#'          - Flank lengths are calculated from GTF data for missing UTRs (with variance correction, check ?calculate_gtf_stats)
#'          - FASTA Nucleotide Sequences for given TRANSCRIPT_REGIONS are fetched from BIOMART/Genome
#'
#'     2) FIND_ORTHOLOGS() -
#'
#' @author Vishvesh Karthik (MDC-Berlin), [vishveshkarthik@gmail.com]
#'
#' @usage
#' (1) EXTRACT_DATA()
#' (2) FIND_TRANSCRIPT_ORTHOLOGS()
#'
#' @seealso [COMPLETE::EXTRACT_DATA()], [COMPLETE::FIND_ORTHOLOGS()]
#' @md
COMPLETE_PIPELINE_DESIGN <- function(){
}
