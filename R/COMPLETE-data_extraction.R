#' Pretty Print tictoc::toc() output
#'
#' This function takes the output of tictoc::toc() and
#' pretty prints it to the console
#'
#' @param clk output of function tictoc::toc()
#' @return (string) A formatted string to print to stdout/console
print_toc <- function(clk){
  return(paste(tictoc::toc.outmsg(clk$tic,clk$toc,clk$msg),"\n\n", sep = ""))
}

#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kept.
#'
#' @param x (data.frame) - Data Frame with columns named as required
#' @param col_name (string) - Name of the column for which the value is to be extracted
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
#' @param org (string) - Name of the organism
#' @param org_ver (string) - Dummy, just for verbosity
#' @param db (string) - Dummy, just for verbosity
#' @return (string) Name of the dataset in biomaRt
#' @export
check_mart_dataset <- function(org, org_ver, db){

  if(isTRUE(is.null(COMPLETE_env$org.meta.list))){
    traceback(3)
    stop("Reload R-COMPLETE, marts for 'ensembl' are not initialized properly\n")
  }

  split_org <- stringi::stri_split(gsub(pattern = "[[:punct:]]|[[:space:]]",x = org, replacement = "_"),fixed="_",simplify = T)
  if(stringi::stri_isempty(split_org[,2])){
    traceback(3)
    stop(paste(org,"name not in proper format, eg danio_rerio\n"))
  }
  mart.dataset <- grep(x = COMPLETE_env$org.meta.list$dataset, pattern=stringr::regex(split_org[2],ignore_case = T),fixed=F, value = T)
  if(length(mart.dataset)==0){
    message(paste("Dataset not found in BIOMART for :",org, org_ver, db,"\nYou can provide the organism in user data\n"))
    # warning(paste("Dataset not found in BIOMART for :",org, org_ver, db,"\nYou can provide the organism in user data\n"))
    return(NULL)
  }else if (length(unique(mart.dataset))==1) {
    return(unique(mart.dataset))
  }else{
    f_name <- paste(tolower(substring(split_org[1],1,1)),split_org[2],sep="")
    mart.dataset <- grep(x = mart.dataset,pattern=paste("^",f_name,sep="") ,ignore.case = T,value = T) #grep(x = mart.dataset,pattern= ,ignore.case = T,value = T)
    if (length(unique(mart.dataset))==1) {
      return(unique(mart.dataset))
    }else{
      traceback(3)
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
#'     gtf_stats <- calculate_stats(org, gtf_data)
#'     gtf_stats <- dplyr::bind_rows(gtf_stats)
#'     names(gtf_stats)[grep(pattern="seqnames",names(gtf_stats))] <- "transcript_id"
#'
#' @param org (string) Name of the organism (Only for DEBUG)
#' @param gtf_data (data.frame) GTF data obtained from biomaRt for an organism
#' @param allow_strand (string) Only allow the specified strand. ("+","-",Default - "" (or) " " (or) "*")
#' @param n_threads (integer) Number of Threads
#' @return (data.frame) Transcript Statistics from GTF data
#' @export
calculate_stats <- function(org, gtf_data, allow_strand="", n_threads=tryCatch(parallel::detectCores(all.tests = T, logical = T), error=function(cond){return(2)})) {

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
    # print(str(gtf_data)) #DEBUG
    # print(head(gtf_data)) #DEBUG
    traceback(3)
    stop(paste("calculate_stats(): Missing columns for", org,"\nHave:", paste(names(gtf_data), collapse=","), "\nRequire :",paste(req_columns,collapse = ","),"\n\n"))
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
      if(any(is.null(gtf_gr), length(gtf_gr)==0, all(is.na(gtf_gr)))){ ##Probably the mrna is not in the required strand so we can safely discard it
        #message(paste("No",strandedness,"strand info (or) region of gtf missing for the gene : ", g_name))
        return()
      }
    }

    transcript_stats <- future.apply::future_lapply(base::split(gtf_x, as.factor(gtf_x$ensembl_transcript_id)), function(x){
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
#' @param p_cmd (string) Command to be executed
#' @param p_args (list) Arguments to be passed to the command
#' @param verbose (bool) Print status messages?
#' @param logfile (string) Redirect output (stdout & stderr) to this file
#' @param params_list (string/list/COMPLETE-options) Output of load_params()
#' @return (processx::process() object) Process ID from processx::new() which can be used for further monitoring of the process
#' @export
add_to_process <- function(p_cmd, p_args=list(), verbose=F, logfile=NULL, params_list){
  
  if (is.null(logfile)) logfile=""
  if (is.null(COMPLETE_env$process_list)) COMPLETE_env$process_list <<- list() # Initialize as list, not vector c()

  # --- FIX: Correctly identifying indices of dead processes ---
  if(length(COMPLETE_env$process_list) > 0){
    # Get INDICES of dead processes
    dead_indices <- which(sapply(COMPLETE_env$process_list, function(x) !x$is_alive()))
    
    if(length(dead_indices) > 0){
       COMPLETE_env$process_list[dead_indices] <- NULL 
    }
  }

  if (verbose) cat(paste("\nAdding process to list...(",length(COMPLETE_env$process_list),"):",p_cmd, " ",paste(p_args, collapse = " "),"\n",sep=""))

  # print(length(COMPLETE_env$process_list) >= params_list$numWorkers) #DEBUG
  # print(length(COMPLETE_env$process_list)) #DEBUG
  # print(params_list$numWorkers) #DEBUG
  # print(ps::ps_num_fds() >= COMPLETE_env$max_file_handles-1) #DEBUG
  # print(ps::ps_num_fds()) #DEBUG
  # print(COMPLETE_env$max_file_handles-1) #DEBUG

  # Check if Queue is Full
  if(length(COMPLETE_env$process_list) >= params_list$numWorkers || ps::ps_num_fds() >= COMPLETE_env$max_file_handles-1){ 
    if (verbose) {
      cat(paste("Process Q Full...Waiting for a process to end(",length(COMPLETE_env$process_list),")\n",sep=""))
    }
    
    # Wait for the oldest process (FIFO) to free up a slot
    # This acts as a blocking call until a slot opens
    if(length(COMPLETE_env$process_list) > 0){
       COMPLETE_env$process_list[[1]]$wait(timeout=-1)
      if (inherits(COMPLETE_env$process_list[[1]], "process")) {
        COMPLETE_env$process_list[[1]]$finalize() # Or explicitly close connections
      }
       # Clean up immediately after waiting
       if(!COMPLETE_env$process_list[[1]]$is_alive()){
         COMPLETE_env$process_list[[1]] <- NULL
       }
    }
  }

  # log_con <- file(logfile, open = "a")
  # Spawn new process
  proc <- processx::process$new(command=p_cmd, args = p_args, cleanup = T, cleanup_tree = T, supervise = TRUE, stdout = logfile, stderr = logfile)
  COMPLETE_env$process_list <- append(COMPLETE_env$process_list, proc)
  # on.exit(close(log_con))
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
#' @param fasta_path (string) Path of folder to check
#' @param org (string) Name of the organism (format important, eg. "danio_rerio")
#' @param genes (vector/string) Filename/Vector of genes to check for
#' @param verbose (bool) Print check result messages?
#' @param params_list (string/list/COMPLETE-options) Output of load_params()
#' @return (bool) TRUE if check passed, FALSE if check failed
check_files <-function(fasta_path,org,genes, verbose=T, params_list){
  if(isTRUE(stringi::stri_isempty(org) || stringi::stri_isempty(fasta_path))){
    stop(paste("[check_files()] Error: org and fasta_path cannot be empty:",org, fasta_path))
  }

  if(length(genes) == 1 && file.exists(genes)){
    genes <- gsub('[[:punct:]]+','_', factor(scan(genes, character(), quiet = T)))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(genes))
  }

tryCatch({
    if(all(dir.exists(fasta_path), length(dir(fasta_path,all.files = F)) > 0, dir.exists(file.path(params_list$OUT_PATH,"genes",org,"")))){ 
      missing_genes <- c()
      available_genes <- c()

      odb_genes <- c()
      if(file.exists(file.path(params_list$OUT_PATH,"genes",org,"odb.list"))){
        odb_genes <- tolower(gsub('[[:punct:]]+','_', factor(scan(file.path(params_list$OUT_PATH,"genes",org,"odb.list"), character(), quiet = T))))
      }
      if(file.exists(file.path(params_list$OUT_PATH,"genes",org,"MISSING_GENES"))){
        missing_genes <- tolower(gsub('[[:punct:]]+','_', factor(scan(file.path(params_list$OUT_PATH,"genes",org,"MISSING_GENES"), character(), quiet = T))))
      }
      if(file.exists(file.path(params_list$OUT_PATH,"genes",org,"AVAILABLE_GENES"))){
        available_genes <- tolower(gsub('[[:punct:]]+','_', factor(scan(file.path(params_list$OUT_PATH,"genes",org,"AVAILABLE_GENES"), character(), quiet = T))))
      }
      files_in_dir <- tolower(unique(gsub('[[:punct:]]+','_',sapply(list.files(fasta_path,no.. = T,recursive = F), FUN=function(x){stringi::stri_split(str = x, fixed='.',simplify = T)[,1]}))))
      missing_genes <- tolower(missing_genes[is.na(match(missing_genes, available_genes))])

      if(all(!is.na(match(intersect(available_genes,genes),files_in_dir))) && all(is.na(match(intersect(available_genes,genes),missing_genes))) ){ 
        if(verbose){
          # print(files_in_dir)
          # print(odb_genes)
          # print(genes)
          # print(length(files_in_dir) - length(intersect(odb_genes,files_in_dir)))
          cat(paste("Org:",org,"Path:", fasta_path,"\n","Org:",org,"Genes in Dir(Matching + ODB):","(",length(files_in_dir) - length(intersect(odb_genes,files_in_dir)),"+",length(odb_genes) - length(setdiff(odb_genes,files_in_dir)),")",", Genes in Dir:",length(files_in_dir),", Unavailable:",length(genes[is.na(match(genes,files_in_dir))]),", User Genes:",length(genes)," : Check PASSED!\n"))
        }

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
  })

}

#' Get GTF data from BIOMART (using biomaRt)
#'
#' This function downloads GTF data from ensembl BIOMART using biomaRt. The attributes downloaded are
#' c("chromosome_name", "strand", "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name",
#' "start_position","end_position","exon_chrom_start","exon_chrom_end","transcript_version","transcript_length",
#' "5_utr_start","5_utr_end","transcript_start","transcription_start_site","transcript_end","cds_start","cds_end","3_utr_start","3_utr_end")
#'
#' @examples
#'     gtf_data <- get_gtf_mart(org, genes)
#'
#' @param org (string) Name of the organism (format important, eg. "danio_rerio")
#' @param org_ver (string) Version of the organism (ONLY for Debug)
#' @param db (string) Query Database (ONLY for Debug)
#' @param gene_list (vector/string) Vector or File with genes to fetch the GTF data for
#' @return (data.frame) GTF Data
#' @export
get_gtf_mart <- function(org,org_ver,db,gene_list){

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

  mart.dataset <- check_mart_dataset(org, org_ver, db)

  # print("Checking mart dataset...COMPLETE") #DEBUG
  if(any(stringi::stri_isempty(mart.dataset),is.null(mart.dataset))){
    return(NULL)
  }
  # print("Connecting to mart...") #DEBUG

  using.mart.data <- mart_connect(biomaRt::useMart,args = list(COMPLETE_env$ENSEMBL_MART, mart.dataset))

  # print("Connecting to mart...COMPLETE") #DEBUG

  # print(list("Listing mart attributes...", list(COMPLETE_env$ENSEMBL_MART, mart.dataset))) #DEBUG

  mart.attributes <- biomaRt::listAttributes(using.mart.data)

  # print("Listing mart attributes...COMPLETE") #DEBUG

  if( !all(!is.na(match(gtf_attributes,mart.attributes$name))) ){
    traceback(3)
    warning(paste("GTF Attributes were not available for :",org,"\n"))
    return(NULL)
  }

  if (length(gene_list) == 1 && file.exists(gene_list)) {
    genes <- factor(scan(gene_list, character(), quiet = T)) 
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
  }

  ##Splitting genes to the recomended number of queries for biomart through biomaRt
  #bm_gtf <- dplyr::bind_rows(lapply(split(genes, ceiling(seq_along(genes)/500)), function(split_genes){
  #  return(mart_connect(biomaRt::getBM,args=list(mart=using.mart.data,attributes=gtf_attributes,uniqueRows=T, useCache=F, filters = c("external_gene_name"), values = split_genes, curl=COMPLETE_env$curl_handle)))
  #}))

  # print("Connecting to mart 2...") #DEBUG

  bm_gtf <- mart_connect(biomaRt::getBM,args=list(mart=using.mart.data,attributes=gtf_attributes,uniqueRows=T, useCache=F, filters = c("external_gene_name"), values = genes)) #curl=COMPLETE_env$curl_handle

  # print("Connecting to mart 2...COMPLETE") #DEBUG

  bm_gtf <- dplyr::arrange(bm_gtf,chromosome_name,start_position)

  bm_gtf$strand[bm_gtf$strand==1] <- "+"
  #bm_gtf$strand[bm_gtf$strand== -1] <- "-"
  bm_gtf$strand[is.na(match(bm_gtf$strand,"+"))] <- "-"
  #bm_gtf$source_name <- "Ensembl"
  return(bm_gtf)
}

#' Get FASTA data from BIOMART (using biomaRt)
#'
#' This function downloads FASTA data from ensembl BIOMART using biomaRt. The Data Frame from calculate_stats() is required.
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
#' @param org (string) Name of the organism (format important, eg. "danio_rerio")
#' @param gtf_stats (data.frame) Data Frame from calculate_stats()
#' @param fasta_path (string) Output path for transcript FASTA files
#' @param params_list (string/list/COMPLETE-options) Output of load_params()
#' @param verbose (bool) Verbosity (Default: T)
#' @return (data.frame) Data Frame with transcripts for which FASTA was obtained
#' @export
fetch_FASTA_mart <- function(org,org_ver,db,gtf_stats, fasta_path, params_list, ...=...){
  # print(c("fetch_FASTA_mart(): ", org,gtf_stats, fasta_path)) #DEBUG
  dot_args <- list(...)
  verbose=T
  if("verbose" %in% names(dot_args)){
    verbose<-dot_args$verbose
  }

  req_columns <- c("gene_name","gene_id","transcript_id","total_exon_len","total_cds_len", "five_len","three_len", "exon_count","cds_count","g.exon_start","g.exon_end","transcript_length","five_flank","three_flank","g.transcript_start","g.transcript_end","chromosome_name","strand")

  
  if(any(grepl(x = class(gtf_stats), pattern = "list",ignore.case = T))){
    gtf_stats <- dplyr::bind_rows(gtf_stats)
  }

  if (any(grepl(pattern="seqnames",names(gtf_stats)))) {
    names(gtf_stats)[grep(pattern="seqnames",names(gtf_stats))] <- "transcript_id"
  }

  if (any(is.na(match(req_columns, names(gtf_stats))))) {
    traceback(3)
    stop(paste("fetch_FASTA_mart(): Missing columns for", org, "\nHave:", paste(names(gtf_stats), collapse=","), "\nRequire :",paste(req_columns,collapse = ","),"\n\n"))
  }

  if (!any(grepl(pattern="safe_gene_name",names(gtf_stats)))) {
    gtf_stats$safe_gene_name <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(gtf_stats$gene_name))
  }

  gtf_stats <- gtf_stats[gtf_stats$gene_name!=0 | gtf_stats$gene_id!=0,] # cleaning up

  fs::dir_create(fasta_path,recurse=T)
  seq_attributes <- c("5utr","coding","3utr") #"external_gene_name"
  mart.dataset <- check_mart_dataset(org,org_ver,db)
  if(any(stringi::stri_isempty(mart.dataset),is.null(mart.dataset))){
    return(NULL)
  }
  using.mart.data <- mart_connect(biomaRt::useMart,args = list(COMPLETE_env$ENSEMBL_MART, mart.dataset))

  if(!all(!is.na(match(seq_attributes, biomaRt::listAttributes(using.mart.data)[,c("name")])))){
    traceback(3)
    stop(paste("Atrributes :",paste(seq_attributes,collapse = ","), ": not availabe for :", org,"\n"))
  }
  check_rows_idx <- unlist(sapply(seq_along(1:nrow(gtf_stats)), FUN=function(x){
    #print(gtf_stats[x,]<=3)
    if(any(na.omit(gtf_stats[x,c("five_len","three_len","five_flank","three_flank")]<=3))){ #any(is.na(gtf_stats[x,])) , not checking for CDS==NA because I was able to obtain coding sequences even when cds_len==NA
      return(x)
    }
  }))

  if (length(check_rows_idx) > 0) {
    invalid_stats <- gtf_stats[check_rows_idx,]
    valid_transcripts <- unique(gtf_stats[-check_rows_idx,]$transcript_id)
  } else {
    invalid_stats <- gtf_stats[0,] # Empty dataframe with same columns
    valid_transcripts <- unique(gtf_stats$transcript_id)
  }
  bm_df <- c()

  if(length(valid_transcripts) > 0){
    mart_lock <- filelock::lock(COMPLETE_env$ensembl_lock, timeout = Inf)
    bm_seq <- parallel::mclapply(seq_attributes, function(x){return(mart_connect(biomaRt::getBM,args=list(mart=using.mart.data,attributes=c("ensembl_transcript_id",x),uniqueRows=T, useCache=F, filters = c("ensembl_transcript_id"), values = valid_transcripts )))},mc.cores = params_list$numWorkers) #,curl=COMPLETE_env$curl_handle
    # Sys.sleep(4.5)
    filelock::unlock(mart_lock)
    bm_seq <- future.apply::future_lapply(bm_seq, function(x){
      if(!any(grepl(x = class(x), fixed = T, pattern = "try-error"))){
        return(x)
      }
    })
    #print(bm_seq) #DEBUG
    bm_df <- purrr::reduce(bm_seq, dplyr::full_join, by = "ensembl_transcript_id")
    unavailable_transcripts <- unique( unlist(apply(bm_df,MARGIN = 1, FUN = function(x){
      row_check <- grepl("unavailable",x=x,ignore.case = T)
      if(all(row_check) || isTRUE(row_check[2]==T)){ #Removing transcripts which do not have any regions or CDS
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

  # bm_seqs <- data.frame()
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
      bm_flanks_cds <- dplyr::tibble(ensembl_transcript_id=flank_values$transcript_id, coding="Sequence unavailable")
      tryCatch({
          bm_flanks_cds <- mart_connect(biomaRt::getSequence,args = list(id=flank_values$transcript_id,type="ensembl_transcript_id",seqType="coding", mart=using.mart.data, useCache = F) )
        }, error=function(e){
          traceback(3)
          message(paste("CDS Error:",org,":", e,"\n\n"))
      })
      
      bm_flanks_five <- dplyr::tibble(ensembl_transcript_id=flank_values$transcript_id, `5utr`="Sequence unavailable", flanking_5utr=F)
      tryCatch({
          bm_flanks_five <- mart_connect(biomaRt::getSequence,args = list(id=flank_values$transcript_id,type="ensembl_transcript_id",seqType="coding_transcript_flank", upstream=flank_values$five_flank , mart=using.mart.data, useCache = F) )
          names(bm_flanks_five)[grep(pattern="coding_transcript_flank",names(bm_flanks_five))] <- "5utr"
          bm_flanks_five <- bm_flanks_five %>% mutate(flanking_5utr=T)
        }, error=function(e){
          message(paste("Upstream Flank Warning:", org,":", e))          
      })
      
      bm_flanks_three <- dplyr::tibble(ensembl_transcript_id=flank_values$transcript_id, `3utr`="Sequence unavailable", flanking_3utr=F)
      tryCatch({
          bm_flanks_three <- mart_connect(biomaRt::getSequence,args = list(id=flank_values$transcript_id,type="ensembl_transcript_id",seqType="coding_transcript_flank", downstream = flank_values$three_flank , mart=using.mart.data, useCache = F))
          names(bm_flanks_three)[grep(pattern="coding_transcript_flank",names(bm_flanks_three))] <- "3utr"
          bm_flanks_three <- bm_flanks_three %>% mutate(flanking_3utr=T)
        }, error=function(e){
          message(paste("Downstream Flank Warning:", org,":", e))
      })

      transcripts_with_data <- intersect(flank_values$transcript_id,unique(c(bm_flanks_five$ensembl_transcript_id,bm_flanks_three$ensembl_transcript_id, bm_flanks_cds$ensembl_transcript_id)))

      bm_seqs <- purrr::reduce(list(bm_flanks_five,bm_flanks_cds,bm_flanks_three), full_join, by="ensembl_transcript_id")

      unavailable_transcripts <- unique( unlist(apply(bm_seqs,MARGIN = 1, FUN = function(x){
        row_check <- grepl("unavailable",x=x,ignore.case = T)
        names(row_check) <- names(bm_seqs)
        #print(row_check)
        if(all(row_check) || isTRUE(row_check["coding"]==T)){ #Removing transcripts which do not have any regions or CDS
          return(as.character(x["ensembl_transcript_id"]))
        }
      })) )
      
      if(length(unavailable_transcripts)>0){
        flank_stats <- flank_stats[!is.na(match(flank_stats$transcript_id,unavailable_transcripts)),]
      }
      
      if(length(transcripts_with_data)>0){
        flank_stats <- flank_stats[!is.na(match(flank_stats$transcript_id,transcripts_with_data)),]
      }
      
      if (nrow(flank_stats) > 0) {
        invalid_stats <- flank_stats
      }
      
    }
    
  }else {
    # If flank_stats is empty, we simply log it and let the function proceed to process valid_transcripts
    message(paste("No short UTRs/Flanks detected for", org, ". Skipping flank sequence retrieval."))
    bm_seqs <- data.frame() 
  }
  
  
  if (!is.null(bm_df) && nrow(bm_df) > 0) {
    # If we have both standard and flanking sequences, join them
  
    if (nrow(bm_seqs) > 0) {
  
      bm_df <- purrr::reduce(list(bm_seqs, bm_df), dplyr::full_join, 
                            by = c("5utr", "ensembl_transcript_id", "coding", "3utr", "flanking_5utr", "flanking_3utr"))
  
    } 
    # If bm_seqs is empty, we just keep bm_df as is and proceed
  } else {
    # Only if BOTH are empty do we stop with an error
    if (nrow(bm_seqs) > 0) {
      bm_df <- bm_seqs
    } else {
      traceback(3)
      stop(paste("Error fetching FASTA for :", org, "\n"))
      # return(NULL)
    }
  }

  names(bm_df)[grep(pattern="ensembl_transcript_id",names(bm_df))] <- "transcript_id"
  names(bm_df)[grep(pattern="coding",names(bm_df))] <- "cds"

  bm_df <- inner_join(gtf_stats[,c("gene_name","safe_gene_name","transcript_id","strand")],bm_df, by = "transcript_id")

  bm_df <- unique(bm_df)

  final_unavailable_transcripts <- unique( unlist(apply(bm_df,MARGIN = 1, FUN = function(x){
    row_check <- grepl("unavailable",x=x,ignore.case = T)
    names(row_check) <- colnames(bm_df)
    #print(row_check)
    if(all(row_check[params_list$TRANSCRIPT_REGIONS]) || isTRUE(row_check["cds"]==T)){ #Removing transcripts which do not have any regions or CDS
      return(as.character(x["transcript_id"]))
    }
  })) )

  final_unavailable_transcripts <- unique( c(unavailable_transcripts, final_unavailable_transcripts) )
  if(length(final_unavailable_transcripts)>0){
    bm_df <- bm_df[which(is.na(match(bm_df$transcript_id,final_unavailable_transcripts))),]
    message(paste("(Some) Data missing for : ",org,": Stored in :",file.path(params_list$OUT_PATH,"genes",org,db,org_ver,"non_coding_biomart.csv",sep=""),". Maybe CDS or all regions are missing for the transcripts. This could happen for non-protein coding transcripts or retained introns"))

    non_coding_data <- unique(gtf_stats[which(!is.na(match(gtf_stats$transcript_id,final_unavailable_transcripts))),c("gene_name","transcript_id")])
    data.table::fwrite(list(non_coding_data),file = file.path(params_list$OUT_PATH,"genes",org,db,org_ver,"non_coding_biomart.csv",sep=""),quote = F,row.names = F,col.names = T,sep = ",",na = "-", nThread = params_list$numWorkers)
  }

  parallel::mclapply(base::split(bm_df[,c("transcript_id","gene_name","safe_gene_name",params_list$TRANSCRIPT_REGIONS,"flanking_5utr","flanking_3utr","strand")],as.factor(bm_df$gene_name)), function(x){
    future.apply::future_lapply(params_list$TRANSCRIPT_REGIONS, function(y){
      if (isTRUE(any(grepl("3utr|5utr",y,ignore.case = T)))) { ##Adding "_FLANK" for transcripts with UTR Flanks
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

#' Extract Transcripts from Genome
#' 
#' Requires a genome and an annotation. This function invokes external SHELL function extract_transcript_regions from fs::path_package("COMPLETE","exec","functions.sh") (just like the piepline for user data) and cannot be monitored
#' 
#' @param genome_path (string) Path to genome
#' @param gtf_path (string) Path to annotation
#' @param gene_list (string/vector) Filename of gene list
#' @param org_name (string) Name of the organism
#' @param accession (string) Accession of the organism (For filenaming)
#' @param taxid (integer) taxid of the organism (Dummy, returned back to the user)
#' @param org_ver (string) Version of the organism (For filenaming)
#' @param db (string) Query DB (For folder structure)
#' @param params_list (string/list/COMPLETE-options) Parameter list from load_params()
#' @param keep_data (bool) Keep fetched data after extraction?
#' @param verbose (bool) Verbose?
#' @param seed (integer) Seed value for future::future()
#' @return (future::future()) Future object that can be resolves to a named vector with organism details
#' @export
extract_transcript_regions <- function(genome_path, gtf_path, gene_list, org_name, accession, taxid, org_ver, db, params_list, keep_data, verbose, seed) {
  # print(c("extract_transcript_regions(): ",genome_path, gtf_path, gene_list, org_name)) #DEBUG
  
  if(!file.exists(genome_path) || !file.exists(gtf_path)){
    traceback(3)
    stop(paste("MISSING:",genome_path,gtf_path, "\n"))
  }

  return(future::future({
      
      #Spawn process
      extraction_proc <- do.call(add_to_process, list(
        p_cmd = COMPLETE_env$SHELL, 
        p_args = c(
          fs::path_package("COMPLETE", "exec", "functions.sh"),
          "extract_transcript_regions", genome_path, gtf_path, gene_list, 
          org_name, org_ver, db, params_list$param_file, 
          COMPLETE_env$parallel, dplyr::if_else(keep_data, "TRUE", "FALSE")
        ), 
        logfile = paste(params_list$TEMP_PATH, "/", org_name,"_",org_ver,"_",db, ".log", sep = ""), 
        params_list = params_list, 
        verbose = verbose
      ))
      
      gtf_stats_path <- file.path(params_list$OUT_PATH,"genes", org_name, "gtf_stats.csv")
      
      # 3. Use tryCatch inside the future to handle crashes
      result <- tryCatch({
        extraction_proc$wait(timeout = -1)
        if (inherits(extraction_proc, "process")) {
          extraction_proc$finalize() # Or explicitly close connections
        }
        status <- extraction_proc$get_exit_status()
        
        if (status == 0) {
          # SUCCESS: Return the character vector
          return(c(org=org_name, accession=accession, taxid=taxid, version=org_ver, db=db, genome=genome_path, gtf=gtf_path, source="r-biomartr"))
        } else {
          if(!file.exists(gtf_stats_path) || file.info(gtf_stats_path)$size <= 0){
            #Organism does not have any genes, append it to unavailable_orgs.txt
            message(paste0("[extract_transcript_regions()]",org_name,": Organism does not have any genes, append it to unavailable_orgs.txt:",status))
            message(paste("UNAVAILABLE:", org_name,org_ver,db, "\n\n"))
            cat(paste(org_name,org_ver,db,accession,sep="\t"), sep="\n", file = file.path(params_list$OUT_PATH,"unavailable_orgs.txt"), append = TRUE)
          }else{
            warning(paste("UNKNOWN ERROR (Check Logfile):", org_name,org_ver,db,accession, "\n\n"))
          }
          return(NULL)
        }
        
      }, error = function(err) {
          message(paste("[extract_transcript_regions()] ERROR OCCURRED:", err))
          
          if(!file.exists(gtf_stats_path) || file.info(gtf_stats_path)$size <= 0){
            #Organism does not have any genes, append it to unavailable_orgs.txt
            message(paste("[extract_transcript_regions()] Error:",org_name,org_ver,db,":",err))
            message(paste("[extract_transcript_regions()] UNAVAILABLE:", org_name,org_ver,db,accession, "\n\n"))
            cat(paste(org_name,org_ver,db,accession,sep="\t"), sep="\n", file = file.path(params_list$OUT_PATH,"unavailable_orgs.txt"), append = TRUE)
          }
        
        message(print_toc(tictoc::toc(quiet = T, log = T)))
        return(NULL) # This becomes the value of the future on failure
      })
      
      return(result)
    }, seed=seed)
    ) 
}
 
#' Rgb - read.gtf implementation
#'
#' Sourcing and using read.gtf script from the deprecated R-Genome Browser (Rgb) package. Credits to the original author(s), Sylvain Mareschal <maressyl@gmail.com>
#'
#' @author Sylvain Mareschal <maressyl@gmail.com>
#' @param file (string) Path to GTF
#' @param attr (string) Method to process GTF attributes. Accepted values are "split", "intact", "skip"
#' @param features (string) Vector with feature names to extract from GTF
#' @param quiet (bool) Print Messages? (Default: F)
#' @return (data.frame) Extracted GTF
#' @export
read.gtf <- function(file, attr=c("split", "intact", "skip"), features=NULL, quiet=FALSE) {
  # Checks
  attr <- match.arg(attr)
  
  if(!isTRUE(quiet)) message("File parsing ... ", appendLF=FALSE)
  
  # print(file) #DEBUG

  # Parsing
  columns <- list("seqname"=character(0), "source"=character(0), "feature"=character(0), "start"=integer(0), "end"=integer(0), "score"=double(0), "strand"=character(0), "frame"=integer(0), "attributes"=character(0))
  content <- scan(file=file, what=columns, sep="\t", dec=".", comment.char="#", na.strings=".", quote="\"", quiet=TRUE)
  names(content) <- names(columns)
  
  # Strand
  content$strand[ is.na(content$strand) ] <- "."
  content$strand[ content$strand == "?" ] <- NA
  content$strand <- factor(content$strand, levels=c("-","+","."))
  
  # As data.frame
  class(content) <- "data.frame"
  rownames(content) <- 1:length(content$seqname)
  
  # Feature filtering
  if(!is.null(features)) content <- content[ content$feature %in% features , , drop=FALSE ]
  
  if(!isTRUE(quiet)) message(nrow(content), " rows processed")
  
  # Attributes
  if(isTRUE(attr == "skip")) {
    # No attribute
    content$attributes <- NULL
  } else if(isTRUE(attr == "split")) {
    
    if(!isTRUE(quiet)) message("Attribute splitting ... ", appendLF=FALSE)
    
    # Split attributes of each row
    att <- strsplit(content$attributes, split=" *; *")
    
    # Vectorize all attributes, keeping a parallel row index for each
    attRows <- rep.int(1:length(att), times=unlist(future.apply::future_lapply(att, length)))
    att <- unlist(att)
    
    # Split all name-value pairs
    regex <- regexpr(pattern="^(?<id>[A-Za-z][A-Za-z0-9_]*).(?<value>.+)$", text=att, perl=TRUE)
    attNames <- substr(att, 1, attr(regex, "capture.length")[,"id"])
    attValues <- substr(att, attr(regex, "capture.start")[,"value"], nchar(att))
    
    if(!isTRUE(quiet)) message(length(attValues), " pairs processed")
    
    if(!isTRUE(quiet)) message("Attribute sorting ... ", appendLF=FALSE)
    
    # Initialize a storage matrix (character)
    allNames <- unique(attNames)
    attMtx <- matrix(as.character(NA), nrow=nrow(content), ncol=length(allNames), dimnames=list(NULL, allNames))
    
    # Fill the storage matrix
    attMtx[ cbind(attRows, match(attNames, allNames)) ] <- attValues
    
    if(!isTRUE(quiet)) message(ncol(attMtx), " tags found")
    
    if(!isTRUE(quiet)) message("Attribute binding ...")
    
    # Convert character matrix to typed data.frame
    attDf <- as.data.frame(attMtx, stringsAsFactors=FALSE)
    for(i in 1:ncol(attDf)) attDf[[i]] <- utils::type.convert(attDf[[i]], as.is=TRUE)
    
    # Append to content
    content$attributes <- NULL
    content <- cbind(content, attDf, stringsAsFactors=FALSE)
  }
  
  if(!isTRUE(quiet)) message("done")
  
  # Return
  return(content)
}

#' Split GTF attributes into induvidual columns
#'
#' @param gtf_data (data.frame) GTF Data
#' @param attribute (string) Attribute name
#' @return (string) Attribute value
get_gtf_attributes <- function(gtf_data,attribute){
  return(unlist(future.apply::future_lapply(strsplit(gtf_data$attributes, split="; "), function(x){
    attr_present <- grepl(x=x,pattern = attribute,ignore.case = T)
    if(any(attr_present)){
      #return(unlist(strsplit(x = x[which(attr_present)], split=" "))[2])
      return(sub(x=unlist(strsplit(x = x[which(attr_present)], split=" "))[2],replacement = "",pattern = "[[:punct:]]$"))
    }else{
      return("")
    }
  })))
}

#' Internal Function - Get FASTA data from BIOMART (using biomartr)
#'
#' This function downloads FASTA data from BIOMART using biomartr. The GENOME and GTF for the organism are downloaded and passed through the shell pipeline for extracting transcripts.
#' This function invokes external SHELL function extract_transcript_regions from fs::path_package("COMPLETE","exec","functions.sh") (just like the piepline for user data) and cannot be monitored
#'
#' @examples
#'     params_list <- load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt"))
#'     fetch_FASTA_biomartr(c(name="danio_rerio",version="106",accession="acc1",taxid="dummy1"),params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'
#' @param org_row (Named vector) Named vector with name of the organism (format important, eg. "danio_rerio") and other details eg, c(name="danio_rerio",version="106",accession="acc1",taxid="dummy1")
#' @param db (string) Any supported DB from biomartr::listGenomes(), Valid Inputs: "ensembl", "genbank", "refseq", "" (Default: "ensembl"). Note: If empty, data is retrieved from both "genbank" and "ensembl"
#' @param params_list (string/list/COMPLETE-options) Output of load_params()
#' @param gene_list (string/vector) Vector or File containing list of genes
#' @param data_types (vector) Data types to fetch, one or all of "gtf", "genome" (Default: c("gtf","genome"))
#' @param keep_data (bool) Keep downloaded genomes and GTF data? (Default: FALSE)
#' @param only_fetch (bool) Only fetch data?, without extraction (only genomes & GTF annotations from biomartr) (Default: FALSE)
#' @param type (string) Filter to organisms. Input of biomartr::listGenomes() ["all", "kingdom", "group", "subgroup"] (Default: "all")
#' @param subset (string/vector) Filter to organisms. Input of biomartr::listGenomes() [Differs between DBs and filter type, check manually - e.g,if type == "group|subgroup" then check biomartr::getGroups()] (Default: NULL)
#' @param verbose (bool) Verbosity (Default: T)
#' @return (Named vector) Organism details on successful extraction, NULL otherwise
fetch_FASTA_biomartr <- function(org_row, db="", params_list, gene_list, data_types=c("gtf","genome"), keep_data=F, ...){ 
  dot_args <- list(...)
  verbose=T
  if("verbose" %in% names(dot_args)){
    verbose<-dot_args$verbose
  }
  only_fetch=F
  if("only_fetch" %in% names(dot_args)){
    only_fetch<-dot_args$only_fetch
  }
  type="all"
  if("type" %in% names(dot_args)){
    type<-dot_args$type
  }
  subset=NULL
  if("subset" %in% names(dot_args)){
    subset<-dot_args$subset
  }
  seed <- 123
  if("seed" %in% names(dot_args)){
    seed <- dot_args$seed
  }
  unavailable_orgs <- data.frame()
  if("unavailable_orgs" %in% names(dot_args)){
    unavailable_orgs <- dot_args$unavailable_orgs
  }
  # print(c("fetch_FASTA_biomartr(): ", org_row, params_list, gene_list)) #DEBUG
  
  if(!any(grepl(x= data_types, pattern = "gtf|genome", ignore.case = T, perl = T))){
    stop("data_types can be one/all of gtf, genome\n")
  }

  # if(stringi::stri_isempty(db)){
  #   return(dplyr::full_join(fetch_FASTA_biomartr(org_row=org_row, db="ensembl", params_list=params_list, gene_list=gene_list, data_types=data_types, keep_data=keep_data, ...=...),
  #   fetch_FASTA_biomartr(org_row=org_row, db="genbank", params_list=params_list, gene_list=gene_list, data_types=data_types, keep_data=keep_data,...=...)))
  # }

  # print(org_row) #DEBUG
  
  if(isTRUE(stringi::stri_isempty(db)) || !isTRUE(any(grepl(x=db, pattern="genbank|refseq|ensembl")))){
    stop("db should be one of genbank|refseq|ensembl.\n")
  }
  
  org_row[["db"]] <- db
  
  # print(org_row) #DEBUG
  if(isTRUE(all(c("name","version","accession","taxid") %in% names(org_row)))){
    org_row[["name"]] <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org_row[["name"]],perl = T, replacement = "_")))
    org_row[["version"]] <- trimws(gsub(pattern = "[[:space:]]", x=org_row[["version"]],perl = T, replacement = "_"))
    org_row[["accession"]] <- trimws(gsub(pattern = "[[:space:]]", x=org_row[["accession"]],perl = T, replacement = "_"))
  }else if(isTRUE(any(c("genbank","refseq") %in% org_row[["db"]]))){
    if("organism_name" %in% names(org_row)){
      org_row[["organism_name"]] <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org_row[["organism_name"]],perl = T, replacement = "_")))
      org_row[["name"]] <- org_row[["organism_name"]]
    }
    if("asm_name" %in% names(org_row)){
      org_row[["asm_name"]] <- trimws(gsub(x = org_row[["asm_name"]], pattern = "[[:space:]]", replacement = "_"))
      org_row[["version"]] <- org_row[["asm_name"]]
    }
    if("assembly_accession" %in% names(org_row)){
      org_row[["assembly_accession"]] <- trimws(gsub(x = org_row[["assembly_accession"]], pattern = "[[:space:]]", replacement = "_"))
      org_row[["accession"]] <- org_row[["assembly_accession"]]
    }
    if("taxid" %in% names(org_row)){
      org_row[["taxid"]] <- trimws(gsub(x = org_row[["taxid"]], pattern = "[[:punct:]]|[[:space:]]", replacement = "_"))
    }
  }else if(isTRUE("ensembl" %in% org_row[["db"]])){
    if("name" %in% names(org_row)){
      org_row[["name"]] <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org_row[["name"]],perl = T, replacement = "_")))
    }
    if("release" %in% names(org_row)){
      org_row[["version"]] <- trimws(gsub(pattern = "[[:space:]]", x=org_row[["release"]],perl = T, replacement = "_"))
    }
    if("accession" %in% names(org_row)){
      org_row[["accession"]] <- trimws(gsub(pattern = "[[:space:]]", x=org_row[["accession"]],perl = T, replacement = "_"))
    }
    if("taxon_id" %in% names(org_row)){
      org_row[["taxid"]] <- trimws(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org_row[["taxon_id"]],perl = T, replacement = "_"))
    }
  }else{
    stop(paste("Org:", org_row,"does not contain name|version|accession|taxid"))
  }

  if(stringi::stri_isempty(org_row[["name"]]) || stringi::stri_isempty(org_row[["version"]]) || stringi::stri_isempty(org_row[["accession"]]) || stringi::stri_isempty(org_row[["taxid"]])){
    # print(org_row[["name"]])
    # print(org_row[["version"]])
    stop("[fetch_FASTA_biomartr()] Input must be a named list/vector with organism name and version. Accepted names: name, accession, taxid, version.\n")
  }

  if(isTRUE(nrow(unavailable_orgs) > 0)){
    is_match <- unavailable_orgs$org_name %in% org_row[["name"]] & 
      unavailable_orgs$org_ver %in% org_row[["version"]] & 
      unavailable_orgs$db %in% org_row[["db"]]
    # print(is_match)
    if(any(is_match)){
      message(paste("UNAVAILABLE:", org_row[["name"]],org_row[["version"]],org_row[["db"]], "\n\n"))
      return(NULL)
    }
  }
  
  tictoc::tic(msg=paste("[biomartr] Processed:",org_row[["name"]], org_row[["db"]], org_row[["version"]], org_row[["accession"]]))
  tmp_gene_list <- c()
  if(any(!is.na(match(COMPLETE_env$org.meta$name,org_row[["name"]])))){    
    fs::dir_create(file.path(params_list$GENOMES_PATH,org_row[["name"]]),recurse = T)
    fs::dir_create(file.path(params_list$ANNOS_PATH,org_row[["name"]]),recurse = T)
    genome_path<-file.path(params_list$GENOMES_PATH,org_row[["name"]],paste(org_row[["db"]],"_",org_row[["version"]],".fa.gz",sep = ""))
    gtf_path<-file.path(params_list$ANNOS_PATH,org_row[["name"]], paste(org_row[["db"]],"_",org_row[["version"]],".gtf.gz",sep = ""))
    gtf_path_tmp<-file.path(params_list$ANNOS_PATH,org_row[["name"]], paste(org_row[["db"]],"_",org_row[["version"]],".gtf",sep = ""))
    gff_path<-file.path(params_list$ANNOS_PATH,org_row[["name"]], paste(org_row[["db"]],"_",org_row[["version"]],".gff.gz",sep = ""))
    tmp_file_path <- tempfile() #paste0(gtf_path,".tmp")
    extracted_tmp_file <- paste0(tmp_file_path,".tmp") #gsub(".gz$", "", tmp_file_path)
    org_fasta_path <- file.path(params_list$FASTA_OUT_PATH ,org_row[["name"]], org_row[["db"]], org_row[["version"]])
    
    # print(paste(genome_path,gtf_path))
    
    if("ensembl" %in% org_row[["db"]]){ 
      
      if(any(grepl(x= data_types, pattern = "gtf", ignore.case = T, perl = T)) && any(!file.exists(gtf_path), file.info(gtf_path)$size <= 20, params_list$CLEAN_EXTRACT)){
        gtf_ori <- biomartr::getGTF(organism = org_row[["name"]], db=org_row[["db"]],path = params_list$ANNOS_PATH, mute_citation=T) 
        if(!is.logical(gtf_ori)){
          fs::file_move(tools::file_path_as_absolute(gtf_ori),gtf_path)
        }
      }
      
      if(any(grepl(x= data_types, pattern = "genome", ignore.case = T, perl = T)) && any(!file.exists(genome_path), file.info(genome_path)$size <= 20, params_list$CLEAN_EXTRACT)){
        genome_ori <- biomartr::getGenome(organism = org_row[["name"]], db=org_row[["db"]],path = params_list$GENOMES_PATH,reference = T,gunzip = F, mute_citation=T) 
        if(!is.logical(genome_ori)){
          fs::file_move(tools::file_path_as_absolute(genome_ori),genome_path)
        }
      }
      
    }else{

      is_match <- unavailable_orgs$org_name %in% org_row[["name"]] & 
        unavailable_orgs$org_ver %in% org_row[["version"]] & 
        unavailable_orgs$db %in% org_row[["db"]] &
        unavailable_orgs$accession %in% org_row[["accession"]]
      # print(is_match)
      if(any(is_match)){
        message(paste("UNAVAILABLE:", org_row[["name"]],org_row[["version"]],org_row[["db"]],org_row[["accession"]], "\n\n"))
        cat(paste(org_row[["name"]],org_row[["version"]],org_row[["db"]],org_row[["accession"]],sep="\t"), sep="\n", file = file.path(params_list$OUT_PATH,"unavailable_orgs.txt"), append = TRUE)
        return(NULL)
      }
      #if db is genbank then we can use the ftp path to fetch the genome and gtf because biomartr cannot seem to auto-fetch
      if(all(any(c("genbank","refseq") %in% org_row[["db"]]),"asm_name" %in% colnames(COMPLETE_env$org.meta))){
        clean_asm_name <- trimws(gsub(pattern = "[[:space:]]", x=org_row[["asm_name"]],perl = T, replacement = "_"))
        if(any(grepl(x= data_types, pattern = "gtf", ignore.case = T, perl = T)) && any(!file.exists(gtf_path), file.info(gtf_path)$size <= 20, params_list$CLEAN_EXTRACT)){
          ncbi_gtf <- paste0(org_row[["ftp_path"]],"/",org_row[["assembly_accession"]],"_",clean_asm_name,"_genomic.gtf.gz")
          ncbi_lock <- filelock::lock(COMPLETE_env$ncbi_lock, timeout = Inf)
          ret_code <- curl::curl_fetch_disk(URLencode(ncbi_gtf), path = tmp_file_path)
          filelock::unlock(ncbi_lock)
          # print(ret_code)
          if(isTRUE(ret_code$status_code!=200)){
            message(paste("Error: GTF RET CODE :",ret_code$status_code,org_row[["name"]],":",ncbi_gtf))
            #Trying GFF because GTF not found
            ncbi_gff <- paste0(org_row[["ftp_path"]],"/",org_row[["assembly_accession"]],"_",clean_asm_name,"_genomic.gff.gz")
            ncbi_lock <- filelock::lock(COMPLETE_env$ncbi_lock, timeout = Inf)
            ret_code <- curl::curl_fetch_disk(URLencode(ncbi_gff), path = tmp_file_path)
            filelock::unlock(ncbi_lock)
            if(isTRUE(ret_code$status_code!=200)){
              fs::file_delete(tmp_file_path)
              #Neither GFF or GTF found
              message(paste("Error",ret_code$status_code,": check URL(s) :",org_row[["name"]],":",ncbi_gtf,ncbi_gff))
              message(print_toc(tictoc::toc(quiet = T, log = T)))
              return(NULL)
            }
            #Got GFF, convert to GTF
            # rtracklayer::export(rtracklayer::import.gff3(gzfile(gff_path,open="r"), format="gff"),con=gzfile(gtf_path, open="w"), format="gtf")
            # fs::file_delete(gff_path)
            
            R.utils::gunzip(tmp_file_path, destname = extracted_tmp_file, remove = T)
            convert_proc <- do.call(add_to_process, list(
              p_cmd = "agat_convert_sp_gff2gtf.pl",
              p_args = c(
                "-v","4","--gtf_version","relax","--gff", extracted_tmp_file, "-o", gtf_path_tmp
              ), 
              logfile = paste(params_list$TEMP_PATH, "/", org_row[["name"]],"_",org_row[["version"]],"_",db, "_gffconversion.log", sep = ""), 
              params_list = params_list, 
              verbose = verbose
            ))
          }else{
            
            R.utils::gunzip(tmp_file_path, destname = extracted_tmp_file, remove = T)
            convert_proc <- do.call(add_to_process, list(
              p_cmd = "agat_convert_sp_gff2gtf.pl",
              p_args = c(
                "-v","4","--gtf_version","relax","--gtf", extracted_tmp_file, "-o", gtf_path_tmp
              ), 
              logfile = paste(params_list$TEMP_PATH, "/", org_row[["name"]],"_",org_row[["version"]],"_",db, "_gtfconversion.log", sep = ""), 
              params_list = params_list, 
              verbose = verbose
            ))
          }
          convert_proc$wait(timeout = -1)
          if (inherits(convert_proc, "process")) {
            convert_proc$finalize() # Or explicitly close connections
          }
          gzfile_con <- gzfile(gtf_path, "w")
          writeLines(readLines(gtf_path_tmp), gzfile_con)
          close(gzfile_con)
          fs::file_delete(gtf_path_tmp)
        }
        
        if(isTRUE(any(grepl(x= data_types, pattern = "genome", ignore.case = T, perl = T)) && any(!file.exists(genome_path), file.info(genome_path)$size <= 20, params_list$CLEAN_EXTRACT))){
          ncbi_genome <- paste0(org_row[["ftp_path"]],"/",org_row[["assembly_accession"]],"_",clean_asm_name,"_genomic.fna.gz")
          ncbi_lock <- filelock::lock(COMPLETE_env$ncbi_lock, timeout = Inf)
          ret_code <- curl::curl_fetch_disk(URLencode(ncbi_genome), path = genome_path)
          # print(ret_code)
          # Sys.sleep(4.5)
          filelock::unlock(ncbi_lock)
          if(isTRUE(ret_code$status_code!=200)){
            message(paste("Error",ret_code$status_code,": check genome URL :",org_row[["name"]],":",ncbi_genome))
            message(print_toc(tictoc::toc(quiet = T, log = T)))
            print(org_row) #DEBUG
            return(NULL)
          }
        }
        
      }else{
        message(paste("Organism not available :", org_row[["accession"]],org_row[["name"]],org_row[["version"]],"in",org_row[["db"]],"\n\n"))
        return(NULL)
      }
    }
    
    if(isTRUE(only_fetch)){
      return(future::future({
        return(c(org=org_row[["name"]],accession=org_row[["accession"]],taxid=org_row[["taxid"]], version=org_row[["version"]], db=org_row[["db"]], genome=genome_path, gtf=gtf_path, source="r-biomartr"))
      }, label=paste0("return fetch_FASTA_biomartr():", paste(org_row[["accession"]],org_row[["name"]], org_row[["version"]],org_row[["db"]],sep=":"))))
    }
    
    if(isTRUE(length(gene_list) == 1 && file.exists(gene_list))){
      genes <- factor(scan(gene_list, character(), quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
      genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
    }else{
      genes <- tolower(as.vector(gene_list))
      tmp_gene_list <- tempfile(pattern="genelist",tmpdir = params_list$TEMP_PATH)
      data.table::fwrite(x = list(gene_list),file = tmp_gene_list ,quote = F,row.names = F,col.names = F, nThread = params_list$numWorkers)
      gene_list <- tmp_gene_list
    }
    
    if(isTRUE(params_list$CLEAN_EXTRACT || !check_files(fasta_path = org_fasta_path,org = org_row[["name"]],genes = genes, verbose = verbose, params_list = params_list))){
      if ( (!is.logical(gtf_path) && !is.logical(genome_path) && file.exists(gtf_path) && file.exists(genome_path) && file.info(gtf_path)$size > 20 && file.info(genome_path)$size > 20)) { 
        fs::dir_create(org_fasta_path, recurse = T)
        cat(paste("Logfile : ",params_list$TEMP_PATH,"/",org_row[["name"]],"_",org_row[["version"]],"_",db,".log\n",sep=""))
        return(extract_transcript_regions(genome_path,gtf_path,gene_list,org_row[["name"]],org_row[["accession"]],org_row[["taxid"]],org_row[["version"]], org_row[["db"]], params_list, keep_data, verbose, seed=seed))
      }else{
        if(!keep_data){
          fs::file_delete(genome_path)
          fs::file_delete(gtf_path)
        }
        ex_genome <- file.path(params_list$GENOMES_PATH,org,paste0("user_",org_row[["version"]],".fa"))
        ex_genome_idx <- file.path(params_list$GENOMES_PATH, org,paste0("user_",org_row[["version"]],".fa.fai"))
        if(file.exists(ex_genome))
          fs::file_delete(ex_genome)
        if(file.exists(ex_genome))
          fs::file_delete(ex_genome_idx)
        if(!is.null(tmp_gene_list)){
          fs::file_delete(tmp_gene_list)   
        }
        message(print_toc(tictoc::toc(quiet = T, log = T)))
        return(NULL)
      }
    }else{
      if(file.exists(genome_path))
        fs::file_delete(genome_path)
      if(file.exists(gtf_path))
        fs::file_delete(gtf_path)
      ex_genome <- file.path(params_list$GENOMES_PATH,paste0(org_row[["name"]],".fa"))
      ex_genome_idx <- file.path(params_list$GENOMES_PATH,paste0(org_row[["name"]],".fa.fai"))
      if(file.exists(ex_genome))
        fs::file_delete(ex_genome)
      if(file.exists(ex_genome_idx))
        fs::file_delete(ex_genome_idx)
      if(!is.null(tmp_gene_list)){
        fs::file_delete(tmp_gene_list)   
      }   
      message(print_toc(tictoc::toc(quiet = T, log = T)))
      return(future::future({
        return(c(org=org_row[["name"]],accession=org_row[["accession"]],taxid=org_row[["taxid"]],version=org_row[["version"]], db=org_row[["db"]], genome=genome_path, gtf=gtf_path, source="r-biomartr"))
      }, label=paste0("final return fetch_FASTA_biomartr():", paste(org_row[["name"]], org_row[["version"]],org_row[["db"]],sep=":"))))
    }
  }else{
    if(!is.null(tmp_gene_list)){
      fs::file_delete(tmp_gene_list)
    }
    message(print_toc(tictoc::toc(quiet = T, log = T)))
    message(paste("Organism not available :",  org_row[["accession"]], org_row[["name"]], org_row[["db"]], org_row[["version"]],"\n\n"))
    return(NULL)
  }
}

#' Internal Function - Get FASTA data
#'
#' Main function to download FASTA data.
#'
#' @note This function wraps around other fetch_FASTA_*() function (Except for fetch_FASTA_user() which calls fetch_FASTA() if a genome or a gtf is not provided)
#'
#' @examples
#'     params_list <- load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt"))
#'     fetch_FASTA(c(name="danio_rerio",version="106",accession="acc1",taxid="dummy1"), params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'
#' @param org_row (Named vector) Named vector with name of the organism (format important, eg. "danio_rerio") and other details eg, c(name="danio_rerio",version="106",accession="acc1",taxid="dummy1")
#' @param db (string) Any supported DB from biomartr::listGenomes(), Valid Inputs: "ensembl", "genbank", "refseq", "" (Default: "ensembl"). Note: If empty, data is retrieved from both "genbank" and "ensembl"
#' @param params_list (string/list/COMPLETE-options) Output of load_params()
#' @param gene_list (string/vector) Vector or File containing list of genes
#' @param keep_data (bool) Keep the downloaded genomes and GTF data? (Default: FALSE)
#' @param only_fetch (bool) Only fetch data, without extraction (only genomes & GTF annotations from biomartr) (Default: FALSE)
#' @param type (string) Filter to organisms. Input of biomartr::listGenomes() ["all", "kingdom", "group", "subgroup"] (Default: "all")
#' @param subset (string/vector) Filter to organisms. Input of biomartr::listGenomes() [Differs between DBs and filter type, check manually - e.g,if type == "group|subgroup" then check biomartr::getGroups()] (Default: NULL)
#' @param verbose (bool) Verbosity (Default: T)
#' @return (Named vector) Organism details on successful extraction, NULL otherwise
#' @export
fetch_FASTA <- function(org_row, db="", params_list, gene_list, keep_data=F, ...) {
  # print(c("fetch_data(): ", org_row)) #DEBUG
  dot_args <- list(...)
  verbose=T
  if("verbose" %in% names(dot_args)){
    verbose<-dot_args$verbose
  }
  only_fetch=F
  if("only_fetch" %in% names(dot_args)){
    only_fetch<-dot_args$only_fetch
  }
  type="all"
  if("type" %in% names(dot_args)){
    type<-dot_args$type
  }
  subset=NULL
  if("subset" %in% names(dot_args)){
    subset<-dot_args$subset
  }
  seed <- 123
  if("seed" %in% names(dot_args)){
    seed <- dot_args$seed
  }
  unavailable_orgs <- data.frame()
  if("unavailable_orgs" %in% names(dot_args)){
    unavailable_orgs <- dot_args$unavailable_orgs
  }
  
  # if(stringi::stri_isempty(db)){
  #   return(dplyr::full_join(future::value(fetch_FASTA(org_row=org_row, db="genbank", params_list=params_list, gene_list=gene_list, keep_data=keep_data, ...=...)),
  #   future::value(fetch_FASTA(org_row=org_row, db="ensembl", params_list=params_list, gene_list=gene_list, keep_data=keep_data, ...=...))))
  # }

  if(isTRUE(stringi::stri_isempty(db)) || !isTRUE(any(grepl(x=db, pattern="genbank|refseq|ensembl")))){
    stop("db should be one of genbank|refseq|ensembl.\n")
  }
  
  org_ver <- ""
  org_name <- ""
  accession <- NA
  taxid <- NA
  
  if(isTRUE(all(c("name","version","accession","taxid") %in% names(org_row)))){
    org_row[["name"]] <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org_row[["name"]],perl = T, replacement = "_")))
    org_name <- org_row[["name"]]
    org_row[["version"]] <- trimws(gsub(pattern = "[[:space:]]", x=org_row[["version"]],perl = T, replacement = "_"))
    org_ver <- org_row[["version"]]
    org_row[["accession"]] <- trimws(gsub(pattern = "[[:space:]]", x=org_row[["accession"]],perl = T, replacement = "_"))
    accession <- org_row[["accession"]]
    taxid <- org_row[["taxid"]]
    org_row[["organism_name"]] <- org_row[["name"]]
    org_row[["asm_name"]] <- org_row[["version"]]
    org_row[["assembly_accession"]] <- org_row[["accession"]]
    org_row[["release"]] <- org_row[["version"]]
    org_row[["taxon_id"]] <- org_row[["taxid"]]
    # if(isTRUE("ftp_path" %in% names(org_row))){
    #   org_row[["ftp_path"]] <- org_row[["ftp_path"]]
    # }
    # print(paste(org_name,org_ver,accession,taxid)) #DEBUG
  }else if(isTRUE(any(c("genbank","refseq") %in% db))){
    if("organism_name" %in% names(org_row)){
      org_row[["organism_name"]] <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org_row[["organism_name"]],perl = T, replacement = "_")))
      org_name <- org_row[["organism_name"]]
    }
    if("asm_name" %in% names(org_row)){
      org_row[["asm_name"]] <- trimws(gsub(x = org_row[["asm_name"]], pattern = "[[:space:]]", replacement = "_"))
      org_ver <- org_row[["asm_name"]]
    }
    if("assembly_accession" %in% names(org_row)){
      org_row[["assembly_accession"]] <- trimws(gsub(x = org_row[["assembly_accession"]], pattern = "[[:space:]]", replacement = "_"))
      accession <- org_row[["assembly_accession"]]
    }
    if("taxid" %in% names(org_row)){
      org_row[["taxid"]] <- trimws(gsub(x = org_row[["taxid"]], pattern = "[[:punct:]]|[[:space:]]", replacement = "_"))
      taxid <- org_row[["taxid"]]
    }
  }else if(isTRUE("ensembl" %in% db)){
    if("name" %in% names(org_row)){
      org_row[["name"]] <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org_row[["name"]],perl = T, replacement = "_")))
      org_name <- org_row[["name"]]
    }
    if("release" %in% names(org_row)){
      org_row[["version"]] <- trimws(gsub(pattern = "[[:space:]]", x=org_row[["release"]],perl = T, replacement = "_"))
      org_ver <- org_row[["version"]]
    }
    if("accession" %in% names(org_row)){
      org_row[["accession"]] <- trimws(gsub(pattern = "[[:space:]]", x=org_row[["accession"]],perl = T, replacement = "_"))
      accession <- org_row[["accession"]]
    }
    if("taxon_id" %in% names(org_row)){
      org_row[["taxid"]] <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org_row[["taxon_id"]],perl = T, replacement = "_")))
      taxid <- org_row[["taxid"]]
    }
  }else{
    stop(paste("Org:", org_row,"does not contain name|version|accession|taxid"))
  }
  
  if(isTRUE(stringi::stri_isempty(org_name) || stringi::stri_isempty(org_ver) || stringi::stri_isempty(accession) || stringi::stri_isempty(taxid))){
    # print(org_name)
    # print(org_ver)
    stop("[fetch_FASTA()] Input must be a named list/vector with organism name and version. Accepted names: name, accession, taxid, version.\n")
  }
  
  org_row[["db"]] <- db
  
  if(isTRUE(any(!is.na(match(COMPLETE_env$org.meta$name,org_name))))){    
    
    if(isTRUE(nrow(unavailable_orgs) > 0)){
      is_match <- unavailable_orgs$org_name %in% org_name & 
        unavailable_orgs$org_ver %in% org_ver & 
        unavailable_orgs$db %in% db &
        unavailable_orgs$accession %in% accession
      # print(is_match)
      if(isTRUE(any(is_match))){
        message(paste("UNAVAILABLE:", org_name,org_ver,db, "\n\n"))
        return(NULL)
      }
    }
  
    tictoc::tic(msg=paste("[fetch_FASTA()] Processed:",org_name, db, org_ver, accession))
  
    org_fasta_path <- file.path(params_list$FASTA_OUT_PATH,org_name, db, org_ver)
    fs::dir_create(org_fasta_path, recurse=T)
    tmp_gene_list <- NULL
  
    if (length(gene_list) == 1 && file.exists(gene_list)) {
      genes <- factor(scan(gene_list, character(),quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
      genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
    }else{
      genes <- tolower(as.vector(gene_list))
      tmp_gene_list <- tempfile(pattern="genelist",tmpdir = params_list$TEMP_PATH)
      data.table::fwrite(x = list(gene_list),file = tmp_gene_list ,quote = F,row.names = F,col.names = F, nThread = params_list$numWorkers)
      gene_list <- tmp_gene_list
    }
  
    #print(check_files(fasta_path = org_fasta_path,org = org_name,genes = genes,params_list = params_list))
  
    if(isTRUE(!params_list$CLEAN_EXTRACT && check_files(fasta_path = org_fasta_path,org = org_name,genes = genes,params_list = params_list, verbose = verbose) && !keep_data)){
      cat(print_toc(tictoc::toc(quiet = T, log = T)))
      return(future::future({
        return(c(org=org_name,accession=accession,taxid=taxid,version=org_ver,db=db, genome="-",gtf="-",source="r-biomaRt"))
        }, label=paste0("return fetch_FASTA():", paste(org_name, org_ver,db,sep=":"))))
    }
  
    if(any(c("genbank","refseq") %in% db)){
      org <- COMPLETE_env$org.meta[which(COMPLETE_env$org.meta$name %in% org_name & COMPLETE_env$org.meta$asm_name %in% org_ver & COMPLETE_env$org.meta$assembly_accession %in% accession),]
      if(isTRUE(is.null(org) || nrow(org) == 0)){
        #possible the version is not found or version info is provided by the user
        #try an easier fallback filter 
        org <- COMPLETE_env$org.meta[which(COMPLETE_env$org.meta$name %in% org_name | COMPLETE_env$org.meta$asm_name %in% org_ver | COMPLETE_env$org.meta$assembly_accession %in% accession),]
      }
    }else{
      org <- COMPLETE_env$org.meta[which(COMPLETE_env$org.meta$name %in% org_name & COMPLETE_env$org.meta$assembly %in% org_ver & COMPLETE_env$org.meta$accession %in% accession),]
      if(isTRUE(is.null(org) || nrow(org) == 0)){
        #possible the version is not found or version info is provided by the user
        #try an easier fallback filter 
        org <- COMPLETE_env$org.meta[which(COMPLETE_env$org.meta$name %in% org_name | COMPLETE_env$org.meta$assembly %in% org_ver | COMPLETE_env$org.meta$accession %in% accession),]
      }
    }
    
    if(isTRUE(nrow(org) <= 0)){
      message(print_toc(tictoc::toc(quiet = T, log = T)))
      stop(paste("[fetch_FASTA()] Error: No organism passed the filters.",org_name, org_ver, db, accession))
    }
    
    org <- org %>% dplyr::mutate(db=db)
    
    if(isTRUE(any(c("genbank","refseq") %in% db))){
      if("assembly_accession" %in% colnames(org)){
        org$accession <- trimws(gsub(x = org$assembly_accession, pattern = "[[:punct:]]|[[:space:]]", replacement = "_"))
      }else{
        stop("[fetch_FASTA_biomartr()] assembly_accession required in column names when using genbank DB.")
      }
      if("organism_name" %in% colnames(org)){
        org$name <- trimws(tolower(gsub(x = org$organism_name, pattern = "[[:punct:]]|[[:space:]]", replacement = "_")))
      }else{
        stop("[fetch_FASTA_biomartr()] organism_name required in column names when using genbank DB.")
      }
      if("asm_name" %in% colnames(org)){
        org$version <- trimws(gsub(x = org$asm_name, pattern = "[[:punct:]]|[[:space:]]", replacement = "_"))
      }else{
        stop("[fetch_FASTA_biomartr()] asm_name required in column names when using genbank DB.")
      }
    }else if(isTRUE("ensembl" %in% db)){
      if("accession" %in% colnames(org)){
        org$accession <- trimws(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org$accession,perl = T, replacement = "_"))
      }else{
        stop("[fetch_FASTA_biomartr()] assembly_accession required in column names when using ensembl DB.")
      }
      if("name" %in% colnames(org)){
        org$name <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org$name,perl = T, replacement = "_")))
      }else{
        stop("[fetch_FASTA_biomartr()] name required in column names when using ensembl DB.")
      }
      if("release" %in% colnames(org)){
        org$version <- trimws(gsub(pattern = "[[:punct:]]|[[:space:]]", x=org$release,perl = T, replacement = "_"))
      }else{
        stop("[fetch_FASTA_biomartr()] release required in column names when using ensembl DB.")
      }
    }

      return(future::future({   
        ret_vals <- future.apply::future_lapply(purrr::transpose(org), function(orgx){
          mart_check_res <- NULL
          if(isTRUE(any(grepl(x=orgx[["db"]],pattern = "ensembl|user",ignore.case = T)))){
            mart_check_res <- check_mart_dataset(orgx[["name"]],orgx[["version"]],orgx[["db"]]);
          }
          # print(paste("HERE",org_row,orgx[["db"]], mart_check_res))
          if(isTRUE(is.null(mart_check_res))){
            tryCatch({
              ret_val <- fetch_FASTA_biomartr(org_row = orgx, db=orgx[["db"]], params_list = params_list, gene_list = gene_list, data_types=c("gtf","genome"), unavailable_orgs=unavailable_orgs, keep_data=keep_data, verbose=verbose, only_fetch=only_fetch, type=type, subset=subset, seed=seed) 
              # print(paste("HERE", ret_val))
              if(isTRUE(!is.null(ret_val))){  
                return(future::value(ret_val))
              }else{
                return(ret_val)
              }
            }, error=function(cond2){
              warning(paste(cond2,"\n"))
              stop(cond2)
            }) 
          }
  
          odb_list <- file.path(params_list$OUT_PATH,"genes",orgx[["name"]],orgx[["db"]],orgx[["version"]],"odb.list")
          
          odb_gene_map <- file.path(params_list$OUT_PATH,"genes",orgx[["name"]],orgx[["db"]],orgx[["version"]],"odb.final_map")
          
          gtf_stats_file <- file.path(params_list$OUT_PATH,"genes",orgx[["name"]],orgx[["db"]],orgx[["version"]],"gtf_stats.csv")
        
          if(isTRUE(params_list$CLEAN_EXTRACT)){
            fs::file_delete(org_fasta_path,recursive = T,force=T,expand = T)
            fs::file_delete(odb_list,force=T,expand = T)
            fs::file_delete(odb_gene_map,force=T,expand = T)
            fs::file_delete(gtf_stats_file,force=T,expand = T)
          }
        
          fs::dir_create(file.path(params_list$OUT_PATH,"genes",orgx[["name"]],orgx[["db"]],orgx[["version"]], sep=""),recurse = T)
        
          odb_list_genes <- c()
          if(isTRUE(COMPLETE_env$USE_ORTHODB)){
            if(isTRUE(any(params_list$CLEAN_EXTRACT, !file.exists(odb_list), file.info(odb_list)$size == 0))){
              proc <- do.call(add_to_process,list(p_cmd = COMPLETE_env$SHELL, p_args = c(fs::path_package("COMPLETE","exec","functions.sh"), "check_OrthoDB",orgx[["name"]], orgx[["version"]], orgx[["db"]], gene_list, odb_list, odb_gene_map,params_list$param_file, COMPLETE_env$SELECT_ALL_GENES), params_list=params_list,verbose = verbose))
              proc$wait(timeout=-1)
              if (inherits(proc, "process")) {
                proc$finalize() # Or explicitly close connections
              }
            }
            if(isTRUE(file.exists(odb_list) && file.info(odb_list)$size > 0)){
              odb_list_genes <- factor(scan(odb_list, character(), quiet = T))
              odb_list_genes <- odb_list_genes[grep("gene",tolower(odb_list_genes), invert = T, fixed = T)]
            }else{
              message(paste("ODB gene list could not be found for : ",orgx[["name"]]))
            }
          }
    
          gtf_data <- c()
          gtf_data <- get_gtf_mart(org = orgx[["name"]],orgx[["version"]],orgx[["db"]], gene_list = unique(c(genes,odb_list_genes)))
          if(isTRUE(nrow(gtf_data) == 0)){
            tryCatch(
            {
              gtf_details <- fetch_FASTA_biomartr(org_row = org_row, db=orgx[["db"]], params_list = params_list, gene_list = genes, data_types="gtf", unavailable_orgs=unavailable_orgs, keep_data=T, verbose=verbose, only_fetch=only_fetch, type=type, subset=subset, seed=seed)#...=...)
              # print(paste("HERE5: GTF DETAILS:"))
              # print(gtf_details) #DEBUG
              gtf_details <- future::value(gtf_details)
              if(isTRUE(length(gtf_details) == 0)){
                traceback(3)
                stop(paste("Error downloading Genome || GTF : ",orgx[["name"]], "\n"))
              }
              #Filter GTF
              # check_OrthoDB $f_org_name $GENE_LIST $OUT_PATH/genes/$f_org_name/odb.list $OUT_PATH/genes/$f_org_name/odb.final_map $param_file
              # time zgrep -i $MODE -f <(printf -- '%s\n' "${eexp_gene[@]}") $ANNO_FILE
              filter_proc <- do.call(add_to_process,list(p_cmd = COMPLETE_env$SHELL, p_args = c(fs::path_package("COMPLETE","exec","functions.sh"),"filter_GTF", gtf_details[["gtf"]], gene_list, orgx[["name"]], orgx[["version"]], gtf_details[["db"]], params_list$param_file, COMPLETE_env$parallel), logfile=paste(params_list$TEMP_PATH,"/",orgx[["name"]],"_",orgx[["version"]],"_",db,".log",sep=""), params_list = params_list,verbose = verbose))
              filter_proc$wait(timeout=-1)
              if (inherits(filter_proc, "process")) {
                filter_proc$finalize() # Or explicitly close connections
              }
              
              if(filter_proc$get_exit_status() > 0){
                if(!keep_data){
                  fs::file_delete(gtf_details[["gtf"]])
                }
                # print(paste0(orgx[["name"]],": Organism does not have any query genes, append it to unavailable_orgs.txt"))
                cat(paste(orgx[["name"]],orgx[["version"]],orgx[["db"]],orgx[["accession"]],sep="\t"), sep="\n", file = file.path(params_list$OUT_PATH,"unavailable_orgs.txt"), append = TRUE)
                message(paste("UNAVAILABLE: Organism does not have any query genes:", orgx[["name"]],orgx[["version"]],orgx[["db"]], "\n\n"))
                return(NULL)
              }
            
              gtf <- COMPLETE::read.gtf(gtf_details[["gtf"]],attr = c("intact"),quiet=T)
            
              gtf$ensembl_gene_id <- get_gtf_attributes(gtf,"gene_id")
            
              gtf$ensembl_transcript_id <- get_gtf_attributes(gtf,"transcript_id")
            
              gtf$external_gene_name <- tolower(get_gtf_attributes(gtf,"gene_name"))
    
              gtf <- gtf %>% rename(transcript_start = start) %>% rename(transcript_end = end)
    
              gtf <- dplyr::bind_rows(purrr::map(gene_list, .f = function(x){
                  matched_idx <- grep(x=gtf$external_gene_name,pattern = x,ignore.case = T)
                  if(isTRUE(length(matched_idx) > 0)){
                      return(gtf[matched_idx,])
                  }
              }))
    
              if(!keep_data){
                fs::file_delete(gtf_details[["genome"]])
                fs::file_delete(gtf_details[["gtf"]])
              }
              gtf_data <- gtf %>% dplyr::distinct()
            }
            , error=function(cond2){
              message(print_toc(tictoc::toc(quiet = T, log = T)))
              traceback(3)
              stop(paste(cond2,"\n"))
            })
          }else{
            gtf_data <- gtf_data %>% dplyr::distinct()
          }
    
          req_columns <- c("external_gene_name",
                        "ensembl_gene_id",
                        "ensembl_transcript_id",
                        "strand",
                        "transcript_start",
                        "transcript_end")
    
          # print(paste("HERE6: Calculate Stats:",orgx[["name"]]))  
          gtf_stats <- calculate_stats(orgx[["name"]], gtf_data,allow_strand = params_list$STRAND, n_threads = params_list$numWorkers)  
          gtf_stats <- dplyr::bind_rows(gtf_stats)
          # print(paste("HERE7: GTF Stats:",orgx[["name"]]))
          # print(str(gtf_stats))
      
          if(isTRUE(nrow(gtf_stats) > 0)){
            names(gtf_stats)[grep(pattern="seqnames",names(gtf_stats))] <- "transcript_id"
            # print(str(gtf_stats))
            # print(paste("fetch_FASTA():",paste(colnames(gtf_stats),collapse=",")))
            gtf_stats$safe_gene_name <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(gtf_stats$gene_name))
            gtf_stats <- gtf_stats[gtf_stats$gene_name!=0 | gtf_stats$gene_id!=0,] # cleaning up
          }
          if(isTRUE(nrow(gtf_stats)==0)){
            message(paste("No genes passed filters for : ",orgx[["name"]]))
            gtf_stats_path <- file.path(params_list$OUT_PATH,"genes", orgx[["name"]], "gtf_stats.csv")
            # print(paste("fetch_FASTA():",gtf_stats_path)) #DEBUG
            # print(!file.exists(gtf_stats_path)) #DEBUG
            # print(file.info(gtf_stats_path)$size) #DEBUG
          
            if(!file.exists(gtf_stats_path) || file.info(gtf_stats_path)$size <= 0){
              #Organism does not have any genes, append it to unavailable_orgs.txt
              message(paste0("[fetch_FASTA()]",orgx[["name"]],": Organism does not have any genes, appending it to unavailable_orgs.txt:"))
              message(paste("UNAVAILABLE:", orgx[["name"]],orgx[["version"]],orgx[["db"]],orgx[["accession"]],"\n\n"))
              cat(paste(orgx[["name"]],orgx[["version"]],orgx[["db"]], orgx[["accession"]],sep="\t"), sep="\n", file = file.path(params_list$OUT_PATH,"unavailable_orgs.txt"), append = TRUE)
            }
            message(print_toc(tictoc::toc(quiet = T, log = T)))
            # print(paste("HERE8: GTF Stats == 0"))
            return(NULL)
          }
        
          gtf_stats <- fetch_FASTA_mart(org = orgx[["name"]],gtf_stats = gtf_stats,fasta_path = org_fasta_path, params_list = params_list, verbose=verbose, only_fetch=only_fetch, type=type, subset=subset) 
          if(isTRUE(is.null(gtf_stats) || nrow(gtf_stats) == 0)){
            gtf_details <- fetch_FASTA_biomartr(org_row = org_row, params_list = params_list, gene_list = genes, unavailable_orgs=unavailable_orgs, only_fetch=T, keep_data=T, verbose = F, seed=seed)
            # print(gtf_details) #DEBUG
            gtf_details <- future::value(gtf_details)
            if(isTRUE(is.null(gtf_details) || length(gtf_details) == 0)){
                traceback(3)
                stop(paste("Error downloading Genome || GTF : ",orgx[["name"]]))
            }
            #Filter GTF
            # check_OrthoDB $f_org_name $GENE_LIST $OUT_PATH/genes/$f_org_name/odb.list $OUT_PATH/genes/$f_org_name/odb.final_map $param_file
            # time zgrep -i $MODE -f <(printf -- '%s\n' "${eexp_gene[@]}") $ANNO_FILE
            filter_proc <- do.call(add_to_process,list(p_cmd = COMPLETE_env$SHELL, p_args = c(fs::path_package("COMPLETE","exec","functions.sh"), "check_OrthoDB",orgx[["name"]],gene_list, odb_list, odb_gene_map,params_list$param_file, COMPLETE_env$SELECT_ALL_GENES), params_list=params_list,verbose = verbose))
            filter_proc$wait(timeout=-1)
            if (inherits(filter_proc, "process")) {
              filter_proc$finalize() # Or explicitly close connections
            }
            gtf <- read.gtf(gtf_details[["gtf"]],attr = c("intact"),quiet=T)
            
            gtf$ensembl_gene_id <- get_gtf_attributes(gtf,"gene_id")
            
            gtf$ensembl_transcript_id <- get_gtf_attributes(gtf,"transcript_id")
            
            gtf$external_gene_name <- tolower(get_gtf_attributes(gtf,"gene_name"))
            
            gtf <- gtf %>% rename(transcript_start = start) %>% rename(transcript_end = end)
            
            gtf <- dplyr::bind_rows(purrr::map(gene_list, .f = function(x){
              matched_idx <- grep(x=gtf$external_gene_name,pattern = x,ignore.case = T)
              if(isTRUE(length(matched_idx) > 0)){
                return(gtf[matched_idx,])
              }
            }))
            
            if(!keep_data){
              fs::file_delete(gtf_details[["genome"]])
              fs::file_delete(gtf_details[["gtf"]])
            }
            gtf_data <- gtf %>% dplyr::distinct()
            gtf_stats <- calculate_stats(orgx[["name"]], gtf_data,allow_strand = params_list$STRAND, n_threads = params_list$numWorkers)  
          }

          if(isTRUE(nrow(gtf_stats)==0 || is.null(gtf_stats))){
            if(!is.null(tmp_gene_list))
              fs::file_delete(tmp_gene_list)   
            warning(paste("[fetch_FASTA()] Error fetching FASTA for : ",orgx[["name"]]))
            message(print_toc(tictoc::toc(quiet = T, log = T)))
            return(NULL)
          }
          #print(paste(org_fasta_path,orgx[["name"]],genes,odb_gene_map,params_list)) #DEBUG
          
          label_FASTA_files(fasta_path = org_fasta_path,org = orgx[["name"]],db=orgx[["db"]], org_ver=orgx[["version"]], gene_list = genes,odb_gene_map = odb_gene_map,params_list = params_list, duplicates.method = "merge")
          names(gtf_stats)[grep(pattern="transcript_length",names(gtf_stats))] <- "transcript_length.annotated"
          gtf_stats <- gtf_stats %>% mutate(org=orgx[["name"]])
          gtf_stats$transcript_length.estimated <- gtf_stats$five_flank + gtf_stats$total_cds_len + gtf_stats$three_flank
          gtf_stats <- gtf_stats[,c("gene_name","gene_id","transcript_id","total_cds_len","five_len","three_len","exon_count","cds_count","five_flank","three_flank","transcript_length.estimated","transcript_length.annotated","org")]
          gtf_stats$gene_name <- tolower(gtf_stats$gene_name)
          #print(gtf_stats_file)
          #save("gtf_stats",file="gtf_stats.RData")
          
          data.table::fwrite(x = gtf_stats,file = gtf_stats_file,quote = F,sep = ",",row.names = F,col.names = T,na = "-", nThread = params_list$numWorkers)
          unlink( unlist(parallel::mclapply(list.dirs(path = org_fasta_path, full.names=TRUE,recursive = F), function(x) {
            if (file.info(x)$size == 0) {
              return(x) #print(x)
            }
          }, mc.cores = params_list$numWorkers, mc.silent = T, mc.preschedule = T)) ,recursive = F, force = T, expand =T)
      
          orgs_files <- tolower(unique(tools::file_path_sans_ext(list.files(org_fasta_path,no.. = T,recursive = F))))
          
          data.table::fwrite(x= list(orgs_files),file=file.path(params_list$OUT_PATH,"genes",orgx[["name"]],orgx[["db"]],orgx[["version"]],"AVAILABLE_GENES"),quote=F,row.names=F,col.names=F,na = "-" , nThread = params_list$numWorkers )
          
          data.table::fwrite(x= list(tolower(unique(genes[is.na(match( tolower(genes), orgs_files ))]))) ,file=file.path(params_list$OUT_PATH,"genes",orgx[["name"]],orgx[["db"]],orgx[["version"]],"MISSING_GENES"),quote=F,row.names=F,col.names=F ,na = "-", nThread = params_list$numWorkers ) 
      
          if(!is.null(tmp_gene_list))
            fs::file_delete(tmp_gene_list)   
          
          rm(gtf)
          rm(gtf_stats)
          rm(gtf_data)
          rm(orgs_files)
          
          # invisible(gc(full = TRUE))
          
          cat(print_toc(tictoc::toc(quiet = T, log = T)))

            
          return(c(org=orgx[["name"]],accession=orgx[["accession"]],taxid=orgx[["taxid"]],version=orgx[["version"]],db=orgx[["db"]], genome=gtf_details[["genome"]],gtf=gtf_details[["gtf"]],source="r-biomaRt")) 
        
        })
        return(unlist(ret_vals))
      }, globals = list(org=org, COMPLETE_env = COMPLETE_env, params_list=params_list, genes=genes, gene_list=gene_list, tmp_gene_list=tmp_gene_list, org_fasta_path=org_fasta_path, org_row=org_row, unavailable_orgs=unavailable_orgs, verbose=verbose, keep_data=keep_data, only_fetch=only_fetch, type=type, subset=subset, seed=seed, check_files=check_files, check_mart_dataset=check_mart_dataset, get_gtf_mart=get_gtf_mart, fetch_FASTA_biomartr=fetch_FASTA_biomartr, add_to_process=add_to_process, read.gtf=read.gtf, get_gtf_attributes=get_gtf_attributes, calculate_stats=calculate_stats, fetch_FASTA_mart=fetch_FASTA_mart, label_FASTA_files=label_FASTA_files, print_toc=print_toc), packages = c("S4Vectors", "IRanges", "dplyr", "purrr", "fs", "tictoc"), label=paste("fetch_FASTA():",paste(org_name,org_ver,db,sep=":")), seed=seed)) #END - future::future
  }else{
    if(!is.null(tmp_gene_list)){
      fs::file_delete(tmp_gene_list)
    }
    message(print_toc(tictoc::toc(quiet = T, log = T)))
    # traceback(3)
    stop(paste("Organism not available :", org_name, db, org_ver,"\n\n"))
  }
}

#' Get FASTA for User data
#'
#' Function to download FASTA data for user data. fetch_FASTA() is invoked if a genome and a gtf are not provided.
#'
#' @note Use "-" in Genome/GTF to check for the organism in BIOMART
#'
#' @examples
#'     params_list <- load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt"))
#'     fetch_FASTA_user(c(org="danio_rerio",genome="http://some.link",gtf="some.file",version="106",accession="acc1",taxid="dummy1"),params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'     fetch_FASTA_user(c(org="danio_rerio",genome="-",gtf="-",version="106",accession="acc1",taxid="dummy1"),params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'     fetch_FASTA_user(c(org="xenopus_laevis",genome="https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz",gtf="https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/XENLA_10.1_GCF.gff3.gz",version="106",accession="acc1",taxid="dummy1"),params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"))
#'
#' @param data (Named vector) Named vector with name of the organism (format important, eg. "danio_rerio"), Genome, GTF and other details eg, c(org="danio_rerio",genome="http://some.link",gtf="some.file",version="106",accession="acc1",taxid="dummy1").
#' @param params_list (string/list/COMPLETE-options) Output of load_params()
#' @param gene_list (string/vector) Filename with the list of genes
#' @param keep_data (bool) Keep downloaded Genomes and GTF data? (Default: FALSE)
#' @param only_fetch (bool) Only fetch data, without extraction (only genomes & GTF annotations from biomartr) (Default: FALSE)
#' @param type (string) Filter to organisms. Input of biomartr::listGenomes() ["all", "kingdom", "group", "subgroup"] (Default: "all")
#' @param subset (string/vector) Filter to organisms. Input of biomartr::listGenomes() [Differs between DBs and filter type, check manually - e.g,if type == "group|subgroup" then check biomartr::getGroups()] (Default: NULL)
#' @param verbose (bool) Verbosity (Default: T)
#' @return (Named vector) Organism details on successful extraction, NULL otherwise
fetch_FASTA_user <- function(data, params_list, gene_list, keep_data=F, ...){
  dot_args <- list(...)
  verbose=T
  if("db" %in% names(dot_args)){
    db<-dot_args$db
  }
  if("verbose" %in% names(dot_args)){
    verbose<-dot_args$verbose
  }
  only_fetch=F
  if("only_fetch" %in% names(dot_args)){
    only_fetch<-dot_args$only_fetch
  }
  type="all"
  if("type" %in% names(dot_args)){
    type<-dot_args$type
  }
  subset=NULL
  if("subset" %in% names(dot_args)){
    subset<-dot_args$subset
  }
  seed <- 123
  if("seed" %in% names(dot_args)){
    seed <- dot_args$seed
  }
  unavailable_orgs <- data.frame()
  if("unavailable_orgs" %in% names(dot_args)){
    unavailable_orgs <- dot_args$unavailable_orgs
  }
  
  genome_path <- c()
  gtf_path <- c()
  #print(data)
  genome <- data[["genome"]]
  gtf <- data[["gtf"]]
  org <- data[["name"]]
  org_ver <- data[["version"]]
  accession <- data[["accession"]]
  taxid <- data[["taxid"]]

  if(isTRUE(COMPLETE_env$SKIP_USER_DATA)){
    return(NULL)
    # return(c(data, source="user-provided", db="user"))
  }

  if(isTRUE(nrow(unavailable_orgs) > 0)){
    is_match <- unavailable_orgs$org_name %in% org & 
      unavailable_orgs$org_ver %in% org_ver & 
      unavailable_orgs$db %in% "user" &
      unavailable_orgs$accession %in% accession
    # print(is_match)
    if(isTRUE(any(is_match))){
      message(paste("UNAVAILABLE:", org,org_ver,"user", "\n\n"))
      return(NULL)
    }
  }
  
  if(isTRUE(length(gene_list) == 1 && file.exists(gene_list))){
    genes <- factor(scan(gene_list, character(),quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
    data.table::fwrite(x = list(gene_list),file = file.path(params_list$TEMP_PATH,"gene_list.txt"),quote = F,row.names = F,col.names = F, nThread = params_list$numWorkers)
    gene_list <- file.path(params_list$TEMP_PATH,"gene_list.txt")
  }

  tictoc::tic(msg = paste("[User] Processed:",org, "user",org_ver, accession))
  org_fasta_path <- file.path(params_list$FASTA_OUT_PATH ,org,"user", org_ver)
  if(isTRUE(!params_list$CLEAN_EXTRACT && check_files(fasta_path = org_fasta_path,org = org,genes = genes, params_list = params_list, verbose = verbose) && !keep_data)){
    message(print_toc(tictoc::toc(quiet = T, log = T)))
    # print(c(org=org, version=org_ver,db=db,genome=genome,gtf=gtf, source="user-provided")) #DEBUG
    return(c(org=org, accession=accession, taxid=taxid, version=org_ver,db="user",genome=genome,gtf=gtf, source="user-provided"))
  }

  if(isTRUE(any(genome == "-", gtf == "-", genome == " ", gtf == " ", stringi::stri_isempty(genome), stringi::stri_isempty(gtf)))){
    # message(paste("Trying marts for:", org))
    tryCatch({
        # print(data) #DEBUG
        ret_val <- fetch_FASTA(org_row = c(name=as.character(org),accession=accession,taxid=taxid,version=org_ver,genome="-",gtf="-"), db=db, params_list = params_list, gene_list = genes, unavailable_orgs=unavailable_orgs, keep_data=keep_data, verbose=verbose, only_fetch=only_fetch, type=type, subset=subset, seed=seed) #verbose=F
        # message(print_toc(tictoc::toc(quiet = T, log = T)))
        # print(paste("HERE4.1: RET_VAL:", org))
        # print(ret_val)
        # print(future::value(ret_val))
        return(ret_val)
      }
      ,error=function(cond){
      message(print_toc(tictoc::toc(quiet = T, log = T)))
      traceback(3)
      stop(paste(cond,"\n"))
      #return(NULL)
    }) 
  }
  
  if(grepl(x = basename(URLdecode(genome)), pattern="gz")){
    genome_path<-file.path(params_list$GENOMES_PATH,org,paste("user_",org_ver,".",tools::file_ext(tools::file_path_sans_ext(basename(URLdecode(genome)))),".gz",sep = ""))  
  }else{
    genome_path<- file.path(params_list$GENOMES_PATH, org, paste("user_",org_ver,".",tools::file_ext(basename(URLdecode(genome))),sep = ""))
  }
  
  if(grepl(x = basename(URLdecode(gtf)), pattern="gz")){
    gtf_path<-file.path(params_list$ANNOS_PATH,org,paste("user_",org_ver,".",tools::file_ext(tools::file_path_sans_ext(basename(URLdecode(gtf)))),".gz",sep = ""))  #".gtf.gz"
  }else{
    gtf_path<-file.path(params_list$ANNOS_PATH,org,paste("user_",org_ver,".",tools::file_ext(basename(URLdecode(gtf))),sep = ""))  #".gtf.gz"
  }

  fs::dir_create(file.path(params_list$GENOMES_PATH,org),recurse = T)
  fs::dir_create(file.path(params_list$ANNOS_PATH,org),recurse = T)

  if(grepl("://|http|ftp|www",genome)){
    if(any(!file.exists(genome_path), !file.info(genome_path)$size > 20, params_list$CLEAN_EXTRACT)){
      ret_code <- curl::curl_fetch_disk(genome, paste(params_list$GENOMES_PATH, "/",basename(URLdecode(genome)),sep=""))
      if(isTRUE(ret_code$status_code!=200)){
        message(paste("Error",ret_code$status_code,": check genome URL :",org,"-",gtf))
        message(print_toc(tictoc::toc(quiet = T, log = T)))
        return(NULL)
      }
      fs::file_move(paste(params_list$GENOMES_PATH, "/",basename(URLdecode(genome)),sep=""),genome_path)
    } 
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
    if(any(!file.exists(gtf_path), !file.info(gtf_path)$size > 20, params_list$CLEAN_EXTRACT)){
      ret_code <- curl::curl_fetch_disk(gtf, paste(params_list$ANNOS_PATH, "/",basename(URLdecode(gtf)),sep=""))
      if(isTRUE(ret_code$status_code!=200)){
        message(paste("Error", ret_code$status_code,": check gtf URL :",org,"-",gtf))
        message(print_toc(tictoc::toc(quiet = T, log = T)))
        return(NULL)
      }
      #zcat -f $ANNO_FILE | gffread - -T -O -E -o - | gzip -c > $params_list$ANNOS_PATH/"$f_org_name".gtf.gz
      if(!grepl(pattern = "gtf",x = basename(URLdecode(gtf)),ignore.case = T)){

      }
      fs::file_move(paste(params_list$ANNOS_PATH, "/",basename(URLdecode(gtf)),sep=""),gtf_path)
    }
  }else{
    if(file.exists(gtf) && file.info(gtf)$size > 0){
      gtf_path <- gtf
    }else{
      message(paste("User GTF not found :",org,"-",gtf))
      message(print_toc(tictoc::toc(quiet = T, log = T)))
      return(NULL)
    }
  }

  #print(paste(genome_path,gtf_path,org))
  ex_genome <- file.path(params_list$GENOMES_PATH,org,paste0("user_",org_ver,".fa"))
  ex_genome_idx <- file.path(params_list$GENOMES_PATH, org,paste0("user_",org_ver,".fa.fai"))

  if(isTRUE(params_list$CLEAN_EXTRACT || !check_files(fasta_path = org_fasta_path,org = org,genes = genes, params_list = params_list, verbose = verbose))){
    if ( (!is.logical(gtf_path) && !is.logical(genome_path)) && !isTRUE(COMPLETE_env$SKIP_USER_DATA)) {
      fs::dir_create(org_fasta_path,recurse = T)
      cat(paste("Logfile : ",params_list$TEMP_PATH,"/",org,"_",org_ver,"_",db,".log\n",sep=""))
      return(extract_transcript_regions(genome_path,gtf_path,gene_list,org,accession,taxid,org_ver,"user", params_list, keep_data, verbose, seed))
    }else{
      if(!keep_data){
          fs::file_delete(genome_path)
          fs::file_delete(gtf_path)
      }
      if(file.exists(ex_genome))
        fs::file_delete(ex_genome)
      if(file.exists(ex_genome))
        fs::file_delete(ex_genome_idx)
      message(print_toc(tictoc::toc(quiet = T, log = T)))
      return(NULL)
    }
  }else{
    if(!keep_data){
          fs::file_delete(genome_path)
          fs::file_delete(gtf_path)
      }
      if(file.exists(ex_genome))
        fs::file_delete(ex_genome)
      if(file.exists(ex_genome))
        fs::file_delete(ex_genome_idx)
    message(print_toc(tictoc::toc(quiet = T, log = T)))
    return(future::future({
      return(c(org=org,accession=accession,taxid=taxid,version=org_ver,db="user",genome=genome,gtf=gtf, source="user-provided"))
      }, label=paste0("final return fetch_FASTA_user():", paste(org, org_ver,"user",sep=":"))))
  }
}

#' Deduplicate FASTA Sequence
#'
#' Merge, Make Unique or Delete FASTA sequences with duplicated names. IF the duplicate sequence names are CDS blocks, merge them. If the duplicate seq names are EXON blocks, make them unique. If duplicate seq names are sequence duplicates, delete them.
#'
#' @param fasta_path (string) Path to FASTA File
#' @param duplicates.method (string) merge/delete/make_unique. How to handle sequences with duplicate names?. For CDS/UTR blocks - merge (Concatenate sequences with duplicate names), for Exon blocks - make_unique (Sequence names/IDs are made unique), delete - for deleting all duplicate sequences (first seq is kept).
#' @param n_threads (integer) Number of Threads (Default=4)
#' @return (data.frame) Data frame with seq_names and seqs
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
#' @param fasta_path (string) Path to FASTA File
#' @param org (string) Name of the organism
#' @param db (string) Query DB
#' @param org_ver (string) Version of the organism
#' @param gene (string/vector) Filename of gene list or Vector of gene names
#' @param odb_clusters (string) ODB Clusters (String collapsed with ',')
#' @param duplicates.method (string) merge/delete/make_unique. How to handle sequences with duplicate names?. For CDS/UTR blocks - merge (Concatenate sequences with duplicate names), for Exon blocks - make_unique (Sequence names/IDs are made unique), delete - for deleting all duplicate sequences (first seq is kept).
#' @param params_list (string/list/COMPLETE-options) Output of load_params()
#' @export
label_seqIDs_lengthen <- function(fasta_path, org, db, org_ver, gene, odb_clusters, duplicates.method, params_list) {

  if(is.null(duplicates.method) || !grepl(pattern = c("merge|delete|make_unique"), x = duplicates.method,ignore.case = T)){
    traceback(3)
    stop(paste("duplicates.method must be one of merge|delete|make_unique.\n"))
  }

  seq_set <- deduplicate_FASTA(fasta_path=fasta_path, duplicates.method=duplicates.method, n_threads=params_list$numWorkers)

  split_seq_names <- stringi::stri_split(str = names(seq_set), fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)
  if( ncol(split_seq_names) == 1 ){ 
    #print(paste(names(seq_set),org,x,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)) #DEBUG
    #print(names(seq_set)) #DEBUG
    names(seq_set) <- paste(names(seq_set),paste(org,db,org_ver,sep=.Platform$file.sep),gene,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)
    Biostrings::writeXStringSet(x = seq_set,filepath = fasta_path,append = F,format = "fasta")
    #return(paste(names(seq_set),org,x,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)) #DEBUG
  }else if( ncol(split_seq_names) > length(COMPLETE_env$FORMAT_ID_INDEX) || ncol(split_seq_names) < length(COMPLETE_env$FORMAT_ID_INDEX) ){
    names(seq_set) <- paste(split_seq_names[,COMPLETE_env$FORMAT_ID_INDEX$TRANSCRIPT_ID] ,paste(org,db,org_ver,sep=.Platform$file.sep),gene,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM) 
    Biostrings::writeXStringSet(x = seq_set,filepath = fasta_path,append = F,format = "fasta")
  }else if(ncol(split_seq_names) != length(COMPLETE_env$FORMAT_ID_INDEX)){
    message(paste("Error : Check sequence ID format of", y))
  }
}

#' Label Sequence IDs of all FASTA Files in a folder (Refer ?COMPLETE_PIPELINE_DESIGN about COMPLETE.format.IDs)
#'
#' Convert short sequence IDs in a FASTA file into longer COMPLETE format IDs
#'
#' @param fasta_path (string) Path of FASTA File
#' @param org (string) Name of the organism
#' @param gene_list (string/vector) Filename of gene list or Vector of gene names
#' @param odb_gene_map (string) Filename of OrthoDB Ortholog Clusters mapped to gene names  (Optional). If not provided, sequences cluster ID/name "ungrouped" is used
#' @param duplicates.method (string) merge/delete/make_unique. How to handle sequences with duplicate names?. For CDS/UTR blocks - merge (Concatenate sequences with duplicate names), for Exon blocks - make_unique (Sequence names/IDs are made unique), delete - for deleting all duplicate sequences (first seq is kept).
#' @param params_list (string/list/COMPLETE-options) Output of load_params()
#' @export
label_FASTA_files <- function(fasta_path,org,db,org_ver,gene_list,odb_gene_map=NULL,params_list, duplicates.method){

  if(is.null(duplicates.method) || !grepl(pattern = c("merge|delete|make_unique"), x = duplicates.method,ignore.case = T)){
    traceback(3)
    stop(paste("duplicates.method must be one of merge|delete|make_unique.\n"))
  }

  tictoc::tic(msg = "Labelling Sequence IDs...")
  if(!is.null(odb_gene_map)){
    if(file.exists(odb_gene_map) && file.info(odb_gene_map)$size > 0){
      odb_gene_map <- read.table(file = odb_gene_map,header = F,quote = "",sep = "\t")
      #local ortho_cluster=$(grep -w $gene_name $odb_clusters | awk -F'\t' '{if (length(c) == 0){c=$1;}else{c=c","$1;}}END{print c}')
    }else{
      message(paste(odb_gene_map,"does not exist!"))
      odb_gene_map <- NULL
    }
  }
  if (length(gene_list) == 1 && file.exists(gene_list)) {
    genes <- factor(scan(gene_list, character(),quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
  }

  furrr::future_map(genes, function(x){ 
    safe_gene <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(x))
    file_path <- paste(fasta_path,safe_gene,sep="")
    if(!is.null(odb_gene_map)){
      odb_clusters <- paste(odb_gene_map[grep(pattern = x,x = odb_gene_map[,2],ignore.case = T,value = F),1],collapse = ",")
    }else{
      odb_clusters <- "ungrouped"
    }
    furrr::future_map(list.files(path = fasta_path,pattern = safe_gene,full.names = T), function(y){ #lapply
      label_seqIDs_lengthen(fasta_path = y, org = org, db=db, org_ver=org_ver, gene = x, odb_clusters = odb_clusters, params_list = params_list, duplicates.method = duplicates.method)

      # Biostrings::writeXStringSet(x = seq_set,filepath = y,append = F,format = "fasta")
    }, .options = furrr::furrr_options(seed=T, scheduling=params_list$numWorkers))

  },.options = furrr::furrr_options(seed=T, scheduling=params_list$numWorkers))#, mc.cores = params_list$numWorkers,mc.silent = T)
  cat(print_toc(tictoc::toc(quiet = T)))
}

#' Internal Function - Merge and Format OrthoDB Flat Files
#'
#' This step is essential for speeding up the extraction process. The gene information (Gene IDs and Gene Names) are merged with Ortholog Groups (Cluster IDs and Gene IDs) and the format is converted into a Tab Delimited file wrapped over a CSV of Gene IDs and Gene Names. This is performed because one Ortholog Group/Cluster can encompass multiple genes (across organisms). The final format looks like this [cluster1][tab][geneID1,geneID2..geneIDN][tab][gene_name1,gene_name2..gene_nameN].
#'
#' @note Required Flat files from OrthoDB are \*_OG2genes.tab.gz,\*_genes.tab.gz and \*_species.tab.gz. Output files are stored in paste(odb_prefix,"_OG2genes_fixed.tab.gz",sep=""), paste(odb_prefix,"_genes_fixed_user.tab.gz",sep="") and paste(odb_prefix,"_OG2genes_fixed_user.tab.gz",sep="")
#'
#' @param params_list (string/list/COMPLETE-options) Output of load_params() that contains prefix to the Flat Files from OrthoDB (Prefix of \*_OG2genes.tab.gz,\*_genes.tab.gz,\*_species.tab.gz) in ORTHODB_PREFIX parameter
#' @param quick.check (bool) Only check if files exist? (TRUE/FALSE). FALSE - When running the pipeline with a new list of genes (Default:T)
#' @param n_threads (integer) Number of threads
#' @param gene_list (string) File name of gene list
#' @return (bool) TRUE if output files exist, FALSE otherwise
#' @export
merge_OG2genes_OrthoDB <- function(params_list,quick.check=T,n_threads=tryCatch(parallel::detectCores(all.tests = T, logical = T), error=function(cond){return(2)}),gene_list){
  tictoc::tic(msg = paste("Checking and Transforming ODB Files..."))
  odb_prefix <- params_list$ORTHODB_PREFIX
  # print(c(gene_list, odb_prefix))
  if(!quick.check){
    processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"merge_OG2genes_OrthoDB",odb_prefix,!quick.check,n_threads,gene_list ) ,spinner = T,stdout = "1",stderr = "2>&1")
  }
  
  if( (file.exists(paste(odb_prefix,"_OG2genes_fixed.tab.gz",sep="")) && file.info(paste(odb_prefix,"_OG2genes_fixed.tab.gz",sep=""))$size > 23) && (file.exists(paste(odb_prefix,"_OG2genes_fixed_user.tab.gz",sep="")) && file.info(paste(odb_prefix,"_OG2genes_fixed_user.tab.gz",sep=""))$size > 23) && (file.exists(paste(odb_prefix,"_genes_fixed_user.tab.gz",sep="")) && file.info(paste(odb_prefix,"_genes_fixed_user.tab.gz",sep=""))$size > 23) && (file.exists(file.path(dirname(odb_prefix),"odb.gene_list")) && file.info(file.path(dirname(odb_prefix),"odb.gene_list"))$size > 0)){
    if(!all(tools::md5sum(file.path(dirname(odb_prefix),"odb.gene_list")) == tools::md5sum(gene_list))){
      cat(paste("New gene list detected:",print_toc(tictoc::toc(quiet = T))))
      fs::file_delete(x = c(paste(params_list$OUT_PATH,"available_orgs.txt"),
                 file.path(params_list$OUT_PATH,"unavailable_orgs.txt"),
                 file.path(params_list$OUT_PATH,"selected_ORGS.txt"),
                 file.path(params_list$OUT_PATH,"ALL_GROUPS.txt"),
                 file.path(params_list$OUT_PATH,"all_gtf_stats.csv")),force = T,expand = T)
      return(FALSE)
    }
    cat("Success:",paste(print_toc(tictoc::toc(quiet = T))))
    return(TRUE)
  }else{
    cat(paste("Check Failed:",print_toc(tictoc::toc(quiet = T))))
    return(FALSE)
  }

}

#' (1) - Extracts Sequences for Protein Coding Transcripts from Organisms
#'
#' This is the main function which calls all the other functions and performs and end-end execution of data extraction part of the pipeline. It requires a filename of a formatted parameter file and a gene list (check the github repo for an example) or fs::path_package("COMPLETE","pkg_data","parameters.txt").
#'
#' @examples
#'     COMPLETE::EXTRACT_DATA(db="ensembl", params_list = fs::path_package("COMPLETE","pkg_data","parameters.txt"), gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"), user_data = fs::path_package("COMPLETE","pkg_data", "user_data.txt"), only.user.data = F )
#'     COMPLETE::EXTRACT_DATA(db="genbank", params_list = fs::path_package("COMPLETE","pkg_data","parameters.txt"), gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"), user_data = NULL, only.user.data = F )
#'     COMPLETE::EXTRACT_DATA(db="ensembl", params_list = fs::path_package("COMPLETE","pkg_data","parameters.txt"), gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"), user_data = NULL, only.user.data = T )
#'     COMPLETE::EXTRACT_DATA(db="genbank", params_list=COMPLETE::load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt")), user_data="files/user_data.txt", gene_list = "files/genelist.txt", keep_data=T,db="genbank", type="subgroup", subset=c("Fishes","Amphibians"), skip_bacteria=T)
#'     COMPLETE::EXTRACT_DATA(params_list=COMPLETE::load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt")), user_data="files/user_data.txt", gene_list = "files/genelist.txt", keep_data=T,db="ensembl", type="kingdom", subset=c("EnsemblVertebrates","EnsemblBacteria"))
#'     
#' @note If samtools/bedtools are not available in $PATH user data (genomes & GTFs) is not processed (unless "-" is used where the organism is looked-up in BIOMART using biomaRt). Files are still downloaded and saved in params_list$GENOMES_PATH and params_list$ANNOS_PATH
#'
#' @param db (string) Any supported DB from biomartr::listGenomes(), Valid Inputs: "ensembl", "genbank", "refseq", "" (Default: "ensembl"). Note: If empty, data is retrieved from both "genbank" and "ensembl"
#' @param params_list (string/list/COMPLETE-options) Filename of a formatted parameter file (check the github repo for an example) or Output of load_params().
#' @param gene_list (string/vector) Vector or File with a list of genes to extract data for(check the github repo for an example)
#' @param user_data (string/data.frame) File name or table with user-specified organisms(genomes,GTFs). File must be in CSV format and should not contain header and column names are not required for the table. If data.frame is provided then it should have all the columns in c("name","accession","taxid","version","genome","gtf"). Check system.file("exec", "pkg_data", "user_data.txt", mustWork = T ,package = "COMPLETE") for an example user-data file.
#' @param only.user.data (bool) ONLY Process user data and not all available organisms? (TRUE/FALSE). (Default: FALSE)
#' @param keep_data (bool) Keep downloaded Genomes and GTF data? (TRUE/FALSE). (Default: FALSE)
#@param rerun_count Number of times to rerun the function to account for extraction/network failures and other errors. Default 3
#' @param only_fetch (bool) Only fetch data, without extraction (only genomes & GTF annotations from biomartr) (Default: FALSE)
#' @param type (string) Filter to organisms. Input of biomartr::listGenomes() ["all", "kingdom", "group", "subgroup"] (Default: "all")
#' @param subset (string/vector) Filter to organisms. Input of biomartr::listGenomes() [Differs between DBs and filter type, check manually - e.g,if type == "group|subgroup" then check biomartr::getGroups()] (Default: NULL)
#' @param skip_bacteria (bool) [biomartr::listGenomes()] Due to its enormous dataset size (> 700MB as of July 2023), the bacterial summary file will not be loaded by default anymore. If users wish to gain insights for the bacterial kingdom they needs to actively specify skip_bacteria = FALSE. When skip_bacteria = FALSE is set then the bacterial summary file will be downloaded. (Default: F - Unless subset contains "bact", T Otherwise)
#' @param verbose (bool) Verbosity (Default: T)
#' @seealso [biomartr::listGenomes()], [COMPLETE::FIND_ORTHOLOGS()]
#' @export
EXTRACT_DATA <- function(db="", params_list, gene_list, user_data=NULL, only.user.data=F, keep_data=F, ...){
  if(isTRUE(stringi::stri_isempty(db) || !any(grepl(x=db, pattern="genbank|refseq|ensembl")))){
    stop("db should be one of genbank|refseq|ensembl.\n")
  }
  
  dot_args <- list(...)
  verbose=T
  if("verbose" %in% names(dot_args)){
    verbose<-dot_args$verbose
  }
  only_fetch=F
  if("only_fetch" %in% names(dot_args)){
    only_fetch<-dot_args$only_fetch
  }
  type="all"
  if("type" %in% names(dot_args)){
    if(isTRUE(any(grepl(x=db, pattern="ensembl")) && any(grepl(x=dot_args$type, pattern="group|subgroup")))){
      message(paste("DB: 'ensembl' support only 'kingdom' filter type...Setting type='all'"))
    }else{
      type<-dot_args$type
    }
  }
  subset=NULL
  if("subset" %in% names(dot_args)){
    if(isTRUE(any(grepl(x=db, pattern="ensembl")) && any(grepl(x=dot_args$type, pattern="group|subgroup")))){
      message(paste("DB: 'ensembl' support only 'kingdom' filter type...Setting subset=NULL"))
    }else{
      subset<-dot_args$subset
    }
  }
  skip_bacteria <- dplyr::if_else(any(grepl(x = subset, pattern = "bact", ignore.case = T)),FALSE,TRUE)
  if("skip_bacteria" %in% names(dot_args)){
    if(isTRUE(dot_args$skip_bacteria))
      skip_bacteria <- dot_args$skip_bacteria
  }
  seed <- 123
  if("seed" %in% names(dot_args)){
    seed <- dot_args$seed
  }
  set.seed(seed)
  
  if(isTRUE(!is.null(subset) && !any(grepl(x=db, pattern="ensembl")))){
    message("Checking subset filter...")

    valid_groups <- biomartr::getGroups(db = db,"all")
    matching_subsets <- which(unlist(purrr::map(subset, function(x){return(grepl(x=valid_groups, pattern=x, fixed=T))})))
    if(length(matching_subsets) != length(subset)){
      stop(paste(valid_groups[matching_subsets],"not in db/type. Subset parameter should be one of :", paste0("biomartr::getGroups(db ='", db,"', 'all')")))
    }
  }

  # if(isTRUE(stringi::stri_isempty(db))){
  #   return(dplyr::full_join(EXTRACT_DATA(db="ensembl", params_list=params_list, gene_list=gene_list, user_data=user_data, only.user.data=only.user.data, keep_data=keep_data, rerun_count=rerun_count, verbose=verbose),
  #   EXTRACT_DATA(db="genbank", params_list=params_list, gene_list=gene_list, user_data=user_data, only.user.data=only.user.data, keep_data=keep_data, rerun_count=rerun_count, ...=...)))
  # }

  if (!curl::has_internet()) {
    traceback(3)
    stop("Check if there is an internet connection\n")
  }
  
  if( only.user.data && is.null(user_data) ){
    traceback(3)
    stop("only.user.data=TRUE but user_data is not provided!\n")
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
      traceback(3)
      stop("ERROR: Parameters file is missing and is required\n")
    }
    loaded_PARAMS <- load_params(params_list)
  }else{
    if(any(grepl(x = class(params_list), pattern = "COMPLETE-options"))){
      loaded_PARAMS <- params_list
    }else{
      traceback(3)
      stop("Error: params_list not valid!\n")
    }
  }
  print(loaded_PARAMS)

  install_parallel()

  print(paste("MAX PROCESSES:",loaded_PARAMS$numWorkers))

  future::plan(future::multisession, workers = loaded_PARAMS$numWorkers)
  
  tictoc::tic.clearlog()
  invisible(gc(full = TRUE))

  if (loaded_PARAMS$CLEAN_EXTRACT) {
    fs::file_delete(x = c(paste(loaded_PARAMS$OUT_PATH,"available_orgs.txt"),
                 file.path(loaded_PARAMS$OUT_PATH,"unavailable_orgs.txt"),
                 file.path(loaded_PARAMS$OUT_PATH,"selected_ORGS.txt"),
                 file.path(loaded_PARAMS$OUT_PATH,"ALL_GROUPS.txt"),
                 file.path(loaded_PARAMS$OUT_PATH,"all_gtf_stats.csv")),force = T,expand = T)
  }

  fs::dir_create(loaded_PARAMS$OUT_PATH, recurse = T)
  unlink(loaded_PARAMS$TEMP_PATH, recursive = T,force = T,expand = T)
  fs::dir_create(loaded_PARAMS$TEMP_PATH, recurse = T)
  fs::dir_create(loaded_PARAMS$FASTA_OUT_PATH, recurse = T)
  fs::dir_create(loaded_PARAMS$GENOMES_PATH, recurse = T)
  fs::dir_create(loaded_PARAMS$ANNOS_PATH, recurse = T)

  if(isTRUE(any(grepl(x=db,pattern = "ensembl", ignore.case = T)))){  
    if(isTRUE(is.null(COMPLETE_env$using.mart))){
      tryCatch({
        COMPLETE_env$using.mart <- mart_connect(biomaRt::useMart,args=list(COMPLETE_env$ENSEMBL_MART)) #For biomaRt
      }, error=function(err){
        stop(paste("(R-COMPLETE) Error initializing marts and datasets for 'ensembl' DB:",err))
      })
    }
    COMPLETE_env$org.meta.list <- mart_connect(biomaRt::listDatasets,args=list(mart=COMPLETE_env$using.mart)) #For biomaRt
  }

  COMPLETE_env$org.meta <- data.frame()
  tryCatch({
    if(isTRUE(!is.null(subset))){
      filtered_subset <- mart_connect(biomartr::listGenomes,args=list(db = db, type = type, subset=subset, details = T, skip_bacteria=skip_bacteria))
      if (isTRUE(!is.null(filtered_subset) && "taxon_id" %in% names(filtered_subset))) {
        filtered_subset <- dplyr::mutate(filtered_subset, taxon_id = as.character(taxon_id))
      }
      COMPLETE_env$org.meta <- filtered_subset %>% dplyr::distinct()
    }
  },error=function(cond){
      # message(paste0("biomartr::getGroups(db =", db,",'all')"))
      stop(paste("(R-COMPLETE) biomartr Error:",cond))
  })
  

  unavailable_orgs <- c()
  all_orgs <- c()
  if(isTRUE(any(c("genbank","refseq") %in% db))){
    all_orgs <- unique(trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=COMPLETE_env$org.meta$organism_name,perl = T, replacement = "_"))))
  }else if(isTRUE("ensembl" %in% db)){
    all_orgs <-  unique(COMPLETE_env$org.meta$name)
  }
  
  if(file.exists(file.path(loaded_PARAMS$OUT_PATH,"unavailable_orgs.txt")) && file.info(file.path(loaded_PARAMS$OUT_PATH,"unavailable_orgs.txt"))$size > 0){
    unavailable_orgs <- read.table(file = file.path(loaded_PARAMS$OUT_PATH,"unavailable_orgs.txt"),header = F,sep = "\t",fill = T,na.strings = "",as.is = T, colClasses = "character") 
    colnames(unavailable_orgs) <- c("org_name","org_ver","db", "accession")
    if(isTRUE(length(setdiff(unavailable_orgs$org_name,all_orgs)) == length(unavailable_orgs$org_name) && length(intersect(all_orgs,unavailable_orgs$org_name)) == 0)){
      #unavailable organisms and all_orgs do not intersect, might be due to new set of orgs (and/or different subset filter)
      message(paste("Unavailable orgs and all orgs do not intersect, might be due to new set of orgs (and/or different subset filter)..."))
      unavailable_orgs <- c()
    }
    orgs_to_fetch <- unique(COMPLETE_env$org.meta[which(is.na(match(all_orgs,unavailable_orgs$org_name))),])
  }else{
    orgs_to_fetch <- unique(COMPLETE_env$org.meta)
  }

  if(isTRUE(any(c("genbank","refseq") %in% db))){
    # print(colnames(orgs_to_fetch)) #DEBUG
    orgs_to_fetch <- orgs_to_fetch %>% dplyr::filter(assembly_level == "Complete Genome")
    orgs_to_fetch$name <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=orgs_to_fetch$organism_name,perl = T, replacement = "_")))
    COMPLETE_env$org.meta$name <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=COMPLETE_env$org.meta$organism_name,perl = T, replacement = "_")))
    org_ver <- orgs_to_fetch$asm_name
    orgs_to_fetch$version <- trimws(gsub(x = org_ver, pattern = "-", replacement = "_"))
    COMPLETE_env$org.meta$version <- trimws(gsub(x = COMPLETE_env$org.meta$asm_name, pattern = "-", replacement = "_"))
  }else if(isTRUE("ensembl" %in% db)){
    orgs_to_fetch$version <- trimws(tolower(gsub(pattern = "[[:punct:]]|[[:space:]]", x=orgs_to_fetch$assembly,perl = T, replacement = "_")))
  }

  if(isTRUE(!is.null(user_data) && is.character(user_data))){
    user_data <- read.csv(user_data,header = F)
  }

  if(isTRUE(!is.null(user_data))){ 
    names(user_data) <- c("name","accession","taxid","version","genome","gtf")
    all_orgs <- c(all_orgs,user_data$name)

    if(isTRUE(!is.null(unavailable_orgs) && length(unavailable_orgs) > 0)){
      user_data <- user_data[which(is.na(match(user_data$name,unavailable_orgs$org_name))),]
    }
    if(isTRUE(!is.null(user_data)  && nrow(user_data) > 0)){
      orgs_to_fetch <- orgs_to_fetch[which(is.na(match(orgs_to_fetch$name,user_data$name))),]
    }

  }else{
    message(paste("User Data not provided or is empty!"))
  }

  tmp_gene_list <- NULL
  if (length(gene_list) == 1 && file.exists(gene_list)) {
    genes <- factor(scan(gene_list, character(), quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
    tmp_gene_list <- tempfile(pattern="genelist",tmpdir = params_list$TEMP_PATH)
    data.table::fwrite(x = list(gene_list),file = tmp_gene_list ,quote = F,row.names = F,col.names = F, nThread = loaded_PARAMS$numWorkers)
    gene_list <- tmp_gene_list
  }

  if(isTRUE(COMPLETE_env$USE_ORTHODB)){
    if(!merge_OG2genes_OrthoDB(params_list = loaded_PARAMS,quick.check = !loaded_PARAMS$CLEAN_EXTRACT,n_threads = loaded_PARAMS$numWorkers,gene_list = gene_list)){
      merge_OG2genes_OrthoDB(params_list = loaded_PARAMS,quick.check = F,n_threads = loaded_PARAMS$numWorkers,gene_list = gene_list)
      }
  }

  cat("Checking and downloading transcripts...\n")

  cat(paste("User Data Logs saved in : ",loaded_PARAMS$TEMP_PATH,"/*.log\n\n", sep=""))

  user_saved_meta <- c()
  saved_meta <- c()

  tictoc::tic(msg = "Total Extraction Time :")

  if(!isTRUE(COMPLETE_env$SKIP_USER_DATA) && !is.null(user_data)){ 
    if(isTRUE(nrow(user_data)!=0 && !is.null(user_data))){
      # print(user_data) #DEBUG
      user_saved_meta <- apply(user_data, MARGIN = 1, function(x){
      tryCatch({
          # print(x) #DEBUG
          return(fetch_FASTA_user(data = x, params_list = loaded_PARAMS, gene_list = genes, unavailable_orgs=unavailable_orgs, keep_data=keep_data, db=db, verbose=verbose, only_fetch=only_fetch, type=type, subset=subset, seed=seed))
          #cat(paste("DONE :", x["org"],"\n"));
        },error=function(cond){
          message(cond)
          return(NULL)
        })
      })
    }
  }else{
    message("User Data is skipped!")
  }

  if( !only.user.data ){ 
    saved_meta_futures <- apply(orgs_to_fetch, MARGIN = 1, FUN = function(x){
    tryCatch({
        # print(x) #DEBUG
        return(fetch_FASTA(org_row = x, db=db, params_list = loaded_PARAMS, gene_list = genes, unavailable_orgs=unavailable_orgs, keep_data=keep_data, verbose=verbose, only_fetch=only_fetch, type=type, subset=subset, seed=seed))
      },error=function(cond){
        message(cond)
        return(NULL)
      })
    })
  }

  progress_bar <- progressr::progressor(steps = length(user_saved_meta) + length(saved_meta_futures)) # Create the progressor
  # print(str(saved_meta_futures)) #DEBUG
  saved_meta  <- lapply(list(user_saved_meta, saved_meta_futures), FUN=function(x){
    if(isTRUE(!is.null(x))){
      org_details <- future::value(x)
      print(paste("HERE4:ORG_DETAILS", length(org_details)))
      print(str(org_details))
      progress_bar()
      if(isFALSE(is.null(org_details) && length(org_details) == 0)){
        lapply(org_details, function(row){
          if(isTRUE(!is.null(row))){
            org_fasta_path <- file.path(params_list$FASTA_OUT_PATH ,row[["org"]],row[["db"]], row[["version"]])
            check_files(fasta_path = org_fasta_path,org = row[["org"]],genes = genes, verbose = verbose, params_list = params_list)
          }
        })
      }
      return(org_details)
    }
  })
  
  # print(COMPLETE_env$process_list) #DEBUG
  
  try(parallel::mclapply(COMPLETE_env$process_list, function(x){ 
    # print(x) #DEBUG
    if(isTRUE(!is.null(x) && x$is_alive())){
      x$wait(timeout=-1)
      if (inherits(x, "process")) {
        x$finalize() # Or explicitly close connections
      }
    }
  }, mc.cores =  loaded_PARAMS$numWorkers))

  final_toc_print <- tictoc::toc(quiet = T, log = T)
  while(!is.null(final_toc_print)){
    cat(print_toc(final_toc_print))
    final_toc_print <- tictoc::toc(quiet = T, log = T)
  }

  #save(saved_meta, file ="saved_meta.RData")

  #print(saved_meta) #DEBUG
  #save(saved_meta, file ="saved_meta.RData") #DEBUG
  if(isTRUE(!is.null(saved_meta) && length(saved_meta) > 0)){
    saved_meta[sapply(saved_meta, is.null)] <- NULL
    saved_meta <- t(saved_meta) #dplyr::bind_rows(t(saved_meta))
    # save(saved_meta, file="saved_meta.rds")
    saved_meta <- unlist(saved_meta, recursive=F)
    saved_meta <- dplyr::bind_rows(saved_meta[!sapply(saved_meta,is.null)])
    data.table::fwrite(x = saved_meta,file = file.path(loaded_PARAMS$OUT_PATH,"org_meta.txt"), quote = F,sep=",", row.names = F,col.names = T,na = "-", nThread = loaded_PARAMS$numWorkers)
  }

  #DELETING EMPTY DIRECTORIES - THESE ORGANISMS COULD NOT BE FETCHED
  parallel::mclapply(list.dirs(path = loaded_PARAMS$FASTA_OUT_PATH, full.names=TRUE,recursive = F), function(x) { 
    fi <- file.info(x)
    if (fi$isdir) {
      f <- list.files(x, all.files=TRUE, recursive=TRUE, full.names=TRUE)
      sz <- sum(file.info(f)$size)

      #as precaution, print to make sure before using unlink(x, TRUE)
      if (sz==0L) unlink(x,recursive = T, force = T, expand =T) #print(x)
    }
  }, mc.cores = loaded_PARAMS$numWorkers, mc.silent = T, mc.preschedule = T)

  tictoc::tic(msg= "Coercing metdata from available organisms ...")
  
  missing_genes_list <- parallel::mclapply(list.files(path = file.path(loaded_PARAMS$OUT_PATH,"genes"),include.dirs=TRUE, full.names=TRUE),function(x){
    if(file.exists(file.path(x,"MISSING_GENES")) && file.info(file.path(x,"MISSING_GENES"))$size > 0 ){
      return(scan(file.path(x,"MISSING_GENES"), character(), quiet = T))
    }
  }, mc.cores =  loaded_PARAMS$numWorkers)
  missing_genes <- tolower(purrr::reduce(missing_genes_list, unique))
  available_genes_list <- parallel::mclapply(list.files(path = file.path(loaded_PARAMS$OUT_PATH,"genes"),include.dirs=TRUE, full.names=TRUE),function(x){
    if(file.exists(file.path(x,"AVAILABLE_GENES")) && file.info(file.path(x,"AVAILABLE_GENES"))$size > 0 ){
      return(scan(file.path(x,"AVAILABLE_GENES"), character(), quiet = T))
    }
  }, mc.cores =  loaded_PARAMS$numWorkers)
  available_genes <- tolower(purrr::reduce(available_genes_list, unique))

  available_orgs <- list.dirs(path= loaded_PARAMS$FASTA_OUT_PATH, full.names = F,recursive = F) #factor(scan(paste(loaded_PARAMS$OUT_PATH,"/available_orgs.txt
  data.table::fwrite(x = list(available_orgs) ,file = file.path(loaded_PARAMS$OUT_PATH,"available_orgs.txt"), quote = F, row.names = F,col.names = F, nThread = loaded_PARAMS$numWorkers)

  if(isTRUE(nrow(unavailable_orgs) > 0)){
    if(isTRUE(length(available_orgs) > 0)){ 
      # print(unavailable_orgs) #DEBUG
      unavailable_orgs <- unavailable_orgs %>% dplyr::filter(org_name %in% all_orgs[which(is.na(match(all_orgs,available_orgs)))])
      cat(paste("Unavailable Organisms :",paste(unique(unavailable_orgs$org_name), collapse = ","),"\n"))
    }else{
      unavailable_orgs <- unavailable_orgs %>% dplyr::filter(org_name %in% all_orgs)
      traceback(3)
      stop("No organisms were available!. Rety with other options or a different gene list.\n")
    }
    data.table::fwrite(x = list(unavailable_orgs),file = file.path(loaded_PARAMS$OUT_PATH,"unavailable_orgs.txt"), quote = F, row.names = F,col.names = F, nThread = loaded_PARAMS$numWorkers, append=T)
  }

  data.table::fwrite(x = list(list.files(path = loaded_PARAMS$FASTA_OUT_PATH,include.dirs=TRUE, full.names=F)),file = file.path(loaded_PARAMS$OUT_PATH,"selected_ORGS.txt"), quote = F, row.names = F, col.names = F, nThread = loaded_PARAMS$numWorkers)
  selected_orgs <-  factor(scan( file.path(loaded_PARAMS$OUT_PATH,"selected_ORGS.txt"), character(), quiet = T))

  all_gtf_stats <- dplyr::bind_rows(parallel::mclapply(list.files(path = file.path(loaded_PARAMS$OUT_PATH,"genes"),include.dirs=TRUE, full.names=TRUE),function(x){
    if(file.exists(file.path(x,"gtf_stats.csv")) && file.info(file.path(x,"gtf_stats.csv"))$size > 0 ){
      tmp_tab <- read.table(file = file.path(x,"gtf_stats.csv"),header = T,sep = ",",fill = T,na.strings = "",as.is = T, colClasses = "character")
      tmp_tab <- tmp_tab[stats::complete.cases(tmp_tab),]
      if(length(tmp_tab) > 0){
        return(tmp_tab)
      }else{
        return(NULL)
      }
    }
  }, mc.cores =  loaded_PARAMS$numWorkers, mc.preschedule = T))
  write.table(x = list(all_gtf_stats),file = file.path(loaded_PARAMS$OUT_PATH,"all_gtf_stats.csv"),sep = ",", quote = F, row.names = F,col.names = T,na = "-") 

  # 1. Get the directory info
  dirs_info <- fs::dir_info(loaded_PARAMS$FASTA_OUT_PATH, type = "directory", recurse = 3)
  # 2. Filter by checking if the number of files inside is > 0
  # Use map_int to count children for each directory
  non_empty_dirs <- dirs_info %>%
    mutate(n_files = map_int(path, ~ length(fs::dir_ls(.x, all = TRUE, no.. = TRUE)))) %>%
    filter(n_files > 0)
  # Optional: Still sort by the directory's own metadata size if needed
  non_empty_dirs <- non_empty_dirs %>% arrange(desc(size))
  data.table::fwrite(x = list(non_empty_dirs$path) ,file = file.path(loaded_PARAMS$OUT_PATH,"dir_paths.txt"), quote = F, row.names = F,col.names = F, nThread = loaded_PARAMS$numWorkers)
  
  # time ./find_orthologs.sh files/selected_ORGS.txt $1 #100 ##This also selects the transcripts
  # time ./align_seqs.sh $1
  # time ./predict_structures.sh $1
  # rm $TEMP_PATH/*
  # Rscript gene_stats.R >> files/stats.txt

  if(!is.null(tmp_gene_list))
    fs::file_delete(tmp_gene_list)   

  rm(unavailable_orgs)
  rm(all_gtf_stats)
  rm(available_orgs)
  rm(missing_genes)
  rm(missing_genes_list)
  rm(saved_meta)
  # rm(COMPLETE_env$process_list)
  rm(user_saved_meta)
  rm(saved_meta_futures)
  rm(tmp_gene_list)
  rm(genes)
  rm(gene_list)
  rm(orgs_to_fetch)
  rm(user_data)
  rm(all_orgs)
  
  cat(print_toc(tictoc::toc(quiet = T, log = T)))

  #tictoc::tic.clearlog()
  cat(paste("Extraction Time Log:\n"))
  cat(paste(tictoc::tic.log(),collapse = "\n"))
  cat(paste("\n\n"))
  
  # tictoc::clear()
  tictoc::tic.clearlog()
  invisible(gc(full = TRUE))
  # if(rerun_count > 0)
  #   EXTRACT_DATA(params_list, gene_list, user_data=user_data, db=db, only.user.data=only.user.data, keep_data=keep_data, rerun_count=rerun_count - 1)

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
#'          - odb12v2_species.tab.gz - Ortho DB organism ids based on NCBI taxonomy ids (mostly species level) 
#'          - odb12v2_genes.tab.gz  -Ortho DB genes with some info 
#'          - odb12v2_OG2genes.tab.gz - OGs to genes correspondence 
#'          (OR)
#'          - odb12v2_OG2genes_fixed.tab.gz - Merged & Transformed ODB file (Done within pipeline)
#'          - odb12v2_genes_fixed_user.tab.gz - Subset of ODB genelist BASED on user gene list (Done within pipeline)
#'          - odb12v2_OG2genes_fixed_user.tab.gz - Merged & Transformed ODB file BASED on user gene list (Done within pipeline)
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
#'               >$transcript_id $transcripID_delimiter $transcript_region ($strand) $seqID_delimiter $org_name/$DB/$org_version $seqID_delimiter $gene_name $seqID_delimiter $ortho_cluster
#'               >SOME_TRANSCRIPT||cds(+)::SOMEORG/DB/VERSION::RANDOMGENE::ORTHOLOG_CLUSTERS
#'               >ENSDART00000193157||cds(+)::danio_rerio/ensembl/115::sulf1::18335at7898,51668at7742,360590at33208
#'
#' + FLOW :
#'
#'     1) EXTRACT_DATA() - 
#'            Extracts the transcript regions for Protein Coding Transcripts (provided in parameters, pipeline requires cds,5utr,3utr)
#'        from BIOMART(ensembl),genbank(ncbi) and/or User provided genomes & GTFs. This functions uses biomaRt/biomartr for extracting data from BIOMART
#'        and BASH function extract_transcript_regions() for user provided data. Extraction priority/flow : User Data > biomaRt > biomartr
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
