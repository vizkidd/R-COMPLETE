#@param sep2 Delimiter 2 of the BLAST File columns. Default - c("","|",""). Check ?data.table::fwrite or ?data.table::fread
#@param transcript_ID_metadata Tab-delimited File with the filenames, indexed transcript IDs and the long transcript IDs.
#@param col.names Name of the columns of the BLAST File
#' Load BLAST Hits into a data.frame
#'
#' Give the path to a BLAST Hits file to load it into a data.frame(BLAST HITs Table). The column names can be provided as col.names. Rows with NAs are automatically removed.
#' 
#' Note : If use.feather is enabled, then the functions returns an arrow::RecordBatchStreamReader object. Call wrap it in an iterators::iter() (like batch_iter <- iter(function(){ batch_reader$read_next_batch() })) and call iterators::nextElem(batch_iter) to get each batch of the BLAST Hits.
#'
#' Optional : You can provide the file/table with indexed Transcsript IDs. The format must be "file"[tab]"long_id"[tab]"index" (without a header). It can be generated with the  index_FASTA_IDs() (check ?index_FASTA_IDs or index_fastaIDs() in fs::path_package("COMPLETE","exec","functions.sh"))
#'
#' @examples
#'     LoadBLASTHits(infile, transcript_ID_metadata=NULL, col.names=NULL)
#'
#' @note Column indices will be used when processing the data. Column names are only for user reference. First column must be the query sequence ID and the second column must be the subject sequence ID
#'
#' @param infile BLAST hits filename (not a connection) (Gzipped files supported)
#' @param sep Delimiter of the BLAST File columns. Default - '\\t'
#' @param header Does the file have a header? . Default - FALSE
#' @param use.feather Should feather API be used to read BLAST Hits? - Default - FALSE
#' @return Data Frame with BLAST Results
#' @export
LoadBLASTHits <- function(infile, sep="\t", header = F,use.feather=F){ # transcript_ID_metadata=NULL, col.names=NULL,  ##,sep2 = c("","|","") # gzipped=F

  # if(!is.null(transcript_ID_metadata) && is.character(transcript_ID_metadata)){
  #   transcript_ID_metadata <- read.table(file = transcript_ID_metadata,header = F,sep="\t",quote = "")
  # }
  if(try(file.exists(infile))&& file.info(infile)$size > 0){ #any(grepl(x = class(infile), pattern = "gzfile|connection", ignore.case = T)) #
    #if(gzipped){
    #  infile <- gzfile(description = infile, open = "r")
    #}
    #blast_results <- read.table(file = infile,header = header,sep=sep,quote = "", blank.lines.skip = T, fill = t,na.strings = NA)
    if(!use.feather){
      blast_results <- iterators::iread.table(file = infile,row.names = NULL,header = header,sep=sep,quote = "", blank.lines.skip = T, fill=T,na.strings="NA") #data.table::fread(file = infile,header = header,sep=sep,quote = "", blank.lines.skip = T, nThread = n_threads)
      return(blast_results)
    }else{
      arrow_lfs <- arrow::LocalFileSystem$create()
      arrow_i_stream <- arrow_lfs$OpenInputStream(infile)
      batch_reader <- arrow::RecordBatchStreamReader$create(arrow_i_stream)
      #blast_results <- arrow::read_feather(file = infile, mmap = T)
      return(batch_reader) #batch_iter <- iter(function())
    }
    
  }else{
    stop(paste("File",infile,"does not exist or size 0"))
  }
  # if(!is.null(col.names)){
  #   if(ncol(blast_results) == length(col.names)){
  #     colnames(blast_results) <- col.names
  #   }else{
  #     stop(paste("Number of columns in BLAST file and length of given column names are not equal!",ncol(blast_results),"!=",length(col.names)))
  #   }
  # }

  # if(!is.null(transcript_ID_metadata)){
  #   transcript_ID_metadata <- as.data.frame(transcript_ID_metadata)
  #   colnames(transcript_ID_metadata) <- c("file","id","index")
  # 
  #   query_long_IDs <- transcript_ID_metadata[match(blast_results[,c(1)],transcript_ID_metadata$index),c("id")]
  #   subject_long_IDs <- transcript_ID_metadata[match(blast_results[,c(2)],transcript_ID_metadata$index),c("id")]
  # 
  #   blast_results[,c(1)] <- query_long_IDs
  #   blast_results[,c(2)] <- subject_long_IDs
  # }
  # 
  # blast_results <- blast_results[apply(blast_results, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),] #REMOVING NAs
  # 
  # return(blast_results)
  return(NULL)
}

#@param sep2 Delimiter 2 of the BLAST File columns. Default - c("","|",""). Check ?data.table::fwrite or ?data.table::fread
#' Create A GRanges Object from BLAST RESULTS
#'
#' Functions accepts a filename or a BLAST table and converts it into a GenomicRanges::GRanges object (which can be used by R-WISARD). The accepted BLAST format is 6 and can be formatted from BLAST format 11 using convert_BLAST_format (shown below). Columns indices can provided with the col.indices option. BLAST table/file must be non-indexed (IDs not shortened) if COMPLETE.format.ids=T
#'
#' @note
#' * Output GRanges Object uses column names different from BLAST column names
#' * Require column indices of "qseqid","sseqid","pident","length","qstart","qend","sstart","send","evalue","bitscore","gaps","frames","qcovhsp","sstrand","qlen","slen"
#' * Convert to compatible BLAST format with convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive")) (in_file must be of BLAST format 11)
#'
#' @examples
#'     GRObject_from_BLAST(blast_input = in_file, COMPLETE.format.ids = F, col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list=NULL)
#'
#' @param blast_input Filename with BLAST hits or a BLAST table (not a connection) (Gzipped files supported)
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN)
#' @param col.indices A Named List with indices of columns Query sequence ID (qseqid), Subject sequence ID (sseqid), E Value (evalue), Query start (qstart), Query end (qend), Subject start (sstart), Subject end (send), Bitscore (bitscore), Query HSP coverage(qcovhsp), Query length (qlen), Subject length(slen), Frames (frames), Percentage Identity (pident), Number of Gaps (gaps), Alignment Length (length), Subject strand (sstrand). Default values are col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17)
#' @param params_list Output of load_params() (Optional - Only give it when using Files generated by R-COMPLETE and if COMPLETE.format.ids=TRUE)
#' @param sep Delimiter of the BLAST File columns. Default - '\\t'
#' @param header Does the file have a header? . Default - FALSE
#' @param use.feather Should feather API be used to read BLAST Hits? - Default - FALSE
#' @return A GRanges object of BLAST hits
#' @export
GRObject_from_BLAST <- function(blast_input, COMPLETE.format.ids=F, col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list=NULL, sep="\t", header=F, use.feather=F){ #sep2=c("","|","")

  # if(!is.null(blast_input) && is.character(blast_input) &&){
  #   blast_input <- LoadBLASTHits(blast_input,sep=sep,header=header)
  # }
  #print(head(blast_input)) #DEBUG
  # if(!is.null(col.names)){
  #   if(ncol(blast_input) != length(col.names)){
  #     stop("ncol(blast_input) != length(col.names)")
  #   }
  #   colnames(blast_input) <- col.names
  # }

  if(!is.null(blast_input) && is.character(blast_input)){ # any(grepl(x = class(blast_input), pattern = "gzfile|connection", ignore.case = T)) ||
    if(file.exists(blast_input) && file.info(blast_input)$size > 0){ #any(grepl(x = class(blast_input), pattern = "gzfile|connection", ignore.case = T)) ||
      blast_input <- LoadBLASTHits(infile = blast_input, sep=sep, header=header, use.feather=use.feather)
    }else{
      stop(paste(blast_input,"does not exist!"))
    }
  }else if(!is.null(blast_input) && !is.character(blast_input)){
    blast_input <- blast_input
    #if(!is.null(col.names)){
    #  colnames(blast_table) <- col.names
    #}
  }else{
    stop(paste("Input does not exist or is not a table!"))
  }

  req_columns <- c("qseqid","sseqid","qstart","qend","sstart","send", "sstrand") #"pident","length", "evalue","bitscore","gaps","frames","qcovhsp","sstrand","qlen","slen"

  if (any(is.na(match(req_columns, names(col.indices))))) {
    stop(paste("Missing columns, Require indices of :",paste(req_columns,collapse = ",")))
  }

  if(COMPLETE.format.ids && !is.null(params_list)){
    blast_input <- blast_input %>% dplyr::mutate(subject_gene=unlist(purrr::map(blast_input[,col.indices[["sseqid"]]],function(x){
      #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      genes <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(genes[,COMPLETE_env$FORMAT_ID_INDEX$GENE])
    })) )

    blast_input <- blast_input %>% dplyr::mutate(query_gene=unlist(purrr::map(blast_input[,col.indices[["qseqid"]]],function(x){
      #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      genes <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(genes[,COMPLETE_env$FORMAT_ID_INDEX$GENE])
    })) )

    blast_input <- blast_input %>% dplyr::mutate(subject_transcript_id=unlist(purrr::map(blast_input[,col.indices[["sseqid"]]],function(x){
      #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      tx_id <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(stringi::stri_split_fixed(tx_id[,COMPLETE_env$FORMAT_ID_INDEX$TRANSCRIPT_ID],pattern = params_list$TRANSCRIPT_ID_DELIM, simplify=T)[,1])
    })) )

    blast_input <- blast_input %>% dplyr::mutate(query_transcript_id=unlist(purrr::map(blast_input[,col.indices[["qseqid"]]],function(x){
      #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      tx_id <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(stringi::stri_split_fixed(tx_id[,COMPLETE_env$FORMAT_ID_INDEX$TRANSCRIPT_ID],pattern = params_list$TRANSCRIPT_ID_DELIM, simplify=T)[,1])
    })) )

    if(!is.null(col.indices[["evalue"]])) blast_input <- blast_input[which(blast_input[,col.indices[["evalue"]]] < params_list$E_VALUE_THRESH),]

    if(nrow(blast_input) == 0){
      return(stop("GRObject_from_BLAST() - No hits passed the E-Value threshold!."))
    }

  }else{
    warning(paste("is.null(params_list)==",is.null(params)," && COMPLETE.format.ids==",COMPLETE.format.ids,". Must be COMPLETE.format.ids==TRUE && is.null(params_list)==FALSE",sep=""))
  }

  ##CHANGE strands otherwise IRanges will not work
  change_qstrand <- which(blast_input[,col.indices[["qstart"]]] > blast_input[,col.indices[["qend"]]])
  if (length(change_qstrand) > 0) {
    # print("changing qstart")
    tmp <- blast_input[change_qstrand,col.indices[["qstart"]]]
    blast_input[change_qstrand,col.indices[["qstart"]]] <- blast_input[change_qstrand,col.indices[["qend"]]]
    blast_input[change_qstrand,col.indices[["qend"]]] <- tmp
    tmp <- NULL
  }

  change_sstrand <- which(blast_input[,col.indices[["sstart"]]] > blast_input[,col.indices[["send"]]])
  if (length(change_sstrand) > 0) {
    # print("changing sstart")
    tmp <- blast_input[change_sstrand,col.indices[["sstart"]]]
    blast_input[change_sstrand,col.indices[["sstart"]]] <- blast_input[change_sstrand,col.indices[["send"]]]
    blast_input[change_sstrand,col.indices[["send"]]] <- tmp
    tmp <- NULL
  }

  no_strand_info <- which(stringi::stri_cmp_eq("N/A",blast_input[,col.indices[["sstrand"]]]))
  blast_input[no_strand_info,col.indices[["sstrand"]]] <- "*"

  #print(blast_input) #DEBUG
  #tmp_gr <- GenomicRanges::GRanges(S4Vectors::Rle(blast_input[,col.indices["sseqid"]]), ranges =  IRanges::IRanges(blast_input[,col.indices["sstart"]], end = blast_input[,col.indices["send"]],strand= blast_input[,col.indices["sstrand"]])) #, names = orths$sseqid))

##NEW CODE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! - Chat GPT
  tmp_df <- data.frame(seq_names = blast_input[, col.indices[["sseqid"]]])
  
  tmp_df <- tmp_df %>%
    dplyr::mutate(sstart = blast_input[, col.indices[["sstart"]]],
           send = blast_input[, col.indices[["send"]]],
           sstrand = blast_input[, col.indices[["sstrand"]]],
           Hsp_num = 1:nrow(blast_input))
  
  if (!is.null(col.indices[["bitscore"]]))
    tmp_df <- tmp_df %>% dplyr::mutate(Hsp_bit.score = blast_input[, col.indices[["bitscore"]]])
  
  if (!is.null(col.indices[["qcovhsp"]]))
    tmp_df <- tmp_df %>% dplyr::mutate(Hsp_score = blast_input[, col.indices[["qcovhsp"]]])
  
  if (!is.null(col.indices[["evalue"]]))
    tmp_df <- tmp_df %>% dplyr::mutate(Hsp_evalue = blast_input[, col.indices[["evalue"]]])
  
  tmp_df <- tmp_df %>%
    dplyr::mutate(subject_HSP_from = blast_input[, col.indices[["sstart"]]],
           subject_HSP_to = blast_input[, col.indices[["send"]]],
           query_id = blast_input[, col.indices[["qseqid"]]])
  
  if (!is.null(col.indices[["qlen"]]))
    tmp_df <- tmp_df %>% dplyr::mutate(query_len = blast_input[, col.indices[["qlen"]]])
  
  if (!is.null(col.indices[["slen"]]))
    tmp_df <- tmp_df %>% dplyr::mutate(subject_len = blast_input[, col.indices[["slen"]]])
  
  tmp_df <- tmp_df %>%
    dplyr::mutate(query_HSP_from = blast_input[, col.indices[["qstart"]]],
           query_HSP_to = blast_input[, col.indices[["qend"]]])
  
  if (!is.null(col.indices[["frames"]])) {
    tmp_df <- tmp_df %>%
      dplyr::mutate(query.frame = furrr::future_map_chr(as.character(blast_input[, col.indices[["frames"]]]), ~as.integer(unlist(stringi::stri_split_fixed(.x, pattern = "/")))[1], .options = furrr::furrr_options(seed = TRUE, scheduling=F)),
             subject.frame = furrr::future_map_chr(as.character(blast_input[, col.indices[["frames"]]]), ~as.integer(unlist(stringi::stri_split_fixed(.x, pattern = "/")))[2], .options = furrr::furrr_options(seed = TRUE, scheduling=F)))
  }
  
  if (COMPLETE.format.ids && !is.null(params_list)) {
    tmp_df <- dplyr::mutate(tmp_df, query_org = unlist(furrr::future_map(blast_input[,col.indices[["qseqid"]]],function(x){
      #org <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      org <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(org[,COMPLETE_env$FORMAT_ID_INDEX$ORG])
    }, .options = furrr::furrr_options(scheduling=F))) )
    tmp_df <- dplyr::mutate(tmp_df, subject_org = unlist(furrr::future_map(blast_input[,col.indices[["sseqid"]]],function(x){
      #org <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      org <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(org[,COMPLETE_env$FORMAT_ID_INDEX$ORG])
    }, .options = furrr::furrr_options(scheduling=F))) )
  }
  
  if (!is.null(col.indices[["pident"]]))
    tmp_df <- tmp_df %>% dplyr::mutate(pidentity = blast_input[, col.indices[["pident"]]])
  
  if (!is.null(col.indices[["gaps"]]))
    tmp_df <- tmp_df %>% dplyr::mutate(Hsp_gaps = blast_input[, col.indices[["gaps"]]])
  
  if (!is.null(col.indices[["length"]]))
    tmp_df <- tmp_df %>% dplyr::mutate(Hsp_align.len = blast_input[, col.indices[["length"]]])
  
  if (COMPLETE.format.ids && !is.null(params_list)) {
    tmp_df <- tmp_df %>%
      dplyr::mutate(subject_gene = blast_input[, c("subject_gene")],
             query_gene = blast_input[, c("query_gene")],
             query_transcript_id = blast_input[, c("query_transcript_id")],
             subject_transcript_id = blast_input[, c("subject_transcript_id")])
  }
  
  df_cols <- names(unlist(col.indices)[order(unlist(col.indices))])
  colnames(tmp_df) <- df_cols
  print(tmp_df) #DEBUG
  tmp_gr <- GenomicRanges::makeGRangesFromDataFrame(df = tmp_df,seqnames.field = "sseqid",start.field = "sstart",end.field = "send",strand.field = "sstrand" ,keep.extra.columns = T,ignore.strand = F)

  return(tmp_gr)
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
    invisible(processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"index_fastaIDs", index_out,path) ,spinner = T,stdout = "",stderr = ""))
  }else{
    stop("Give path and output file for index")
  }
}

#' Shorten FASTA IDs
#'
#' WARNING - This will convert all the long-format FASTA IDs into short-format(indices) in the files for downstream compatibility. This can be reverted back with Lengthen FASTA IDs
#'
#' @param transcript_metadata Transcript Metadata file (from index_FASTA_IDs())
#' @param fasta_files Vector of FASTA files to shorten IDs (Optional)
#' @export
shorten_FASTA_IDs <- function(transcript_metadata, fasta_files=c()){
  warning("WARNING - This will convert all the long-format FASTA IDs into short-format(indices) for downstream compatibility. This can be reverted back with lengthen_FASTA_IDs()")
  if (!is.null(transcript_metadata)){
    install_parallel()
    invisible(processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"shorten_fastaIDs", transcript_metadata, COMPLETE_env$parallel, paste(fasta_files,collapse = " ")) ,spinner = T,stdout = "",stderr = ""))
  }else{
    stop("Give transcript metadata file (can be generated with index_FASTA_IDs())")
  }
}

#' Lengthen FASTA IDs
#'
#' WARNING - This will convert all the short-format(indices) into long-format FASTA IDs in the files.
#'
#' @param transcript_metadata Transcript Metadata file (from index_FASTA_IDs())
#' @param fasta_files Vector of FASTA files (Optional)
#' @export
lengthen_FASTA_IDs <- function(transcript_metadata, fasta_files=c()){
  warning("WARNING - This will convert all the short-format(indices) into long-format FASTA IDs in the files.")
  if (!is.null(transcript_metadata)){
    install_parallel()
    invisible(processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"lengthen_fastaIDs", transcript_metadata, COMPLETE_env$parallel, paste(fasta_files,collapse = " ")) ,spinner = T,stdout = "",stderr = ""))
  }else{
    stop("Give transcript metadata file (can be generated with index_FASTA_IDs())")
  }
}

#' Convert between BLAST output formats
#'
#' Convert between formats using the blast_formatter program
#'
#' @examples
#'     convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'
#' @note blast_formatter must be present in $PATH and be accessible with Sys.which("blast_formatter"). (Optional) Input file must be of the BLAST format 11 to convert it between formats
#'
#' @param infile Input BLAST Hits File
#' @param outfile Output BLAST Hits File
#' @param outformat Format to convert to (Default 6)
#' @param cols BLAST Columns to output. Check BLAST Output formats for more details
#' @param conversion_prg_dir Path to directory of blast_formatter if not found in $PATH
#' @param verbose Print DEBUG Messages?
#' @export
convert_BLAST_format <- function(infile, outfile,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"),conversion_prg_dir = dirname(Sys.which("blast_formatter")), verbose=F){
  if(stringi::stri_isempty(conversion_prg_dir)){
    stop(paste("convert_BLAST_format() - blast_formatter not found in $PATH or",conversion_prg_dir ,"..cannot continue!"))
  }


  if(verbose){
    cmd_verbose <- ""
  }else{
    cmd_verbose <- NULL
  }

  invisible(processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"convert_BLAST_format",conversion_prg_dir,infile,outfile,outformat,paste(cols,collapse = " ") ) ,spinner = T,stdout = cmd_verbose,stderr = cmd_verbose))

}

#' Make BLAST DB
#'
#' @param fasta_file Path to FASTA File
#' @param blast_bin Path to directory of BLAST+ programs makeblastdb and blastdb_path
#' @param clean_extract Remove existing BLAST DB?
#' @param verbose Print DEBUG Messages?
#' @export
make_BLAST_DB <- function(fasta_file,blast_bin=dirname(Sys.which("tblastx")),clean_extract=T, verbose=F){
  if(file.exists(fasta_file)){
    if(verbose){
      cmd_verbose <- ""
    }else{
      cmd_verbose <- NULL
    }

    db_check <- processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"check_DB",fasta_file,blast_bin), stdout = "")

    if(is.null(db_check$stdout) || clean_extract){
      invisible(processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"make_BLAST_DB",fasta_file, blast_bin) ,spinner = T,stdout = cmd_verbose,stderr = cmd_verbose))
    }
  }else{
    if(verbose){ stop(paste("make_BLAST_DB() - Missing File - ",fasta_file)) }
    #return(fasta_file)
  }
}

#@param return_f_callback Function to execute on returning data. Set NULL for discarding returning data
#@param ... Named List. Parameters for return_f_callback(list(BLAST_hits,BLAST_file,run_name,seed),...)
#' Perform Nucleotide BLAST
#'
#' BLAST between two FASTA files. BLAST output is saved temporarily in BLAST format 11, internally converted to format 6 with convert_BLAST_format() and BLAST Hits are returned as a GRanges Object. If blast_out is provided then the intermediate BLAST Format 11 output can be found in the same directory with the same name and *.blast11 extension. If blast_DB_dir is provided then the Subject FASTA is copied into this directory and then BLASTed.
#'
#' @note ASSUMES Nucleotide FASTA sequences. Files are copied to tempdir() and split by blast.sequence.limit before blasting. Not checking/Not working for Protein/Peptide sequences and FASTQ files. Columns of GRanges BLAST 6 output are c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive")
#'
#' @examples
#'     run_BLAST(query_path = "query.fasta",subject_path = "subject.fasta",blast_DB_dir = "files/blastdb", blast_program="tblastx", run_name = "blast_positive",blast_options = "-strand plus", gzip.output=F, use.feather=T)
#'     run_BLAST(query_path = Biostrings::DNAString(gsub("[\r\n]", "", "ATTGTCGAAGTTGTCGCTCGAGAGGCGGGAGTTTACCGACACTTTTCCTCAGAAGTTTAC
#'     CGTGAAGCTGACCGGAGAACGGCGGAGTCTGTGCTGAATCTGCCATCATGTCCAGGCGGA
#'     ")), subject_path=Biostrings::DNAString(gsub("[\r\n]", "", "TTACCGACACTTTTCCTCAGAAGTTTACCGTGAAGCTGACCGGAGAACGGCGGAGTCTGT
#'     GCTGAATCTGCCATCATGTCCAGGCGGAGCTCACTCATAGTGCCGATGAAGAGTATTGAG
#'     ")))
#'
#' @param query_path Path to Query FASTA or Biostrings::DNAString|AAString|RNAString
#' @param subject_path Path to Subject FASTA or or Biostrings::DNAString|AAString|RNAString
#' @param blast_DB_dir Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed. Default is NULL
#' @param blast_out Path to BLAST output file, Default BLAST FORMAT is 11. It is converted internally to BLAST Format 6 and returned as a GRanges Object. Default is tempdir(). Can also be NULL
#' @param blast_program Give path to the BLAST program. eg, Default - Sys.which("tblastx") if tblastx is in SHELL $PATH.
#' @param run_name Name of the BLAST run. Only for logging (Optional)
#' @param blast_options Extra Options to be passed to the BLAST program
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE (Default) otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (Optional)
#' @param params_list Output from load_params() (Optional)
#' @param blast.sequence.limit Maximum number of sequences allowed in each BLAST file. Default - 0. If the query FASTA sequences > blast.sequence.limit, the sequences are split into multiple files and BLASTed. Use 0 to not split & copy the files into temporary files
#' @param n_threads Number of threads. Default - 8
#' @param clean_extract Delete Output file if exists? (Optional). Default - F
#' @param verbose Print DEBUG Messages?
#' @param gzip.output Should the output files be Gzipped? Default - FALSE
#' @param seed Seed value for reproducibility. Default - 123
#' @param return_data Logical. Should GRanges Object of BLAST Hits be returned? Default - TRUE
#' @param use.feather Should feather API be used to read BLAST Hits? - Default - TRUE
#' @param header Does the BLAST format output have header? (Some do, check it, by default fmt 6 is used which does not have a header) - FALSE
#' @return BLAST Hits as GRanges Object
#' @export
run_BLAST <- function(query_path, subject_path,blast_DB_dir = file.path(tempdir(),run_name), blast_out=NULL, blast_program=Sys.which("tblastx"), run_name="BLAST",blast_options="", COMPLETE.format.ids=F,params_list=NULL, blast.sequence.limit=0,n_threads=8, gzip.output=F, clean_extract=F,verbose=F, seed=123, return_data=TRUE, use.feather=T, header=F){ #return_f_callback=return,... #keep.output.files=T
set.seed(seed)
  tictoc::tic(msg = paste(run_name,":"))

  if(verbose){
    cmd_verbose <- ""
  }else{
    cmd_verbose <- NULL
  }

  if(use.feather && gzip.output){
    warning("run_BLAST: Both use.feather && gzip.output cannot be TRUE. Setting use.feather=FALSE")
    use.feather=F
  }
  
  # if(!is.null(return_f_callback) && class(return_f_callback)=="function"){
  #   callback_args <- rlang::dots_list(..., .named = T)
  # }

  if(is.null(blast_out) && gzip.output==T){ #keep.output.files==F
    warning("is.null(blast_out) && gzip.output==T. Output will be saved temporarily.")
    gzip.output=FALSE
  }else if(grepl(pattern = "gz", x = blast_out,ignore.case = T)){
    gzip.output=TRUE
  }

  if(!is.null(blast_out)){
    if(any(all(file.exists(sprintf("%s.%s", blast_out, "gz")), file.info(sprintf("%s.%s", blast_out, "gz"))$size >0), all(file.exists(blast_out), file.info(blast_out)$size >0, !clean_extract))){
      if(file.exists(sprintf("%s.%s", blast_out, "gz")) && file.info(sprintf("%s.%s", blast_out, "gz"))$size >0) {
        unlink(x = blast_out,recursive = F,force = T,expand = T)
        blast_out <- sprintf("%s.%s", blast_out, "gz")
        }

      #return(blast_GR)
      if(verbose){
        message(paste(blast_out,"exists!"))
      }
      cat(print_toc(tictoc::toc(quiet = T, log = T)))

      # if(is.null(return_f_callback) && class(return_f_callback)=="function"){
      #   if(is.null(unlist(callback_args, recursive = F,use.names = T))){
      #     return_f_callback(list(BLAST_hits=blast_GR,BLAST_file=blast_outfile,run_name=run_name,seed=seed))
      #   }else{
      #     do.call(what = return_f_callback, args = c(list(BLAST_hits=blast_GR,BLAST_file=blast_out,run_name=run_name,seed=seed),callback_args))
      #   }
      #
      # }else{
      #   return(NULL)
      # }
      if(return_data){
        blast_GR <- GRObject_from_BLAST(blast_input = blast_out,COMPLETE.format.ids = COMPLETE.format.ids,col.indices = c(qseqid = 1, sseqid = 2, evalue = 11, qstart = 7, qend = 8, sstart = 9, send = 10, bitscore = 12, qcovhsp = 16, qlen = 18, slen = 19, frames = 15, pident = 3, gaps = 14, length = 4, sstrand = 17), params_list = params_list, use.feather = use.feather)
        #blast_GR <- LoadBLASTHits(blast_out, , use.feather = use.feather)
        return(list(BLAST_hits=blast_GR,BLAST_file=blast_out,run_name=run_name))
        }else{
        return(NULL)
      }
    }else{
      unlink(c(blast_out,sprintf("%s.%s", blast_out, "gz")), recursive = F,force = T,expand = T)
    }
  }
  unlink(blast_DB_dir,force = T,recursive = T,expand = T)
  dir.create(blast_DB_dir,showWarnings = F)
  tryCatch({
    fasta_in_paths <- parallel::mclapply(list(query_path,subject_path),function(x){
      if(any(grepl(x = class(x),pattern = "DNAString|AAString|RNAString",ignore.case = T,fixed = F))){
        if(blast.sequence.limit > 0){
          num_seqs <- length(x)
          x <- split(Biostrings::DNAStringSet(x), ceiling(seq_along(x) / blast.sequence.limit))
          num_files <- length(x)
          if(verbose){
            cat(paste("Sequence Count :", num_seqs,", File Count :", num_files))
          }
        }
        tmp_fasta_file <- tempfile(pattern = paste("tmp.",seq(1:length(x)),".",sep=""),fileext = ".fa", tmpdir = blast_DB_dir) #tempfile(pattern = "tmp",fileext = ".FASTA")
        return( unlist( parallel::mclapply(seq_along(x), function(y){
          Biostrings::writeXStringSet(x = Biostrings::DNAStringSet(x[[y]]),filepath = tmp_fasta_file[[y]], append = F,format = "fasta")
          return(tmp_fasta_file[[y]])
        }, mc.preschedule = T,mc.cores = n_threads, mc.set.seed = seed) ) )
        #return(tmp_fasta_file)
      }else{
        if(blast.sequence.limit > 0){
          x <- Biostrings::readDNAStringSet(x)
          num_seqs <- length(x)
          x <- split(Biostrings::DNAStringSet(x), ceiling(seq_along(x) / blast.sequence.limit))
          num_files <- length(x)
          if(verbose){
            cat(paste("\nSequence Count :", num_seqs,", File Count :", num_files))
          }
          tmp_fasta_file <- tempfile(pattern = paste("tmp.",seq(1:length(x)),".",sep=""),fileext = ".fa",tmpdir = blast_DB_dir) #tempfile(pattern = "tmp",fileext = ".FASTA")
          return( unlist( parallel::mclapply(seq_along(x), function(y){
            Biostrings::writeXStringSet(x = Biostrings::DNAStringSet(x[[y]]),filepath = tmp_fasta_file[[y]], append = F,format = "fasta")
            return(tmp_fasta_file[[y]])
          }, mc.preschedule = T,mc.cores = n_threads, mc.silent = !verbose, mc.set.seed = seed) ))
        }else{
          return(x)
        }
      }
    }, mc.preschedule = T,mc.cores = n_threads, mc.silent = !verbose, mc.set.seed = seed)
  }

  ,error=function(cond){
    stop(cond)
  } )

  query_path <- fasta_in_paths[[1]]
  subject_path <- fasta_in_paths[[2]]

  tryCatch({
    #if (stringi::stri_isempty(COMPLETE_env$parallel)) {
    #  stop("Problem with GNU parallel installation. Reload R-COMPLETE")
    #}

    if(stringi::stri_isempty(blast_program)){
      stop("BLAST+ not found in $PATH. Provide blast_program")
    }else{
      BLAST_BIN <- dirname(blast_program)
    }

    #subject_path <- tools::file_path_as_absolute(subject_path)
    #if (!is.null(blast_DB_dir)) {
      #dir.create(blast_DB_dir,showWarnings = F,recursive = T)
      #query_DB <- tools::file_path_as_absolute(paste(blast_DB_dir,"/",basename(query_path),sep=""))
      subject_DB <- file.path(blast_DB_dir,basename(subject_path)) #tools::file_path_as_absolute()
      #if(!stringi::stri_cmp_eq(query_path,query_DB)){
      #  file.copy(query_path,query_DB,overwrite = T)
      #}
      if(!all(stringi::stri_cmp_eq(subject_path,subject_DB))){
        file.copy(subject_path, subject_DB, overwrite = T)
      #invisible(tryCatch(file.copy(subject_path,subject_DB,overwrite = T), error=function(cond){
      #  #stop(cond)
      #}))
      }
      #query_path <- query_DB
      #if(verbose){
      # print(subject_DB)
      #}
      subject_path <- subject_DB
    #}

    blast_outfile <- blast_out
    #if(is.null(blast_out)){
    if(!is.null(params_list)){
      blast_out <- tempfile(pattern="blast_out", tmpdir = params_list$TEMP_PATH)
    }else{
      blast_out <- tempfile(pattern="blast_out", tmpdir = tempdir())
    }
    #}
    final_blast_out <- paste(blast_out,seq(1:length(subject_path)),"arrow",sep=".")
    #blast_out <- paste(tools::file_path_sans_ext(final_blast_out),".blast11",sep="")

    if(verbose){
      print(query_path)
      print(subject_path)
      #print(blast_out)
      print(final_blast_out)
    }

    #print(length(query_path))
    #print(length(subject_path))

    #print(final_blast_out)

    #MAKE BLAST DB of FASTA files
    #Only subject fasta files needs to be a BLAST DB
    #processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"make_BLAST_DB",query_path, dirname(blast_program)) ,spinner = T,stdout = "",stderr = "")
    invisible( furrr::future_map(subject_path,.f = function(x){
      make_BLAST_DB(fasta_file=x,blast_bin=BLAST_BIN, verbose=verbose, clean_extract = T)
    }, .options = furrr::furrr_options(seed = seed, scheduling=F)) ) #n_threads

    path_combinations <- unique(tidyr::crossing(query_path,subject_path))
    #print(path_combinations) #DEBUG
    
    #tmp_list <- list()#
    arrow_lfs <- arrow::LocalFileSystem$create()
    arrow_ao_stream <- arrow_lfs$OpenAppendStream(blast_outfile) #arrow::FileOutputStream$create(blast_outfile)
    furrr::future_map(.x = seq_along(1:nrow(path_combinations)), .f = function(i){
    #parallel::mclapply(seq_along(1:nrow(path_combinations)), function(i){
      q_x <- as.character(path_combinations[i,1])
      s_y <- as.character(path_combinations[i,2])
      
      #processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"do_BLAST",COMPLETE_env$parallel,run_name,q_x,s_y,final_blast_out[i],blast_program,n_threads,blast_options) ,spinner = T,stdout = "",stderr = cmd_verbose) #stdout = cmd_verbose
      tmp_proc <- processx::process$new(command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"do_BLAST",COMPLETE_env$parallel,run_name,q_x,s_y,final_blast_out[i],blast_program,n_threads,blast_options),stdout = "|", stderr = cmd_verbose,supervise = T) 
      invisible(tmp_proc$poll_io(timeout=5)) ##DONT POLL without timeout when there is no output (or) when using linux fifo
      
      if(all(tmp_proc$is_alive(), tmp_proc$has_output_connection())){ # , !tmp_proc$is_incomplete_output()
      #  #tmp_proc$wait(timeout=5)
      #  #tmp_list <<- append(tmp_list,read.table(textConnection(processx::processx_conn_read_lines(tmp_proc$get_output_connection()),open = "r")))
      #  #tmp_list <- append(tmp_list,read.table(fifo(description = paste(final_blast_out[i],"_pipe",sep=""), open = "r",blocking = F))) #textConnection(readLines(tmp_proc$get_output_connection(),open = "r"))
      #tmp_read <- arrow::BufferReader$create(read.table(fifo(description = paste(final_blast_out[i],"_pipe",sep=""), open = "r",blocking = F)))
      #read.table(readLines(fifo(description = paste(final_blast_out[i],"_pipe",sep=""), open = "r",blocking = T))) %>% arrow::write_dataset(path=paste(final_blast_out[i],"_arrow",sep=""), format = "feather")
        #print(paste(final_blast_out[i],"_arrow",sep=""))
        
        file_pipe <- fifo(description = final_blast_out[i], open = "r", blocking = T) #paste(final_blast_out[i],"_pipe",sep="")
        out_tbl <- arrow::Table$create(read.table( file_pipe , row.names = NULL, header = header, fill = T, na.strings = "NA"))
      arrow::write_feather(x= out_tbl, sink =arrow_ao_stream, compression = "lz4")  #fill = T, na.strings = "NA"
      
      }
      #if(verbose){
      #  message(paste("Converting",blast_out[i],"to BLAST Format 6 :", final_blast_out[i]))
      #}
      #convert_BLAST_format(infile = blast_out[i],outfile = final_blast_out[i],conversion_prg_dir = BLAST_BIN, verbose = verbose ) 

    }, .options = furrr::furrr_options(seed = seed, scheduling=F))#, mc.cores = n_threads,mc.silent = !verbose,mc.preschedule = T , mc.set.seed = seed)#, .options = furrr::furrr_options(seed = TRUE, scheduling=n_threads))
    #tmp_list <- dplyr::bind_rows(tmp_list)
    
    #}
    arrow_ao_stream$close()
   
    # tmp_o_stream <- arrow::FileOutputStream$create(blast_outfile)
    #purrr::map(final_blast_out,.f = function(b_out){
    #  arrow::write_feather(x=as.data.frame(arrow::read_feather(paste(b_out,"_arrow",sep=""))), sink = tmp_o_stream,compression = "lz4")  
    #  unlink(x = c(b_out,paste(b_out,"_arrow",sep="")),force = T,expand = T)
    #})
    #tmp_o_stream$close()
    
    #print("here2")

    # if(is.null(blast_outfile)){
    #   if(!is.null(params_list)){
    #     blast_outfile <- tempfile(pattern="blast_out", tmpdir = params_list$TEMP_PATH)
    #   }else{
    #     blast_outfile <- tempfile(pattern="blast_out", tmpdir = tempdir())
    #   }
    # }
    #print(blast_outfile) #DEBUG
    #print(final_blast_out) #DEBUG
    #processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"cat_files",blast_outfile, paste(final_blast_out,collapse = " ") ) ,spinner = T,stdout = cmd_verbose,stderr = cmd_verbose)

    unlink(x = c(final_blast_out,blast_out, blast_DB_dir),recursive = T,force = T,expand = T)

    if(return_data){ 
    blast_GR <- GRObject_from_BLAST(blast_input = blast_outfile,COMPLETE.format.ids = COMPLETE.format.ids,col.indices = c(qseqid = 1, sseqid = 2, evalue = 11, qstart = 7, qend = 8, sstart = 9, send = 10, bitscore = 12, qcovhsp = 16, qlen = 18, slen = 19, frames = 15, pident = 3, gaps = 14, length = 4, sstrand = 17), params_list = params_list, use.feather = use.feather)
    #blast_GR <- LoadBLASTHits(blast_outfile, use.feather = use.feather)
    if(gzip.output){
         gzip_outfile <- sprintf("%s.%s", blast_outfile, "gz")
      R.utils::compressFile(filename=blast_outfile, destname=gzip_outfile, ext="gz", temporary=FALSE, skip=TRUE, overwrite=FALSE, remove=TRUE, FUN=gzfile)
      blast_outfile <- gzip_outfile
    }#else{
    # unlink(x = blast_outfile,recursive = T,force = T,expand = T)
    #}
    cat(print_toc(tictoc::toc(quiet = T, log = T)))
    return(blast_GR)
    }else{
      if(gzip.output){
        gzip_outfile <- sprintf("%s.%s", blast_outfile, "gz")
        R.utils::compressFile(filename=blast_outfile, destname=gzip_outfile, ext="gz", temporary=FALSE, skip=TRUE, overwrite=FALSE, remove=TRUE, FUN=gzfile)
        blast_outfile <- gzip_outfile
      }
      cat(print_toc(tictoc::toc(quiet = T, log = T)))
      return(NULL)
    }
    
  },error=function(cond){
    stop(cond)
  })
}

#@param return_f_callback Function to execute on returning data. Set NULL for discarding returning data - in run_BLAST()
#@param ... Parameters for return_f_callback() - in run_BLAST()
#' Execute one2one BLAST
#'
#' Executes One-to-One BLAST between two lists of organisms/genes/clusters. The BLAST Hits are
#'
#' @examples one2one_BLAST(first_set,second_set,blast_DB_dir=blast_DB_dir,blast_program,output_dir=output_dir, blast_options=blast_options, input_prefix_path=input_prefix_path, params_list=params_list,COMPLETE.format.ids=COMPLETE.format.ids, verbose=verbose, use.feather=use.feather)
#'
#' @note  ASSUMES Nucleotide sequences. Not checking/Not working for Protein/Peptide sequences.
#'
#' @param first_list Vector of PATHS to FASTA
#' @param second_list Vector of PATHS to FASTA
#' @param file_ext File extension of input files. eg- "cds" or "fa"
#' @param run_name Name of the BLAST run and name of the output file
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program. Default is tempdir(). Can also be NULL
#' @param blast_DB_dir Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed
#' @param blast_options Extra Options to be passed to the BLAST program
#' @param output_dir Path to BLAST output
#' @param input_prefix_path If input lists/vectors are filenames, then provide input folder to prefix path
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE (Default) otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (Optional)
#' @param gzip.output Should the output files be Gzipped? Default - FALSE
#' @param blast.sequence.limit Maximum number of sequences allowed in each BLAST file. Default - 0. If the query FASTA sequences > blast.sequence.limit, the sequences are split into multiple files and BLASTed. Use 0 to not split & copy the files into temporary files
#' @param params_list Output of load_params() (Optional)
#' @param clean_extract Delete Output file if exists? (Optional). Default - F
#' @param n_threads Number of Threads
#' @param verbose Print DEBUG Messages?
#' @param seed Seed value for reproducibility. Default - 123
#' @param return_data Logical. Should GRanges Object of BLAST Hits be returned? Default - TRUE
#' @param use.feather Should feather API be used to read BLAST Hits? - Default - FALSE
#' @return A GRObject of BLAST Hits
#' @export
one2one_BLAST <- function(first_list,second_list, file_ext="fa", run_name="BLAST.one2one" ,blast_DB_dir=tempdir(),blast_program,output_dir="./", blast_options="", input_prefix_path=NULL, params_list=NULL,COMPLETE.format.ids=F, gzip.output=F,blast.sequence.limit=0, clean_extract=F, n_threads=8, verbose=F, seed=123,return_data=TRUE, use.feather=FALSE){ #return_f_callback=return,...
set.seed(seed)
  #tictoc::tic(msg = "one2one_BLAST :")

  if(use.feather && gzip.output){
    warning("one2one_BLAST: Both use.feather && gzip.output cannot be TRUE. Setting use.feather=FALSE")
    use.feather=F
  }
  
  if(!is.null(input_prefix_path)){
    first_list <- paste(input_prefix_path,"/",first_list,".",file_ext,sep="")
    #first_list <- grep(pattern = first_list, x= grep(pattern = file_ext_pattern, x = processx::process$new(command = Sys.which("find"), args = c(input_prefix_path,"-type","f"), stderr = NULL, stdout = "|")$read_all_output_lines() , ignore.case = T,value = T), ignore.case = T,value = T) #list.files(path = input_prefix_path,pattern = first_list, full.names = T,recursive = T,include.dirs = F,ignore.case = T)
    second_list <- paste(input_prefix_path,"/",second_list,".",file_ext,sep="")
    #second_list <- grep(pattern = second_list, x= grep(pattern = file_ext_pattern, x = processx::process$new(command = Sys.which("find"), args = c(input_prefix_path,"-type","f"), stderr = NULL, stdout = "|")$read_all_output_lines() , ignore.case = T,value = T), ignore.case = T,value = T)  #list.files(path = input_prefix_path,pattern = second_list, full.names = T,recursive = T,include.dirs = F,ignore.case = T)
  }

  list_combos <- unique(tidyr::crossing(first_list[order(first_list)], second_list[order(second_list)]))

  #return_data <- furrr::future_map2(.x=first_list[order(first_list)], .y=second_list[order(second_list)], .f=function(x,y){
  return_data <- parallel::mclapply(seq_along(1:nrow(list_combos)), function(idx){
    x <- toString(list_combos[idx,1])
    y <- toString(list_combos[idx,2])
    if (file.exists(x) && file.exists(y) && file.info(x)$size > 0 && file.info(y)$size > 0) {
      # run_name1 <- tools::file_path_sans_ext(BiocGenerics::basename(x)) #tools::file_path_as_absolute()
      # run_name2 <- tools::file_path_sans_ext(BiocGenerics::basename(y)) #tools::file_path_as_absolute()
      # run_name <- paste( run_name1,run_name2,"all2all" ,sep=".")
      if(verbose){
        cat(paste(run_name,"\n",sep = ""))
        print(x)
        print(y)
      }
      out_file <- paste( output_dir, "/",run_name,sep="")

      #if(file.exists(out_file) && file.info(out_file)$size > 0){
      #  if(clean_extract){
      #    unlink(x = out_file,recursive = F,force = T,expand = T)
      #  }#else{
      #  #  stop(paste("Output file :",out_file,"exists!"))
      #  #}
      #}

      # return(tryCatch({
      #   return_BLAST <- run_BLAST(query_path = x,subject_path = y,blast_DB_dir = blast_DB_dir, blast_program=blast_program, blast_out = out_file, run_name = run_name,blast_options = blast_options,COMPLETE.format.ids = COMPLETE.format.ids,params_list = params_list, clean_extract=clean_extract,blast.sequence.limit = blast.sequence.limit, n_threads=n_threads, verbose = verbose, gzip.output=gzip.output)
      #   return(return_BLAST)
      #   }, error=function(cond){
      #   message(cond)
      #   return(NULL)
      # }))
      return_BLAST <- NULL
      try(return_BLAST <- run_BLAST(query_path = x,subject_path = y,blast_DB_dir = blast_DB_dir, blast_program=blast_program, blast_out = out_file, run_name = run_name,blast_options = blast_options,COMPLETE.format.ids = COMPLETE.format.ids,params_list = params_list, clean_extract=clean_extract,blast.sequence.limit = blast.sequence.limit, n_threads=n_threads, verbose = verbose, gzip.output=gzip.output, seed=seed,return_data=return_data, use.feather=use.feather)) #, return_f_callback=return_f_callback,...=...
      if(return_data){
        return(return_BLAST)
      }else{
        return(NULL)
      }

    }
  }, mc.cores = n_threads, mc.silent = !verbose) #, .options = furrr::furrr_options(seed = seed, scheduling=F)) # #params_list$numWorkers

  #cat(print_toc(tictoc::toc(quiet = T, log = T)))

  return(return_data)

}

#@param return_f_callback Function to execute on returning data. Set NULL for discarding returning data - in run_BLAST()
#@param ... Parameters for return_f_callback() - in run_BLAST()
#' Execute all2all BLAST
#'
#' Executes All-to-All BLAST between two lists of organisms/genes/clusters. Output BLAST files are bi-directional and are stored in the format filename1.filename2.all2all under output_dir. (All-to-All is simply Many-to-Many association)
#'
#' @examples
#'  all2all_BLAST(first_list = list.files(path="fasta",pattern = "\*.fa",all.files = T,full.names = T,include.dirs = F,recursive = F), second_list = list.files(path="fasta",pattern = "\*.fa",all.files = T,full.names = T,include.dirs = F,recursive = F),blast_DB_dir = "files/blastdb",blast_program = "tblastx",output_dir = "files/all2all",blast_options = "-strand plus",input_prefix_path = NULL,params_list=NULL, use.feather=T,gzip.output=F)
#'
#' @note ASSUMES Nucleotide sequences. Not checking/Not working for Protein/Peptide sequences.
#'
#' @param first_list Vector of PATHS to FASTA
#' @param second_list Vector of PATHS to FASTA
#' @param file_ext File extension of input files. eg- "cds" or "fa"
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program. Default is tempdir(). Can also be NULL
#' @param blast_DB_dir Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed
#' @param blast_options Extra Options to be passed to the BLAST program
#' @param output_dir Path to BLAST output
#' @param input_prefix_path If input lists/vectors are filenames, then provide input folder to prefix path
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE (Default) otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (Optional)
#' @param gzip.output Should the output files be Gzipped? Default - FALSE
#' @param blast.sequence.limit Maximum number of sequences allowed in each BLAST file. Default - 0. If the query FASTA sequences > blast.sequence.limit, the sequences are split into multiple files and BLASTed. Use 0 to not split & copy the files into temporary files
#' @param clean_extract Delete Output file if exists? (Optional). Default - F
#' @param params_list Output of load_params() (Optional)
#' @param n_threads Number of Threads
#' @param verbose Print DEBUG Messages?
#' @param seed Seed value for reproducibility. Default - 123
#' @param return_data Logical. Should GRanges Object of BLAST Hits be returned? Default - TRUE
#' @param use.feather Should feather API be used to read BLAST Hits? - Default - TRUE
#' @return Named List of GRObjects (of BLAST Hits)
#' @export
all2all_BLAST <- function(first_list,second_list,file_ext="fa",blast_DB_dir=tempdir(),blast_program,output_dir="./", blast_options="", input_prefix_path=NULL, params_list=NULL,COMPLETE.format.ids=F, gzip.output=F, blast.sequence.limit=0, clean_extract=F, n_threads=8, verbose=F, seed=123,return_data=TRUE, use.feather=T){ #,return_f_callback=return,...

  # if(is.null(params_list)){
  #   tryCatch(numWorkers <- parallel::detectCores(all.tests = T, logical = T), error=function(){numWorkers <- 2})
  # }else{
  #   #warning("Parameter file not loaded with load_params || COMPLETE.format.ids==FALSE")
  #   numWorkers <- params_list$numWorkers
  # }

  if(use.feather && gzip.output){
    warning("all2all_BLAST: Both use.feather && gzip.output cannot be TRUE. Setting use.feather=FALSE")
    use.feather=F
  }
  
  if(use.feather){
    f_ext <- "arrows"
  }else{
    f_ext <- "gz"  
  }
  
  #tictoc::tic(msg = "all2all_BLAST :")
set.seed(seed)
  if(verbose){
    cat(paste("All2All BLAST Started...","\n",sep = ""))
    print(paste(first_list,collapse = ","))
    print(paste(second_list,collapse = ","))
  }

  dir.create(path = output_dir,recursive = T,showWarnings = F)

  list_combinations <- unique(tidyr::crossing(first_list,second_list))
  #return_data <- furrr::future_map2(.x=list_combinations$first_list, .y=list_combinations$second_list, .f=function(first_set,second_set){
  return_data <- parallel::mclapply(seq_along(1:nrow(list_combinations)), function(idx){
    first_set <- list_combinations[idx,1]
    second_set <- list_combinations[idx,2]
    fw_dir <- NULL
    bk_dir <- NULL
    tryCatch(fw_dir <- unlist(one2one_BLAST(first_list = first_set,second_list = second_set, file_ext=file_ext,run_name = paste(first_set,second_set,"fw.all2all",f_ext,sep="."),blast_DB_dir=blast_DB_dir,blast_program = blast_program,output_dir=output_dir, blast_options=blast_options, input_prefix_path=input_prefix_path, params_list=params_list,COMPLETE.format.ids=COMPLETE.format.ids, gzip.output=gzip.output, clean_extract=clean_extract, n_threads=1, verbose=verbose, seed=seed,return_data=return_data, use.feather=use.feather)),
             error=function(cond){ if(verbose){message(cond)}}) #return_f_callback=return_f_callback,...=...
    tryCatch(bk_dir <- unlist(one2one_BLAST(first_list = second_set,second_list = first_set,file_ext=file_ext, run_name = paste(first_set,second_set,"bk.all2all",f_ext,sep="."),blast_DB_dir=blast_DB_dir,blast_program = blast_program,output_dir=output_dir, blast_options=blast_options, input_prefix_path=input_prefix_path, params_list=params_list,COMPLETE.format.ids=COMPLETE.format.ids, gzip.output=gzip.output, clean_extract=clean_extract, n_threads=1, verbose=verbose, seed=seed,return_data=return_data, use.feather=use.feather)),
             error=function(cond){ if(verbose){message(cond)}}) #,return_f_callback=return_f_callback,...=...
    if(return_data){
      return(list(fw_dir=fw_dir,bk_dir=bk_dir))
    }else{
      return(NULL)
    }

    #return(NULL)
    #  }, mc.cores = floor(sqrt(numWorkers)) )
    #}, mc.cores = floor(sqrt(numWorkers)) )

  }, mc.cores = n_threads, mc.preschedule=T ) #, .options = furrr::furrr_options(seed = seed, scheduling=F)) #

  #cat(print_toc(tictoc::toc(quiet = T, log = T)))
  return(return_data)
}

