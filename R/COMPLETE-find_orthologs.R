#' Load BLAST Hits into a data.frame
#'
#' Give the path to a BLAST Hits file to load it into a data.frame(BLAST HITs Table). The column names can be provided as col.names. Rows with NAs are automatically removed
#'
#' Optional : You can provide the file/table with indexed Transcsript IDs. The format must be "file"[tab]"long_id"[tab]"index" (without a header). It can be generated with the  index_FASTA_IDs() (check ?index_FASTA_IDs or index_fastaIDs() in system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"))
#'
#' @examples
#'     LoadBLASTHits(infile, transcript_ID_metadata=NULL, col.names=NULL)
#'
#' @note Column indices will be used when processing the data. Column names are only for user reference. First column must be the query sequence ID and the second column must be the subject sequence ID
#'
#' @param infile BLAST hits filename
#' @param transcript_ID_metadata Tab-delimited File with the filenames, indexed transcript IDs and the long transcript IDs.
#' @param col.names Name of the columns of the BLAST File
#' @param sep Delimiter of the BLAST File columns. Default - '\t'
#' @param header Does the file have a header? . Default - FALSE
#' @return Data Frame with BLAST Results
#' @export
LoadBLASTHits <- function(infile, transcript_ID_metadata=NULL, col.names=NULL, sep="\t", header = F){

  if(!is.null(transcript_ID_metadata) && is.character(transcript_ID_metadata)){
    transcript_ID_metadata <- read.table(file = transcript_ID_metadata,header = F,sep="\t",quote = "")
  }
  if(file.exists(infile) && file.info(infile)$size > 0){
    blast_results <- read.table(file = infile,header = header,sep=sep,quote = "", blank.lines.skip = T)
  }else{
    stop(paste("File",infile,"does not exist"))
  }
  if(!is.null(col.names)){
    if(ncol(blast_results) == length(col.names)){
      colnames(blast_results) <- col.names
    }else{
      stop(paste("Number of columns in BLAST file and length of given column names are not equal!",ncol(blast_results),"!=",length(col.names)))
    }
  }

  if(!is.null(transcript_ID_metadata)){
    transcript_ID_metadata <- as.data.frame(transcript_ID_metadata)
    colnames(transcript_ID_metadata) <- c("file","id","index")

    query_long_IDs <- transcript_ID_metadata[match(blast_results[,c(1)],transcript_ID_metadata$index),c("id")]
    subject_long_IDs <- transcript_ID_metadata[match(blast_results[,c(2)],transcript_ID_metadata$index),c("id")]

    blast_results[,c(1)] <- query_long_IDs
    blast_results[,c(2)] <- subject_long_IDs
  }

  blast_results <- blast_results[apply(blast_results, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),] #REMOVING NAs

  return(blast_results)
}

#' Index Table Column
#'
#' Convert longer FASTA IDs into short indices. BLAST table IDs need to be short for downstream processing (eg, graph package does not accept long names). The IDs are stored as an extra column and returned along with the BLAST table.
#'
#' @note WARNING : Do not index it more than once (although you may do index it as much as you like)
#'
#' @examples
#'  #Indexing columns 1 & 2
#' in1_data <- index_BLAST_table(in1_data,2,offset = 0)
#' in1_data <- index_BLAST_table(in1_data,1,offset = length(levels(factor( in1_data[,1] ))))
#' #TO INDEX TWO TABLES
#' in1_data <- index_BLAST_table(in1_data,2,offset = 0)
#' in1_data <- index_BLAST_table(in1_data,1,offset = length(levels(factor( in2_data[,1] ))) + length(levels(factor( in1_data[,2] ))))
#' in2_data <- index_BLAST_table(in2_data,2,offset = length(levels(factor( in1_data[,1] ))) + length(levels(factor( in2_data[,2] ))))
#' in2_data <- index_BLAST_table(in2_data,1,offset = length(levels(factor( in2_data[,1] ))) + length(levels(factor( in1_data[,2] ))))
#'
#' @param in_table Input table
#' @param index_col Column index of BLAST table to index (shorten IDs)
#' @param offset Offset value to add to the indices
#' @return BLAST table with indexed IDs from index_col (original index_col is attached to the table)
#' @export
index_table_column <- function(in_table, index_col, offset=0){
  old_name = colnames(in_table)[index_col]
  in_table <- dplyr::mutate(in_table, in_table[,index_col])

  seq_ids <- levels(factor(in_table[,index_col]))
  indexed_ids <- data.frame(index=paste("i",seq(1:length(seq_ids))+offset,sep=""))
  rownames(indexed_ids) <- seq_ids
  #in_table[,index_col] <- paste("i",seq(1:length(seq_ids))+offset,sep="")
  in_table[,index_col] <- sapply(in_table[,index_col], function(x){
    return(indexed_ids[x,"index"])
  })

  colnames(in_table)[index_col] <- paste("indexed",old_name,sep="_")
  colnames(in_table)[ncol(in_table)] <- old_name
  #print(blast_table)
  return(in_table)
}

#' Index BLAST Table(s)
#'
#' Convert longer FASTA IDs into short indices. BLAST table IDs need to be short for downstream processing (eg, graph package does not accept long names). The IDs are stored as an extra column and returned along with the BLAST table.
#'
#' @note WARNING : Do not index it more than once (although you may do index it as much as you like)
#'
#' @examples
#'  in_data <- index_BLAST_tables(blast_tables = list(data1,data2),query_index_cols = 1, subject_index_cols = 2)
#'
#' @param blast_tables List of BLAST tables
#' @param query_index_cols List of Column indices of query sequence IDs from each table
#' @param subject_index_cols List of Column indices of subject sequence IDs from each table
#' @param offset Offset value to add to the indices
#' @return BLAST table with indexed IDs from index_col (original index_col is attached to the table)
#' @export
index_BLAST_tables <- function(blast_tables, query_index_cols, subject_index_cols, offset=0){

  if(length(query_index_cols) != length(blast_tables) && length(query_index_cols) == 1){
    query_index_cols <- lapply(1:length(blast_tables),function(x){
      return(query_index_cols)
    })
  }else{
    stop("length(query_index_cols) != length(blast_tables)")
  }

  if(length(subject_index_cols) != length(blast_tables) && length(subject_index_cols) == 1){
    subject_index_cols <- lapply(1:length(blast_tables),function(x){
      return(subject_index_cols)
    })
  }else{
    stop("length(subject_index_cols) != length(blast_tables)")
  }

  old_query_cols = lapply(seq_along(blast_tables), function(x){
    colnames(blast_tables[[x]])[query_index_cols[[x]]]
  })  #list(colnames(blast_table1)[query_index_cols[[1]]],colnames(blast_table2)[query_index_cols[[2]]])
  old_subject_cols = lapply(seq_along(blast_tables), function(x){
    colnames(blast_tables[[x]])[subject_index_cols[[x]]]
  }) #list(colnames(blast_table1)[subject_index_cols[[1]]],colnames(blast_table2)[subject_index_cols[[2]]])

  all_ids <- unique( unlist( lapply(seq_along(blast_tables), function(x){
    return(levels(c(factor(blast_tables[[x]][,query_index_cols[[x]]]),factor(blast_tables[[x]][,subject_index_cols[[x]]]))))
  }),recursive = T,use.names = F )  ) #levels(unique(c(factor(blast_table1[,query_index_cols[[1]]]),factor(blast_table2[,query_index_cols[[2]]]),factor(blast_table1[,subject_index_cols[[1]]]),factor(blast_table2[,subject_index_cols[[2]]]))))
  #print(head(all_ids)) #DEBUG
  indexed_ids <- data.frame(id=all_ids,index=paste("i",seq(1:length(all_ids))+offset,sep=""))
  #print(head(indexed_ids)) #DEBUG
  blast_tables <- parallel::mclapply(seq_along(blast_tables), function(x){
    blast_tables[[x]] <- blast_tables[[x]] %>% mutate(blast_tables[[x]][,query_index_cols[[x]]])
    blast_tables[[x]][,query_index_cols[[x]]] <- indexed_ids[match(blast_tables[[x]][,query_index_cols[[x]]],indexed_ids$id), c("index")]
    colnames(blast_tables[[x]])[query_index_cols[[x]]] <- paste("indexed",old_query_cols[[x]],sep="_")
    colnames(blast_tables[[x]])[ncol(blast_tables[[x]])] <- old_query_cols[[x]]

    blast_tables[[x]] <- blast_tables[[x]] %>% mutate(blast_tables[[x]][,subject_index_cols[[x]]])
    blast_tables[[x]][,subject_index_cols[[x]]] <- indexed_ids[match(blast_tables[[x]][,subject_index_cols[[x]]],indexed_ids$id), c("index")]
    colnames(blast_tables[[x]])[subject_index_cols[[x]]] <- paste("indexed",old_subject_cols[[x]],sep="_")
    colnames(blast_tables[[x]])[ncol(blast_tables[[x]])] <- old_subject_cols[[x]]
    return(blast_tables[[x]])
  },mc.silent = T,mc.cores = COMPLETE$numWorkers)

  return(blast_tables)

}

#' De-Index BLAST Table Column
#'
#' Convert indexed FASTA IDs into original IDs. The indexed ID column is removed and replaced with original ID column
#'
#' @note Indexed Column name must be prefixed with "indexed_" followed by original column name, eg "indexed_query"
#'
#' @examples
#' in1_data <- index_BLAST_table(in1_data,"sseqid",offset = 0)
#' in1_data <- index_BLAST_table(in1_data,"qseqid",offset = length(levels(factor( in1_data[,"qseqid"] ))))
#' #TO INDEX TWO TABLES
#' in1_data <- index_BLAST_table(in1_data,"sseqid",offset = 0)
#' in1_data <- index_BLAST_table(in1_data,"qseqid",offset = length(levels(factor( in2_data[,"qseqid"] ))) + length(levels(factor( in1_data[,"sseqid"] ))))
#' in2_data <- index_BLAST_table(in2_data,"sseqid",offset = length(levels(factor( in1_data[,"qseqid"] ))) + length(levels(factor( in2_data[,"sseqid"] ))))
#' in2_data <- index_BLAST_table(in2_data,"qseqid",offset = length(levels(factor( in2_data[,"qseqid"] ))) + length(levels(factor( in1_data[,"sseqid"] ))))
#'
#' @param blast_table BLAST table
#' @param index_col Column index of BLAST table to de-index
#' @return BLAST table with original IDs from index_col (index_col is removed from the table and replaced with original column)
#' @export
deindex_BLAST_table <- function(blast_table, index_col){
  old_name = colnames(blast_table)[index_col]
  new_name = stringi::stri_split_fixed(old_name, pattern = "_", n=2,simplify = T)[,2]

  blast_table[,old_name] <- blast_table[,new_name]
  blast_table[,new_name] <- NULL
  colnames(blast_table)[which(!is.na(match(colnames(blast_table), old_name)))] <- new_name
  return(blast_table)
}

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
#' @param blast_input Filename with BLAST hits or a BLAST table
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN)
#' @param col.indices A Named List with indices of columns Query sequence ID (qseqid), Subject sequence ID (sseqid), E Value (evalue), Query start (qstart), Query end (qend), Subject start (sstart), Subject end (send), Bitscore (bitscore), Query HSP coverage(qcovhsp), Query length (qlen), Subject length(slen), Frames (frames), Percentage Identity (pident), Number of Gaps (gaps), Alignment Length (length), Subject strand (sstrand). Default values are col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17)
#' @param params_list Output of load_params() (Optional - Only give it when using Files generated by R-COMPLETE and if COMPLETE.format.ids=TRUE)
#' @return A GRanges object of BLAST hits
#' @export
GRObject_from_BLAST <- function(blast_input, COMPLETE.format.ids=F, col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list=NULL){

  if(!is.null(blast_input) && is.character(blast_input)){
    blast_input <- LoadBLASTHits(blast_input)
  }

  # if(!is.null(col.names)){
  #   if(ncol(blast_input) != length(col.names)){
  #     stop("ncol(blast_input) != length(col.names)")
  #   }
  #   colnames(blast_input) <- col.names
  # }

  req_columns <- c("qseqid","sseqid","pident","length","qstart","qend","sstart","send","evalue","bitscore","gaps","frames","qcovhsp","sstrand","qlen","slen")

  if (any(is.na(match(req_columns, names(col.indices))))) {
    stop(paste("Missing columns, Require indices of :",paste(req_columns,collapse = ",")))
  }

  if(COMPLETE.format.ids && !is.null(params_list)){
    blast_input <- blast_input %>% mutate(subject_gene=unlist(purrr::map(blast_input[,col.indices[["sseqid"]]],function(x){
      #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      genes <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(genes[,COMPLETE$FORMAT_ID_INDEX$GENE])
    })) )

    blast_input <- blast_input %>% mutate(query_gene=unlist(purrr::map(blast_input[,col.indices[["qseqid"]]],function(x){
      #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      genes <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(genes[,COMPLETE$FORMAT_ID_INDEX$GENE])
    })) )

    blast_input <- blast_input[which(blast_input[,col.indices[["evalue"]]] < params_list$E_VALUE_THRESH),]

  }#else{
  #  stop("Parameter file not loaded with load_params(), params_list is NULL & COMPLETE.format.ids==FALSE")
  #}

  ##CHANGE strands
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

  #tmp_gr <- GenomicRanges::GRanges(S4Vectors::Rle(blast_input[,col.indices["sseqid"]]), ranges =  IRanges::IRanges(blast_input[,col.indices["sstart"]], end = blast_input[,col.indices["send"]],strand= blast_input[,col.indices["sstrand"]])) #, names = orths$sseqid))

  tmp_df <- data.frame(seq_names=blast_input[,col.indices[["sseqid"]]])

  #tmp_df <- mutate(tmp_df, seqnames=S4Vectors::Rle(blast_input[,col.indices["sseqid"]]))
  #tmp_df <- mutate(tmp_df, ranges=IRanges::IRanges(start = blast_input[,col.indices["sstart"]], end = blast_input[,col.indices["send"]],strand= blast_input[,col.indices["sstrand"]]))
  tmp_df <- mutate(tmp_df, sstart = blast_input[,col.indices[["sstart"]]])
  tmp_df <- mutate(tmp_df, send = blast_input[,col.indices[["send"]]])
  tmp_df <- mutate(tmp_df, sstrand= blast_input[,col.indices[["sstrand"]]])
  tmp_df <- mutate(tmp_df, Hsp_num=c(1:nrow(blast_input))) #seq(1,nrow(tmp_blast),1)
  tmp_df <- mutate(tmp_df, Hsp_bit.score = blast_input[,col.indices[["bitscore"]]])
  tmp_df <- mutate(tmp_df, Hsp_score = blast_input[,col.indices[["qcovhsp"]]])
  tmp_df <- mutate(tmp_df, Hsp_evalue = blast_input[,col.indices[["evalue"]]])
  tmp_df <- mutate(tmp_df, subject_HSP_from = blast_input[,col.indices[["sstart"]]])
  tmp_df <- mutate(tmp_df, subject_HSP_to = blast_input[,col.indices[["send"]]])
  tmp_df <- mutate(tmp_df, query_id = blast_input[,col.indices[["qseqid"]]])
  tmp_df <- mutate(tmp_df, query_len = blast_input[,col.indices[["qlen"]]])
  tmp_df <- mutate(tmp_df, subject_len = blast_input[,col.indices[["slen"]]])
  tmp_df <- mutate(tmp_df, query_HSP_from = blast_input[,col.indices[["qstart"]]])
  tmp_df <- mutate(tmp_df, query_HSP_to = blast_input[,col.indices[["qend"]]])
  tmp_df <- mutate(tmp_df, Hsp_query.frame =  unlist(purrr::map(blast_input[,col.indices[["frames"]]],function(x){
    frames <- as.integer(unlist(stringi::stri_split_fixed(x,pattern = "/")))
    return(frames[1])
  })) )
  tmp_df <- mutate(tmp_df, Hsp_hit.frame =  unlist(purrr::map(blast_input[,col.indices[["frames"]]],function(x){
    frames <- as.integer(unlist(stringi::stri_split_fixed(x,pattern = "/")))
    return(frames[2])
  })) )
  if(COMPLETE.format.ids && !is.null(params_list)){
    tmp_df <- mutate(tmp_df, query_org = unlist(purrr::map(blast_input[,col.indices[["qseqid"]]],function(x){
      #org <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      org <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(org[,COMPLETE$FORMAT_ID_INDEX$ORG])
    })) )
    tmp_df <- mutate(tmp_df, subject_org = unlist(purrr::map(blast_input[,col.indices[["sseqid"]]],function(x){
      #org <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      org <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
      return(org[,COMPLETE$FORMAT_ID_INDEX$ORG])
    })) )
  }else{
    warning("Parameter file not loaded with load_params & COMPLETE.format.ids==FALSE")
  }

  tmp_df <- mutate(tmp_df, pidentity = blast_input[,col.indices[["pident"]]])
  #gr$Hsp_positive <-
  tmp_df <- mutate(tmp_df, Hsp_gaps = blast_input[,col.indices[["gaps"]]])
  tmp_df <- mutate(tmp_df, Hsp_align.len =  blast_input[,col.indices[["length"]]])
  #print(head(blast_input)) #DEBUG
  if(COMPLETE.format.ids  && !is.null(params_list)){
    tmp_df <- mutate(tmp_df, subject_gene =  blast_input[,c("subject_gene")])
    tmp_df <- mutate(tmp_df, query_gene =  blast_input[,c("query_gene")])
  }

  tmp_gr <- GenomicRanges::makeGRangesFromDataFrame(df = tmp_df,seqnames.field = "seq_names",start.field = "sstart",end.field = "send",strand.field = "sstrand" ,keep.extra.columns = T,ignore.strand = F)

  return(tmp_gr)
}

#' Select highest scoring interval of non-overlapping HSPs from Bi-Directional BLAST Hits
#'
#' Select the highest scoring pairs (HSPs) which give the maximum coverage over the BLAST alignments of each transcript (without overlaps/minimal overlaps). These HSPs will then be used to find Transcript level orthologs across gene orthologs across organisms. This function only accepts bi-directional (Query <-> Subject) BLAST Hits formatted with GRObject_from_BLAST()
#'
#' @examples
#'
#' blast_GO <- GRObject_from_BLAST(blast_input = in_file, COMPLETE.format.ids = T, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
#' wis_GO <- run_WISARD(blast_hits = blast_GO,score_col = "Hsp_score",COMPLETE.format.ids = T,params_list=NULL) #score_col=16
#' melt_wisard_list(wis_GO)
#'
#' @param blast_hits GRanges Object of BLAST Hits (Query -> Subject) (from GRObject_from_BLAST())
#' @param score_col Column in the GRanges Object used for scoring intervals in WISARD, Default would be "Hsp_score" (qcovhsp (Query Coverage HSP))
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, Default - FALSE otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN)
#' @param params_list Output of load_params() (Optional - Only give it when using Files generated by R-COMPLETE and if COMPLETE.format.ids=TRUE)
#' @param verbose Print Output?
#' @return A GRanges object of BLAST hits
#' @export
run_WISARD <- function(blast_hits, score_col, COMPLETE.format.ids=F,params_list=NULL, verbose=F){

  if(!all(grepl(pattern="GRanges",x = class(blast_hits),ignore.case = T))){
    stop("This function only accepts formatted GRanges Object from GRObject_from_BLAST()")
  }

  if(is.null(params_list)){
    tryCatch(numWorkers <- parallel::detectCores(all.tests = T, logical = T), error=function(){numWorkers <- 2})
  }else{
    #warning("Parameter file not loaded with load_params || COMPLETE.format.ids==FALSE")
    numWorkers <- params_list$numWorkers
  }
  gr <- blast_hits
  query_vector <- unique(gr$query_id)
  #print(S4Vectors::runValue(GenomicRanges::seqnames(gr)))
  subject_vector <- unique(S4Vectors::runValue(GenomicRanges::seqnames(gr)))
  child_results <- BiocGenerics::unlist(parallel::mcmapply(FUN=function(query_id,subject_ids, gr,score_col){
    result_list <- c()
    result_list <- parallel::mclapply(subject_ids,function(sub_id){
      #if(verbose){
      #  result <- wisard::get_wis(gr[unique(intersect(which(gr$query_id == query_id),which(S4Vectors::runValue(GenomicRanges::seqnames(gr)) == sub_id))),],max_score = score_col,overlap = 0)
      #}else{
      #  invisible(result <- wisard::get_wis(gr[unique(intersect(which(gr$query_id == query_id),which(S4Vectors::runValue(GenomicRanges::seqnames(gr)) == sub_id))),],max_score = score_col,overlap = 0))
      #}
      result <- suppressMessages(invisible(wisard::get_wis(gr[unique(intersect(which(gr$query_id == query_id),which(S4Vectors::runValue(GenomicRanges::seqnames(gr)) == sub_id))),],max_score = score_col,overlap = 0)))
      if(result$max_score > 0 && length(result$alignments) > 0){
        #result_list <- c(result_list, result)
        return(result)
        #print(result_list)
      }
      #print(result_list)

    },mc.cores = numWorkers)

    return(result_list)
  }, query_vector, MoreArgs=list(subject_ids=subject_vector,gr=gr,score_col=score_col) ,mc.cores = numWorkers, SIMPLIFY = FALSE,USE.NAMES = F),recursive = F, use.names = T)

  all_results <- list()
  if(COMPLETE.format.ids && !is.null(params_list)){
    for(child in child_results) {
      if(verbose){
        print(child)
      }
      g_name <- NULL
      g_name <- unique(unlist(stringi::stri_split_fixed(child$alignments$subject_gene,pattern=params_list$SEQUENCE_ID_DELIM,n = 1,tokens_only = T)))
      if(verbose){
        print(g_name)
      }
      #all_results<- c(all_results, child)
      if(!is.null(g_name)){
        if(is.null(all_results) || length(all_results)==0){
          all_results[[g_name]] <- child
        }else{
          all_results[[g_name]] <- merge_wisard_lists(all_results[[g_name]], child)
        }
        # names(all_results[[g_name]]) <- g_name
      }
    }
  }else{
    all_results <- child_results #purrr::reduce(child_results, merge_wisard_table)
  }

  return(all_results)

}

#' Internal Function - To Merge WISARD OUTPUTS Lists
#'
#' @param x WIZARD output List x
#' @param y WIZARD output List y
#' @return A Merged List of WISARD Outputs
#' @export
merge_wisard_lists <- function(x, y){
  new_list <- list()
  new_list$alignments <- c(x$alignments, y$alignments)
  #new_list$alignments$max_score <- c(rep(x$max_score, length(x$alignments)),rep(y$max_score, length(y$alignments)))
  new_list$max_score <- sum(x$max_score, y$max_score)
  #print(new_list)
  return(unlist(new_list))
}

#' Internal Function - Convert/Melt WISARD List to Data.Frame
#'
#' @param x WIZARD output List x
#' @return Data Frame
#' @export
melt_wisard_list <- function(x){
  wis_df <- dplyr::bind_rows(lapply(x, function(wis_obj){
    tmp_df <- data.frame(wis_obj[["alignments"]])
    tmp_df <- mutate(tmp_df, max_score=wis_obj[["max_score"]])
    return(tmp_df)
  }))
  names(wis_df)[grep(pattern="seqnames",names(wis_df))] <- "subject_id"
  if(any(grepl(x = names(wis_df), pattern = "ID", fixed = T))){
    wis_df[,"ID"] <- NULL
  }
  return(wis_df)
}

#' Internal Function - To calculate gene conservation score (GSC)
#'
#'This function is only designed to work with long ID format followed by the R-COMPLETE pipeline and assumes that the BLAST table is filtered by WISARD. The BLAST File/Table must be of the format 6. BLAST formats can be converted between each other using. By default, Percentage Sequence Identity is used as a score metric, the function expects BLAST Hits from tblastx (because the sequence is translated before it is BLASTed and the Percentage Identity is from the protein sequence). Other columns can be used as a metric with the score_col parameter
#'
#' convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'
#' @param blast_table BLAST table
#' @param gene BLAST table
#' @param score_col Column name or number to be used as a score metric. Percentage Sequence Identity is used by default (column 3)
#' @param available_orgs Vector of available organisms (in the params_list$FASTA_OUT_PATH)
#' @param num_transcripts Number of Transcripts available for the gene. This can be obtained from the blastdb of a gene
#' @param REF_ORGS Vector of reference organisms
#' @return Named Vector with an double value between minimum and maximum values in the score column (0-100 in case of percentage identity) which is the gene conservation score for the gene
calculate_gene_conservation <- function(blast_table, gene, score_col=3,available_orgs, num_transcripts, REF_ORGS){

  if(is.character(blast_table)){
    blast_table <- LoadBLASTHits(blast_table)
  }

  identities_out_path <- paste(params_list$OUT_PATH,"/idents",sep="")
  thres_out <-  paste(params_list$OUT_PATH,"/gene_thresholds.txt",sep="")
  #aorgs_path <- available_orgs
  #sorgs_path <- REF_ORGS
  oinfo_out_path <-  paste(params_list$OUT_PATH,"/orgs",sep="")

  dir.create(identities_out_path,showWarnings = F,recursive = T)
  dir.create(oinfo_out_path)

  blast_table <- blast_table %>% mutate(from_org=unlist(lapply(stringi::stri_split_fixed(blast_table[,1],delimiter,n=4,tokens_only = T),function(x){return(x[COMPLETE$FORMAT_ID_INDEX$ORG])})))
  blast_table <- blast_table %>% mutate(to_org=unlist(lapply(stringi::stri_split_fixed(blast_table[,2],delimiter,n=4,tokens_only = T),function(x){return(x[COMPLETE$FORMAT_ID_INDEX$ORG])})))

  blast_table[,1] <- factor(blast_table[,1])
  blast_table[,2] <- factor(blast_table[,2])
  blast_table$from_org <- factor(blast_table$from_org)
  blast_table$to_org <- factor(blast_table$to_org)

  #REF_ORGS <- factor(scan(params_list$REF_ORGS_FILE, character()))

  min_seq_identity <- data.frame(matrix(numeric(), nrow = length(REF_ORGS), ncol = length(available_orgs)))
  max_seq_identity <- data.frame(matrix(numeric(), nrow = length(REF_ORGS), ncol = length(available_orgs)))
  rownames(min_seq_identity) <- REF_ORGS
  rownames(max_seq_identity) <- REF_ORGS
  colnames(min_seq_identity) <- available_orgs
  colnames(max_seq_identity) <- available_orgs

  min_seq_identity[is.na(min_seq_identity)] <- 101
  max_seq_identity[is.na(max_seq_identity)] <- -1

  invisible(lapply(REF_ORGS, function(s_org){
    lapply(available_orgs,function(a_org){
      if(s_org!=a_org){
        min_seq_identity[s_org,a_org] <- min(blast_table[which(blast_table$"from_org" == s_org & blast_table$"to_org" == a_org),score_col] )
        max_seq_identity[s_org,a_org] <- max(blast_table[which(blast_table$"from_org" == s_org & blast_table$"to_org" == a_org),score_col] )
      }
    })
  }))

  min_seq_identity <- apply(min_seq_identity, MARGIN = 2 , function(x) replace(x, is.infinite(x), 0))
  max_seq_identity <- apply(max_seq_identity, MARGIN = 2 , function(x) replace(x, is.infinite(x), 0))

  if(identical(class(max_seq_identity),"numeric") && identical(class(min_seq_identity),"numeric")){  #If only one organism is selected then we have to adjust datatypes
    max_seq_identity <- t(as.data.frame(max_seq_identity))
    rownames(max_seq_identity) <- REF_ORGS
    min_seq_identity <- t(as.data.frame(min_seq_identity))
    rownames(min_seq_identity) <- REF_ORGS
  }

  if(nrow(max_seq_identity)>1){
    png(filename =paste(file.path(params_list$PLOT_PATH,gene),"max_ident.png",sep = "."),width = 15, height = 15 , units = "in", res = 100)
    heatmap.2(max_seq_identity, trace = "row", tracecol = "black", margins = c(12,12), scale = "none",cexRow=1.2,cexCol=1.2)
    dev.off()
  }
  if(nrow(min_seq_identity)>1){
    png(filename =paste(file.path(params_list$PLOT_PATH,gene),"min_ident.png",sep = "."),width = 15, height = 15 , units = "in", res = 100)
    heatmap.2(min_seq_identity, trace = "row", tracecol = "black", margins = c(12,12),scale = "none",cexRow=1.2,cexCol=1.2)
    dev.off()
  }

  discardable_orgs <- available_orgs[which(is.na(match(available_orgs, levels(blast_table$to_org))))] # We have no sequences of the gene for these orgs
  saveable_orgs <- available_orgs[which(!is.na(match(available_orgs, levels(blast_table$to_org))))]

  nref_orgs <- c()
  nref_orgs <- BiocGenerics::unlist(lapply(REF_ORGS,function(s_org){
    if(max_seq_identity[which(rownames(max_seq_identity)==s_org),which(colnames(max_seq_identity)==s_org)]==0){
      warning(paste(gene,": is not present in : ", s_org))
      #nref_orgs <- c(nref_orgs, s_org)
      return(s_org)
    }
  }))

  if(any(is.na(match(REF_ORGS, levels(blast_table$from_org))))){ ## Are all reference species blasted
    warning(paste(gene,": Either gene missing/Please rerun BLAST for ", gene," (",blast_results_path,")"))
    #quit(status = 2)
  }

  if(any(is.na(match(available_orgs, levels(blast_table$to_org))))){ ## Do blast results contain all the available species?
    warning(paste(gene, ": Species not blasted against: ")) # because gene not available
    print(discardable_orgs)
  }

  if(any(is.na(match(REF_ORGS, levels(blast_table$to_org))))){
    warning(paste("WARNING: ",gene, " not present in all reference species:"))
    print(nref_orgs)
  } ##Between reference species

  max_seq_identity[,discardable_orgs] <- NA
  max_seq_identity[nref_orgs,] <- NA
  min_seq_identity[,discardable_orgs] <- NA
  min_seq_identity[nref_orgs,] <- NA
  max_seq_identity <- as.data.frame(max_seq_identity)
  min_seq_identity <- as.data.frame(min_seq_identity)

  semi_ortho_orgs <- unique(colnames(max_seq_identity)[which(max_seq_identity == 0, arr.ind = T)[,"col"]])

  semi_orths <- max_seq_identity[, c(semi_ortho_orgs)]

  for (a_org in names(semi_orths)) {
    if(all(is.na(match(semi_orths[, (names(semi_orths) %in% a_org)],0)))){
      discardable_orgs <- c(discardable_orgs, a_org)
      semi_orths <- semi_orths[, !(names(semi_orths) %in% a_org)]
      semi_ortho_orgs <- semi_ortho_orgs[!semi_ortho_orgs %in% a_org]
    }
  }
  warning("These organisms are semi-orthologous (ie, they are not orthologous to all reference species):")
  print(semi_ortho_orgs)
  print(semi_orths)
  write.table(semi_orths, file = file.path(identities_out_path,paste(gene,".semiorgs",sep = "")),quote = F, row.names = T, col.names = T)

  ##Correct isoform representation, for all reference species
  ##As in, we do not care about the similarity of the reference sequence and it's isoforms. So we make it 0
  ##Having these values also gives rise to ridiculous scores for poorly conserved genes (camk2g scores - before correction = 100, after correction = 48)
  #min_seq_identity[which(min_seq_identity == 100, arr.ind = T)] <- 0
  #min_seq_identity[,(names(min_seq_identity) %in% orgs.selected)] <- 0
  for (s_org in REF_ORGS) {
    min_seq_identity[s_org,s_org] <- 0
  }
  if(length(semi_ortho_orgs) > 0){
    max_seq_identity <- max_seq_identity[, !(names(max_seq_identity) %in% semi_ortho_orgs)]
    min_seq_identity <- min_seq_identity[,!(names(min_seq_identity) %in% semi_ortho_orgs)]
  }
  if(length(discardable_orgs) > 0){
    max_seq_identity <- max_seq_identity[,!(names(max_seq_identity) %in% discardable_orgs)]
    min_seq_identity <- min_seq_identity[,!(names(min_seq_identity) %in% discardable_orgs)]
  }

  min_gene_conservation <- max(min(max_seq_identity,na.rm = T),max(min_seq_identity,na.rm = T))
  #print(min_gene_conservation)
  #print(max_seq_identity)
  low_orthology_orgs <- unique(names(max_seq_identity)[which(apply(max_seq_identity,MARGIN=c(1,2),function(x){ x < min_gene_conservation }), arr.ind = T)[,"col"]])
  if(!identical(low_orthology_orgs, character(0))){
    warning(paste("These organisms did not pass the minimum threshold(",min_gene_conservation,"):"))
    print(all_of(low_orthology_orgs))
  }

  saveable_orgs <- anti_join(as.data.frame(list(saveable_orgs), col.names=c("orgs")),as.data.frame(list(semi_ortho_orgs), col.names=c("orgs")))
  saveable_orgs <- as.vector(t(saveable_orgs))
  if(length(discardable_orgs) > 0){
    fwrite(list(discardable_orgs), file = file.path(oinfo_out_path,paste(gene,".dorgs",sep = "")))
  }
  if(length(saveable_orgs) > 0){
    fwrite(list(saveable_orgs), file = file.path(oinfo_out_path,paste(gene,".sorgs",sep = "")))
  }
  if(length(low_orthology_orgs) > 0){
    fwrite(list(low_orthology_orgs), file = file.path(oinfo_out_path,paste(gene,".lorgs",sep = "")))
  }
  if(length(nref_orgs) > 0){
    fwrite(list(nref_orgs), file = file.path(oinfo_out_path,paste(gene,".nref",sep = "")))
  }
  write.table(min_seq_identity, file = file.path(identities_out_path,paste(gene,".min",sep = "")), quote = F, row.names = T, col.names = T)
  write.table(max_seq_identity, file = file.path(identities_out_path,paste(gene,".max",sep = "")), quote = F, row.names = T, col.names = T)

  #fileConn<-file(thres_out,open = "at")
  #writeLines(paste(gene, min_gene_conservation,length(nref_orgs)==0,num_transcripts,sep = ","), fileConn)
  #close(fileConn)

  names(min_gene_conservation) <- gene
  return(min_gene_conservation)
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
#' @param conversion_prg Path to blast_formatter if not found in $PATH
#' @export
convert_BLAST_format <- function(infile, outfile,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"),conversion_prg = Sys.which("blast_formatter")){
  if(stringi::stri_isempty(conversion_prg)){
    stop("blast_formatter not found in $PATH..cannot continue!")
  }

  processx::run( command = COMPLETE$SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"convert_BLAST_format",conversion_prg,infile,outfile,outformat,paste(cols,collapse = " ") ) ,spinner = T,stdout = "",stderr = "")

}

#' Perform Nucleotide BLAST
#'
#' BLAST between two FASTA files. BLAST output is saved temporarily in BLAST format 11, internally converted to format 6 with convert_BLAST_format() and BLAST Hits are returned as a GRanges Object. If blast_out is provided then the intermediate BLAST Format 11 output can be found in the same directory with the same name and *.blast11 extension. If blast_DB_dir is provided then the Subject FASTA is copied into this directory and then BLASTed.
#'
#' @note ASSUMES Nucleotide FASTA sequences. Not checking/Not working for Protein/Peptide sequences and FASTQ files. Columns of GRanges BLAST 6 output are c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive")
#'
#' @examples
#'     run_BLAST(query_path = "query.fasta",subject_path = "subject.fasta",blast_DB_dir = "files/blastdb", blast_program="tblastx", run_name = "blast_positive",blast_options = "-strand plus")
#'
#' @param query_path Path to Query FASTA
#' @param subject_path Path to Subject FASTA
#' @param blast_DB_dir Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed. Default is tempdir(). Can also be NULL
#' @param blast_out Path to BLAST output file, Default BLAST FORMAT is 11. It is converted internally to BLAST Format 6 and returned as a GRanges Object
#' @param blast_program Give path to the BLAST program. eg, Sys.which("tblastx") if tblastx is in SHELL $PATH.
#' @param run_name Name of the BLAST run. Only for logging (Optional)
#' @param blast_options Extra Options to be passed to the BLAST program
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE (Default) otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (Optional)
#' @param params_list Output from load_params() (Optional)
#' @param keep.output.files TRUE(Default)/FALSE - Keep Output and BLAST DB Files? (Optional)
#' @param verbose Print DEBUG Messages?
#' @return BLAST Hits as GRanges Object
#' @export
run_BLAST <- function(query_path, subject_path,blast_DB_dir = tempdir(), blast_out, blast_program=COMPLETE$BLAST_BIN, run_name="BLAST",blast_options="", COMPLETE.format.ids=F,params_list=NULL, keep.output.files=T, verbose=T){

  tryCatch({
    if (stringi::stri_isempty(COMPLETE$parallel)) {
      stop("Problem with GNU parallel installation. Reload R-COMPLETE")
    }

    if(stringi::stri_isempty(blast_program)){
      stop("BLAST+ not found in $PATH. Provide blast_program")
    }else{
      BLAST_BIN <- dirname(blast_program)
    }

    #subject_path <- tools::file_path_as_absolute(subject_path)
    if (!is.null(blast_DB_dir)) {
      dir.create(blast_DB_dir,showWarnings = F,recursive = T)
      #query_DB <- tools::file_path_as_absolute(paste(blast_DB_dir,"/",basename(query_path),sep=""))
      subject_DB <- paste(blast_DB_dir,"/",basename(subject_path),sep="") #tools::file_path_as_absolute()
      #if(!stringi::stri_cmp_eq(query_path,query_DB)){
      #  file.copy(query_path,query_DB,overwrite = T)
      #}
      #if(!stringi::stri_cmp_eq(subject_path,subject_DB)){
      tryCatch(file.copy(subject_path,subject_DB,overwrite = T), error=function(cond){stop(cond)})
      #}
      #query_path <- query_DB
      if(verbose){
        print(subject_DB)
      }
      subject_path <- subject_DB
    }

    if(is.null(blast_out)){
      if(is.null(params_list)){
        blast_out <- tempfile(pattern="blast_out", tmpdir = params_list$TEMP_PATH)
      }else{
        blast_out <- tempfile(pattern="blast_out", tmpdir = tempdir())
      }
    }
    final_blast_out <- blast_out
    blast_out <- paste(tools::file_path_sans_ext(blast_out),".blast11",sep="")

    if(verbose){
      print(query_path)
      print(subject_path)
      print(blast_out)
    }

    #MAKE BLAST DB of FASTA files
    #Only subject fasta files needs to be a BLAST DB
    #processx::run( command = COMPLETE$SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"make_BLAST_db",query_path, dirname(blast_program)) ,spinner = T,stdout = "",stderr = "")
    processx::run( command = COMPLETE$SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"make_BLAST_db",subject_path, BLAST_BIN ) ,spinner = T,stdout = "",stderr = "")

    processx::run( command = COMPLETE$SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"do_BLAST",COMPLETE$parallel,run_name,query_path,subject_path,blast_out,blast_program,blast_options) ,spinner = T,stdout = "",stderr = "")

    message(paste("Converting",blast_out,"to BLAST Format 6 :", final_blast_out))

    convert_BLAST_format(infile = blast_out,outfile = final_blast_out,conversion_prg = tools::file_path_as_absolute(paste(BLAST_BIN,"/blast_formatter",sep="")) )
    blast_GR <- GRObject_from_BLAST(blast_input = final_blast_out,COMPLETE.format.ids = COMPLETE.format.ids,col.indices = c(qseqid = 1, sseqid = 2, evalue = 11, qstart = 7, qend = 8, sstart = 9, send = 10, bitscore = 12, qcovhsp = 16, qlen = 18, slen = 19, frames = 15, pident = 3, gaps = 14, length = 4, sstrand = 17), params_list = params_list)

    if(verbose){
      print(head(blast_GR))
    }

    if(!keep.output.files){
      unlink(x = c(final_blast_out,blast_out), recursive = T,force = T,expand = T)
      if(!is.null(blast_DB_dir)){
        unlink(x = list.files(blast_DB_dir,pattern = paste(basename(subject_path),".",sep=""),full.names = T ), recursive = T,force = T,expand = T)
      }else{
        unlink(x = list.files(dirname(subject_path),pattern = paste(basename(subject_path),".",sep=""),full.names = T ), recursive = T,force = T,expand = T)
      }

    }

    return(blast_GR)
  },error=function(cond){
    stop(cond)
  })
}

#' Execute one2one BLAST
#'
#' Executes One-to-One BLAST between two lists of organisms/genes/clusters.
#'
#' @examples one2one_BLAST(first_set,second_set,blast_DB_dir=blast_DB_dir,blast_program,output_dir=output_dir, blast_options=blast_options, input_prefix_path=input_prefix_path, params_list=params_list,COMPLETE.format.ids=COMPLETE.format.ids, keep.output.files=keep.output.files, verbose=verbose)
#'
#' @note  ASSUMES Nucleotide sequences. Not checking/Not working for Protein/Peptide sequences.
#'
#' @param first_list Vector of PATHS to FASTA
#' @param second_list Vector of PATHS to FASTA
#' @param run_name Name of the BLAST run and name of the output file
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program. Default is tempdir(). Can also be NULL
#' @param blast_DB_dir Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed
#' @param blast_options Extra Options to be passed to the BLAST program
#' @param output_dir Path to BLAST output
#' @param input_prefix_path If input lists/vectors are filenames, then provide input folder to prefix path
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE (Default) otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (Optional)
#' @param keep.output.files TRUE(Default)/FALSE - Keep Output and BLAST DB Files? (Optional)
#' @param params_list Output of load_params() (Optional)
#' @param verbose Print DEBUG Messages?
#' @export
one2one_BLAST <- function(first_list,second_list, run_name="BLAST.one2one" ,blast_DB_dir=tempdir(),blast_program,output_dir="./", blast_options="", input_prefix_path=NULL, params_list=NULL,COMPLETE.format.ids=F, keep.output.files=T, verbose=F){

  if(!is.null(input_prefix_path)){
    #first_list <- paste(input_prefix_path,"/",first_list,sep="")
    first_list <- list.files(path = input_prefix_path,pattern = first_list, full.names = T,recursive = T,include.dirs = F,ignore.case = T)
    #second_list <- paste(input_prefix_path,"/",second_list,sep="")
    second_list <- list.files(path = input_prefix_path,pattern = second_list, full.names = T,recursive = T,include.dirs = F,ignore.case = T)
  }
  purrr::map2(first_list[order(first_list)], second_list[order(second_list)], function(x,y){
    if (file.exists(x) && file.exists(y) && file.info(x)$size > 0 && file.info(y)$size > 0) {
      # run_name1 <- tools::file_path_sans_ext(BiocGenerics::basename(x)) #tools::file_path_as_absolute()
      # run_name2 <- tools::file_path_sans_ext(BiocGenerics::basename(y)) #tools::file_path_as_absolute()
      # run_name <- paste( run_name1,run_name2,"all2all" ,sep=".")
      if(verbose){
        cat(paste(run_name,"\n",sep = ""))
        print(x)
        print(y)
      }
      out_file <- paste( output_dir, run_name,sep="")

      if(file.exists(out_file) && file.info(out_file)$size > 0){
        if(params_list$CLEAN_EXTRACT){
          unlink(x = out_file,recursive = F,force = T,expand = T)
        }else{
          stop(paste("Output file :",out_file,"exists!"))
        }
      }

      tryCatch(run_BLAST(query_path = x,subject_path = y,blast_DB_dir = blast_DB_dir, blast_program=blast_program, blast_out = out_file, run_name = run_name,blast_options = blast_options,COMPLETE.format.ids = COMPLETE.format.ids,params_list = params_list,keep.output.files = keep.output.files, verbose = verbose), error=function(cond){
        message(cond)
      })
    }
  })

}

#' Execute all2all BLAST
#'
#' Executes All-to-All BLAST between two lists of organisms/genes/clusters. Output BLAST files are stored in the format filename1.filename2.all2all under output_dir. (All-to-All is simply Many-to-Many association)
#'
#' @examples
#'  all2all_BLAST(first_list = list.files(path="fasta",pattern = "\*.fa",all.files = T,full.names = T,include.dirs = F,recursive = F), second_list = list.files(path="fasta",pattern = "\*.fa",all.files = T,full.names = T,include.dirs = F,recursive = F),blast_DB_dir = "files/blastdb",blast_program = "tblastx",output_dir = "files/all2all",blast_options = "-strand plus",input_prefix_path = NULL,params_list=NULL)
#'
#' @note ASSUMES Nucleotide sequences. Not checking/Not working for Protein/Peptide sequences.
#'
#' @param first_list Vector of PATHS to FASTA
#' @param second_list Vector of PATHS to FASTA
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program. Default is tempdir(). Can also be NULL
#' @param blast_DB_dir Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed
#' @param blast_options Extra Options to be passed to the BLAST program
#' @param output_dir Path to BLAST output
#' @param input_prefix_path If input lists/vectors are filenames, then provide input folder to prefix path
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE (Default) otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (Optional)
#' @param keep.output.files TRUE(Default)/FALSE - Keep Output and BLAST DB Files? (Optional)
#' @param params_list Output of load_params() (Optional)
#' @param verbose Print DEBUG Messages?
#' @export
all2all_BLAST <- function(first_list,second_list,blast_DB_dir=tempdir(),blast_program,output_dir="./", blast_options="", input_prefix_path=NULL, params_list=NULL,COMPLETE.format.ids=F, keep.output.files=T, verbose=F){

  if(is.null(params_list)){
    tryCatch(numWorkers <- parallel::detectCores(all.tests = T, logical = T), error=function(){numWorkers <- 2})
  }else{
    #warning("Parameter file not loaded with load_params || COMPLETE.format.ids==FALSE")
    numWorkers <- params_list$numWorkers
  }

  if(verbose){
    cat(paste("All2All BLAST Started...","\n",sep = ""))
    print(paste(first_list,collapse = ","))
    print(paste(second_list,collapse = ","))
  }

  dir.create(path = output_dir,recursive = T,showWarnings = F)

  list_combinations <- unique(tidyr::crossing(first_list,second_list))
  purrr::map2(list_combinations$first_list, list_combinations$second_list, function(first_set,second_set){
    #parallel::mclapply(first_list, function(first_set){
    #  parallel::mclapply(second_list,function(second_set){

    try(one2one_BLAST(first_set,second_set, run_name = paste(first_set,second_set,"all2all",sep="."),blast_DB_dir=blast_DB_dir,blast_program = blast_program,output_dir=output_dir, blast_options=blast_options, input_prefix_path=input_prefix_path, params_list=params_list,COMPLETE.format.ids=COMPLETE.format.ids, keep.output.files=keep.output.files, verbose=verbose))

    #return(NULL)
    #  }, mc.cores = floor(sqrt(numWorkers)) )
    #}, mc.cores = floor(sqrt(numWorkers)) )
  })
}

## @param transcript_region_lengths A data.frame (with one column for regional lengths and BLAST/Transcript IDs/Transcript Names as row.names) or a Named Vector or a Named List. Assumed to have short IDs
#' Calculate HSP Coverage
#'
#' This function calculates the coverage of HSPs (sum(HSP Alignment lengths)/mRNA CDS Length) between each Query->Subject Hits. Hits are filtered based on the minimum coverage filter value.
#'
#' @param blast_table Filename or BLAST Table with Query->Subject Hits
#' @param col.indices A Named List with indices of columns Query sequence ID (qseqid) and Subject sequence ID (sseqid), Query Length (query_len), Subject Length (subject_len) and Alignment Length (HSP alignment length in our case)(align_len). Eg col.indices=list(qseqid=12,sseqid=1,query_len=13,subject_len=14,align_len=23)
#' @param group Name of the group/BLAST Run. Default -  "ungrouped"
#' @param run.mode "both" or "coverage_distance" (Default) or "coverage_filter". "coverage_distance" - Hits are filtered based on distance between bi-directional minimum HSP coverages (coverage_distance <= min_coverage_filter). This option selects more BLAST hits and should be used when the coverage values are very low (and the BLAST Hits/sequences are distant). "coverage_filter" - Filters Hits based on minimum coverage of HSPs from either direction. Use this option when the coverage values are high (and the BLAST Hits/sequences are closely related). "both" - Uses both "coverage_distance" and "coverage_filter" and is very strict
#' @param min_coverage_filter Minimum HSP Coverage value to filter out Hits (Default - 0.5)
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, Default - FALSE otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (ONLY FOR blast_table,transcript_region_lengths are assumed to have short IDs)
#' @param params_list Output of load_params()
#' @param sep Delimiter for the input Files. Only valid if blast_table is a file. Default - '\t'
#' @param header Does the input files have header?. Only valid if blast_table is a file. Default - FALSE
#' @param verbose Print Output Messages?
#' @return BLAST table with Hits which pass min_coverage_filter
#' @export
calculate_HSP_coverage <- function(blast_table,col.indices, group="ungrouped",run.mode="coverage_distance",min_coverage_filter=0.5, COMPLETE.format.ids=F,params_list, sep="\t", header=F, verbose=F){ #transcript_region_lengths

  if(!grepl(pattern ="coverage_distance|coverage_filter|both",ignore.case = T,x = run.mode) || is.null(run.mode)){
    stop("run.mode must be either 'both' or 'coverage_distance' or 'coverage_filter'")
  }

  if(!is.null(blast_table) && is.character(blast_table) && file.exists(blast_table)){
    blast_table <- LoadBLASTHits(infile = blast_table, sep=sep, header=header)
  }else if(!is.null(blast_table) && !is.character(blast_table)){
    blast_table <- blast_table
    #if(!is.null(col.names)){
    #  colnames(blast_table) <- col.names
    #}
  }else{
    stop(paste(blast_table,"does not exist!"))
  }

  # if(is.vector(transcript_region_lengths)){
  #   region_lengths <- as.data.frame(x=transcript_region_lengths,row.names = names(transcript_region_lengths)) #vector
  # }else if(inherits(transcript_region_lengths, "list")){ #is.list(transcript_region_lengths)
  #   region_lengths <- as.data.frame(x=unlist(transcript_region_lengths,use.names = T,recursive = T),row.names = names(transcript_region_lengths)) # list
  # }else if(inherits(transcript_region_lengths, "data.frame") && ncol(transcript_region_lengths) == 1){ #is.data.frame(transcript_region_lengths)
  #   region_lengths <- transcript_region_lengths
  # }else{
  #   stop("transcript_region_lengths must be a Data.Frame (with one column for (CDS/UTR) lengths and BLAST IDs/Transcript IDs as row.names) or a Named Vector or a Named List")
  # }

  query_ids <- unique(blast_table[,col.indices[["qseqid"]]])
  subject_ids <- unique(blast_table[,col.indices[["sseqid"]]])

  #print(region_lengths) #DEBUG
  #print(paste(query_ids,collapse = ",")) #DEBUG
  #print(paste(subject_ids,collapse = ",")) #DEBUG

  id_combinations <- unique(tidyr::crossing(query_ids,subject_ids))
  passed_coverage <- purrr::map2(id_combinations$query_ids, id_combinations$subject_ids, function(x,y){
    if(COMPLETE.format.ids){
      x_short <- stringi::stri_split(str = stringi::stri_split(str = x,fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)[,1], fixed = params_list$TRANSCRIPT_ID_DELIM, simplify=T)[,COMPLETE$FORMAT_ID_INDEX$TRANSCRIPT_ID]
      y_short <- stringi::stri_split(str = stringi::stri_split(str = y,fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)[,1], fixed = params_list$TRANSCRIPT_ID_DELIM, simplify=T)[,COMPLETE$FORMAT_ID_INDEX$TRANSCRIPT_ID]
    }else{
      x_short <- x
      y_short <- y
    }
    if(x_short!=y_short){
      #q_length <- as.numeric(region_lengths[x_short,])
      #s_length <- as.numeric(region_lengths[y_short,])

      #if(verbose) print(paste(c(x,x_short,q_length,y,y_short,s_length),collapse = ":"))

      query_hits <- blast_table[which(!is.na(match(blast_table[,col.indices[["qseqid"]]],x))),]
      subject_hits <- blast_table[which(!is.na(match(blast_table[,col.indices[["sseqid"]]],y))),]

      ##s_overlaps <- dissolve_GR_Overlaps(subject_result)
      #s_align_length <- sum(subject_hits[,col.indices[["send"]]] - subject_hits[,col.indices[["sstart"]]])
      #q_align_length <- sum(query_hits[,col.indices[["qend"]]] - query_hits[,col.indices[["qstart"]]])
      #cov_q <- q_align_length/q_length
      #cov_s <- s_align_length/s_length

      #print(subject_hits[,col.indices[["subject_len"]]])
      #print(subject_hits[,col.indices[["Hsp_align.len"]]])
      #print(subject_hits[,col.indices[["subject_len"]]] / subject_hits[,col.indices[["Hsp_align.len"]]] ) #DEBUG
      #print(query_hits[,col.indices[["query_len"]]])
      #print(query_hits[,col.indices[["Hsp_align.len"]]])
      #print(query_hits[,col.indices[["query_len"]]] / query_hits[,col.indices[["Hsp_align.len"]]] ) #DEBUG
      #print(subject_hits[,col.indices[["align_len"]]])

      cov_s <- 1 / ( sum(subject_hits[,col.indices[["align_len"]]]) / unique(subject_hits[,col.indices[["subject_len"]]]) )  #sum(subject_hits[,col.indices[["send"]]] - subject_hits[,col.indices[["sstart"]]])
      cov_q <- 1 / ( sum(query_hits[,col.indices[["align_len"]]]) / unique(query_hits[,col.indices[["query_len"]]]) ) #sum(query_hits[,col.indices[["qend"]]] - query_hits[,col.indices[["qstart"]]])

      coverage_distance= 1 - (cov_q/cov_s) #1- (1/(cov_q/cov_s)) #sqrt((1-1/(cov_q/cov_s))^2)

      #if(verbose) print(paste(paste(x,"(",cov_q,")",sep = ""),paste(y,"(",cov_s,")",sep = ""),sep="->"),)
      #if(verbose) print(paste(x,"[","q_align_len/q_CDS_length:",q_align_length,"/",q_length,"]",sep=""))
      #if(verbose) print(paste(y,"[","s_align_len/s_CDS_length:",s_align_length,"/",s_length,"]",sep=""))
      if(verbose) print(paste(x,"->",y,"[","cov_q/cov_s:",cov_q,"/",cov_s,"]=",cov_q/cov_s," Coverage Distance:",coverage_distance,sep=""))
      if(verbose) print(coverage_distance)
      #if(cov_q >= min_coverage_filter && cov_s >= min_coverage_filter){
      #  return(data.frame(query=x,subject=y))
      #}

      if(grepl(pattern ="both",ignore.case = T,x = run.mode)){
        if(coverage_distance <= min_coverage_filter && cov_q >= min_coverage_filter && cov_s >= min_coverage_filter ){
          return(data.frame(query=x,subject=y,coverage_distance=coverage_distance,cov_q=cov_q,cov_s=cov_s, min=min(cov_q,cov_s),max=max(cov_q,cov_s),group=group))
        }
      }else if(grepl(pattern ="coverage_distance",ignore.case = T,x = run.mode)){
        #IF Min_Coverage Distance <= min_coverage_filter, return (min coverage distance=0 is full/same coverage between directions, min coverage distance=1 is no coverage)
        if(coverage_distance <= min_coverage_filter ){ #&& cov_q >= min_coverage_filter && cov_s >= min_coverage_filter
          return(data.frame(query=x,subject=y,coverage_distance=coverage_distance,cov_q=cov_q,cov_s=cov_s, min=min(cov_q,cov_s),max=max(cov_q,cov_s),group=group))
        }
      }else if(grepl(pattern ="coverage_filter",ignore.case = T,x = run.mode)){
        if(cov_q >= min_coverage_filter && cov_s >= min_coverage_filter){
          return(data.frame(query=x,subject=y,coverage_distance=coverage_distance,cov_q=cov_q,cov_s=cov_s, min=min(cov_q,cov_s),max=max(cov_q,cov_s),group=group))
        }
      }

    }
  })

  #tmp_passed_coverage <<- passed_coverage
  passed_coverage <- dplyr::bind_rows(passed_coverage[!sapply(passed_coverage, is.null)])
  blast_table <- blast_table[!is.na(match(blast_table[,col.indices[["qseqid"]]],unique(c(passed_coverage$query,passed_coverage$subject)))),]
  blast_table <- blast_table[!is.na(match(blast_table[,col.indices[["sseqid"]]],unique(c(passed_coverage$query,passed_coverage$subject)))),]
  return(list(coverage=passed_coverage, blast_table=blast_table))
}

#all_gtf_stats Coerced GTF stats from all the organisms. Found in paste(params_list$OUT_PATH,"/all_gtf_stats.csv",sep=""). This file is generated by EXTRACT_DATA() and is required for calculating HSP Coverage
#' Internal Function - Transcript Ortholog Extraction Function for R-COMPLETE pipeline
#'
#' This function calls the Transcript Ortholog Extraction pipeline which is used to reduce the pool of genes (step 1), reduce the pool of organisms and create sets of organisms (step 2), find transcript level orthologs (step 3). It takes only one argument which is the path/name of the BLAST program to use and refers to the values from the parameters file for other variables.
#'
#'  * Step 1 - Genes which are available in all the reference organisms are chosen
#'  * Step 2 - A Per-Gene Conservation Score (GSC) is calculated from the availability of a gene across organisms (literally the count of organisms which have the gene, normalized to 1 relative to other genes). Genes which have GSC score below GENE_DROP_THRESHOLD (parameter) are dropped (GENE_DROP_THRESHOLD=0 does not omit any genes). Sets of organisms are created based on the available genes after GSC filtering. I can suggest reference organisms based on which ones have the maximum number of genes
#'  * Step 3 - Two way BLAST followed by HSP selection with WISARD and Two way RBH are performed. Only transcripts which are bi-directional best hits are kept for further analysis (RBH from both the directions, not RBH in itself is bi-directionaly from the point of the QUERY, We can also do an RBH from the context of the SUBJECT to verify if it did not pas RBH by chance (even though it is very unlikely))
#'
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program. BLAST options are taken from params_list$BLAST_OPTIONS
#' @param params_list Output of load_params()
#' @param clusters_left Vector of file names (Set/Subset) in input_dir (OG Clusters/Genes) to BLAST clusters_right with (Can be same as clusters_right)
#' @param clusters_right Vector of file names (Set/Subset) in input_dir (OG Clusters/Genes) to BLAST clusters_left with (Can be same as clusters_left)
#' @param input_dir Give the directory with FASTA files (to BLAST between them using blast_program)
#' @param output_dir Directory for saving output files (\*.out, \*.all2all, \*.wis_out,\*rbh_out)
#' @export
extract_transcript_orthologs <- function(blast_program, params_list, clusters_left, clusters_right, input_dir,output_dir){ #all_gtf_stats

  if(!any(grepl(x = class(params_list), pattern = "COMPLETE-options"))){
    stop("Error: params_list is not a COMPLETE-options class. Use load_params()")
  }

  #unique_lengths <- unique(all_gtf_stats[,c("transcript_id","total_cds_len")])
  #tx_CDS_lengths <- data.frame(length=unique_lengths$total_cds_len, row.names = unique_lengths$transcript_id)

  blast_options <- params_list$BLAST_OPTIONS
  #blast_DB_dir <- params_list$BLAST_DB_PATH

  all2all_BLAST(first_list = clusters_left, second_list = clusters_right,blast_program = blast_program,output_dir =output_dir,blast_options = blast_options,input_prefix_path = input_dir, params_list = params_list, COMPLETE.format.ids = T, keep.output.files = T ) #blast_DB_dir = blast_DB_dir #second_list = grep(clusters_left,clusters_right,ignore.case = T,invert = T,value=T)
  all2all_BLAST(first_list = clusters_right, second_list = clusters_left,blast_program = blast_program,output_dir = output_dir,blast_options = blast_options,input_prefix_path = input_dir, params_list = params_list, COMPLETE.format.ids = T, keep.output.files = T ) #blast_DB_dir = blast_DB_dir #first_list = grep(clusters_left,clusters_right,ignore.case = T,invert = T,value=T)

  # parallel::mclapply(list.files(path = output_dir,pattern = "*.all2all", ignore.case = T,full.names = T),function(in_file){
  #   out_file <- paste(output_dir,tools::file_path_sans_ext(BiocGenerics::basename(in_file)),".out",sep="")
  #   convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
  # }, mc.cores = params_list$numWorkers )

  # parallel::mclapply(list.files(path = output_dir,pattern = "*.out", ignore.case = T,full.names = T),function(in_file){
  parallel::mclapply(list.files(path = output_dir,pattern = c(paste(clusters_left,clusters_right,"all2all",sep="."),paste(clusters_right,clusters_left,"all2all",sep=".")), ignore.case = T,full.names = T),function(in_file){ #c(paste(clusters_left,clusters_right,"all2all",sep="."),paste(clusters_right,clusters_left,"all2all",sep=".")) #list.files(path = output_dir,pattern = "*.all2all", ignore.case = T,full.names = T)
    out_file <- paste(output_dir,tools::file_path_sans_ext(BiocGenerics::basename(in_file)),".wis_out",sep="")
    try(
      if(!file.exists(out_file)){
        blast_GO <- GRObject_from_BLAST(blast_input = in_file, COMPLETE.format.ids = T, col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list = params_list)
        wis_GO <- invisible(run_WISARD(blast_hits = blast_GO,score_col = "Hsp_score",COMPLETE.format.ids = T, params_list = params_list)) #score_col=16)
        wis_GO <- melt_wisard_list(wis_GO)
        write.table(x = wis_GO,file = out_file,quote = F,col.names = T,row.names = F,sep = "\t")
      }else{
        stop(paste(out_file,"exists!"))
      } )
  }, mc.cores = params_list$numWorkers )

  #save("wisard_results", file="files/all2all/wisard_results.RData")
  #load("files/all2all/wisard_results.RData")

  ##RUN RBH
  group_combinations <- unique(tidyr::crossing(clusters_left,clusters_right))
  purrr::map2(group_combinations$clusters_left,group_combinations$clusters_right ,function(query,subject){
    #lapply(clusters_left, function(query){
    #  parallel::mclapply(clusters_right,function(subject){ #grep(clusters_left,clusters_right,ignore.case = T,invert = T,value=T)
    in1 <- paste(output_dir,query,".",subject,".wis_out",sep="")
    in2 <- paste(output_dir,subject,".",query,".wis_out",sep="")
    out1 <- paste(output_dir,query,".",subject,".rbh_out",sep="")
    out2 <- paste(output_dir,subject,".",query,".rbh_out",sep="")
    #print(paste(in1,in2)) #DEBUG
    try(
      if (!file.exists(out1) && !file.exists(out2) && file.exists(in1) && file.exists(in2) && file.info(in1)$size > 0 && file.info(in2)$size > 0 ) {
        #RBH(in1 = in1, in2 = in2, index.tables = T, col.indices = list(qseqid=12,sseqid=1,weight.col=22),col.names = c("subject_id","start","end","width","strand","Hsp_num","Hsp_bit.score","Hsp_score","Hsp_evalue","Hsp_query.from","Hsp_query.to","query_id","query_len","subject_len","Hsp_hit.from","Hsp_hit.to","Hsp_query.frame","Hsp_hit.frame","Hsp_pidentity","Hsp_gaps","Hsp_align.len","max_score"))
        #print(paste(out1,out2))
        RBH_out <- RBH(in1 = in1, in2 = in2, index.tables = T, col.indices = list(qseqid=12,sseqid=1), header = T,n_threads = params_list$numWorkers) #,weight.col=c(22,8) ), unique.hit.weights = T, process.weights.func = max)

        write.table(x = RBH_out$in1,file = out1,quote = F,col.names = T,row.names = F, sep = "\t")
        write.table(x = RBH_out$in2,file = out2,quote = F,col.names = T,row.names = F, sep = "\t")
        #save(RBH_out, tx_CDS_lengths, file = "tmp.RData")
        #calculate_HSP_coverage(RBH_out$in1,transcript_region_lengths = tx_CDS_lengths, col.indices=list(qseqid=12,sseqid=1,qstart=10,qend=11,sstart=2,send=3), COMPLETE.format.ids = T,params_list = params_list)
        #calculate_HSP_coverage(RBH_out$in2,transcript_region_lengths = tx_CDS_lengths, col.indices=list(qseqid=12,sseqid=1,qstart=10,qend=11,sstart=2,send=3), COMPLETE.format.ids = T,params_list = params_list)
      }else{
        stop(paste(out1,"exists!"))
      })
    # }, mc.cores = params_list$numWorkers )
    #})
  })

  final_blast_tables <- purrr::map2(group_combinations$clusters_left,group_combinations$clusters_right ,function(query,subject){
    in_data <- paste(output_dir,query,".",subject,".wis_out",sep="")
    out_data <- paste(output_dir,query,".",subject,".final_out",sep="")
    out_cov_data <- paste(output_dir,query,".",subject,".min_cov",sep="")
    try(if(file.exists(in_data) && !file.exists(out_data) && !file.exists(out_cov_data)){
      final_blast_table <- calculate_HSP_coverage(in_data,col.indices=list(qseqid=12,sseqid=1,query_len=13,subject_len=14,align_len=21), group=paste(query,subject,sep="."),COMPLETE.format.ids = T,params_list = params_list)
      write.table(x = final_blast_table$blast_table,file = out_data,quote = F,col.names = T,row.names = F, sep = "\t")
      write.table(x = final_blast_table$coverage,file = out_cov_data,quote = F,col.names = T,row.names = F, sep = "\t")
      return(final_blast_table)
    }else{
      stop(paste(out_data,"&",out_cov_data,"exists! OR",in_data, "does not exist!"))
    })
  })

  save(final_blast_tables, file="final_blast_tables.RData")

  ##PLOT_CODE
  tmp2 <- c()
  tmp3 <- c()
  tmp4 <- c()
  tmp5 <- c()
  tmp6 <- c()
  plot_filename <- paste(ref_org,org,sep = "-")
  tmp3 <- sapply(HSP_fw, USE.NAMES = T ,function(x){
    return(data.frame(query=x$query,subject=x$subject,min_cov=x$min_cov,max_cov=x$max_cov,cov_q=x$cov_q,cov_s=x$cov_s,same_CDS_count=as.logical(x$same_CDS_count),hit_from=x$hit_from,hit_to=x$hit_to,query_len=x$query_len,subject_len=x$subject_len,q_align_len=x$q_align_len,s_align_len=x$s_align_len,query_from=x$query_from,query_to=x$query_to, pident=x$pident))
  })
  tmp4 <- sapply(HSP_bk, USE.NAMES = T ,function(x){
    return(data.frame(query=x$query,subject=x$subject,min_cov=x$min_cov,max_cov=x$max_cov,cov_q=x$cov_q,cov_s=x$cov_s,same_CDS_count=as.logical(x$same_CDS_count),hit_from=x$hit_from,hit_to=x$hit_to,query_len=x$query_len,subject_len=x$subject_len,q_align_len=x$q_align_len,s_align_len=x$s_align_len,query_from=x$query_from,query_to=x$query_to, pident=x$pident))
  })
  #print(str(tmp3))
  if(length(tmp3)>0){ ##IF length(tmp3) == 0 then there are no genes matching between organisms
    for(i in 1:ncol(tmp3)){
      #print(data.frame(tmp3[,i],colnames(tmp3)[i]))
      #print(colnames(tmp3)[i])
      #print(data.frame(t(tmp3[[i]]),colnames(tmp3)[i]))
      tmp2 <- rbind(data.frame(tmp3[,i],gene=colnames(tmp3)[i]),tmp2)
    }
    tmp2 <- tmp2 %>% mutate(ref_org=rep(ref_org)) %>% mutate(org=rep(org)) %>% mutate(direction=rep("forward"))
    if(length(tmp4)>0){
      for(i in 1:ncol(tmp4)){
        #print(data.frame(tmp3[,i],colnames(tmp3)[i]))
        #print(colnames(tmp3)[i])
        #print(data.frame(t(tmp3[[i]]),colnames(tmp3)[i]))
        tmp5 <- rbind(data.frame(tmp4[,i],gene=colnames(tmp4)[i]),tmp5)
      }
      tmp5 <- tmp5 %>% mutate(ref_org=rep(org)) %>% mutate(org=rep(ref_org)) %>% mutate(direction=rep("backward"))
      #tmp6 <- rbind(tmp2,tmp5)
      if(mean(tmp2$min_cov)>=mean(tmp5$min_cov)){
        tmp6 <- tmp2
      }else{
        tmp6 <- tmp5
      }}else{
        tmp6 <- tmp2
      }
    #print(tmp6[which(tmp6$min_cov > 1),])
    ##Remove outliers
    tmp6 <- tmp6[!is.na(tmp6$min_cov),]
    tmp6 <- tmp6[!is.na(tmp6$same_CDS_count),]
    #tmp6 <- tmp6[which(tmp6$min_cov <= 1),] #tmp6[!which(tmp6$min_cov > 1),]
    #tmp6 <- tmp6[which(tmp6$max_cov <= 1),]
    tmp6 <- unique(tmp6)

    #print(head(tmp6))
    #ggplot(tmp6, aes(x=min_cov*100, fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion") + facet_wrap(~direction)
    if(nrow(tmp6)>0){
      #print(head(tmp6))

      #print(tmp2[which(is.na(tmp2)),])
      #DENSITY PLOT FOR ALL GENES COLOURED BY SAME_CDS_COUNT
      density_plot <- ggplot(tmp6, aes(x=min_cov*100, fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion") + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))
      density_plot_max <- ggplot(tmp6, aes(x=max_cov*100, fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Max.Coverage of all transcripts(pairwise)") + ylab("Proportion") + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))

      #GROUPED BY GENE
      genewise_plot <- ggplot(tmp6, aes(x=min_cov*100, group=gene ,fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion") + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))
      genewise_plot_max <- ggplot(tmp6, aes(x=max_cov*100, group=gene ,fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Max.Coverage of all transcripts(pairwise)") + ylab("Proportion") + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))

      # facet_wrap(~gene, drop=FALSE, scales=c("free"))
      total_pages <- n_pages(genewise_plot + facet_wrap_paginate(facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y")))
      total_pages_max <- n_pages(genewise_plot_max + facet_wrap_paginate(facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y")))
      #GROUPED BY SAME_CDS_COUNT
      cds_count_plot <- ggplot(tmp6, aes(x=min_cov*100 ,group=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity') + xlab("Same_CDS_Count") + ylab("Proportion") + facet_wrap(~same_CDS_count, drop=FALSE, scales=c("free")) + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))
      cds_count_plot_max <- ggplot(tmp6, aes(x=max_cov*100 ,group=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity') + xlab("Same_CDS_Count") + ylab("Proportion") + facet_wrap(~same_CDS_count, drop=FALSE, scales=c("free")) + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))

      ##Do a boxplot with facet_wrap(~ref organisms)
      #print(head(tmp6))
      if(file.exists(paste(plot_out_path,paste(plot_filename,".pdf",sep=""),sep="/"))){
        file.remove(paste(plot_out_path,paste(plot_filename,".pdf",sep=""),sep="/"))
      }

      pdf(file =paste(plot_out_path,paste(plot_filename,".pdf",sep=""),sep="/"),title = plot_filename)
      plot.new()
      text(.5, .5, "MINIMUM_COVERAGE")
      print(density_plot)
      #print(total_pages)
      if(total_pages > 0){
        for(i in 1:total_pages){
          try(print(genewise_plot + facet_wrap_paginate(facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y"), page=i))) #facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,scales=c("free_y"),
        }
      }
      #print(genewise_plot)
      print(cds_count_plot)
      plot.new()
      text(.5, .5, "MAXIMUM_COVERAGE")
      print(density_plot_max)
      #print(total_pages)
      if(total_pages_max > 0){
        for(i in 1:total_pages_max){
          try(print(genewise_plot_max + facet_wrap_paginate(facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y"), page=i))) #facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,scales=c("free_y"),
        }
      }
      #print(genewise_plot)
      print(cds_count_plot_max)
      plot.new()
      text(.5, .5, "BLAST Coverage Plots")
      for(row_q in unique(tmp6$query)){
        for(row_s in unique(tmp6$subject)){
          subset_data <- tmp6[which(tmp6$query==row_q & tmp6$subject==row_s),]
          #print(subset_data)
          if(nrow(subset_data)>0){
            subset_data <- subset_data %>% mutate(group=1:nrow(subset_data))
            line_coords <- pivot_longer(subset_data[,c("query","subject","hit_from","hit_to","query_from","query_to","pident","min_cov","group","max_cov")], cols = c("hit_from","hit_to","query_from","query_to") ,
                                        names_to = "direction", values_to = "coords")
            line_coords$direction[which(line_coords$direction=="hit_from" | line_coords$direction=="hit_to")] <- "hit"
            line_coords$direction[which(line_coords$direction=="query_from" | line_coords$direction=="query_to")] <- "query"
            data_rows=nrow(subset_data)
            #print(data_rows)
            line_coords <- data.frame(query=rep(unique(line_coords$query),2*data_rows),subject=rep(unique(line_coords$subject),2*data_rows),from=line_coords$coords[which(line_coords$direction=="query")],to=line_coords$coords[which(line_coords$direction=="hit")], groups=rep(1:data_rows,each=2),pident=rep(subset_data$pident,each=2),min_cov=rep(subset_data$min_cov,each=2),max_cov=rep(subset_data$max_cov,each=2))
            line_coords$groups <- factor(line_coords$groups)
            line_coords$min_cov <- as.numeric(line_coords$min_cov)
            line_coords$to <- as.numeric(line_coords$to)
            line_coords$from <- as.numeric(line_coords$from)
            plot_labels <- as.vector(t(apply(subset_data[,c("min_cov","max_cov")], MARGIN = c(1,2),FUN=function(x){
              return(round(x*100,2))
            })))
            try(print(ggplot(line_coords,aes(x=from,y=to,group=groups,color=pident)) + ylim(0,unique(subset_data$query_len)) + xlim(0,unique(subset_data$subject_len)) +
                        geom_line(na.rm = T) + geom_point(na.rm = T) + geom_text(aes(label=plot_labels)) + ylab(paste(unique(line_coords$query),"(",unique(subset_data$query_len),")")) + xlab(paste(unique(line_coords$subject),"(",unique(subset_data$subject_len),")"))))
            HSP <- rbind(HSP,subset_data)
          }
        }
      }

      dev.off()
    }
  }
  ###PLOT-CODE OVER
  #create organism sets

  #calculate cluster occupancy - number of genes per cluster && number of organisms per cluster

  # ##ITERATION 2 - two way RBH
  #
  # #cat $reference_ORGS files/oneway/SET > files/oneway/set.tmp
  #
  # #time all2all_refblast $reference_ORGS $fasta_path files/all2all/all2all.genelist files/all2all_final $blastdb_path $region $region tblastx files/oneway/set.tmp
  #
  # #time Rscript wisard.R files/oneway/set.tmp files/all2all_final files/gtf_stats.csv
  #
  # #readarray set_orgs < files/oneway/set.tmp
  #
  # #time parallel -j $((${#ref_orgs[@]}*${#set_orgs[@]})) "twoway_RBH files/all2all_final {1} {2} $PY3_PATH $RBH_SCRIPT $SAME_GENE" ::: ${ref_orgs[@]} ::: ${set_orgs[@]}


}
#' Internal Function - Select groups which contain any of the reference organisms
#'
#' This function removes groups from params_list$GROUPS_PATH which does not contain any/all of the reference organisms (based on params_list$SELECT_REF_ORG_GROUPS_METHOD)
#'
#' @param params_list Output of load_params()
select_ref_org_groups <- function(params_list){
  run_status <- processx::run( command = COMPLETE$SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"select_ref_org_groups",params_list$REF_ORGS_FILE, params_list$GROUPS_PATH,COMPLETE$parallel, params_list$SELECT_REF_ORG_GROUPS_METHOD, params_list$numWorkers) ,spinner = T,stdout = "",stderr = "")
  if(run_status$status>0){
    stop("Error in select_ref_org_groups()")
  }else{
    return(unique((tools::file_path_sans_ext(list.files(path=params_list$GROUPS_PATH,full.names = F,recursive = F,include.dirs = F)))))
  }
}

# #' Group FASTA Sequences into Clusters (Ortholog Clusters)
# #'
# #' Groups the Sequences in FASTA files in a folder(recursively) into Clusters and save the Cluster FASTA files in params_list$GROUPS_PATH. Wrapper function for group_FASTA(). This function is used with run.mode="cluster" in FiND_TRANSCRIPT_ORTHOLOGS(. FASTA IDs are required to be in R-COMPLETE's long format (?COMPLETE_PIPELINE_DESIGN)
# #'
# #' @param params_list Output of load_params()
# #' @return Vector of all available OG Clusters
# #' @export
# group_FASTA_clusters <- function(params_list){
#   tictoc::tic(msg = paste("Grouping FASTA into Clusters ..."))
#   group_FASTA(params_list, id.col.index=4)
#   cat(print_toc(tictoc::toc(quiet = T)))
#   message(paste("Ortholog Clusters are stored in :", params_list$GROUPS_PATH))
#
#   all_clusters_list <- parallel::mclapply(list.files(path = paste(params_list$OUT_PATH,"/genes/",sep=""),include.dirs=TRUE, full.names=TRUE),function(x){
#     if(file.exists(paste(x,"/ORG_CLUSTERS.4",sep="")) && file.info(paste(x,"/ORG_CLUSTERS.4",sep=""))$size > 0 ){
#       return(scan(paste(x,"/ORG_CLUSTERS.4",sep=""), character(), quiet = T))
#     }
#   }, mc.cores =  params_list$numWorkers, mc.preschedule = T)
#   all_clusters <- purrr::reduce(all_clusters_list, unique)
#   return(all_clusters)
# }
#
# #' Group FASTA Sequences into Genes (Gene Names)
# #'
# #' Groups the Sequences in FASTA files in a folder(recursively) into Gene Names and save the Gene-FASTA files in params_list$GROUPS_PATH. Wrapper function for group_FASTA(). This function is used with run.mode="gene" in FIND_TRANSCRIPT_ORTHOLOGS(). FASTA IDs are required to be in R-COMPLETE's long format (?COMPLETE_PIPELINE_DESIGN)
# #'
# #' @param params_list Output of load_params()
# #' @return Vector of all available Genes
# #' @export
# group_FASTA_genes <- function(params_list){
#   tictoc::tic(msg = paste("Grouping FASTA into Genes ..."))
#   group_FASTA(params_list, id.col.index=3)
#   cat(print_toc(tictoc::toc(quiet = T)))
#   message(paste("Gene Clusters are stored in :", params_list$GROUPS_PATH))
#
#   all_genes_list <- parallel::mclapply(list.files(path = paste(params_list$OUT_PATH,"/genes/",sep=""),include.dirs=TRUE, full.names=TRUE),function(x){
#     if(file.exists(paste(x,"/ORG_CLUSTERS.3",sep="")) && file.info(paste(x,"/ORG_CLUSTERS.3",sep=""))$size > 0 ){
#       return(scan(paste(x,"/ORG_CLUSTERS.3",sep=""), character(), quiet = T))
#     }
#   }, mc.cores =  params_list$numWorkers, mc.preschedule = T)
#   all_genes <- purrr::reduce(all_genes_list, unique)
#   return(all_genes)
# }

#' Group FASTA Sequences Based on a column index of COMPLETE.format.ids (?COMPLETE_PIPELINE_DESIGN)
#'
#' group_FASTA_genes() for run.more="gene" & group_FASTA_clusters() for run.mode="cluster" wrap this function. run.mode is used in FIND_TRANSCRIPT_ORTHOLOGS(). This function write groupings (clusters/genes) of each organism (org_name) into paste(params_list$OUT_PATH,"/genes/",org_name,"/ORG_CLUSTERS.", id.col.index,sep="")
#'
#' @param params_list Output of load_params()
#' @param id.col.index The index of Column of COMPLETE.format.ids (?COMPLETE_PIPELINE_DESIGN) to groups sequences into. Use id.col.index=1 for grouping sequences based on Transcript IDs, id.col.index=2 to group sequences into Organisms, id.col.index=3 for grouping sequences based on Gene Names and id.col.index=4 to groups sequences into Ortholog Clusters. Check COMPLETE$FORMAT_ID_INDEX for indices
#' @param verbose Print DEBUG Messages?
#' @export
group_FASTA <- function(params_list, id.col.index, verbose=F){
  if(!any(grepl(x = class(params_list), pattern = "COMPLETE-options"))){
    stop("Error: params_list is not a COMPLETE-options class. Use load_params()")
  }

  id.col.index <- as.numeric(id.col.index)
  grouping_by <- names(COMPLETE$FORMAT_ID_INDEX[id.col.index])
  tictoc::tic(msg = paste("Grouping FASTA into", grouping_by,"..."))
  unlink(x = params_list$GROUPS_PATH,recursive = T,force = T,expand = T)
  dir.create(path = params_list$GROUPS_PATH,showWarnings = F,recursive = T)
  fasta_files <- list.files(path = params_list$FASTA_OUT_PATH,all.files = T,full.names = T,recursive = T,include.dirs = F)
  parallel::mclapply(fasta_files, function(x){
    fasta_recs <- Biostrings::readDNAStringSet(filepath = x,use.names = T, format = "fasta")
    #print(names(fasta_recs)) #DEBUG
    split_recs <- stringi::stri_split(str = names(fasta_recs), fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)
    if(ncol(split_recs)==length(COMPLETE$FORMAT_ID_INDEX)){ ##CHECKING IF FASTA IDs are COMPLETE.format.ids
      org_name <- unique( split_recs[, COMPLETE$FORMAT_ID_INDEX$ORG] )
      #print(org_name) #DEBUG
      all_clusters <- unique(unlist(purrr::map2(seq_along(fasta_recs),names(fasta_recs), function(rec_num, rec_name){
        split_rec <- stringi::stri_split(str = rec_name, fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)
        rec_clusters <- unique(stringi::stri_split(str = split_rec[, id.col.index], fixed = ",", simplify=T)) #ncol(split_rec)
        #print(rec_clusters) #DEBUG
        lapply(rec_clusters, function(each_cluster){
          #print(paste(params_list$GROUPS_PATH,"/",each_cluster,".",tools::file_ext(x),sep = ""))  #DEBUG
          Biostrings::writeXStringSet(x = fasta_recs[rec_num],filepath = paste(params_list$GROUPS_PATH,"/",each_cluster,".",tools::file_ext(x),sep = ""), append = T,format = "fasta")
        })
        return(rec_clusters)
      })))
      write.table(x = all_clusters,file = paste(params_list$OUT_PATH,"/genes/",org_name,"/ORG_CLUSTERS.", grouping_by,sep=""),quote = F,row.names = F,col.names = F)
    }else{ #DEBUG
      #print(x) #DEBUG
      if(verbose){
        message(paste(x," : does not have COMPLETE.format.ids"))
      }
    } #DEBUG
  }, mc.cores = params_list$numWorkers,mc.preschedule = T,mc.silent = !verbose)
  cat(print_toc(tictoc::toc(quiet = T)))

  # all_groups_list <- parallel::mclapply(list.files(path = paste(params_list$OUT_PATH,"/genes/",sep=""),include.dirs=TRUE, full.names=TRUE),function(x){
  #   if(file.exists(paste(x,"/ORG_CLUSTERS.",grouping_by,sep="")) && file.info(paste(x,"/ORG_CLUSTERS.",grouping_by,sep=""))$size > 0 ){
  #     return(scan(paste(x,"/ORG_CLUSTERS.",grouping_by,sep=""), character(), quiet = T))
  #   }
  # }, mc.cores =  params_list$numWorkers, mc.preschedule = T, mc.silent = !verbose)
  all_groups <- unique((tools::file_path_sans_ext(list.files(path=params_list$GROUPS_PATH,full.names = F,recursive = F,include.dirs = F)))) #unique(unlist(all_groups_list,recursive = T)) #purrr::reduce(all_groups_list, union)
  write.table(x = all_groups,file = paste(params_list$OUT_PATH,"/ALL_GROUPS.txt",sep=""), quote = F, row.names = F,col.names = F,na = "-")
  return(all_groups)
}

#' (2) - Find Transcript Orthologs
#'
#' This function can be executed after COMPLETE::EXTRACT_DATA() and is the continuation of R-COMPLETE pipeline
#'
#' This is the main function which calls all the other functions and performs and end-end execution of finding transcript level orthologs. It runs the iterative Transcript Ortholog Extraction pipeline which is used to reduce the pool of genes, reduce the pool of organisms, find transcript level orthologs (check ?extract_transcript_orthologs)
#'
#' @note ONLY USE THIS FUNCTION WHEN RUNNING THE PIPELINE OF R-COMPLETE. Use other helper function to work with custom BLAST files not generated by this R package. run.mode="gene" is NOT RECOMMENDED because the sequences are grouped based on gene names
#'
#' @param params_list Filename of a formatted parameter file (check the github repo for an example) or Output of load_params().
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program. BLAST options are taken from params_list. Default is Sys.which("tblastx")
#' @param gene_list Vector or File with a list of genes to extract data for(check the github repo for an example).
#' @param run.mode A value from COMPLETE$FORMAT_ID_INDEX. Default - COMPLETE$FORMAT_ID_INDEX$CLUSTERS. Find transcript orthologs in the level of Orgs, Genes or Ortholog Clusters. Genes have more tight orthology and fewer transcript orthologs which may be very similar. Ortholog Clusters are a level higher than Genes (Because an Ortholog Cluster can have more than one gene) and have highest number of transcript orthologs with a lot of dissimilarity. run.mode=COMPLETE$FORMAT_ID_INDEX$GENE is NOT RECOMMENDED because the sequences are grouped based on gene names, while run.mode=COMPLETE$FORMAT_ID_INDEX$CLUSTERS groups sequences based on protein identity.
#' @export
FIND_TRANSCRIPT_ORTHOLOGS <- function(params_list, blast_program=Sys.which("tblastx"), gene_list, run.mode=COMPLETE$FORMAT_ID_INDEX$CLUSTERS){
  set.seed(123)

  if(is.na(match(COMPLETE$FORMAT_ID_INDEX$CLUSTERS,COMPLETE$FORMAT_ID_INDEX)) || is.null(run.mode)) {
    stop(paste("run.mode must be one of COMPLETE$FORMAT_ID_INDEX"))
  }

  if(!grepl(x=Sys.info()["sysname"],pattern="linux",ignore.case = T)){
    stop("Pipeline only supports Linux (and bash) :(")
  }

  if(is.character(params_list)){
    if(!file.exists(params_list) || file.info(params_list)$size < 0){
      stop("ERROR: Parameters file is missing and is required\n")
    }
    loaded_PARAMS <- load_params(param_file = params_list)
  }else{
    if(any(grepl(x = class(params_list), pattern = "COMPLETE-options"))){
      loaded_PARAMS <- params_list
    }else{
      stop("Error: params_list not valid. Use load_params()")
    }
  }
  print(loaded_PARAMS)

  #loaded_PARAMS$gene_list <- gene_list
  #loaded_PARAMS$genes <- factor(scan(gene_list, character())) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
  #loaded_PARAMS$genes <- genes[grep("gene",tolower(genes), invert = T, fixed = T)]

  all2all_out <- paste(loaded_PARAMS$OUT_PATH,"/all2all/",sep = "")
  all2allfinal_out <- paste(loaded_PARAMS$OUT_PATH,"/all2all_final/",sep = "")

  #print(paste("MAX PROCESSES:",loaded_PARAMS$numWorkers))

  if(loaded_PARAMS$CLEAN_EXTRACT){
    #unlink(paste(loaded_PARAMS$OUT_PATH,"/oneway",sep = ""), recursive = T,force = T,expand = T)
    unlink(all2all_out, recursive = T,force = T,expand = T)
    unlink(all2allfinal_out, recursive = T,force = T,expand = T)
    unlink(paste(loaded_PARAMS$OUT_PATH,"/gene_thresholds.txt",sep=""), recursive = T,force = T,expand = T)
    unlink(paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), recursive = T,force = T,expand = T)
  }

  #dir.create(paste(loaded_PARAMS$OUT_PATH,"/oneway",sep = ""),showWarnings = F, recursive = T)
  dir.create(all2all_out,showWarnings = F, recursive = T)
  dir.create(all2allfinal_out,showWarnings = F, recursive = T)

  if (file.exists(blast_program)) {

    install_parallel()

    REF_ORGS <- factor(scan(loaded_PARAMS$REF_ORGS_FILE, character(), quiet = T))

    # tryCatch({
    #   all_gtf_stats <- read.table(file = paste(loaded_PARAMS$OUT_PATH,"/all_gtf_stats.csv",sep=""),sep = " ",header = T,quote = "",fill = T,na.strings = "-")
    # }, error= function(cond){
    #   stop(paste("Coerced GTF stats file", paste(loaded_PARAMS$OUT_PATH,"/all_gtf_stats.csv",sep=""),"not found/invalid, Please rerun EXTRACT_DATA()"))
    # })

    available_genes_list <- parallel::mclapply(paste(loaded_PARAMS$OUT_PATH,"/genes/",REF_ORGS,sep=""),function(x){
      if(file.exists(paste(x,"/AVAILABLE_GENES",sep="")) && file.info(paste(x,"/AVAILABLE_GENES",sep=""))$size > 0 ){
        return(scan(paste(x,"/AVAILABLE_GENES",sep=""), character(),quiet = T))
      }
    }, mc.cores =  loaded_PARAMS$numWorkers, mc.preschedule = T)
    available_genes <- unique(purrr::reduce(available_genes_list, union))

    #find which clusters they belong to
    # if(grepl(pattern = "cluster", x=run.mode)) {
    #   # odb_map_list <- parallel::mclapply(paste(loaded_PARAMS$OUT_PATH,"/genes/",REF_ORGS,sep=""),function(x){
    #   #   if(file.exists(paste(x,"/odb.final_map",sep="")) && file.info(paste(x,"/odb.final_map",sep=""))$size > 0 ){
    #   #     return(read.table(paste(x,"/odb.final_map",sep=""),header = F,sep = "\t",quote = "", col.names = c("cluster","genes")))
    #   #   }
    #   # }, mc.cores =  loaded_PARAMS$numWorkers, mc.preschedule = T)
    #   # odb_map <- unique(purrr::reduce(odb_map_list, inner_join, by = c("cluster","genes")))
    #   # available_clusters <- unique(BiocGenerics::unlist(lapply(available_genes, function(gene){
    #   #   return(odb_map[grep(pattern=gene,x = odb_map$genes,ignore.case = T, value = F), c("cluster")])
    #   # })))
    #
    #   #if(length(dir(loaded_PARAMS$GROUPS_PATH)==0)){
    #GROUPING FASTA SEQUENCES
    #print(run.mode) #DEBUG
    available_clusters <- c()
    tryCatch({
      all_clusters <- scan(paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), character(),quiet = T)
      clusters_in_dir <- unique((tools::file_path_sans_ext(list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = F,recursive = F,include.dirs = F)))) #unique(stringi::stri_split(str = list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = F,recursive = F,include.dirs = F), simplify=T, fixed = ".")[,1])
      if(length(which(!is.na(match(all_clusters,clusters_in_dir)))) != length(which(!is.na(match(clusters_in_dir,all_clusters)))) || length(clusters_in_dir) < 1 || length(all_clusters) < 1){
        stop("Regrouping clusters...")
      }
    },error=function(cond){
      message(cond)
      all_clusters <- group_FASTA(params_list = loaded_PARAMS, id.col.index = as.numeric(run.mode))
      #print(all_clusters) #DEBUG
      write.table(x = all_clusters,file = paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), quote = F, row.names = F,col.names = F,na = "-")
    }, finally = {
      #print(all_clusters) #DEBUG
      #unlink(x = grep(x = list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F) , pattern = "cds", ignore.case = T,invert = T, value = T), recursive = F,force = T,expand = T)
      non_cds_file_list <- grep(x = list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F) , pattern = "cds", ignore.case = T,invert = T, value = T)
      dir.create(path = file.path(loaded_PARAMS$GROUPS_PATH,"groups_noncds"), showWarnings = F,recursive = T)
      file.rename(non_cds_file_list,paste(loaded_PARAMS$GROUPS_PATH,"/groups_noncds/",basename(non_cds_file_list),sep = ""))
      available_clusters <- all_clusters
    })
    #}

    #available_clusters <- scan(paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), character(),quiet = T) #unique((tools::file_path_sans_ext(list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = F,recursive = F,include.dirs = F))))
    if(length(available_clusters)==0){
      #available_clusters <- "ungrouped"
      stop("No clusters were found!. Try other values for run.mode")
    }

    #STEP 1 - ONLY for run.mode="cluster" - Place ungrouped sequences into groups (all2allblast BLAST ungrouped cluster againts all clusters)
    #all2allblast and then wisard and then RBH for grouping ungrouped clusters
    if(any(grepl(pattern = "ungrouped",x = available_clusters,ignore.case = T))){
      message("STEP 1 - Placing ungrouped sequences into groups\n")
      tictoc::tic(msg = "Placing ungrouped sequences into groups...")
      parallel::mclapply(available_clusters, function(x){
        extract_transcript_orthologs(blast_program = blast_program, params_list = loaded_PARAMS,clusters_left = "ungrouped",clusters_right = x,input_dir = loaded_PARAMS$GROUPS_PATH,output_dir = all2all_out)
      }, mc.cores = loaded_PARAMS$numWorkers)
      cat(print_toc(tictoc::toc(quiet = T)))
    }else{
      message("STEP 1 - Skipped because all sequences are grouped or run.mode != COMPLETE$FORMAT_ID_INDEX$CLUSTERS\n")
    }
    # }else{
    #   message("STEP 1 - Skipped because run.mode='gene'\n")
    #   all_genes <- group_FASTA(params_list = loaded_PARAMS, id.col.index = COMPLETE$FORMAT_ID_INDEX$GENE)
    #   available_clusters <- list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F)
    # }

    if(loaded_PARAMS$SELECT_REF_ORG_GROUPS){
      all_clusters <- select_ref_org_groups(loaded_PARAMS)
      write.table(x = all_clusters,file = paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), quote = F, row.names = F,col.names = F,na = "-")
    }
    available_clusters <- all_clusters

    #STEP 2 - Select transcript level orthologs with minimum coverage between clusters/genes
    message(paste("STEP 2 - Select transcript level orthologs between orgs/genes/clusters, based on minimum coverage\n"))
    tictoc::tic(msg = "Extracting Transcript Orthologs...")
    parallel::mclapply(available_clusters, function(x){
      extract_transcript_orthologs(blast_program = blast_program, params_list = loaded_PARAMS,clusters_left = x,clusters_right = x,input_dir = loaded_PARAMS$GROUPS_PATH,output_dir = all2allfinal_out)
    }, mc.cores = loaded_PARAMS$numWorkers)
    cat(print_toc(tictoc::toc(quiet = T)))

    #calculate gene conservation - calculate_gene_conservation.R - probably not needed
    ##Maybe write one for cluster conservation/coverage across organisms

  }else{
    stop(paste(blast_program," NOT found."))
  }

}

# cppFunction('RObject Multi_Map(NumericVector x, NumericVector y, Function fill_RBH_mat) {
#   int x_len = x.size(), y_len = y.size();
#   RObject out(x_len,y_len);
#
#   for (int i = 0; i < x_len; i++) {
#     for (int j = 0; j < y_len; j++) {
#       out[i,j] = Rcpp::List::create(fill_RBH_mat(x[i], y[j]));
#     }
#   }
#   return out;
# }')

#' Find Reciprocal Blast Hits (RBH)
#'
#' Find RBH between BLAST results of different organisms/genes/transcripts (FASTA/FASTQ). The BLAST results must be of the format 6 and can be converted from BLAST format 11 with convert_BLAST_format().
#' The command with the required column names are given below.
#'
#' convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'
#' Optional : You can provide the file/table with indexed Transcsript IDs. The format must be "file"[tab]"long_id"[tab]"index" (without a header). It can be generated with the  index_FASTA_IDs() (check ?index_FASTA_IDs or index_fastaIDs() in system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"))
#'
#' @note Both the input files are expected to have a header, the same number of columns with matching column order (and names). Give only indices for col.indices. Order of execution is unique.hit.weights followed by process.weights.func (if any/all these options are set), i.e Unique Weights are chosen for each hit (if unique.hit.weights=T) and then weights are processed using process.weights.func (if process.weights.func is set). If the function fails, try indexing the tables with index.tables=T. Weight Columns in col.indices are optional, when not given all the reciprocal hits are returned
#'
#' @examples
#'    convert_BLAST_format(in_file1,outfile = out_file1,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'    convert_BLAST_format(in_file2,outfile = out_file2,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'    blast_GO1 <- GRObject_from_BLAST(blast_input = out_file1, COMPLETE.format.ids = T, col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
#'    blast_GO2 <- GRObject_from_BLAST(blast_input = out_file2, COMPLETE.format.ids = T, col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
#'    wis1 <- run_WISARD(blast_hits = blast_GO1,score_col = "Hsp_score",COMPLETE.format.ids = T) #score_col=16
#'    wis2 <- run_WISARD(blast_hits = blast_GO2,score_col = "Hsp_score",COMPLETE.format.ids = T) #score_col=16
#'    wis1 <- melt_wisard_list(wis1)
#'    wis2 <- melt_wisard_list(wis2)
#'    RBH(in1 = wis1, in2 = wis2, index.tables = T, col.indices = list(qseqid=12,sseqid=1,weight.col=c(8,22)),col.names = c("subject_id","start","end","width","strand","Hsp_num","Hsp_bit.score","Hsp_score","Hsp_evalue","Hsp_query.from","Hsp_query.to","query_id","query_len","subject_len","Hsp_hit.from","Hsp_hit.to","Hsp_query.frame","Hsp_hit.frame","Hsp_pidentity","Hsp_gaps","Hsp_align.len","max_score"),unique.hit.weights = T, process.weights.func = max)
#'
#' @param in1 Input (query->subject) BLAST/WISARD hits table/filename
#' @param in2 Input (query<-subject) BLAST/WISARD hits table/filename
#' @param sep Delimiter for the input Files. Only valid if in1 and in2 are files. Default - '\t'
#' @param header Does the input files have header?. Only valid if in1 and in2 are files. Default - FALSE
#' @param transcript_ID_metadata Tab-delimited File with the filenames, indexed transcript IDs and the long transcript IDs.
#' @param col.names Columns names for the BLAST/WISARD tables/files
#' @param col.indices A Named List with indices of columns Query sequence ID (qseqid), Subject sequence ID (sseqid), and columns to be used as edge weights (eg, "Hsp_score","max_score" etc). Eg col.indices=list(qseqid=1,sseqid=2,weight.col=c(8,22)) OR col.indices=list(qseqid=1,sseqid=2)
#' @param unique.hit.weights Should only the unique Weights be taken for all Query->Subject Hits? (TRUE/FALSE (Default)). Only valid if weight.col is given in col.indices
#' @param process.weights.func Pass a function name to process the weights (eg, sum/max/min etc) (Default - max). Only valid if weight.col is given in col.indices
#' @param index.tables Should the IDs in the tables be indexed? (TRUE (Default) if COMPLETE.format.ids/Long BLAST Sequence IDs are used)
#' @param n_threads Number of Threads (Optional)
#' @return Named List list(in1,in2) with the selected Hits which are RBHs between in1 and in2 data
#' @export
RBH <- function(in1,in2,sep="\t",header=F, transcript_ID_metadata=NULL, col.names=NULL, index.tables=T,col.indices, unique.hit.weights=F, process.weights.func=max, n_threads=tryCatch(parallel::detectCores(all.tests = T, logical = T), error=function(cond){return(2)})){

  #print(col.indices)

  if(!is.null(transcript_ID_metadata) && is.character(transcript_ID_metadata)){
    transcript_ID_metadata <- read.table(file = transcript_ID_metadata,header = F,sep="\t",quote = "")
  }

  # if(!is.null(in1) && is.character(in1)){
  #   in1_data <- LoadBLASTHits(infile = in1, transcript_ID_metadata = transcript_ID_metadata, col.names = col.names)
  # }else if(!is.null(in1)){
  #   in1_data <- in1
  #   if(!is.null(col.names)){
  #     colnames(in1_data) <- col.names
  #   }
  # }
  # if(!is.null(in2) && is.character(in2)){
  #   in2_data <- LoadBLASTHits(infile = in2, transcript_ID_metadata = transcript_ID_metadata, col.names = col.names)
  # }else if(!is.null(in2)){
  #   in2_data <- in2
  #   if(!is.null(col.names)){
  #     colnames(in2_data) <- col.names
  #   }
  # }

  in_data <- purrr::map(c(in1,in2), function(in_file){
    if(!is.null(in_file) && is.character(in_file) && file.exists(in_file)){
      in_file_data <- LoadBLASTHits(infile = in_file, transcript_ID_metadata = transcript_ID_metadata, col.names = col.names, sep=sep, header=header)
    }else if(!is.null(in_file)){
      in_file_data <- in_file
      if(!is.null(col.names)){
        colnames(in_file_data) <- col.names
      }
    }
  })

  #print(head(in_data)) #DEBUG

  # #Sanity checks for weight.col - Converted to numeric indices
  if(!is.list(col.indices)){
    stop("col.indices must be a Named List!")
  }
  if(!all(stringi::stri_cmp_eq(colnames(in_data[[1]])[col.indices[["weight.col"]]],colnames(in_data[[2]])[col.indices[["weight.col"]]]))){
    stop("Column indices of weight.col and order of columns must match between data!")
  }

  if(index.tables){
    # in_data[[1]] <- index_BLAST_table(in_data[[1]],col.indices[["sseqid"]],offset = 0)
    # in_data[[1]] <- index_BLAST_table(in_data[[1]],col.indices[["qseqid"]],offset = length(levels(factor( in_data[[2]][,col.indices[["qseqid"]]] ))) + length(levels(factor( in_data[[1]][,col.indices[["sseqid"]]] ))))
    # in_data[[2]] <- index_BLAST_table(in_data[[2]],col.indices[["sseqid"]],offset = length(levels(factor( in_data[[1]][,col.indices[["qseqid"]]] ))) + length(levels(factor( in_data[[2]][,col.indices[["sseqid"]]] ))))
    # in_data[[2]] <- index_BLAST_table(in_data[[2]],col.indices[["qseqid"]],offset = length(levels(factor( in_data[[2]][,col.indices[["qseqid"]]] ))) + length(levels(factor( in_data[[1]][,col.indices[["sseqid"]]] ))))
    in_data <- index_BLAST_tables(blast_tables = in_data,query_index_cols = col.indices[["qseqid"]], subject_index_cols = col.indices[["sseqid"]])
  }

  #tmp_in_data <<- in_data

  in1_g <- data.frame(from=in_data[[1]][,col.indices[["qseqid"]]], to=in_data[[1]][,col.indices[["sseqid"]]], stringsAsFactors = T)
  in1_g <- in1_g[apply(in1_g, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),]

  in2_g <- data.frame(from=in_data[[2]][,col.indices[["sseqid"]]], to=in_data[[2]][,col.indices[["qseqid"]]], stringsAsFactors = T)
  in2_g <- in2_g[apply(in2_g, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),]

  #print(in1_g) #DEBUG
  #print(in2_g) #DEBUG

  #tmp_in1_g <<- in1_g
  #tmp_in2_g <<- in2_g

  in1_valid_hits <- unique(intersect(in1_g$from,in2_g$to),intersect(in1_g$to,in2_g$from)) #purrr::reduce(list(in1_g$from,in2_g$to), intersect)
  in2_valid_hits <- unique(intersect(in2_g$from,in1_g$to),intersect(in2_g$to,in1_g$from)) #purrr::reduce(list(in2_g$from,in1_g$to), intersect)

  #print(head(in1_valid_hits)) #DEBUG
  #print(head(in2_valid_hits)) #DEBUG

  # in1_g <- in1_g[which(!is.na(match(in1_g$from,in1_valid_hits))),]
  # in2_g <- in2_g[which(!is.na(match(in2_g$from,in2_valid_hits))),]
  # in1_g <- in1_g[which(!is.na(match(in1_g$to,in2_valid_hits))),]
  # in2_g <- in2_g[which(!is.na(match(in2_g$to,in1_valid_hits))),]

  #MAKE ADJACENCY MATRIX for the graph
  if(length(in1_valid_hits) >  0 && length(in2_valid_hits) > 0){
    ##adj_mat <- matrix(rep(list(), (length(in1_valid_hits) + length(in2_valid_hits)) ^ 2 ), byrow=TRUE, ncol=length(in1_valid_hits) + length(in2_valid_hits))
    adj_mat <- array(rep(list(), (length(in1_valid_hits) + length(in2_valid_hits)) ^ 2 ), dim = c(length(in1_valid_hits) + length(in2_valid_hits),length(in1_valid_hits) + length(in2_valid_hits)), dimnames = list(c(in1_valid_hits,in2_valid_hits),c(in1_valid_hits,in2_valid_hits)))
    #adj_mat <- 1 - diag(length(in1_valid_hits) + length(in2_valid_hits))
    #adj_mat[adj_mat == 1] <- list()
    #dimnames(adj_mat) <- list(c(in1_valid_hits,in2_valid_hits),c(in1_valid_hits,in2_valid_hits))
  }else{
    stop("No hits are reciprocal")
  }
  #rownames(adj_mat) <- c(in1_valid_hits,in2_valid_hits)
  #colnames(adj_mat) <- c(in1_valid_hits,in2_valid_hits)

  #adj_mat[in1_valid_hits,in1_valid_hits] <- 0
  #adj_mat[in2_valid_hits,in2_valid_hits] <- 0

  #print(c(in1_valid_hits,in2_valid_hits))
  #print(adj_mat)

  if(!is.null(col.indices[["weight.col"]])){
    hit_combinations <- unique(tidyr::crossing(in1_valid_hits,in2_valid_hits))
    weight_mats <- parallel::mclapply(col.indices[["weight.col"]], function(idx){
      weight_mat <- adj_mat
      weight_list <- dplyr::bind_rows( unlist(
        purrr::map2(hit_combinations$in1_valid_hits,hit_combinations$in2_valid_hits, function(x,y){
          #parallel::mclapply(in1_valid_hits, function(x){
          #return( parallel::mclapply(in2_valid_hits, function(y){
          if(x!=y){
            df1_rows <- intersect(which(!is.na(match(in_data[[1]][,col.indices[["qseqid"]]], x))), which(!is.na(match(in_data[[1]][,col.indices[["sseqid"]]], y))))
            df2_rows <- intersect(which(!is.na(match(in_data[[2]][,col.indices[["sseqid"]]], y))), which(!is.na(match(in_data[[2]][,col.indices[["qseqid"]]], x))))
            #print(df1_rows)
            #print(df2_rows)
            in1_weights <- in_data[[1]][df1_rows,idx]
            in2_weights <- in_data[[2]][df2_rows,idx]
            if(unique.hit.weights){
              in1_weights <- unique(in1_weights)
              in2_weights <- unique(in2_weights)
            }
            if(!is.null(process.weights.func)){
              in1_weights <- process.weights.func(in1_weights)
              in2_weights <- process.weights.func(in2_weights)
            }

            #weight_mat[x,y] <- list(in1_weights)
            #weight_mat[y,x] <- list(in2_weights)
            return(list(row=x,col=y,weight=list(in1_weights)))
          }
          #}, mc.cores = n_threads,mc.silent = T) )
          #},mc.cores = n_threads,mc.silent = T), recursive = F,use.names = T) )
          #print(head(weight_list))
          #return(weight_list)
        }  ) ) )
      weight_mat[weight_list$row,weight_list$col] <- weight_list$weight
      weight_mat[weight_list$col,weight_list$row] <- weight_list$weight
      return(weight_mat)
    } ,mc.cores = n_threads,mc.silent = T)

    #tmp_weight_mats <<- weight_mats

    RBH_list <- parallel::mclapply(weight_mats, function(weight_mat){
      RBH_hits <- parallel::mclapply(c(in1_valid_hits,in2_valid_hits), function(q_hit){
        #print(q_hit) #DEBUG
        s_hit <- names(which.max(unlist(weight_mat[q_hit,], recursive = T,use.names = T)))
        back_q_hit <- names(which.max(unlist(weight_mat[s_hit,], recursive = T,use.names = T)))
        if(all(stringi::stri_cmp_eq(back_q_hit,q_hit))){
          return(data.frame(q_hit=q_hit,s_hit=s_hit))
        }
      },mc.cores = n_threads,mc.silent = T)
      return(dplyr::bind_rows(RBH_hits))
    } ,mc.cores = n_threads,mc.silent = T)

    #tmp_RBH_list <<- RBH_list

    if(length(weight_mats) > 1){
      if(all(unlist(purrr::reduce(RBH_list,all_equal, ignore_row_order = TRUE)))){
        RBH_final_list <- RBH_list
      }else{
        stop("multiple unidentical weight_mats case. take only identical rows. write code") ##use identical()
      }
    }else{
      RBH_final_list <- RBH_list
    }

    RBH_final_df <- unique(dplyr::bind_rows(RBH_final_list))

    #which(!is.na(match(in_data[[1]][,col.indices[["qseqid"]]], RBH_final_df$q_hit)))
    in1_RBH_rows <- which(!is.na( match(in_data[[1]][,col.indices[["sseqid"]]], RBH_final_df$s_hit) )) #c(RBH_final_df$q_hit, RBH_final_df$s_hit) #unique( which(!is.na( match(in_data[[1]][,col.indices[["qseqid"]]], c(RBH_final_df$q_hit, RBH_final_df$s_hit)) )), which(!is.na( match(in_data[[1]][,col.indices[["sseqid"]]], c(RBH_final_df$q_hit, RBH_final_df$s_hit)) )) )
    in2_RBH_rows <- which(!is.na( match(in_data[[2]][,col.indices[["qseqid"]]], RBH_final_df$q_hit) )) #c(RBH_final_df$q_hit, RBH_final_df$s_hit) #unique( which(!is.na( match(in_data[[2]][,col.indices[["qseqid"]]], c(RBH_final_df$q_hit, RBH_final_df$s_hit)) )), which(!is.na( match(in_data[[2]][,col.indices[["sseqid"]]], c(RBH_final_df$q_hit, RBH_final_df$s_hit)) )) )
    #print(RBH_final_df)
    #print(in1_RBH_rows)
    #print(in2_RBH_rows)
  }else{
    in1_RBH_rows <- which(!is.na( match(in_data[[1]][,col.indices[["sseqid"]]], c(in1_valid_hits,in2_valid_hits) ) ))
    in2_RBH_rows <- which(!is.na( match(in_data[[2]][,col.indices[["qseqid"]]], c(in1_valid_hits,in2_valid_hits) ) ))
  }
  if(index.tables){
    in_data[[1]] <- deindex_BLAST_table(in_data[[1]], col.indices[["qseqid"]])
    in_data[[1]] <- deindex_BLAST_table(in_data[[1]], col.indices[["sseqid"]])
    in_data[[2]] <- deindex_BLAST_table(in_data[[2]], col.indices[["qseqid"]])
    in_data[[2]] <- deindex_BLAST_table(in_data[[2]], col.indices[["sseqid"]])
  }

  return(list(in1=in_data[[1]][in1_RBH_rows,],in2=in_data[[2]][in2_RBH_rows,]))
}

