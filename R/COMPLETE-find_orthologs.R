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
#' @return Data Frame with BLAST Results
#' @export
LoadBLASTHits <- function(infile, transcript_ID_metadata=NULL, col.names=NULL){

  if(!is.null(transcript_ID_metadata) && is.character(transcript_ID_metadata)){
    transcript_ID_metadata <- read.table(file = transcript_ID_metadata,header = F,sep="\t",quote = "")
  }
  if(file.exists(infile) && file.info(infile)$size > 0){
    blast_results <- read.table(file = infile,header = F,sep="\t",quote = "", blank.lines.skip = T)
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

#' Index BLAST Table Column
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
#' @param blast_table BLAST table
#' @param index_col Column index of BLAST table to index (shorten IDs)
#' @param offset Offset value to add to the indices
#' @return BLAST table with indexed IDs from index_col (original index_col is attached to the table)
#' @export
index_BLAST_table <- function(blast_table, index_col, offset=0){
  old_name = colnames(blast_table)[index_col]
  blast_table <- dplyr::mutate(blast_table, blast_table[,index_col])

  seq_ids <- levels(factor(blast_table[,index_col]))
  indexed_ids <- data.frame(index=paste("i",seq(1:length(seq_ids))+offset,sep=""))
  rownames(indexed_ids) <- seq_ids
  #blast_table[,index_col] <- paste("i",seq(1:length(seq_ids))+offset,sep="")
  blast_table[,index_col] <- sapply(blast_table[,index_col], function(x){
    return(indexed_ids[x,"index"])
  })

  colnames(blast_table)[index_col] <- paste("indexed",old_name,sep="_")
  colnames(blast_table)[ncol(blast_table)] <- old_name
  #print(blast_table)
  return(blast_table)
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
  new_name = stringi::stri_split_fixed(old_name, pattern = "_", n=2,simplify = T)[2]

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
#'     GRObject_from_BLAST(blast_input = in_file, COMPLETE.format.ids = F, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list=NULL)
#'
#' @param blast_input Filename with BLAST hits or a BLAST table
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN)
#' @param col.indices A Named Vector with indices of columns Query sequence ID (qseqid), Subject sequence ID (sseqid), E Value (evalue), Query start (qstart), Query end (qend), Subject start (sstart), Subject end (send), Bitscore (bitscore), Query HSP coverage(qcovhsp), Query length (qlen), Subject length(slen), Frames (frames), Percentage Identity (pident), Number of Gaps (gaps), Alignment Length (length), Subject strand (sstrand). Default values are col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17)
#' @param params_list Output of load_params() (Optional - Only give it when using Files generated by R-COMPLETE and if COMPLETE.format.ids=TRUE)
#' @return A GRanges object of BLAST hits
#' @export
GRObject_from_BLAST <- function(blast_input, COMPLETE.format.ids=F, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list=NULL){

  if(!is.null(blast_input) && is.character(blast_input)){
    blast_input <- LoadBLASTHits(blast_input,col.names = col.names)
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
    blast_input$subject_gene <- unlist(purrr::map(blast_input[,col.indices["sseqid"]],function(x){
      genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      return(genes[2])
    }))

    blast_input$query_gene <- unlist(purrr::map(blast_input[,col.indices["qseqid"]],function(x){
      genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      return(genes[2])
    }))

    blast_input <- blast_input[which(blast_input[,col.indices["evalue"]] < params_list$E_VALUE_THRESH),]

  }#else{
  #  stop("Parameter file not loaded with load_params(), params_list is NULL & COMPLETE.format.ids==FALSE")
  #}

  ##CHANGE strands
  change_qstrand <- which(blast_input[,col.indices["qstart"]] > blast_input[,col.indices["qend"]])
  if (length(change_qstrand) > 0) {
    # print("changing qstart")
    tmp <- blast_input[change_qstrand,col.indices["qstart"]]
    blast_input[change_qstrand,col.indices["qstart"]] <- blast_input[change_qstrand,col.indices["qend"]]
    blast_input[change_qstrand,col.indices["qend"]] <- tmp
    tmp <- NULL
  }

  change_sstrand <- which(blast_input[,col.indices["sstart"]] > blast_input[,col.indices["send"]])
  if (length(change_sstrand) > 0) {
    # print("changing sstart")
    tmp <- blast_input[change_sstrand,col.indices["sstart"]]
    blast_input[change_sstrand,col.indices["sstart"]] <- blast_input[change_sstrand,col.indices["send"]]
    blast_input[change_sstrand,col.indices["send"]] <- tmp
    tmp <- NULL
  }

  no_strand_info <- which(stringi::stri_cmp_eq("N/A",blast_input[,col.indices["sstrand"]]))
  blast_input[no_strand_info,col.indices["sstrand"]] <- "*"

  #tmp_gr <- GenomicRanges::GRanges(S4Vectors::Rle(blast_input[,col.indices["sseqid"]]), ranges =  IRanges::IRanges(blast_input[,col.indices["sstart"]], end = blast_input[,col.indices["send"]],strand= blast_input[,col.indices["sstrand"]])) #, names = orths$sseqid))

  tmp_df <- data.frame(seq_names=blast_input[,col.indices["sseqid"]])

  #tmp_df <- mutate(tmp_df, seqnames=S4Vectors::Rle(blast_input[,col.indices["sseqid"]]))
  #tmp_df <- mutate(tmp_df, ranges=IRanges::IRanges(start = blast_input[,col.indices["sstart"]], end = blast_input[,col.indices["send"]],strand= blast_input[,col.indices["sstrand"]]))
  tmp_df <- mutate(tmp_df, sstart = blast_input[,col.indices["sstart"]])
  tmp_df <- mutate(tmp_df, send = blast_input[,col.indices["send"]])
  tmp_df <- mutate(tmp_df, sstrand= blast_input[,col.indices["sstrand"]])
  tmp_df <- mutate(tmp_df, Hsp_num=c(1:nrow(blast_input))) #seq(1,nrow(tmp_blast),1)
  tmp_df <- mutate(tmp_df, Hsp_bit.score = blast_input[,col.indices["bitscore"]])
  tmp_df <- mutate(tmp_df, Hsp_score = blast_input[,col.indices["qcovhsp"]])
  tmp_df <- mutate(tmp_df, Hsp_evalue = blast_input[,col.indices["evalue"]])
  tmp_df <- mutate(tmp_df, Hsp_query.from = blast_input[,col.indices["qstart"]])
  tmp_df <- mutate(tmp_df, Hsp_query.to = blast_input[,col.indices["qend"]])
  tmp_df <- mutate(tmp_df, query_id = blast_input[,col.indices["qseqid"]])
  tmp_df <- mutate(tmp_df, query_len = blast_input[,col.indices["qlen"]])
  tmp_df <- mutate(tmp_df, subject_len = blast_input[,col.indices["slen"]])
  tmp_df <- mutate(tmp_df, Hsp_hit.from = blast_input[,col.indices["sstart"]])
  tmp_df <- mutate(tmp_df, Hsp_hit.to = blast_input[,col.indices["send"]])
  tmp_df <- mutate(tmp_df, Hsp_query.frame =  unlist(purrr::map(blast_input[,col.indices["frames"]],function(x){
    frames <- as.integer(unlist(stringi::stri_split_fixed(x,pattern = "/")))
    return(frames[1])
  })) )
  tmp_df <- mutate(tmp_df, Hsp_hit.frame =  unlist(purrr::map(blast_input[,col.indices["frames"]],function(x){
    frames <- as.integer(unlist(stringi::stri_split_fixed(x,pattern = "/")))
    return(frames[2])
  })) )
  if(COMPLETE.format.ids && !is.null(params_list)){
    tmp_df <- mutate(tmp_df, query_org = unlist(purrr::map(blast_input[,col.indices["qseqid"]],function(x){
      org <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      return(org[3])
    })) )
    tmp_df <- mutate(tmp_df, subject_org = unlist(purrr::map(blast_input[,col.indices["sseqid"]],function(x){
      org <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      return(org[3])
    })) )
  }else{
    warning("Parameter file not loaded with load_params & COMPLETE.format.ids==FALSE")
  }

  tmp_df <- mutate(tmp_df, Hsp_pidentity = blast_input[,col.indices["pident"]])
  #gr$Hsp_positive <-
  tmp_df <- mutate(tmp_df, Hsp_gaps = blast_input[,col.indices["gaps"]])
  tmp_df <- mutate(tmp_df, Hsp_align.len =  blast_input[,col.indices["length"]])
  if(COMPLETE.format.ids){
    tmp_df <- mutate(tmp_df, subject_gene =  blast_input[,c("subject_gene")])
    tmp_df <- mutate(tmp_df, query_gene =  blast_input[,c("query_gene")])
  }

  tmp_gr <- GenomicRanges::makeGRangesFromDataFrame(df = tmp_df,seqnames.field = "seq_names", start.field = "sstart",end.field = "send",strand.field = "sstrand" ,keep.extra.columns = T,ignore.strand = F)

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
#' @return A GRanges object of BLAST hits
#' @export
run_WISARD <- function(blast_hits, score_col, COMPLETE.format.ids=F,params_list=NULL){

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
      result <- wisard::get_wis(gr[unique(intersect(which(gr$query_id == query_id),which(S4Vectors::runValue(GenomicRanges::seqnames(gr)) == sub_id))),],max_score = score_col,overlap = 0)
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
      #print(child)
      g_name <- NULL
      g_name <- unique(unlist(stri_split_fixed(child$alignments$subject_gene,pattern=params_list$SEQUENCE_ID_DELIM,n = 1,tokens_only = T)))
      # print(g_name)
      #all_results<- c(all_results, child)
      if(!is.null(g_name)){
        if(is.null(all_results)){
          all_results[[g_name]] <- child
        }else{
          all_results[[g_name]] <- merge_wis(all_results[[g_name]], child)
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

  blast_table <-  mutate(blast_table,from_org=unlist(lapply(stringi::stri_split_fixed(blast_table[,1],delimiter,n=4,tokens_only = T),function(x){return(x[3])})))
  blast_table <-  mutate(blast_table,to_org=unlist(lapply(stringi::stri_split_fixed(blast_table[,2],delimiter,n=4,tokens_only = T),function(x){return(x[3])})))

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
#' @param blast_DB_dir Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed
#' @param blast_out Path to BLAST output file, Default BLAST FORMAT is 11. It is converted internally to BLAST Format 6 and returned as a GRanges Object
#' @param blast_program Give path to the BLAST program. eg, Sys.which("tblastx") if tblastx is in SHELL $PATH
#' @param run_name Name of the BLAST run. Only for logging (Optional)
#' @param blast_options Extra Options to be passed to the BLAST program
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE (Default) otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (Optional)
#' @param params_list Output from load_params() (Optional)
#' @param keep.output.files TRUE(Default)/FALSE - Keep Output and BLAST DB Files? (Optional)
#' @return BLAST Hits as GRanges Object
#' @export
run_BLAST <- function(query_path, subject_path,blast_DB_dir = NULL, blast_out, blast_program=COMPLETE$BLAST_BIN, run_name="BLAST",blast_options="", COMPLETE.format.ids=F,params_list=NULL, keep.output.files=T){

  tryCatch({
  if (stringi::stri_isempty(COMPLETE$parallel)) {
    stop("Problem with GNU parallel installation. Reload R-COMPLETE")
  }

    if(stringi::stri_isempty(blast_program)){
      stop("BLAST+ not found in $PATH. Provide blast_program")
    }else{
      BLAST_BIN <- dirname(blast_program)
    }

  subject_path <- tools::file_path_as_absolute(subject_path)
  if (!is.null(blast_DB_dir)) {
    dir.create(blast_DB_dir,showWarnings = F,recursive = T)
    #query_DB <- tools::file_path_as_absolute(paste(blast_DB_dir,"/",basename(query_path),sep=""))
    subject_DB <- tools::file_path_as_absolute(paste(blast_DB_dir,"/",basename(subject_path),sep=""))
    #if(!stringi::stri_cmp_eq(query_path,query_DB)){
    #  file.copy(query_path,query_DB,overwrite = T)
    #}
    if(!stringi::stri_cmp_eq(subject_path,subject_DB)){
      file.copy(subject_path,subject_DB,overwrite = T)
    }
    #query_path <- query_DB
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

  #MAKE BLAST DB of FASTA files
  #Only subject fasta files needs to be a BLAST DB
  #processx::run( command = COMPLETE$SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"make_BLAST_db",query_path, dirname(blast_program)) ,spinner = T,stdout = "",stderr = "")
  processx::run( command = COMPLETE$SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"make_BLAST_db",subject_path, BLAST_BIN ) ,spinner = T,stdout = "",stderr = "")

  processx::run( command = COMPLETE$SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"do_BLAST",COMPLETE$parallel,run_name,query_path,subject_path,blast_out,blast_program,blast_options) ,spinner = T,stdout = "",stderr = "")

  convert_BLAST_format(infile = blast_out,outfile = final_blast_out,conversion_prg = tools::file_path_as_absolute(paste(BLAST_BIN,"/blast_formatter",sep="")) )
  blast_GR <- GRObject_from_BLAST(blast_input = final_blast_out,COMPLETE.format.ids = COMPLETE.format.ids,col.indices = c(qseqid = 1, sseqid = 2, evalue = 11, qstart = 7, qend = 8, sstart = 9, send = 10, bitscore = 12, qcovhsp = 16, qlen = 18, slen = 19, frames = 15, pident = 3, gaps = 14, length = 4, sstrand = 17), params_list = params_list)

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
    return(cond)
  })
}

#' Execute all2all BLAST
#'
#' Executes All-to-All BLAST between two lists of organisms/genes/clusters. Output BLAST files are stored in the format filename1.filename2.all2all under output_dir.
#'
#' @examples
#'  all2all_BLAST(first_list = list.files(path="fasta",pattern = "\*.fa",all.files = T,full.names = T,include.dirs = F,recursive = F), second_list = list.files(path="fasta",pattern = "\*.fa",all.files = T,full.names = T,include.dirs = F,recursive = F),blast_DB_dir = "files/blastdb",blast_program = "tblastx",output_dir = "files/all2all",blast_options = "-strand plus",input_prefix_path = NULL,params_list=NULL)
#'
#' @note ASSUMES Nucleotide sequences. Not checking/Not working for Protein/Peptide sequences.
#'
#' @param first_list Vector of PATHS to FASTA
#' @param second_list Vector of PATHS to FASTA
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program
#' @param blast_DB_dir Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed
#' @param blast_options Extra Options to be passed to the BLAST program
#' @param output_dir Path to BLAST output
#' @param input_prefix_path If input lists/vectors are filenames, then provide input folder to prefix path
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE (Default) otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (Optional)
#' @param keep.output.files TRUE(Default)/FALSE - Keep Output and BLAST DB Files? (Optional)
#' @param params_list Output of load_params() (Optional)
#' @export
all2all_BLAST <- function(first_list,second_list,blast_DB_dir=NULL,blast_program,output_dir="./", blast_options="", input_prefix_path=NULL, params_list=NULL,COMPLETE.format.ids=F, keep.output.files=T){

  if(is.null(params_list)){
    tryCatch(numWorkers <- parallel::detectCores(all.tests = T, logical = T), error=function(){numWorkers <- 2})
  }else{
    #warning("Parameter file not loaded with load_params || COMPLETE.format.ids==FALSE")
    numWorkers <- params_list$numWorkers
  }

  dir.create(path = output_dir,recursive = T,showWarnings = F)

  mclapply(first_list, function(first_set){
    mclapply(second_list,function(second_set){
      if(!is.null(input_prefix_path)){
        first_set <- paste(input_prefix_path,"/",first_set,sep="")
        second_set <- paste(input_prefix_path,"/",second_set,sep="")
      }
      if (file.exists(first_set) && file.exists(second_set) && file.info(first_set)$size > 0 && file.info(second_set)$size > 0) {
        run_name1 <- tools::file_path_as_absolute(tools::file_path_sans_ext(BiocGenerics::basename(first_set)))
        run_name2 <- tools::file_path_as_absolute(tools::file_path_sans_ext(BiocGenerics::basename(second_set)))
        run_name <- paste( run_name1,run_name2,"all2all" ,sep=".")
        out_file <- paste( output_dir, run_name,sep="")
        run_BLAST(query_path = first_set,subject_path = second_set,blast_DB_dir = blast_DB_dir, blast_program=blast_program, blast_out = out_file, run_name = run_name,blast_options = blast_options,COMPLETE.format.ids = COMPLETE.format.ids,params_list = params_list,keep.output.files = keep.output.files)
      }
    }, mc.cores = floor(sqrt(numWorkers)) )
  }, mc.cores = floor(sqrt(numWorkers)) )

}

#' Calculate HSP Coverage
#'
#' This function calculates the coverage of HSPs (sum(HSP Alignment lengths)/mRNA CDS Length) between each Query->Subject Hits. Hits are filtered based on the minimum coverage filter value.
#'
#' @param blast_table BLAST Table with Query->Subject Hits
#' @param col.indices A Named List with indices of columns Query sequence ID (qseqid) and Subject sequence ID (sseqid). Eg col.indices=list(qseqid=1,sseqid=2)
#' @param min_coverage_filter Minimum HSP Coverage value to filter out Hits
#' @param params_list Output of load_params()
#' @export
calculate_HSP_coverage <- function(blast_table,col.indices,min_coverage_filter=0.5,params_list){
  purrr::map2(blast_table[,col.indices["qseqid"]], blast_table[,col.indices["sseqid"]], function(x,y){

  })
}

#' Internal Function - Transcript Ortholog Extraction Function for R-COMPLETE pipeline
#'
#' This function calls the Transcript Ortholog Extraction pipeline which is used to reduce the pool of genes (step 1), reduce the pool of organisms and create sets of organisms (step 2), find transcript level orthologs (step 3). It takes only one argument which is the path/name of the BLAST program to use and refers to the values from the parameters file for other variables.
#'
#'  * Step 1 - Genes which are available in all the reference organisms are chosen
#'  * Step 2 - A Per-Gene Conservation Score (GSC) is calculated from the availability of a gene across organisms (literally the count of organisms which have the gene, normalized to 1 relative to other genes). Genes which have GSC score below GENE_DROP_THRESHOLD (parameter) are dropped (GENE_DROP_THRESHOLD=0 does not omit any genes). Sets of organisms are created based on the available genes after GSC filtering. I can suggest reference organisms based on which ones have the maximum number of genes
#'  * Step 3 - Two way BLAST followed by HSP selection with WISARD and Two way RBH are performed. Only transcripts which are bi-directional best hits are kept for further analysis (RBH from both the directions, not RBH in itself is bi-directionaly from the point of the QUERY, We can also do an RBH from the context of the SUBJECT to verify if it did not pas RBH by chance (even though it is very unlikely))
#'
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program
#' @param params_list Output of load_params()
transcript_ortholog_extraction <- function(blast_program, params_list){

  install_parallel()

  blast_options <- params_list$BLAST_OPTIONS
  blast_DB_dir <- params_list$BLAST_DB_PATH
  REF_ORGS <- factor(scan(params_list$REF_ORGS_FILE, character(), quiet = T))
  all2all_out <- paste(loaded_PARAMS$OUT_PATH,"/all2all/",sep = "")

  #STEP 1
  #select genes which are available in all the reference organisms
  available_genes_list <- parallel::mclapply(paste(params_list$OUT_PATH,"/genes/",REF_ORGS,sep=""),function(x){
    if(file.exists(paste(x,"/AVAILABLE_GENES",sep="")) && file.info(paste(x,"/AVAILABLE_GENES",sep=""))$size > 0 ){
      return(scan(paste(x,"/AVAILABLE_GENES",sep=""), character()))
    }
  }, mc.cores =  params_list$numWorkers)
  available_genes <- unique(purrr::reduce(available_genes_list, union))

  #find which clusters they belong to
  odb_map_list <- parallel::mclapply(paste(params_list$OUT_PATH,"/genes/",REF_ORGS,sep=""),function(x){
    if(file.exists(paste(x,"/odb.final_map",sep="")) && file.info(paste(x,"/odb.final_map",sep=""))$size > 0 ){
      return(read.table(paste(x,"/odb.final_map",sep=""),header = F,sep = "\t",quote = "", col.names = c("cluster","genes")))
    }
  }, mc.cores =  params_list$numWorkers)
  odb_map <- unique(purrr::reduce(odb_map_list, inner_join, by = c("cluster","genes")))
  available_clusters <- unique(BiocGenerics::unlist(lapply(available_genes, function(gene){
    return(odb_map[grep(pattern=gene,x = odb_map$genes,ignore.case = T), c("cluster")])
  })))

  if(length(dir(params_list$GROUPS_PATH)==0)){
    group_FASTA_clusters(params_list$FASTA_OUT_PATH)
  }

  ##Place ungrouped sequences into groups (all2allblast BLAST ungrouped cluster againts all clusters)
  #all2allblast and then wisard and then RBH for grouping ungrouped clusters
  all2all_BLAST(first_list = "ungrouped", second_list = grep("ungrouped",available_clusters,ignore.case = T,invert = T,value=T),blast_DB_dir = blast_DB_dir,blast_program = blast_program,output_dir =all2all_out,blast_options = blast_options,input_prefix_path = params_list$GROUPS_PATH, params_list = params_list, COMPLETE.format.ids = T, keep.output.files = T ) #paste(params_list$TEMP_PATH,"/","all2all/",sep="")
  all2all_BLAST(first_list = grep("ungrouped",available_clusters,ignore.case = T,invert = T,value=T), second_list = "ungrouped",blast_DB_dir = blast_DB_dir,blast_program = blast_program,output_dir = all2all_out,blast_options = blast_options,input_prefix_path = params_list$GROUPS_PAT, params_list = params_list, COMPLETE.format.ids = T, keep.output.files = T ) #paste(params_list$TEMP_PATH,"/","all2all/",sep="")

  mclapply(list.files(path = all2all_out,pattern = "*.all2all", ignore.case = T,full.names = T),function(in_file){
    out_file <- paste(all2all_out,tools::file_path_sans_ext(BiocGenerics::basename(in_file)),".out",sep="")
    convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
  }, mc.cores = params_list$numWorkers )

   mclapply(list.files(path = all2all_out,pattern = "*.out", ignore.case = T,full.names = T),function(in_file){
    out_file <- paste(all2all_out,tools::file_path_sans_ext(BiocGenerics::basename(in_file)),".wis_out",sep="")
    blast_GO <- GRObject_from_BLAST(blast_input = in_file, COMPLETE.format.ids = T, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list = params_list)
    wis_GO <- run_WISARD(blast_hits = blast_GO,score_col = "Hsp_score",COMPLETE.format.ids = T, params_list = params_list) #score_col=16
    wis_GO <- melt_wisard_list(wis_GO)
    write.table(x = wis_GO,file = out_file,quote = F,col.names = T,row.names = F)
  }, mc.cores = params_list$numWorkers )

  #save("wisard_results", file="files/all2all/wisard_results.RData")
  #load("files/all2all/wisard_results.RData")

   ##RUN RBH
   lapply("ungrouped", function(query){
     mclapply(grep("ungrouped",available_clusters,ignore.case = T,invert = T,value=T),function(subject){
       in1 <- paste(all2all_out,query,"-",subject,".wis_out",sep="")
       in2 <- paste(all2all_out,subject,"-",query,".wis_out",sep="")
       if (file.exists(in1) && file.exists(in2) && file.info(in1)$size > 0 && file.info(in2)$size > 0 ) {
          #RBH(in1 = in1, in2 = in2, index.tables = T, col.indices = c(qseqid=12,sseqid=1,weight.col=22),col.names = c("subject_id","start","end","width","strand","Hsp_num","Hsp_bit.score","Hsp_score","Hsp_evalue","Hsp_query.from","Hsp_query.to","query_id","query_len","subject_len","Hsp_hit.from","Hsp_hit.to","Hsp_query.frame","Hsp_hit.frame","Hsp_pidentity","Hsp_gaps","Hsp_align.len","max_score"))
         RBH_out <- RBH(in1 = in1, in2 = in2, index.tables = T, col.indices = list(qseqid=12,sseqid=1,weight.col=c(22,8) ), unique.hit.weights = T, process.weights.func = max)
         out1 <- paste(all2all_out,query,"-",subject,".rbh_out",sep="")
         out2 <- paste(all2all_out,subject,"-",query,".rbh_out",sep="")
         write.table(x = RBH_out$in1,file = out1,quote = F,col.names = T,row.names = F)
         write.table(x = RBH_out$in2,file = out2,quote = F,col.names = T,row.names = F)
       }
     }, mc.cores =numWorkers )
   })

  #calculate gene conservation - calculate_gene_conservation.R - probably not needed
   ##Maybe write one for cluster conservation/coverage across organisms



  #create organism sets

  #calculate cluster occupancy - number of genes per cluster && number of organisms per cluster

  ##ITERATION 2 - two way RBH

  #cat $reference_ORGS files/oneway/SET > files/oneway/set.tmp

  #time all2all_refblast $reference_ORGS $fasta_path files/all2all/all2all.genelist files/all2all_final $blastdb_path $region $region tblastx files/oneway/set.tmp

  #time Rscript wisard.R files/oneway/set.tmp files/all2all_final files/gtf_stats.csv

  #readarray set_orgs < files/oneway/set.tmp

  #time parallel -j $((${#ref_orgs[@]}*${#set_orgs[@]})) "twoway_RBH files/all2all_final {1} {2} $PY3_PATH $RBH_SCRIPT $SAME_GENE" ::: ${ref_orgs[@]} ::: ${set_orgs[@]}

}



#' (2) - Find Transcript Orthologs
#'
#' This function can be executed after COMPLETE::EXTRACT_DATA() and is the continuation of R-COMPLETE pipeline
#'
#' This is the main function which calls all the other functions and performs and end-end execution of finding transcript level orthologs. It runs the iterative Transcript Ortholog Extraction pipeline which is used to reduce the pool of genes, reduce the pool of organisms, find transcript level orthologs (check ?transcript_ortholog_extraction)
#'
#' @note ONLY USE THIS FUNCTION WHEN RUNNING THE PIPELINE OF R-COMPLETE. Use other helper function to work with custom BLAST files not generated by this R package
#'
#' @param param_file Filename of a formatted parameter file (check the github repo for an example) or Output of load_params().
#' @param gene_list Vector or File with a list of genes to extract data for(check the github repo for an example).
#' @export
FIND_TRANSCRIPT_ORTHOLOGS <- function(param_file, gene_list){
  set.seed(123)

  if(!grepl(x=Sys.info()["sysname"],pattern="linux",ignore.case = T)){
    stop("Pipeline only supports Linux (and bash) :(")
  }

  if(is.character(params_list)){
    if(!file.exists(params_list) || file.info(params_list)$size < 0){
      stop("ERROR: Parameters file is missing and is required\n")
    }
    loaded_PARAMS <- load_params(params_list)
  }else{
    loaded_PARAMS <- params_list
  }
  print(loaded_PARAMS)

  #loaded_PARAMS$gene_list <- gene_list
  #loaded_PARAMS$genes <- factor(scan(gene_list, character())) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
  #loaded_PARAMS$genes <- genes[grep("gene",tolower(genes), invert = T, fixed = T)]

  print(paste("MAX PROCESSES:",loaded_PARAMS$numWorkers))

  if(CLEAN_EXTRACT){
    unlink(paste(loaded_PARAMS$OUT_PATH,"/oneway",sep = ""), recursive = T,force = T,expand = T)
    unlink(paste(loaded_PARAMS$OUT_PATH,"/all2all", sep=""), recursive = T,force = T,expand = T)
    unlink(paste(loaded_PARAMS$OUT_PATH,"/all2all_final", sep = ""), recursive = T,force = T,expand = T)
    unlink(paste(loaded_PARAMS$OUT_PATH,"/gene_thresholds.txt",sep=""), recursive = T,force = T,expand = T)
  }

  dir.create(paste(loaded_PARAMS$OUT_PATH,"/oneway",sep = ""),showWarnings = F, recursive = T)
  dir.create(paste(loaded_PARAMS$OUT_PATH,"/all2all",sep = ""),showWarnings = F, recursive = T)
  dir.create(paste(loaded_PARAMS$OUT_PATH,"/all2all_final",sep = ""),showWarnings = F, recursive = T)

  if (!stringi::stri_isempty(blast_program)) {
    transcript_ortholog_extraction(blast_program = Sys.which("tblastx"), params_list = loaded_PARAMS)
  }else{
    stop(paste(blast_program," NOT found. Is BLAST+ installed or in $PATH?"))
  }


}

#' Find Reciprocal Blast Hits (RBH)
#'
#' Find RBH between BLAST results of different organisms/genes/transcripts (FASTA/FASTQ). The BLAST results must be of the format 6 and can be converted from BLAST format 11 with convert_BLAST_format().
#' The command with the required column names are given below.
#'
#' convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'
#' Optional : You can provide the file/table with indexed Transcsript IDs. The format must be "file"[tab]"long_id"[tab]"index" (without a header). It can be generated with the  index_FASTA_IDs() (check ?index_FASTA_IDs or index_fastaIDs() in system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"))
#'
#' @note Both the input files are expected to have a header, the same number of columns with matching column order (and names). Give only indices for col.indices. Order of execution is unique.hit.weights followed by process.weights.func (if any/all these options are set), i.e Unique Weights are chosen for each hit (if unique.hit.weights=T) and then weights are processed using process.weights.func (if process.weights.func is set)
#'
#' @examples
#'    convert_BLAST_format(in_file1,outfile = out_file1,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'    convert_BLAST_format(in_file2,outfile = out_file2,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'    blast_GO1 <- GRObject_from_BLAST(blast_input = out_file1, COMPLETE.format.ids = T, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
#'    blast_GO2 <- GRObject_from_BLAST(blast_input = out_file2, COMPLETE.format.ids = T, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
#'    wis1 <- run_WISARD(blast_hits = blast_GO1,score_col = "Hsp_score",COMPLETE.format.ids = T) #score_col=16
#'    wis2 <- run_WISARD(blast_hits = blast_GO2,score_col = "Hsp_score",COMPLETE.format.ids = T) #score_col=16
#'    wis1 <- melt_wisard_list(wis1)
#'    wis2 <- melt_wisard_list(wis2)
#'    RBH(in1 = wis1, in2 = wis2, index.tables = T, col.indices = list(qseqid=12,sseqid=1,weight.col=c(8,22)),col.names = c("subject_id","start","end","width","strand","Hsp_num","Hsp_bit.score","Hsp_score","Hsp_evalue","Hsp_query.from","Hsp_query.to","query_id","query_len","subject_len","Hsp_hit.from","Hsp_hit.to","Hsp_query.frame","Hsp_hit.frame","Hsp_pidentity","Hsp_gaps","Hsp_align.len","max_score"))
#'
#' @param in1 Input (query->subject) BLAST/WISARD hits table/filename
#' @param in2 Input (query<-subject) BLAST/WISARD hits table/filename
#' @param transcript_ID_metadata Tab-delimited File with the filenames, indexed transcript IDs and the long transcript IDs.
#' @param col.names Columns names for the BLAST/WISARD tables/files
#' @param col.indices A Named List with indices of columns Query sequence ID (qseqid), Subject sequence ID (sseqid), and columns to be used as edge weights (eg, "Hsp_score","max_score" etc). Eg col.indices=list(qseqid=1,sseqid=2,weight.col=c(8,22))
#' @param unique.hit.weights Should only the unique Weights be taken for all Query->Subject Hits? (TRUE/FALSE (Default))
#' @param process.weights.func Pass a function name to process the weights (eg, sum/max/min etc) (Default - max)
#' @param index.tables Should the IDs in the tables be indexed? (TRUE (Default) if COMPLETE.format.ids/Long BLAST Sequence IDs are used)
#' @param n_threads Number of Threads (Optional)
#' @return Named List list(in1,in2) with the selected Hits which are RBHs between in1 and in2 data
#' @export
RBH <- function(in1,in2, transcript_ID_metadata=NULL, col.names=NULL, index.tables=T,col.indices, unique.hit.weights=F, process.weights.func=max, n_threads=tryCatch(parallel::detectCores(all.tests = T, logical = T), error=function(cond){return(2)})){

  #print(col.indices)

  if(!is.null(transcript_ID_metadata) && is.character(transcript_ID_metadata)){
    transcript_ID_metadata <- read.table(file = transcript_ID_metadata,header = F,sep="\t",quote = "")
  }

  if(!is.null(in1) && is.character(in1)){
    in1_data <- LoadBLASTHits(infile = in1, transcript_ID_metadata = transcript_ID_metadata, col.names = col.names)
  }else if(!is.null(in1)){
    in1_data <- in1
    if(!is.null(col.names)){
      colnames(in1_data) <- col.names
    }
  }
  if(!is.null(in2) && is.character(in2)){
    in2_data <- LoadBLASTHits(infile = in2, transcript_ID_metadata = transcript_ID_metadata, col.names = col.names)
  }else if(!is.null(in2)){
    in2_data <- in2
    if(!is.null(col.names)){
      colnames(in2_data) <- col.names
    }
  }

    if(index.tables){
      in1_data <- index_BLAST_table(in1_data,col.indices[["sseqid"]],offset = 0)
      in1_data <- index_BLAST_table(in1_data,col.indices[["qseqid"]],offset = length(levels(factor( in2_data[,col.indices[["qseqid"]]] ))) + length(levels(factor( in1_data[,col.indices[["sseqid"]]] ))))
      in2_data <- index_BLAST_table(in2_data,col.indices[["sseqid"]],offset = length(levels(factor( in1_data[,col.indices[["qseqid"]]] ))) + length(levels(factor( in2_data[,col.indices[["sseqid"]]] ))))
      in2_data <- index_BLAST_table(in2_data,col.indices[["qseqid"]],offset = length(levels(factor( in2_data[,col.indices[["qseqid"]]] ))) + length(levels(factor( in1_data[,col.indices[["sseqid"]]] ))))
    }

  # #Sanity checks for weight.col - Converted to numeric indices
  if(!is.list(col.indices)){
    stop("col.indices must be a Named List!")
  }
  if(!all(stringi::stri_cmp_eq(colnames(in1_data)[col.indices[["weight.col"]]],colnames(in2_data)[col.indices[["weight.col"]]]))){
     stop("Column indices of weight.col and order of columns must match between data!")
  }

  in1_g <- data.frame(from=in1_data[,col.indices[["qseqid"]]], to=in1_data[,col.indices[["sseqid"]]], stringsAsFactors = T)
  in1_g <- in1_g[apply(in1_g, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),]

  in2_g <- data.frame(from=in2_data[,col.indices[["sseqid"]]], to=in2_data[,col.indices[["qseqid"]]], stringsAsFactors = T)
  in2_g <- in2_g[apply(in2_g, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),]

  #print(in1_g) #DEBUG
  #print(in2_g) #DEBUG

  in1_valid_hits <- unique(intersect(in1_g$from,in2_g$to),intersect(in1_g$to,in2_g$from)) #purrr::reduce(list(in1_g$from,in2_g$to), intersect)
  in2_valid_hits <- unique(intersect(in2_g$from,in1_g$to),intersect(in2_g$to,in1_g$from)) #purrr::reduce(list(in2_g$from,in1_g$to), intersect)

  # in1_g <- in1_g[which(!is.na(match(in1_g$from,in1_valid_hits))),]
  # in2_g <- in2_g[which(!is.na(match(in2_g$from,in2_valid_hits))),]
  # in1_g <- in1_g[which(!is.na(match(in1_g$to,in2_valid_hits))),]
  # in2_g <- in2_g[which(!is.na(match(in2_g$to,in1_valid_hits))),]

  #MAKE ADJACENCY MATRIX for the graph
  if(length(in1_valid_hits) >  0 && length(in2_valid_hits) > 0){
    #adj_mat <- matrix(rep(list(), (length(in1_valid_hits) + length(in2_valid_hits)) ^ 2 ), byrow=TRUE, ncol=length(in1_valid_hits) + length(in2_valid_hits))
    adj_mat <- array(rep(list(), (length(in1_valid_hits) + length(in2_valid_hits)) ^ 2 ), dim = c(length(in1_valid_hits) + length(in2_valid_hits),length(in1_valid_hits) + length(in2_valid_hits)), dimnames = list(c(in1_valid_hits,in2_valid_hits),c(in1_valid_hits,in2_valid_hits)))
  }else{
    stop("No hits are reciprocal")
  }
  #rownames(adj_mat) <- c(in1_valid_hits,in2_valid_hits)
  #colnames(adj_mat) <- c(in1_valid_hits,in2_valid_hits)

  adj_mat[in1_valid_hits,in1_valid_hits] <- 0
  adj_mat[in2_valid_hits,in2_valid_hits] <- 0

  #print(c(in1_valid_hits,in2_valid_hits))
  #print(adj_mat)
  weight_mats <- unlist(parallel::mclapply(col.indices[["weight.col"]], function(idx){

    return(purrr::map2(in1_valid_hits,in2_valid_hits, function(x,y){
      weight_mat <- adj_mat
      df1_rows <- intersect(which(!is.na(match(in1_data[,col.indices[["qseqid"]]], x))), which(!is.na(match(in1_data[,col.indices[["sseqid"]]], y))))
      df2_rows <- intersect(which(!is.na(match(in2_data[,col.indices[["sseqid"]]], y))), which(!is.na(match(in2_data[,col.indices[["qseqid"]]], x))))
      #print(df1_rows)
      #print(df2_rows)
      in1_weights <- in1_data[df1_rows,idx]
      in2_weights <- in2_data[df2_rows,idx]
      if(unique.hit.weights){
        in1_weights <- unique(in1_weights)
        in2_weights <- unique(in2_weights)
      }
      if(!is.null(process.weights.func)){
        in1_weights <- process.weights.func(in1_weights)
        in2_weights <- process.weights.func(in2_weights)
      }

      weight_mat[x,y] <- list(in1_weights)
      weight_mat[y,x] <- list(in2_weights)
      return(weight_mat)
    }) )

  } ,mc.cores = n_threads,mc.silent = T) ,recursive = F,use.names = T)

  RBH_list <- parallel::mclapply(weight_mats, function(weight_mat){
    RBH_hits <- lapply(c(in1_valid_hits,in2_valid_hits), function(q_hit){
      print(q_hit)
      s_hit <- names(which.max(weight_mat[q_hit,]))
      back_q_hit <- names(which.max(weight_mat[s_hit,]))
      if(all(stringi::stri_cmp_eq(back_q_hit,q_hit))){
        return(data.frame(q_hit=q_hit,s_hit=s_hit))
      }
    })
    return(dplyr::bind_rows(RBH_hits))
  } ,mc.cores = n_threads,mc.silent = T)

  if(length(weight_mats) > 1){
    if(purrr::reduce(RBH_list,all_equal, ignore_row_order = TRUE)){
      RBH_final_list <- RBH_list
    }else{
        print("take only identical rows. write code") ##use identical()
    }
  }else{
    RBH_final_list <- RBH_list
  }

  RBH_final_df <- unique(dplyr::bind_rows(RBH_final_list))

  #which(!is.na(match(in1_data[,col.indices[["qseqid"]]], RBH_final_df$q_hit)))
  in1_RBH_rows <- which(!is.na( match(in1_data[,col.indices[["sseqid"]]], RBH_final_df$s_hit) )) #c(RBH_final_df$q_hit, RBH_final_df$s_hit) #unique( which(!is.na( match(in1_data[,col.indices[["qseqid"]]], c(RBH_final_df$q_hit, RBH_final_df$s_hit)) )), which(!is.na( match(in1_data[,col.indices[["sseqid"]]], c(RBH_final_df$q_hit, RBH_final_df$s_hit)) )) )
  in2_RBH_rows <- which(!is.na( match(in2_data[,col.indices[["qseqid"]]], RBH_final_df$q_hit) )) #c(RBH_final_df$q_hit, RBH_final_df$s_hit) #unique( which(!is.na( match(in2_data[,col.indices[["qseqid"]]], c(RBH_final_df$q_hit, RBH_final_df$s_hit)) )), which(!is.na( match(in2_data[,col.indices[["sseqid"]]], c(RBH_final_df$q_hit, RBH_final_df$s_hit)) )) )
  #print(RBH_final_df)
  #print(in1_RBH_rows)
  #print(in2_RBH_rows)

  if(index.tables){
    in1_data <- deindex_BLAST_table(in1_data, col.indices[["qseqid"]])
    in1_data <- deindex_BLAST_table(in1_data, col.indices[["sseqid"]])
    in2_data <- deindex_BLAST_table(in2_data, col.indices[["qseqid"]])
    in2_data <- deindex_BLAST_table(in2_data, col.indices[["sseqid"]])
  }

  return(list(in1=in1_data[in1_RBH_rows,],in2=in2_data[in2_RBH_rows,]))
}

