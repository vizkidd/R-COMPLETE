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
#' in1_data <- index_BLAST_table(in1_data,1,offset = 0)
#' in1_data <- index_BLAST_table(in1_data,2,offset = nrow(in1_data))
#'
#' @param blast_table BLAST table
#' @param index_col Column index of BLAST table to index (shorten IDs)
#' @param offset Offset value to add to the indices
#' @return BLAST table with indexed IDs from index_col (original index_col is attached to the table)
#' @export
index_BLAST_table <- function(blast_table, index_col, offset=0){
  old_name = colnames(blast_table)[index_col]
  blast_table <- dplyr::mutate(blast_table, blast_table[,index_col])

  blast_table[,index_col] <- paste("i",seq(1:nrow(blast_table))+offset,sep="")

  colnames(blast_table)[index_col] <- paste("indexed",old_name,sep="_")
  colnames(blast_table)[ncol(blast_table)] <- old_name
  return(blast_table)
}

#' Create A GRanges Object from BLAST RESULTS
#'
#' Functions accepts a filename or a BLAST table and converts it into a GenomicRanges::GRanges object (which can be used by R-WISARD). Give col.names when passing BLAST filenames. The accepted BLAST format is 6 and can be formatted from BLAST format 11 using convert_BLAST_format (shown below). Columns indices can provided with the col.indices option. length(col.names) and ncol(blast_input) are assumed to be the same. BLAST table/file must be non-indexed (IDs not shortened) if COMPLETE.format.ids=T
#'
#' @note Convert to compatible BLAST format with convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive")) (in_file must be of BLAST format 11)
#'
#' @examples
#'     GRObject_from_BLAST(blast_input = in_file,col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"), COMPLETE.format.ids = F, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list=NULL)
#'
#' @param blast_input Filename with BLAST hits or a BLAST table
#' @param col.names Columns which for the BLAST tables/files
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE otherwise)
#' @param col.indices A Named Vector with indices of columns Query sequence ID (qseqid), Subject sequence ID (sseqid), E Value (evalue), Query start (qstart), Query end (qend), Subject start (sstart), Subject end (send), Bitscore (bitscore), Query HSP coverage(qcovhsp), Query length (qlen), Subject length(slen), Frames (frames), Percentage Identity (pident), Number of Gaps (gaps), Alignment Length (length), Subject strand (sstrand). Default values are col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17)
#' @param params_list Output of load_params() (Optional - Only give it when using Files generated by R-COMPLETE and if COMPLETE.format.ids=TRUE)
#' @return A GRanges object of BLAST hits
#' @export
GRObject_from_BLAST <- function(blast_input, col.names=NULL, COMPLETE.format.ids=F, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list=NULL){

  if(!is.null(blast_input) && is.character(blast_input)){
    blast_input <- LoadBLASTHits(blast_input,col.names = col.names)
  }
  #tmp_blast <- blast_output
  #col_names <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand", "qlen", "slen", "qseq","sseq","nident","positive")
  #names(tmp_blast) <- col_names

  if(!is.null(col.names)){
    colnames(blast_input) <- col.names
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

  }else{
    stop("Parameter file not loaded with load_params(), params_list is NULL & COMPLETE.format.ids==FALSE")
  }

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
    tmp <- blast_input[,col.indices["sstart"]]
    blast_input[,col.indices["sstart"]] <- blast_input[,col.indices["send"]]
    blast_input[,col.indices["send"]] <- tmp
    tmp <- NULL
  }

  tmp_gr <- GenomicRanges::GRanges(S4Vectors::Rle(blast_input[,col.indices["sseqid"]]), ranges =  IRanges::IRanges(blast_input[,col.indices["sstart"]], end = blast_input[,col.indices["send"]],strand= blast_input[,col.indices["sstrand"]])) #, names = orths$sseqid))

  tmp_gr$Hsp_num <- c(1:nrow(blast_input)) #seq(1,nrow(tmp_blast),1)
  tmp_gr$Hsp_bit.score <- blast_input[,col.indices["bitscore"]]
  tmp_gr$Hsp_score <- blast_input[,col.indices["qcovhsp"]]
  tmp_gr$Hsp_evalue <- blast_input[,col.indices["evalue"]]
  tmp_gr$Hsp_query.from <- blast_input[,col.indices["qstart"]]
  tmp_gr$Hsp_query.to <- blast_input[,col.indices["qend"]]
  tmp_gr$query_id <- blast_input[,col.indices["qseqid"]]
  tmp_gr$query_len <- blast_input[,col.indices["qlen"]]
  tmp_gr$subject_len <- blast_input[,col.indices["slen"]]
  tmp_gr$Hsp_hit.from <- blast_input[,col.indices["sstart"]]
  tmp_gr$Hsp_hit.to <- blast_input[,col.indices["send"]]
  tmp_gr$Hsp_query.frame <-  unlist(purrr::map(blast_input[,col.indices["frames"]],function(x){
    frames <- as.integer(unlist(stringi::stri_split_fixed(x,pattern = "/")))
    return(frames[1])
  }))
  tmp_gr$Hsp_hit.frame <-  unlist(purrr::map(blast_input[,col.indices["frames"]],function(x){
    frames <- as.integer(unlist(stringi::stri_split_fixed(x,pattern = "/")))
    return(frames[2])
  }))
  if(COMPLETE.format.ids && !is.null(params_list)){
    tmp_gr$query_org <- unlist(purrr::map(blast_input[,col.indices["qseqid"]],function(x){
      org <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      return(org[3])
    }))
    tmp_gr$subject_org <- unlist(purrr::map(blast_input[,col.indices["sseqid"]],function(x){
      org <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
      return(org[3])
    }))
  }else{
    warning("Parameter file not loaded with load_params & COMPLETE.format.ids==FALSE")
  }

  tmp_gr$Hsp_pidentity <- blast_input[,col.indices["pident"]]
  #gr$Hsp_positive <-
  tmp_gr$Hsp_gaps <- blast_input[,col.indices["gaps"]]
  tmp_gr$Hsp_align.len <-  blast_input[,col.indices["length"]]
  if(COMPLETE.format.ids){
    tmp_gr$subject_gene <-  blast_input[,c("subject_gene")]
    tmp_gr$query_gene <-  blast_input[,c("query_gene")]
  }
  return(tmp_gr)
}

#' Select highest scoring interval of non-overlapping HSPs from Bi-Directional BLAST Hits
#'
#' Select the highest scoring pairs (HSPs) which give the maximum coverage over the BLAST alignments of each transcript (without overlaps/minimal overlaps). These HSPs will then be used to find Transcript level orthologs across gene orthologs across organisms. This function only accepts bi-directional (Query <-> Subject) BLAST Hits formatted with GRObject_from_BLAST()
#'
#' @examples
#'
#' blast_GO <- GRObject_from_BLAST(blast_input = in_file,col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"), COMPLETE.format.ids = T, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
#' wis_GO <- run_WISARD(blast_hits = blast_GO,score_col = "qcovhsp",COMPLETE.format.ids = T,params_list=NULL) #score_col=16
#'
#' @param blast_hits GRanges Object of BLAST Hits (Query -> Subject)
#' @param score_col Column in the GRanges Object used for scoring intervals in WISARD, Default would be qcovhsp (Query Coverage HSP)
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, Default - FALSE otherwise)
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

#' Internal Function - To Merge WISARD OUTPUTS
#'
#' @param x WIZARD output x
#' @param y WIZARD output y
#' @return A Merged WISARD Object
merge_wisard_table <- function(x, y){
  new_list <- list()
  new_list$alignments <- c(x$alignments, y$alignments)
  #new_list$alignments$max_score <- c(rep(x$max_score, length(x$alignments)),rep(y$max_score, length(y$alignments)))
  new_list$max_score <- sum(x$max_score, y$max_score)
  #print(new_list)
  return(unlist(new_list))
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
#' @export
convert_BLAST_format <- function(infile, outfile,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive")){
  if(stringi::stri_isempty(conversion_prg)){
    stop("blast_formatter not found in $PATH..cannot continue!")
  }

  processx::run( command = SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"convert_BLAST_format",infile,outfile,outformat,paste(cols,collapse = " ") ) ,spinner = T,stdout = "",stderr = "")

}

#' Perform Nucleotide BLAST
#'
#' BLAST between two organisms/genes/clusters.
#'
#' @note ASSUMES Nucleotide sequences. Not checking/Not working for Protein/Peptide sequences.
#'
#' @examples
#'     run_BLAST(query_path = "query.fasta",subject_path = "subject.fasta",blast_DB_dir = "files/blastdb", blast_program="tblastx", blast_out = "blast.out", run_name = "blast_positive",blast_options = "-strand plus")
#'
#' @param query_path Path to Query FASTA
#' @param subject_path Path to Subject FASTA
#' @param blast_DB_dir Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed
#' @param blast_out Path to BLAST output file, Default BLAST FORMAT is 11
#' @param run_name Name of the BLAST run
#' @param blast_options Extra Options to be passed to the BLAST program
#' @export
run_BLAST <- function(query_path, subject_path,blast_DB_dir = NULL, blast_out, blast_program, run_name="BLAST",blast_options=""){

  install_parallel()

  if (!is.null(blast_DB_dir)) {
    dir.create(blast_DB_dir,showWarnings = F,recursive = T)
    query_DB <- paste(blast_DB_dir,"/",query_path)
    subject_DB <- paste(blast_DB_dir,"/",subject_path)
    file.copy(query_path,query_DB,overwrite = T)
    file.copy(subject_path,subject_DB,overwrite = T)
    query_path <- query_DB
    subject_path <- subject_DB
  }

  #MAKE BLAST DB of FASTA files
  processx::run( command = SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"make_BLAST_db",query_path) ,spinner = T,stdout = "",stderr = "")
  processx::run( command = SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"make_BLAST_db",subject_path) ,spinner = T,stdout = "",stderr = "")

  processx::run( command = SHELL ,args=c(system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"),"do_BLAST",run_name,query_path,subject_path,blast_out,blast_program,blast_options) ,spinner = T,stdout = "",stderr = "")

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
#' @param params_list Output of load_params() (Optional)
#' @export
all2all_BLAST <- function(first_list,second_list,blast_DB_dir=NULL,blast_program,output_dir="./", blast_options="", input_prefix_path=NULL, params_list=NULL){

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
        run_BLAST(query_path = first_set,subject_path = second_set,blast_DB_dir = blast_DB_dir, blast_program=blast_program, blast_out = out_file, run_name = run_name,blast_options = blast_options)
      }
    }, mc.cores = floor(sqrt(numWorkers)) )
  }, mc.cores = floor(sqrt(numWorkers)) )

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

  #all2allblast and then wisard and then RBH for grouping ungrouped clusters
  all2all_BLAST(first_list = "ungrouped", second_list = grep("ungrouped",available_clusters,ignore.case = T,invert = T,value=T),blast_DB_dir = blast_DB_dir,blast_program = blast_program,output_dir =all2all_out,blast_options = blast_options,input_prefix_path = params_list$GROUPS_PATH, params_list = params_list ) #paste(params_list$TEMP_PATH,"/","all2all/",sep="")
  all2all_BLAST(first_list = grep("ungrouped",available_clusters,ignore.case = T,invert = T,value=T), second_list = "ungrouped",blast_DB_dir = blast_DB_dir,blast_program = blast_program,output_dir = all2all_out,blast_options = blast_options,input_prefix_path = params_list$GROUPS_PAT, params_list = params_list ) #paste(params_list$TEMP_PATH,"/","all2all/",sep="")

  mclapply(list.files(path = all2all_out,pattern = "*.all2all", ignore.case = T,full.names = T),function(in_file){
    out_file <- paste(all2all_out,tools::file_path_sans_ext(BiocGenerics::basename(in_file)),".out",sep="")
    convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
  }, mc.cores = params_list$numWorkers )

   mclapply(list.files(path = all2all_out,pattern = "*.out", ignore.case = T,full.names = T),function(in_file){
    out_file <- paste(all2all_out,tools::file_path_sans_ext(BiocGenerics::basename(in_file)),".wis_out",sep="")
    blast_GO <- GRObject_from_BLAST(blast_input = in_file,col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"), COMPLETE.format.ids = T, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list = params_list)
    wis_GO <- run_WISARD(blast_hits = blast_GO,score_col = "qcovhsp",COMPLETE.format.ids = T, params_list = params_list) #score_col=16
    #names(wis_GO) <- c(tools::file_path_sans_ext(BiocGenerics::basename(in_file)))
    #return(wis_GO)
    write.table(x = data.frame(wis_GO),file = out_file,quote = F,col.names = T,row.names = F)
  }, mc.cores = params_list$numWorkers )

  #save("wisard_results", file="files/all2all/wisard_results.RData")
  #load("files/all2all/wisard_results.RData")

   ##RUN RBH
   lapply("ungrouped", function(query){
     mclapply(grep("ungrouped",available_clusters,ignore.case = T,invert = T,value=T),function(subject){
       in1 <- paste(all2all_out,query,"-",subject,".wis_out",sep="")
       in2 <- paste(all2all_out,subject,"-",query,".wis_out",sep="")
       if (file.exists(in1) && file.exists(in2) && file.info(in1)$size > 0 && file.info(in2)$size > 0 ) {
          RBH(in1 = in1, in2 = in2, index.tables = T, weight.col = "qcovhsp",col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
       }
     }, mc.cores =numWorkers )
   })

  #calculate gene conservation - calculate_gene_conservation.R - probably not needed
   ##Maybe write one for cluster conservation/coverage across organisms

  #create organism sets

  ##Place ungrouped sequences into groups (oneway BLAST ungrouped cluster againts all genes)

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

  loaded_PARAMS$gene_list <- gene_list
  loaded_PARAMS$genes <- factor(scan(gene_list, character())) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
  loaded_PARAMS$genes <- genes[grep("gene",tolower(genes), invert = T, fixed = T)]

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
#' @note Both the input files are expected to have the same number of columns with matching column order (and names)
#'
#' @examples
#'    convert_BLAST_format(in_file1,outfile = out_file1,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'    convert_BLAST_format(in_file2,outfile = out_file2,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'    blast_GO1 <- GRObject_from_BLAST(blast_input = out_file1,col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"), COMPLETE.format.ids = T, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
#'    blast_GO2 <- GRObject_from_BLAST(blast_input = out_file2,col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"), COMPLETE.format.ids = T, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
#'    wis1 <- run_WISARD(blast_hits = blast_GO1,score_col = "qcovhsp",COMPLETE.format.ids = T) #score_col=16
#'    wis2 <- run_WISARD(blast_hits = blast_GO2,score_col = "qcovhsp",COMPLETE.format.ids = T) #score_col=16
#'    RBH(in1 = wis1, in2 = wis2, index.tables = T, weight.col = "qcovhsp",col.names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'
#' @param in1 Input (query->subject) BLAST hits table/filename
#' @param in2 Input (query<-subject) BLAST hits table/filename
#' @param transcript_ID_metadata Tab-delimited File with the filenames, indexed transcript IDs and the long transcript IDs.
#' @param col.names Columns names for the BLAST tables/files
#' @param weight.col Column index which must be used as edge weight (eg, qcovhsp - HSP coverage scores, bitscore etc)
#' @param sum.hit.weights Should the Weights be summed for each Query->Subject Hit? (TRUE/FALSE)
#' @param index.tables Should the IDs in the tables be indexed?
#' @return Named Vector of the organism details
#' @export
RBH <- function(in1,in2, transcript_ID_metadata=NULL, col.names=NULL, index.tables=F,weight.col, sum.hit.weights=T){

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

  if(sum.hit.weights){

  }

    if(index.tables){
      in1_data <- index_BLAST_table(in1_data,1,offset = 0)
      in1_data <- index_BLAST_table(in1_data,2,offset = nrow(in1_data))
      in2_data <- index_BLAST_table(in2_data,1,offset = nrow(in2_data))
      in2_data <- index_BLAST_table(in2_data,2,offset = 0)
    }

  in1_edge_list <- unlist(lapply(weight.cols, function(x){
    in1_g <- data.frame(from=in1_data[,1], to=in1_data[,2], weight=in1_data[,x], stringsAsFactors = T)
    in1_g <- in1_g[apply(in1_g, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),]
    tmp <- list(in1_g)
    names(tmp) <- colnames(in1_data)[x]
    return(tmp)
  }),recursive = F,use.names = T)

  in2_g <- data.frame(from=in2_data[,1], to=in2_data[,2], weight=in2_data[,weight.col], stringsAsFactors = T)


  in2_g <- in2_g[apply(in2_g, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),]

  in2_edge_list <- unlist(lapply(weight.cols, function(x){
    in2_g <- data.frame(from=in2_data[,1], to=in2_data[,2], weight=in2_data[,x], stringsAsFactors = T)
    in2_g <- in2_g[apply(in2_g, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),]
    tmp <- list(in2_g)
    names(tmp) <- colnames(in2_data)[x]
    return(tmp)
  }),recursive = F,use.names = T)

  #in1_valid_hits_from <- purrr::reduce(list(in1_g$from,in2_g$to), intersect)
  #in2_valid_hits_from <-purrr::reduce(list(in2_g$from,in1_g$to), intersect)
  #in1_valid_hits_to <- purrr::reduce(list(in1_g$to,in2_g$from), intersect)
  #in2_valid_hits_to <-purrr::reduce(list(in2_g$to,in1_g$from), intersect)

  in1_valid_hits <- unique(intersect(in1_g$from,in2_g$to),intersect(in1_g$to,in2_g$from)) #purrr::reduce(list(in1_g$from,in2_g$to), intersect)
  in2_valid_hits <- unique(intersect(in2_g$from,in1_g$to),intersect(in2_g$to,in1_g$from)) #purrr::reduce(list(in2_g$from,in1_g$to), intersect)

  #in1_g <- in1_g[which(!is.na(match(in1_g$from,intersect(in1_valid_hits_from,in2_valid_hits_to)))),]
  #in2_g <- in2_g[which(!is.na(match(in2_g$from,intersect(in2_valid_hits_from,in1_valid_hits_to)))),]
  #in1_g <- in1_g[which(!is.na(match(in1_g$to,intersect(in1_valid_hits_to,in2_valid_hits_from)))),]
  #in2_g <- in2_g[which(!is.na(match(in2_g$to,intersect(in2_valid_hits_to,in1_valid_hits_from)))),]

  in1_g <- in1_g[which(!is.na(match(in1_g$from,in1_valid_hits))),]
  in2_g <- in2_g[which(!is.na(match(in2_g$from,in2_valid_hits))),]
  in1_g <- in1_g[which(!is.na(match(in1_g$to,in2_valid_hits))),]
  in2_g <- in2_g[which(!is.na(match(in2_g$to,in1_valid_hits))),]

  in_g <- graph::MultiGraph(list(in1_g=in1_g,in2_g=in2_g), ignore_dup_edges = T)
  in_g <- graph::graphIntersect(in_g,in_g)

  in_g_list <- graph::extractGraphAM(in_g)

}

