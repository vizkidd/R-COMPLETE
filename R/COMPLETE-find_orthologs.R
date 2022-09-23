#' Load BLAST Hits into a data.frame
#'
#' Give the path to a BLAST Hits file to load it into a data.frame(BLAST HITs Table). The column names can be provided as col.names. Rows with NAs are automatically removed
#'
#' Note : Column indices will be used when processing the data. Column names are only for user reference. First column must be the query sequence ID and the second column must be the subject sequence ID
#'
#' Optional : You can provide the file/table with indexed Transcsript IDs. The format must be "file"\t"long_id"\t"index" (without a header). It can be generated with the  index_FASTA_IDs() (check ?index_FASTA_IDs or index_fastaIDs() in system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"))
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
#' WARNING : Do not index it more than once (although you may do index it as much as you like)
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
#' Functions accepts a filename or a BLAST table and converts it into a GenomicRanges::GRanges object (which can be used by R-WISARD). Give col.names when passing BLAST filenames. The accepted BLAST format is 6 and can be formatted from BLAST format 11 using blast_formatter (shown below). Columns indices can provided with the col.indices option. length(col.names) and ncol(blast_input) are assumed to be the same. BLAST table/file must be non-indexed (IDs not shortened) if COMPLETE.format.ids=T
#'
#' blast_formatter -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive"
#'
#' @param blast_input Filename with BLAST hits or a BLAST table
#' @param col.names Columns which for the BLAST tables/files
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE otherwise)
#' @param col.indices A Named Vector with indices of columns Query sequence ID (qseqid), Subject sequence ID (sseqid), E Value (evalue), Query start (qstart), Query end (qend), Subject start (sstart), Subject end (send), Bitscore (bitscore), Query HSP coverage(qcovhsp), Query length (qlen), Subject length(slen), Frames (frames), Percentage Identity (pident), Number of Gaps (gaps), Alignment Length (length), Subject strand (sstrand). Default values are col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17)
#' @return A GRanges object of BLAST hits
#' @export
GRObject_from_BLAST <- function(blast_input, col.names=NULL, COMPLETE.format.ids=F, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17)){

  if(!is.null(blast_input) && is.character(blast_input)){
    blast_input <- LoadBLASTHits(blast_input,col.names = col.names)
  }
  #tmp_blast <- blast_output
  #col_names <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand", "qlen", "slen", "qseq","sseq","nident","positive")
  #names(tmp_blast) <- col_names

  if(!is.null(col.names)){
    colnames(blast_input) <- col.names
  }

  if(COMPLETE.format.ids && check_params_loaded()){
    blast_input$subject_gene <- unlist(purrr::map(blast_input[,col.indices["sseqid"]],function(x){
      genes <- unlist(stringi::stri_split_fixed(x,pattern = SEQUENCE_ID_DELIM))
      return(genes[2])
    }))

    blast_input$query_gene <- unlist(purrr::map(blast_input[,col.indices["qseqid"]],function(x){
      genes <- unlist(stringi::stri_split_fixed(x,pattern = SEQUENCE_ID_DELIM))
      return(genes[2])
    }))
  }else{
    warning("Parameter file not loaded with load_params & COMPLETE.format.ids==FALSE")
  }

  blast_input <- blast_input[which(blast_input[,col.indices["evalue"]] < E_VALUE_THRESH),]

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
  if(COMPLETE.format.ids && check_params_loaded()){
    tmp_gr$query_org <- unlist(purrr::map(blast_input[,col.indices["qseqid"]],function(x){
      org <- unlist(stringi::stri_split_fixed(x,pattern = SEQUENCE_ID_DELIM))
      return(org[3])
    }))
    tmp_gr$subject_org <- unlist(purrr::map(blast_input[,col.indices["sseqid"]],function(x){
      org <- unlist(stringi::stri_split_fixed(x,pattern = SEQUENCE_ID_DELIM))
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
#' @param blast_hits GRanges Object of BLAST Hits (Query -> Subject)
#' @param score_col Column in the GRanges Object used for scoring intervals in WISARD
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE otherwise)
#' @return A GRanges object of BLAST hits
#' @export
run_WISARD <- function(blast_hits, score_col, COMPLETE.format.ids=F){

  if(!all(grepl(pattern="GRanges",x = class(blast_hits),ignore.case = T))){
    stop("This function only accepts formatted GRanges Object from GRObject_from_BLAST()")
  }

  if(!COMPLETE.format.ids || !check_params_loaded()){
    tryCatch(numWorkers <<- parallel::detectCores(all.tests = T, logical = T), error=function(){numWorkers <<- 2})
  }else{
    warning("Parameter file not loaded with load_params || COMPLETE.format.ids==FALSE")
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
  if(COMPLETE.format.ids && check_params_loaded()){
    for(child in child_results) {
      #print(child)
      g_name <- NULL
      g_name <- unique(unlist(stri_split_fixed(child$alignments$subject_gene,pattern=SEQUENCE_ID_DELIM,n = 1,tokens_only = T)))
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

#'Transcript Ortholog Extraction Function for R-COMPLETE pipeline
#'
#' This function calls the Transcript Ortholog Extraction pipeline which is used to reduce the pool of genes (step 1), reduce the pool of organisms and create sets of organisms (step 2), find transcript level orthologs (step 3). It takes only one argument which is the path/name of the BLAST program to use and refers to the values from the parameters file for other variables. Only exporting code for visibility purposes
#'
#'  * Step 1 - Genes which are available in all the reference organisms are chosen
#'  * Step 2 - A Per-Gene Conservation Score (GSC) is calculated from the availability of a gene across organisms (literally the count of organisms which have the gene, normalized to 1 relative to other genes). Genes which have GSC score below GENE_DROP_THRESHOLD (parameter) are dropped (GENE_DROP_THRESHOLD=0 does not omit any genes). Sets of organisms are created based on the available genes after GSC filtering. I can suggest reference organisms based on which ones have the maximum number of genes
#'  * Step 3 - Two way BLAST followed by HSP selection with WISARD and Two way RBH are performed. Only transcripts which are bi-directional best hits are kept for further analysis (RBH from both the directions, not RBH in itself is bi-directionaly from the point of the QUERY, We can also do an RBH from the context of the SUBJECT to verify if it did not pas RBH by chance (even though it is very unlikely))
#'
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program
#' @export
transcript_ortholog_extraction <- function(blast_program){

  if(!check_params_loaded()){
    stop()
  }

  #STEP 1
  #select genes which are available in all the reference organisms
  available_genes_list <- parallel::mclapply(paste(paste(OUT_PATH,"/files/genes/",sep=""),REF_ORGS,sep=""),function(x){
    if(file.exists(paste(x,"/AVAILABLE_GENES",sep="")) && file.info(paste(x,"/AVAILABLE_GENES",sep=""))$size > 0 ){
      return(scan(paste(x,"/AVAILABLE_GENES",sep=""), character()))
    }
  }, mc.cores =  numWorkers)
  available_genes <- unique(purrr::reduce(available_genes_list, union))

  #find which clusters they belong to
  odb_map_list <- parallel::mclapply(paste(paste(OUT_PATH,"/files/genes/",sep=""),REF_ORGS,sep=""),function(x){
    if(file.exists(paste(x,"/odb.final_map",sep="")) && file.info(paste(x,"/odb.final_map",sep=""))$size > 0 ){
      return(read.table(paste(x,"/odb.final_map",sep=""),header = F,sep = "\t",quote = "", col.names = c("cluster","genes")))
    }
  }, mc.cores =  numWorkers)
  odb_map <- unique(purrr::reduce(odb_map_list, inner_join, by = c("cluster","genes")))
  available_clusters <- unique(BiocGenerics::unlist(lapply(available_genes, function(gene){
    return(odb_map[grep(pattern=gene,x = odb_map$genes,ignore.case = T), c("cluster")])
  })))

  if(length(dir(GROUPS_PATH)==0)){
    group_FASTA_clusters(FASTA_OUT_PATH)
  }

  #calculate gene conservation - calculate_gene_conservation.R

  #create organism sets

  ##STEP 2 - Place ungrouped sequences into groups (oneway BLAST ungrouped cluster againts all genes)

  ##ITERATION 2 - two way RBH

  #cat $reference_ORGS files/oneway/SET > files/oneway/set.tmp

  #time all2all_refblast $reference_ORGS $fasta_path files/all2all/all2all.genelist files/all2all_final $blastdb_path $region $region tblastx files/oneway/set.tmp

  #time Rscript wisard.R files/oneway/set.tmp files/all2all_final files/gtf_stats.csv

  #readarray set_orgs < files/oneway/set.tmp

  #time parallel -j $((${#ref_orgs[@]}*${#set_orgs[@]})) "twoway_RBH files/all2all_final {1} {2} $PY3_PATH $RBH_SCRIPT $SAME_GENE" ::: ${ref_orgs[@]} ::: ${set_orgs[@]}

}

#' START HERE - Find Transcript Orthologs
#'
#' This function can be executed after EXTRACT_DATA()
#'
#' This is the main function which calls all the other functions and performs and end-end execution of finding transcript level orthologs. It runs the iterative Transcript Ortholog Extraction pipeline which is used to reduce the pool of genes, reduce the pool of organisms, find transcript level orthologs (check ?transcript_ortholog_extraction)
#'
#' It requires a filename of a formatted parameter file and a gene list (check the github repo for an example).
#'
#' ONLY USE THIS FUNCTION WHEN RUNNING THE PIPELINE OF COMPLETE. Use other helper function to work with custom BLAST files not generated by this R package
#'
#' @param param_file Filename of a formatted parameter file (check the github repo for an example)
#' @param gene_list File with a list of genes to extract data for(check the github repo for an example)
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program (or Sys.which("tblastx"))
#' @export
FIND_TRANSCRIPT_ORTHOLOGS <- function(param_file, gene_list, blast_program=Sys.which("tblastx")){
  set.seed(123)

  if(!grepl(x=Sys.info()["sysname"],pattern="linux",ignore.case = T)){
    stop("Pipeline only supports Linux (and bash) :(")
  }

  process_list <<- c()

  if(!file.exists(param_file) || file.info(param_file)$size < 0){
    stop("ERROR: parameters.txt is missing and is required\n")
  }

  param_file <<- param_file
  param_table <<- load_params(param_file)
  print(param_table)

  if (grepl(pattern = "bash",ignore.case = T,x = Sys.getenv("SHELL"))) {
    SHELL <<- Sys.getenv("SHELL")
    print(paste("SHELL :",SHELL))
  }else if(file.exists("/bin/bash")){
    SHELL <<- "/bin/bash"
    print(paste("SHELL :",SHELL))
  }else{
    stop(paste("SHELL (",SHELL,") : bash not available, or not in $PATH or SHELL=/bin/bash not set"))
  }

  gene_list <<- gene_list
  genes <<- factor(scan(gene_list, character())) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
  genes <<- genes[grep("gene",tolower(genes), invert = T, fixed = T)]

  REF_ORGS <<- factor(scan(REF_ORGS_FILE, character()))

  print(paste("MAX PROCESSES:",numWorkers))

  if(CLEAN_EXTRACT){
    unlink("files/oneway", recursive = T,force = T,expand = T)
    unlink("files/all2all", recursive = T,force = T,expand = T)
    unlink("files/all2all_final", recursive = T,force = T,expand = T)
    unlink("files/gene_thresholds.txt", recursive = T,force = T,expand = T)
  }

  dir.create("files/oneway",showWarnings = F, recursive = T)
  dir.create("files/all2all",showWarnings = F, recursive = T)
  dir.create("files/all2all_final",showWarnings = F, recursive = T)

  if (!stringi::stri_isempty(blast_program)) {
    transcript_ortholog_extraction(blast_program)
  }else{
    stop(paste(blast_program," NOT found. Is BLAST+ installed or in $PATH?"))
  }


}

#' Find Reciprocal Blast Hits (RBH)
#'
#' Find RBH between BLAST results of different organisms/genes/transcripts (FASTA/FASTQ). The BLAST results must be of the format 6 and can be converted from BLAST format 11 with blast_formatter.
#' The command with the required column names are given below.
#'
#' blast_formatter -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive"
#'
#' \textbf{NOTE : Both the input files are expected to have the same number of columns with matching column order (and names)}
#'
#' Optional : You can provide the file/table with indexed Transcsript IDs. The format must be "file"\t"long_id"\t"index" (without a header). It can be generated with the  index_FASTA_IDs() (check ?index_FASTA_IDs or index_fastaIDs() in system.file("exec", "functions.sh", mustWork = T ,package = "COMPLETE"))
#'
#' @param in1 Input (query->subject) BLAST hits table/filename
#' @param in2 Input (query<-subject) BLAST hits table/filename
#' @param transcript_ID_metadata Tab-delimited File with the filenames, indexed transcript IDs and the long transcript IDs.
#' @param col.names Columns which for the BLAST tables/files
#' @param weight.cols Column indices which must be used as edge weights (eg, mincovhsp - HSP coverage scores, bitscore etc)
#' @param index.tables Should the IDs in the tables be indexed?
#' @return Named Vector of the organism details
#' @export
RBH <- function(in1,in2, transcript_ID_metadata=NULL, col.names=NULL, index.tables=F,weight.col){

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

