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
    tmp_tbl <- blast_tables[[x]]
    tmp_tbl <- tmp_tbl %>% mutate(tmp_tbl[,query_index_cols[[x]]])
    tmp_tbl[,query_index_cols[[x]]] <- indexed_ids[match(tmp_tbl[,query_index_cols[[x]]],indexed_ids$id), c("index")]
    colnames(tmp_tbl)[query_index_cols[[x]]] <- paste("indexed",old_query_cols[[x]],sep="_")
    colnames(tmp_tbl)[ncol(tmp_tbl)] <- old_query_cols[[x]]

    tmp_tbl <- tmp_tbl %>% mutate(tmp_tbl[,subject_index_cols[[x]]])
    tmp_tbl[,subject_index_cols[[x]]] <- indexed_ids[match(tmp_tbl[,subject_index_cols[[x]]],indexed_ids$id), c("index")]
    colnames(tmp_tbl)[subject_index_cols[[x]]] <- paste("indexed",old_subject_cols[[x]],sep="_")
    colnames(tmp_tbl)[ncol(tmp_tbl)] <- old_subject_cols[[x]]
    return(tmp_tbl)
  },mc.silent = T,mc.cores = COMPLETE_vars$numWorkers)

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

#' Select highest scoring interval of non-overlapping HSPs from Bi-Directional BLAST Hits
#'
#' Select the highest scoring pairs (HSPs) which give the maximum coverage over the BLAST alignments of each transcript (without overlaps/minimal overlaps). These HSPs will then be used to find Transcript level orthologs across gene orthologs across organisms. This function only accepts bi-directional (Query <-> Subject) BLAST Hits formatted with GRObject_from_BLAST()
#'
#' @note The BLAST Hits assumed to be on the same frame (1/1,2/2,3/3). Hits across frames are NOT supported (1/2,1/3,...)
#'
#' @examples
#'
#' blast_GO <- GRObject_from_BLAST(blast_input = in_file, COMPLETE.format.ids = T, col.indices=c(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
#' wis_GO <- run_WISARD(blast_hits = blast_GO,score_col = "Hsp_score",COMPLETE.format.ids = T,params_list=NULL) #score_col=16
#' melt_wisard_list(wis_GO)
#'
#' @param blast_hits GRanges Object of BLAST Hits (Query -> Subject) (from GRObject_from_BLAST())
#' @param score_col Column in the GRanges Object used for scoring intervals in WISARD, Default would be "Hsp_score" (qcovhsp (Query Coverage HSP))
#' @param n_threads Number of Threads.
#' @param verbose Print Output?
#' @return A GRanges object of BLAST hits
#' @export
run_WISARD <- function(blast_hits, score_col,n_threads=4, verbose=F){

  if(!all(grepl(pattern="GRanges",x = class(blast_hits),ignore.case = T))){
    stop("This function only accepts formatted GRanges Object from GRObject_from_BLAST()")
  }

  gr <- blast_hits
  if(!any(grepl(pattern = "query_id",x = colnames(as.data.frame(gr)),ignore.case = F))){
    stop("Column Name : query_id is missing from the GRangesObject and is required. Use GRObject_from_BLAST()")
  }
  query_vector <- unique(gr$query_id)
  #print(S4Vectors::runValue(GenomicRanges::seqnames(gr)))
  subject_vector <- unique(S4Vectors::runValue(GenomicRanges::seqnames(gr)))
  bhits_combo <- unique(tidyr::crossing(query_vector,subject_vector))

  bhits_combo <- bhits_combo %>% mutate(num_factors=1:nrow(bhits_combo))

  gr <- dplyr::inner_join(as.data.frame(gr),bhits_combo, by=c("seqnames"="subject_vector","query_id"="query_vector"))

  gr <- split(gr,f = gr$num_factors)

  child_results <- parallel::mclapply(gr,function(x){
    result <- suppressMessages(invisible(wisard::get_wis(GenomicRanges::GRanges(x),max_score = score_col,overlap = 0)))
    #print("here3")
    if(result$max_score > 0 && length(result$alignments) > 0){
      #result_list <- c(result_list, result)
      return(result)
      #print(result_list)
    }
    #print(result_list)

  }, mc.cores = n_threads, mc.silent = !verbose) #!verbose

  # all_results <- list()
  # if(COMPLETE.format.ids && !is.null(params_list)){
  #   parallel::mclapply(child_results,function(child){
  #     #  print(child) #DEBUG
  #     g_name <- NULL
  #     g_name <- unique(child$alignments$subject_gene) #unique(unlist(stringi::stri_split_fixed(child$alignments$subject_gene,pattern=params_list$SEQUENCE_ID_DELIM,n = 1,tokens_only = T)))
  #     if(verbose){
  #       print(g_name)
  #     }
  #     #all_results<- c(all_results, child)
  #     if(!is.null(g_name)){
  #       if(is.null(all_results) || length(all_results)==0){
  #         all_results[[g_name]] <- child
  #       }else{
  #         all_results[[g_name]] <- merge_wisard_lists(all_results[[g_name]], child)
  #       }
  #       # names(all_results[[g_name]]) <- g_name
  #     }
  #   }, mc.preschedule = T, mc.cores = numWorkers, mc.silent = !verbose)
  # }else{
  #   all_results <- child_results #purrr::reduce(child_results, merge_wisard_table)
  # }

  #return(all_results)
  return(child_results)

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

  blast_table <- blast_table %>% mutate(from_org=unlist(lapply(stringi::stri_split_fixed(blast_table[,1],delimiter,n=4,tokens_only = T),function(x){return(x[COMPLETE_vars$FORMAT_ID_INDEX$ORG])})))
  blast_table <- blast_table %>% mutate(to_org=unlist(lapply(stringi::stri_split_fixed(blast_table[,2],delimiter,n=4,tokens_only = T),function(x){return(x[COMPLETE_vars$FORMAT_ID_INDEX$ORG])})))

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
  data.table::fwrite(list(semi_orths), file = file.path(identities_out_path,paste(gene,".semiorgs",sep = "")),quote = F, row.names = T, col.names = T, nThread = COMPLETE_vars$numWorkers)

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
    data.table::fwrite(list(discardable_orgs), file = file.path(oinfo_out_path,paste(gene,".dorgs",sep = "")))
  }
  if(length(saveable_orgs) > 0){
    data.table::fwrite(list(saveable_orgs), file = file.path(oinfo_out_path,paste(gene,".sorgs",sep = "")))
  }
  if(length(low_orthology_orgs) > 0){
    data.table::fwrite(list(low_orthology_orgs), file = file.path(oinfo_out_path,paste(gene,".lorgs",sep = "")))
  }
  if(length(nref_orgs) > 0){
    data.table::fwrite(list(nref_orgs), file = file.path(oinfo_out_path,paste(gene,".nref",sep = "")))
  }
  data.table::fwrite(min_seq_identity, file = file.path(identities_out_path,paste(gene,".min",sep = "")), quote = F, row.names = T, col.names = T, nThread = COMPLETE_vars$numWorkers)
  data.table::fwrite(max_seq_identity, file = file.path(identities_out_path,paste(gene,".max",sep = "")), quote = F, row.names = T, col.names = T, nThread = COMPLETE_vars$numWorkers)

  #fileConn<-file(thres_out,open = "at")
  #writeLines(paste(gene, min_gene_conservation,length(nref_orgs)==0,num_transcripts,sep = ","), fileConn)
  #close(fileConn)

  names(min_gene_conservation) <- gene
  return(min_gene_conservation)
}

#@param sep2 Delimiter 2 of the BLAST File columns. Default - c("","|",""). Check ?data.table::fwrite or ?data.table::fread
## @param transcript_region_lengths A data.frame (with one column for regional lengths and BLAST/Transcript IDs/Transcript Names as row.names) or a Named Vector or a Named List. Assumed to have short IDs
#' Calculate HSP Coverage of Bi-Directional BLAST Tables
#'
#' This function calculates the coverage of HSPs (sum(HSP Alignment lengths)/mRNA CDS Length) between each Query<->Subject Hits (Bi-directional). Hits are filtered based on run.mode option of this function.
#'
#' @param fw_blast_table Filename or BLAST Table with Query->Subject Hits (Forward)
#' @param bk_blast_table Filename or BLAST Table with Query<-Subject Hits (Backward)
#' @param col.indices A Named List with indices of columns Query sequence ID (qseqid) and Subject sequence ID (sseqid), Query Length (query_len), Subject Length (subject_len) and Alignment Length (HSP alignment length in our case)(align_len). Eg col.indices=list(qseqid=12,sseqid=1,query_len=13,subject_len=14,align_len=23)
#' @param group Name of the group/BLAST Run. Default -  "ungrouped"
#' @param run.mode "both" or "coverage_distance" (Default) or "coverage_filter" or "no_filter". "coverage_distance" - Hits are filtered based on distance between bi-directional minimum HSP coverages (coverage_distance <= min_coverage_filter). This option selects more BLAST hits and should be used when the coverage values are very low (and the BLAST Hits/sequences are distant). "coverage_filter" - Filters Hits based on minimum coverage of HSPs from either direction. Use this option when the coverage values are high (and the BLAST Hits/sequences are closely related). "both" - Uses both "coverage_distance" and "coverage_filter" and is very strict. "no_filter" - Only calculates HSP coverages and does not filter any Hits
#' @param min_coverage_filter Minimum HSP Coverage value to filter out Hits (Default - 0.5)
#' @param COMPLETE.format.ids Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, Default - FALSE otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (ONLY FOR blast_table,transcript_region_lengths are assumed to have short IDs)
#' @param params_list Output of load_params()
#' @param sep Delimiter for the input Files. Only valid if blast_table is a file. Default - '\t'
#' @param header Does the input files have header?. Only valid if blast_table is a file. Default - FALSE
#' @param n_threads Number of threads
#' @param verbose Print Output Messages?
#' @param seed Seed Value
#' @return BLAST table with Hits which pass min_coverage_filter
#' @export
calculate_HSP_coverage <- function(fw_blast_table,bk_blast_table,col.indices, group="ungrouped",run.mode="coverage_distance",min_coverage_filter=0.5, COMPLETE.format.ids=F,params_list, sep="\t", header=F, n_threads=8, verbose=T, seed=123){ #,sep2 = c("","|","") #transcript_region_lengths
  set.seed(seed)
  if(!grepl(pattern ="coverage_distance|coverage_filter|both|no_filter",ignore.case = T,x = run.mode) || is.null(run.mode)){
    stop("run.mode must be either 'both' or 'coverage_distance' or 'coverage_filter' or 'no_filter'")
  }

  blast_table <- parallel::mclapply(list(fw_blast_table,bk_blast_table), function(in_file){
    if(!is.null(in_file) && is.character(in_file)){
      if(file.exists(in_file)){
        return(LoadBLASTHits(infile = in_file, sep=sep, header=header))
      }else{
        message(paste(in_file,"does not exist!"))
        return(NULL)
      }
    }else if(!is.null(in_file) && !is.character(in_file)){
      return(in_file)
      #if(!is.null(col.names)){
      #  colnames(blast_table) <- col.names
      #}
    }else{
      message(paste("Input does not exist or is not a table!"))
      return(NULL)
    }
  }, mc.cores = 2, mc.silent = F)
  if (any(sapply(blast_table,is.null))) {
    return(NULL)
  }
  #print(head(blast_table)) #DEBUG
  # if(is.vector(transcript_region_lengths)){
  #   region_lengths <- as.data.frame(x=transcript_region_lengths,row.names = names(transcript_region_lengths)) #vector
  # }else if(inherits(transcript_region_lengths, "list")){ #is.list(transcript_region_lengths)
  #   region_lengths <- as.data.frame(x=unlist(transcript_region_lengths,use.names = T,recursive = T),row.names = names(transcript_region_lengths)) # list
  # }else if(inherits(transcript_region_lengths, "data.frame") && ncol(transcript_region_lengths) == 1){ #is.data.frame(transcript_region_lengths)
  #   region_lengths <- transcript_region_lengths
  # }else{
  #   stop("transcript_region_lengths must be a Data.Frame (with one column for (CDS/UTR) lengths and BLAST IDs/Transcript IDs as row.names) or a Named Vector or a Named List")
  # }

  q_ids <- unique(unlist(c(blast_table[[1]][,col.indices[["qseqid"]]],blast_table[[2]][,col.indices[["seqid"]]])))
  s_ids <- unique(unlist(c(blast_table[[1]][,col.indices[["sseqid"]]],blast_table[[2]][,col.indices[["qseqid"]]])))

  #print(region_lengths) #DEBUG
  #print(paste(fw_ids,collapse = ",")) #DEBUG
  #print(paste(bk_ids,collapse = ",")) #DEBUG

  id_combinations <- unique(tidyr::crossing(q_ids,s_ids)) #id_combinations <- data.frame(q_ids=q_ids,s_ids=s_ids)
  id_combinations <- id_combinations %>% mutate(num_factors=1:nrow(id_combinations))
  #id_combinations$num_factors <- as.factor(id_combinations$num_factors)
  id_combinations_rev <- unique(tidyr::crossing(s_ids,q_ids))
  #id_combinations[match(id_combinations$q_ids,id_combinations$s_ids),c("num_factors")] <- id_combinations[match(id_combinations$s_ids,id_combinations$q_ids),c("num_factors")]
  id_combinations <- dplyr::full_join(id_combinations,id_combinations_rev, by=c("q_ids"="s_ids","s_ids"="q_ids"))
  id_combinations$num_factors <- unlist(apply(id_combinations,MARGIN = 1,FUN=function(row){
    tmp_fct <- row[c("num_factors")]
    if(is.na(tmp_fct)){
      tmp_fct <- id_combinations[id_combinations$q_ids==row[c("s_ids")] & id_combinations$s_ids==row[c("q_ids")],c("num_factors")]
    }
    return(tmp_fct)
  }))



  blast_table[[1]] <- blast_table[[1]] %>% mutate(direction="forward")
  blast_table[[2]] <- blast_table[[2]] %>% mutate(direction="reverse")

  blast_table <- dplyr::full_join(blast_table[[1]],blast_table[[2]])

  # split_btbl <- parallel::mclapply(blast_table,function(b_tbl){
  #   q_col <- colnames(b_tbl)[col.indices[["qseqid"]]]
  #   s_col <- colnames(b_tbl)[col.indices[["sseqid"]]]
  #   id_combinations <- rename_at(id_combinations, "q_ids", ~ q_col)
  #   id_combinations <- rename_at(id_combinations, "s_ids", ~ s_col)
  #   b_tbl$num_factors <- NULL
  #   ret_btbl <- dplyr::inner_join(b_tbl,id_combinations, by = c(s_col, q_col)) #[,c(col.indices[["qseqid"]],col.indices[["sseqid"]])]
  #   ret_btbl <- split(ret_btbl,f = ret_btbl$num_factors)
  #   return(ret_btbl)
  # }, mc.cores = n_threads,mc.silent = F,mc.set.seed = seed)

  q_col <- colnames(blast_table)[col.indices[["qseqid"]]]
  s_col <- colnames(blast_table)[col.indices[["sseqid"]]]
  id_combinations <- rename_at(id_combinations, "q_ids", ~ q_col)
  id_combinations <- rename_at(id_combinations, "s_ids", ~ s_col)
  blast_table$num_factors <- NULL
  split_btbl <- dplyr::full_join(blast_table,id_combinations, by = c(s_col, q_col)) #[,c(col.indices[["qseqid"]],col.indices[["sseqid"]])]
  split_btbl <- split(split_btbl,f = split_btbl$num_factors)


  passed_coverage <- lapply(split_btbl, function(hit_dir){ #parallel::mclapply #furrr::future_map
    parallel::mclapply(split(hit_dir,f = factor(hit_dir$direction)), function(hit_split){#parallel::mclapply #furrr::future_map

      if(is.null(hit_split) || nrow(hit_split)==0){
        return(NULL)
      }

      hits_GO <- GRObject_from_BLAST(blast_input = hit_split,COMPLETE.format.ids = T,col.indices=list(qseqid=12,sseqid=1,qstart=15,qend=16,sstart=2,send=3,sstrand=5),params_list = params_list)
      hits_ovlps <- as.data.frame(GenomicRanges::findOverlaps(hits_GO,hits_GO,type = c("within")))
      hits_rle <- Rle(values = hits_ovlps[,1], lengths = hits_ovlps[,2])
      if(length(hits_rle) > 1 && try(any(duplicated(runLength(hits_rle))))){
        hit_split <- hit_split[-runValue(hits_rle)[-which.min(runLength(hits_rle))],]
      }

      q_length <- hit_split[,col.indices[["query_len"]]] #sum(unique()) # mean()
      s_length <- hit_split[,col.indices[["subject_len"]]] #sum(unique()) # mean()

      #fw_align_length <- mean(hit_split[,col.indices[["align_len"]]]) #sum
      #bk_align_length <- mean(hit_split[,col.indices[["align_len"]]]) #sum
      dir_align_length <- hit_split[,col.indices[["align_len"]]] #sum

      #print(paste(fw_hits[,col.indices[["align_len"]]], q_length) )#DEBUG
      #print(paste(bk_hits[,col.indices[["align_len"]]], s_length) )#DEBUG
      raw_cov_q <-  sum(dir_align_length / q_length)
      raw_cov_s <- sum(dir_align_length / s_length)

      #print(paste(q_length,s_length,sep="/")) #DEBUG
      #print(paste(raw_cov_fw,raw_cov_bk,sep="/")) #DEBUG

      if(length(raw_cov_q) != 0 && length(raw_cov_s) != 0){
        # if(raw_cov_fw > 1){
        #   cov_q <- 1 / raw_cov_fw  #sum(subject_hits[,col.indices[["send"]]] - subject_hits[,col.indices[["sstart"]]])
        # }else{
        #   cov_q <- raw_cov_fw
        # }
        # if(raw_cov_bk > 1){
        #   cov_s <-  1 / raw_cov_bk #sum(query_hits[,col.indices[["qend"]]] - query_hits[,col.indices[["qstart"]]])
        # }else{
        #   cov_s <- raw_cov_bk
        # }

        ##raw_cov_fw[which(raw_cov_fw > 1)] <- 1/raw_cov_fw[which(raw_cov_fw > 1)]
        ##raw_cov_bk[which(raw_cov_bk > 1)] <- 1/raw_cov_bk[which(raw_cov_bk > 1)]
        scaled_covs <- as.vector(scale(c(raw_cov_q,raw_cov_s),center = F))
        cov_q <- scaled_covs[1] #raw_cov_fw
        cov_s <- scaled_covs[2] #raw_cov_bk
        #cov_q <- raw_cov_q
        #cov_s <- raw_cov_s

        coverage_distance= abs(1-(cov_q/cov_s)) #1-((cov_q/cov_s)/100) #cov_q/cov_s #sqrt((1-(cov_q/cov_s))^2) #sqrt((1-1/(cov_q/cov_s))^2) #1 - (cov_q/cov_s) #1- (1/(cov_q/cov_s))

        if(length(coverage_distance) > 0 && length(cov_q) > 0  && length(cov_s) > 0 ){
          #if(verbose) print(paste(paste(x,"(",cov_q,")",sep = ""),paste(y,"(",cov_s,")",sep = ""),sep="->"),)
          #if(verbose) print(paste(x,"[","q_align_len/q_CDS_length:",q_align_length,"/",q_length,"]",sep=""))
          #if(verbose) print(paste(y,"[","s_align_len/s_CDS_length:",s_align_length,"/",s_length,"]",sep=""))
          if(verbose) print(paste(unique(hit_split[,col.indices[["qseqid"]]]),"->",unique(hit_split[,col.indices[["sseqid"]]]),"[","(cov_q, cov_s): (",cov_q,", ",cov_s,")]",sep=""))
          if(verbose) print(paste("Coverage Distance:",coverage_distance,sep=""))
          if(verbose) print(paste(raw_cov_q,raw_cov_s,sep="/"))  #DEBUG
          if(verbose) print(paste("Alignment Length: ",sum(dir_align_length)))
          #if(cov_q >= min_coverage_filter && cov_s >= min_coverage_filter){
          #  return(data.frame(query=x,subject=y))
          #}

          return_data <- list(BLAST_hits=hit_split,"query_id"=unique(hit_split[,col.indices[["qseqid"]]]),"subject_id"=unique(hit_split[,col.indices[["sseqid"]]]),"coverage_distance"=coverage_distance,"raw_cov_q"=raw_cov_q,"raw_cov_s"=raw_cov_s,"min_raw_cov"=min(raw_cov_q,raw_cov_s),"max_raw_cov"=max(raw_cov_q,raw_cov_s),"cov_q"=cov_q,"cov_s"=cov_s, "min"=min(cov_q,cov_s),"max"=max(cov_q,cov_s),"group"=group,"align_length"=sum(dir_align_length),"query_len"=unique(q_length),"subject_len"=unique(s_length),"direction"=unique(hit_split$direction)) #"raw_cov_s"=raw_cov_s,"raw_cov_q"=raw_cov_q #data.frame("query_id"=c(fw_x),"subject_id"=c(bk_y),"coverage_distance"=coverage_distance,"raw_cov_fw"=raw_cov_fw,"raw_cov_bk"=raw_cov_bk,"cov_q"=cov_q,"cov_s"=cov_s, "min"=min(cov_q,cov_s),"max"=max(cov_q,cov_s),"group"=group)

          if(grepl(pattern ="both",ignore.case = T,x = run.mode)){
            if(coverage_distance <= min_coverage_filter && cov_q >= min_coverage_filter && cov_s >= min_coverage_filter ){
              return(return_data)
            }
          }else if(grepl(pattern ="coverage_distance",ignore.case = T,x = run.mode)){
            #IF Min_Coverage Distance <= min_coverage_filter, return (min coverage distance=0 is full/same coverage between directions, min coverage distance=1 is no coverage)
            if(coverage_distance <= min_coverage_filter ){ #&& cov_q >= min_coverage_filter && cov_s >= min_coverage_filter
              return(return_data)
            }
          }else if(grepl(pattern ="coverage_filter",ignore.case = T,x = run.mode)){
            if(cov_q >= min_coverage_filter && cov_s >= min_coverage_filter){
              return(return_data)
            }
          }else{
            return(return_data)
          }

        }
      }

    }, mc.cores = ceiling(((n_threads - 2)+1)/2) ,mc.set.seed = seed,mc.silent = !verbose) #) #, .options = furrr::furrr_options(seed = TRUE, scheduling=F)) #, mc.cores = ceiling(((n_threads - 2)+1)/2) ,mc.set.seed = seed,mc.silent = F)
  }) #, mc.cores = 2,mc.set.seed = seed,mc.silent = F) #) # , .options = furrr::furrr_options(seed = TRUE, scheduling=F))#, mc.cores = 2,mc.set.seed = seed,mc.silent = F)

  #passed_coverage <- furrr::future_map2(.x=id_combinations$fw_ids, .y=id_combinations$bk_ids, .f=function(fw_x,bk_y){
  # passed_coverage <- parallel::mclapply(seq_along(1:nrow(id_combinations)),function(idx){
  #   #print(data.frame(query=x,subject=y))
  #   fw_x <- id_combinations[idx,1]
  #   bk_y <- id_combinations[idx,2]
  #   if(COMPLETE.format.ids){
  #     tx_x <- stringi::stri_split(str = fw_x,fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)[,1]
  #     x_short <- stringi::stri_split(str = tx_x, fixed = params_list$TRANSCRIPT_ID_DELIM, simplify=T)[,COMPLETE_vars$FORMAT_ID_INDEX$TRANSCRIPT_ID]
  #     tx_y <- stringi::stri_split(str = bk_y,fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)[,1]
  #     y_short <- stringi::stri_split(str = tx_y, fixed = params_list$TRANSCRIPT_ID_DELIM, simplify=T)[,COMPLETE_vars$FORMAT_ID_INDEX$TRANSCRIPT_ID]
  #   }else{
  #     x_short <- unlist(fw_x)
  #     y_short <- unlist(bk_y)
  #   }
  #   if(x_short!=y_short && fw_x != bk_y){
  #     #if(verbose) print(paste(c(x,x_short,q_length,y,y_short,s_length),collapse = ":"))
  #
  #     #fw_hits <- blast_table[[1]][which(!is.na(match(blast_table[[1]][,col.indices[["qseqid"]]],x_short)) & !is.na(match(blast_table[[1]][,col.indices[["sseqid"]]],y_short))),]
  #     #bk_hits <- blast_table[[2]][which(!is.na(match(blast_table[[2]][,col.indices[["qseqid"]]],y_short)) & !is.na(match(blast_table[[2]][,col.indices[["sseqid"]]],x_short))),]
  #     fw_hits <- blast_table[[1]][intersect(grep(pattern = x_short,x = blast_table[[1]][,col.indices[["qseqid"]]], fixed = T), grep(pattern = y_short,x = blast_table[[1]][,col.indices[["sseqid"]]], fixed = T)),]
  #     bk_hits <- blast_table[[2]][intersect(grep(pattern = y_short,x = blast_table[[2]][,col.indices[["qseqid"]]], fixed = T), grep(pattern = x_short,x = blast_table[[2]][,col.indices[["sseqid"]]], fixed = T)),]
  #
  #     if(is.null(fw_hits) || is.null(bk_hits) || nrow(fw_hits) == 0 || nrow(bk_hits) == 0){
  #       #print(paste(fw_x,bk_y,sep="///"))
  #       return(NULL)
  #     }
  #
  #     ##ABSORB hits overlapped by other hits
  #     hits_list <- lapply(list(fw_hits,bk_hits), function(hits_dir){
  #       hits_GO <- GRObject_from_BLAST(blast_input = hits_dir,COMPLETE.format.ids = T,col.indices=list(qseqid=12,sseqid=1,qstart=15,qend=16,sstart=2,send=3,sstrand=5),params_list = params_list)
  #       hits_ovlps <- as.data.frame(findOverlaps(hits_GO,hits_GO,type = c("within")))
  #       hits_rle <- Rle(values = hits_ovlps[,1], lengths = hits_ovlps[,2])
  #       if(length(hits_rle) > 1){
  #         hits_dir <- hits_dir[-runValue(hits_rle)[-which.min(runLength(hits_rle))],]
  #       }
  #       return(hits_dir)
  #     })
  #
  #     fw_hits <- hits_list[[1]]
  #     bk_hits <- hits_list[[2]]
  #
  #     #print("---------------------") #DEBUG
  #     #print(head(fw_hits)) #DEBUG
  #     #print(head(bk_hits)) #DEBUG
  #     #print("---------------------") #DEBUG
  #
  #     ##s_overlaps <- dissolve_GR_Overlaps(subject_result)
  #     #s_align_length <- sum(subject_hits[,col.indices[["send"]]] - subject_hits[,col.indices[["sstart"]]])
  #     #q_align_length <- sum(query_hits[,col.indices[["qend"]]] - query_hits[,col.indices[["qstart"]]])
  #     #cov_q <- q_align_length/q_length
  #     #cov_s <- s_align_length/s_length
  #
  #     #print(subject_hits[,col.indices[["subject_len"]]])
  #     #print(subject_hits[,col.indices[["Hsp_align.len"]]])
  #     #print(subject_hits[,col.indices[["subject_len"]]] / subject_hits[,col.indices[["Hsp_align.len"]]] ) #DEBUG
  #     #print(query_hits[,col.indices[["query_len"]]])
  #     #print(query_hits[,col.indices[["Hsp_align.len"]]])
  #     #print(query_hits[,col.indices[["query_len"]]] / query_hits[,col.indices[["Hsp_align.len"]]] ) #DEBUG
  #     #print(subject_hits[,col.indices[["align_len"]]])
  #
  #     #q_length <- as.numeric(region_lengths[x_short,])
  #     #s_length <- as.numeric(region_lengths[y_short,])
  #
  #     fw_s_length <- sum(fw_hits[,col.indices[["subject_len"]]]) #sum(unique()) # mean()
  #     bk_s_length <- sum(bk_hits[,col.indices[["subject_len"]]]) #sum(unique()) # mean()
  #
  #     fw_align_length <- sum(fw_hits[,col.indices[["align_len"]]])
  #     bk_align_length <- sum(bk_hits[,col.indices[["align_len"]]])
  #
  #     #print(paste(fw_hits[,col.indices[["align_len"]]], q_length) )#DEBUG
  #     #print(paste(bk_hits[,col.indices[["align_len"]]], s_length) )#DEBUG
  #     raw_cov_fw <-  fw_align_length / fw_s_length
  #     raw_cov_bk <- bk_align_length / bk_s_length
  #
  #     #print(paste(q_length,s_length,sep="/")) #DEBUG
  #     #print(paste(raw_cov_fw,raw_cov_bk,sep="/")) #DEBUG
  #
  #     if(length(raw_cov_fw) != 0 && length(raw_cov_bk) != 0){
  #       # if(raw_cov_fw > 1){
  #       #   cov_q <- 1 / raw_cov_fw  #sum(subject_hits[,col.indices[["send"]]] - subject_hits[,col.indices[["sstart"]]])
  #       # }else{
  #       #   cov_q <- raw_cov_fw
  #       # }
  #       # if(raw_cov_bk > 1){
  #       #   cov_s <-  1 / raw_cov_bk #sum(query_hits[,col.indices[["qend"]]] - query_hits[,col.indices[["qstart"]]])
  #       # }else{
  #       #   cov_s <- raw_cov_bk
  #       # }
  #
  #       #raw_cov_fw[which(raw_cov_fw > 1)] <- 1/raw_cov_fw[which(raw_cov_fw > 1)]
  #       #raw_cov_bk[which(raw_cov_bk > 1)] <- 1/raw_cov_bk[which(raw_cov_bk > 1)]
  #       scaled_covs <- as.vector(scale(c(raw_cov_fw,raw_cov_bk),center = F))
  #       cov_q <- scaled_covs[1] #raw_cov_fw
  #       cov_s <- scaled_covs[2] #raw_cov_bk
  #
  #       coverage_distance= 1-((cov_q/cov_s)/100) #cov_q/cov_s #sqrt((1-(cov_q/cov_s))^2) #sqrt((1-1/(cov_q/cov_s))^2) #1 - (cov_q/cov_s) #1- (1/(cov_q/cov_s))
  #
  #       if(length(coverage_distance) > 0 && length(cov_q) > 0  && length(cov_s) > 0 ){
  #         #if(verbose) print(paste(paste(x,"(",cov_q,")",sep = ""),paste(y,"(",cov_s,")",sep = ""),sep="->"),)
  #         #if(verbose) print(paste(x,"[","q_align_len/q_CDS_length:",q_align_length,"/",q_length,"]",sep=""))
  #         #if(verbose) print(paste(y,"[","s_align_len/s_CDS_length:",s_align_length,"/",s_length,"]",sep=""))
  #         if(verbose) print(paste(fw_x,"->",bk_y,"[","(cov_q, cov_s): (",cov_q,", ",cov_s,")]",sep=""))
  #         if(verbose) print(paste("Coverage Distance:",coverage_distance,sep=""))
  #         if(verbose) print(paste("raw_cov_fw, raw_cov_bk:",raw_cov_fw,",",raw_cov_bk,sep=""))
  #         #if(verbose) print(paste(q_length,s_length,sep="/"))  #DEBUG
  #         if(verbose) print(paste(sum(fw_hits[,col.indices[["align_len"]]]),sum(bk_hits[,col.indices[["align_len"]]]),sep="/"))
  #         #if(cov_q >= min_coverage_filter && cov_s >= min_coverage_filter){
  #         #  return(data.frame(query=x,subject=y))
  #         #}
  #
  #         return_data <- data.frame("query_id"=c(fw_x),"subject_id"=c(bk_y),"coverage_distance"=coverage_distance,"raw_cov_fw"=raw_cov_fw,"raw_cov_bk"=raw_cov_bk,"cov_q"=cov_q,"cov_s"=cov_s, "min"=min(cov_q,cov_s),"max"=max(cov_q,cov_s),"group"=group) #"query_len"=q_length, "subject_len"=s_length
  #         #return_data <- data.frame(c(x),c(y),q_length,s_length,coverage_distance,raw_cov_fw,raw_cov_bk,cov_q,cov_s, min(cov_q,cov_s),max(cov_q,cov_s),group)
  #         #colnames(return_data) <- c("query","subject","query_len", "subject_len","coverage_distance","raw_cov_fw","raw_cov_bk","cov_q","cov_s", "min","max","group")
  #
  #         if(grepl(pattern ="both",ignore.case = T,x = run.mode)){
  #           if(coverage_distance <= min_coverage_filter && cov_q >= min_coverage_filter && cov_s >= min_coverage_filter ){
  #             return(return_data)
  #           }
  #         }else if(grepl(pattern ="coverage_distance",ignore.case = T,x = run.mode)){
  #           #IF Min_Coverage Distance <= min_coverage_filter, return (min coverage distance=0 is full/same coverage between directions, min coverage distance=1 is no coverage)
  #           if(coverage_distance <= min_coverage_filter ){ #&& cov_q >= min_coverage_filter && cov_s >= min_coverage_filter
  #             return(return_data)
  #           }
  #         }else if(grepl(pattern ="coverage_filter",ignore.case = T,x = run.mode)){
  #           if(cov_q >= min_coverage_filter && cov_s >= min_coverage_filter){
  #             return(return_data)
  #           }
  #         }else{
  #           return(return_data)
  #         }
  #       }
  #     }else{
  #       return(NULL)
  #     }
  #   }else{return(NULL)}
  # }, .options = furrr::furrr_options(seed = TRUE, scheduling=params_list$numWorkers))

  #tmp_passed_coverage <<- passed_coverage
  #passed_coverage <- dplyr::bind_rows(passed_coverage[!sapply(passed_coverage, is.null)])
  #print(head(passed_coverage))
  #if(nrow(passed_coverage) > 0){
    #blast_table <- lapply(blast_table, function(b_table){
    #  b_table <- b_table[!is.na(match(b_table[,col.indices[["qseqid"]]],unique(c(passed_coverage$query_id,passed_coverage$subject_id)))),]
    #  b_table <- b_table[!is.na(match(b_table[,col.indices[["sseqid"]]],unique(c(passed_coverage$query_id,passed_coverage$subject_id)))),]
    #})
    return(list( blast_table=unique(dplyr::bind_rows(lapply(passed_coverage, function(x){return(dplyr::bind_rows(lapply(x,function(y){
      #print(y)
      if(!is.null(y["BLAST_hits"])) {return(y["BLAST_hits"])}
    })))})) ),
    coverage=unique(dplyr::bind_rows(lapply(passed_coverage, function(x){return(dplyr::bind_rows(lapply(x,function(y){
      if(!is.null(y["BLAST_hits"])) {
        y["BLAST_hits"] <- NULL
        return(data.frame(y))
      }
    })))})))) )
  #}else{
  #  message(paste("No hits passed coverage filter :", group))
  #  return(NULL)
  #}

}

# @param all_gtf_stats Coerced GTF stats from all the organisms. Found in paste(params_list$OUT_PATH,"/all_gtf_stats.csv",sep=""). This file is generated by EXTRACT_DATA() and is required for calculating HSP Coverage
#' Internal Function - Transcript Ortholog Extraction Function for R-COMPLETE pipeline
#'
#' This function calls the Transcript Ortholog Extraction pipeline which is used to reduce the pool of genes (step 1), reduce the pool of organisms and create sets of organisms (step 2), find transcript level orthologs (step 3). It takes only one argument which is the path/name of the BLAST program to use and refers to the values from the parameters file for other variables.
#'
#'  * Step 1 - Genes which are available in all the reference organisms are chosen
#'  * Step 2 - A Per-Gene Conservation Score (GSC) is calculated from the availability of a gene across organisms (literally the count of organisms which have the gene, normalized to 1 relative to other genes). Genes which have GSC score below GENE_DROP_THRESHOLD (parameter) are dropped (GENE_DROP_THRESHOLD=0 does not omit any genes). Sets of organisms are created based on the available genes after GSC filtering. I can suggest reference organisms based on which ones have the maximum number of genes
#'  * Step 3 - Two way BLAST followed by HSP selection with WISARD and Two way RBH are performed. Only transcripts which are bi-directional best hits are kept for further analysis (RBH from both the directions, not RBH in itself is bi-directionaly from the point of the QUERY, We can also do an RBH from the context of the SUBJECT to verify if it did not pas RBH by chance (even though it is very unlikely))
#'
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program.
#' @param blast_options BLAST options for the BLAST program
#' @param transcript_region Transcript region (.cds, .3utr, .5utr, .exon, and other gtf regions(not tested)) (eg, ".cds" or ".fa"). Essentially this would be the file extensions of the input FASTA file(names)
#' @param run.mode "both" or "coverage_distance" (Default) or "coverage_filter" or "no_filter". "coverage_distance" - Hits are filtered based on distance between bi-directional minimum HSP coverages (coverage_distance <= min_coverage_filter). This option selects more BLAST hits and should be used when the coverage values are very low (and the BLAST Hits/sequences are distant). "coverage_filter" - Filters Hits based on minimum coverage of HSPs from either direction. Use this option when the coverage values are high (and the BLAST Hits/sequences are closely related). "both" - Uses both "coverage_distance" and "coverage_filter" and is very strict. "no_filter" - Only calculates HSP coverages and does not filter any Hits (Argument used for calculate_HSP_coverage())
#' @param params_list Output of load_params()
#' @param clusters_left Vector of file names (Set/Subset) in input_dir (OG Clusters/Genes) to BLAST clusters_right with (Can be same as clusters_right). This would be the FASTA file name(s) (without file extension)
#' @param clusters_right Vector of file names (Set/Subset) in input_dir (OG Clusters/Genes) to BLAST clusters_left with (Can be same as clusters_left). This would be the FASTA file name(s) (without file extension)
#' @param input_dir Give the directory with FASTA files (to BLAST between them using blast_program)
#' @param output_dir Directory for saving output files (\*.out, \*.all2all, \*.wis_out,\*.rbh_out)
#' @param clean_extract Delete Output file if exists? (Optional). Default - F
#' @param verbose Print DEBUG Messages?. Default - F
#' @export
extract_transcript_orthologs <- function(blast_program, blast_options, transcript_region=".cds",run.mode="both", params_list, clusters_left, clusters_right, input_dir,output_dir,clean_extract=F,verbose=F, seed=123){ #all_gtf_stats
  set.seed(seed)
  if(stringi::stri_isempty(blast_program)){
    stop("extract_transcript_orthologs() - BLAST+ not found in $PATH. Provide blast_program")
  }
  else{
    BLAST_BIN <- dirname(blast_program)
  }

  if(!any(grepl(x = class(params_list), pattern = "COMPLETE-options"))){
    stop("Error: params_list is not a COMPLETE-options class. Use load_params()")
  }

  # unique_lengths <- unique(all_gtf_stats[,c("transcript_id","total_cds_len")])
  # tx_CDS_lengths <- data.frame(length=unique_lengths$total_cds_len, row.names = unique_lengths$transcript_id)

  tictoc::tic(msg=paste("Extracting Transctipt Orthologs :",paste(clusters_left,"<->",clusters_right,sep=""),":"))

  unlink(x = list.files(path = tempdir(check = T), pattern="*blast*", ignore.case = T, full.names = T), recursive = F, expand = T)
  #blast_DB_dir <- params_list$BLAST_DB_PATH

  # file.copy(paste(loaded_PARAMS$GROUPS_PATH,"/ungrouped.cds",sep=""), paste(tempdir(),"/ungrouped.cds",sep=""))
  # make_BLAST_DB(fasta_file=paste(tempdir(),"/ungrouped.cds",sep=""), blast_bin=BLAST_BIN,clean_extract = loaded_PARAMS$CLEAN_EXTRACT, verbose= verbose)

  file.copy(paste(input_dir,"/",clusters_right,".",transcript_region,sep=""), paste(tempdir(),"/",clusters_right,".",transcript_region,sep=""), overwrite = F)
  make_BLAST_DB(fasta_file= paste(tempdir(),"/",clusters_right,".",transcript_region,sep=""), blast_bin=BLAST_BIN,clean_extract = params_list$CLEAN_EXTRACT, verbose= verbose)

  all2all_BLAST(first_list = clusters_left, second_list = clusters_right, file_ext=transcript_region, blast_program = blast_program,output_dir =output_dir,blast_options = blast_options, blast.sequence.limit = 5000, input_prefix_path = input_dir, params_list = params_list, COMPLETE.format.ids = T, clean_extract=clean_extract, n_threads=params_list$numWorkers, verbose = verbose, seed=seed,return_data=FALSE ) # return_f_callback = NULL #all2all_GRObjects <-  #blast_DB_dir = blast_DB_dir #second_list = grep(clusters_left,clusters_right,ignore.case = T,invert = T,value=T)
  #all2all_BLAST(first_list = clusters_right, second_list = clusters_left,blast_program = blast_program,output_dir = output_dir,blast_options = blast_options,input_prefix_path = input_dir, params_list = params_list, COMPLETE.format.ids = T, keep.output.files = T ) #blast_DB_dir = blast_DB_dir #first_list = grep(clusters_left,clusters_right,ignore.case = T,invert = T,value=T)
  #save(all2all_GRObjects, file="all2all_GRObjects.RData") #DEBUG
  # parallel::mclapply(list.files(path = output_dir,pattern = "*.all2all", ignore.case = T,full.names = T),function(in_file){
  #   out_file <- paste(output_dir,tools::file_path_sans_ext(BiocGenerics::basename(in_file)),".out",sep="")
  #   convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
  # }, mc.cores = params_list$numWorkers )

  ## parallel::mclapply(list.files(path = output_dir,pattern = "*.out", ignore.case = T,full.names = T),function(in_file){
  #parallel::mclapply(list.files(path = output_dir,pattern = c(paste(clusters_left,clusters_right,"all2all",sep="."),paste(clusters_right,clusters_left,"all2all",sep=".")), ignore.case = T,full.names = T),function(in_file){ #c(paste(clusters_left,clusters_right,"all2all",sep="."),paste(clusters_right,clusters_left,"all2all",sep=".")) #list.files(path = output_dir,pattern = "*.all2all", ignore.case = T,full.names = T)

  group_combinations <- unique(tidyr::crossing(clusters_left,clusters_right))
  furrr::future_map2(.x=group_combinations$clusters_left,.y=group_combinations$clusters_right ,.f=function(c_left,c_right){
    in_file1 <- paste(output_dir,"/",c_left,".",c_right,".fw.all2all.gz",sep="")
    in_file2 <- paste(output_dir,"/",c_left,".",c_right,".bk.all2all.gz",sep="")
    out_file1 <- paste(output_dir,"/",tools::file_path_sans_ext(BiocGenerics::basename(in_file1)),".wis_out.gz",sep="")
    out_file2 <- paste(output_dir,"/",tools::file_path_sans_ext(BiocGenerics::basename(in_file2)),".wis_out.gz",sep="")
    furrr::future_map2(.x=c(in_file1,in_file2),.y=c(out_file1,out_file2),.f=function(in_file,out_file){
      #try(
      if(!file.exists(out_file) && file.exists(in_file) && file.info(in_file)$size >0 || clean_extract){
        blast_GO <- GRObject_from_BLAST(blast_input = in_file, COMPLETE.format.ids = T, col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17), params_list = params_list)

        #SELECT ONLY FRAMES 1/1
        blast_GO <- blast_GO[blast_GO$query.frame==1]
        blast_GO <- blast_GO[blast_GO$subject.frame==1]

        wis_GO <- invisible(run_WISARD(blast_hits = blast_GO,score_col = "Hsp_score",n_threads = ceiling(params_list$numWorkers/2), verbose = verbose)) #score_col=16)
        wis_GO <- melt_wisard_list(wis_GO)
        wis_GO$num_factors <- NULL
        if(verbose) {print(out_file)}
        #factor_cols <- which(sapply(seq_along(1:ncol(wis_GO)),function(x){is.factor(wis_GO[,x])}))
        #wis_GO[,factor_cols] <- data.frame(lapply(wis_GO[,factor_cols], as.character),stringsAsFactors=FALSE)
        #wis_GO <- data.frame(wis_GO, stringsAsFactors = FALSE)
        #data.table::fwrite(x = as.matrix(wis_GO),file = out_file,quote = F,col.names = T,row.names = F,sep = "\t", nThread = params_list$numWorkers, compress = "auto", verbose = verbose)
        write.table(x = wis_GO,file = gzfile(out_file,open = "w"),quote = F,col.names = T,row.names = F,sep = "\t")
      }
    },.options = furrr::furrr_options( seed = seed, scheduling=F)) #params_list$numWorkers
  },.options = furrr::furrr_options( seed = seed, scheduling=F)) #params_list$numWorkers

  # all2all_wis_out <- parallel::mclapply(all2all_GRObjects, function(all2all_blast_hit){
  #   ret_hit <- parallel::mclapply(all2all_blast_hit, function(blast_hit){
  #
  #     if(is.null(blast_hit)){
  #       return(NULL)
  #     }
  #     blast_hit <- unlist(blast_hit, recursive = T)
  #     blast_GO <- blast_hit$BLAST_hits
  #
  #     if(is.null(blast_GO) || length(blast_GO)==0){
  #       return(NULL)
  #     }
  #
  #     out_file <- paste(dirname(blast_hit$BLAST_file),"/",tools::file_path_sans_ext(blast_hit$run_name),".wis_out",sep="")
  #     blast_GO <- blast_GO[blast_GO$query.frame==1]
  #     blast_GO <- blast_GO[blast_GO$subject.frame==1]
  #
  #     wis_GO <- invisible(run_WISARD(blast_hits = blast_GO,score_col = "Hsp_score",n_threads = params_list$numWorkers, verbose = verbose)) #score_col=16)
  #     wis_GO <- melt_wisard_list(wis_GO)
  #     wis_GO$num_factors <- NULL
  #
  #     if(verbose) {print(out_file)}
  #     #factor_cols <- which(sapply(seq_along(1:ncol(wis_GO)),function(x){is.factor(wis_GO[,x])}))
  #     #wis_GO[,factor_cols] <- data.frame(lapply(wis_GO[,factor_cols], as.character),stringsAsFactors=FALSE)
  #     #data.table::fwrite(x = list(wis_GO),file = out_file,quote = F,col.names = T,row.names = F,sep = "\t", nThread = params_list$numWorkers)
  #     return(list(wis_out=wis_GO,run_name=blast_hit$run_name,BLAST_file=blast_hit$BLAST_file),seed=blast_hit$seed)
  #
  #   }, mc.cores = 2,mc.silent = !verbose, mc.set.seed = seed)
  #   return(ret_hit)
  #
  # }, mc.cores = ceiling(params_list$numWorkers/2),mc.silent = !verbose, mc.set.seed = seed)
  # #rm(all2all_GRObjects)
  # save(all2all_wis_out, file="wisard_results.RData")
  #load("files/all2all/wisard_results.RData")

  ##RUN RBH
  furrr::future_map2(.x=group_combinations$clusters_left,.y=group_combinations$clusters_right ,.f=function(query,subject){
    #lapply(clusters_left, function(query){
    #  parallel::mclapply(clusters_right,function(subject){ #grep(clusters_left,clusters_right,ignore.case = T,invert = T,value=T)
    in1 <- paste(output_dir,"/",query,".",subject,".fw.all2all.wis_out.gz",sep="")
    in2 <- paste(output_dir,"/",query,".",subject,".bk.all2all.wis_out.gz",sep="")
    out1 <- paste(output_dir,"/",query,".",subject,".fw.all2all.rbh_out.gz",sep="")
    out2 <- paste(output_dir,"/",query,".",subject,".bk.all2all.rbh_out.gz",sep="")
    #print(paste(in1,in2)) #DEBUG
    #print(paste(out1,out2)) #DEBUG
    try(
      if ((!file.exists(out1) && !file.exists(out2) || clean_extract) && file.exists(in1) && file.exists(in2) && file.info(in1)$size > 0 && file.info(in2)$size > 0) { #!file.exists(out1) && !file.exists(out2) &&
        RBH_out <- NULL
        try(RBH_out <- RBH(in1 = in1, in2 = in2, index.tables = T, col.indices = list(qseqid=12,sseqid=1), header = T,n_threads = params_list$numWorkers)) #list(qseqid=12,sseqid=1,weight.col=22) #col.names = c("subject_id","start","end","width","strand","Hsp_num","Hsp_bit.score","Hsp_score","Hsp_evalue","Hsp_query.from","Hsp_query.to","query_id","query_len","subject_len","Hsp_hit.from","Hsp_hit.to","Hsp_query.frame","Hsp_hit.frame","Hsp_pidentity","Hsp_gaps","Hsp_align.len","max_score")
        #print(paste(out1,out2))
        #RBH_out <- RBH(in1 = blast_hit[[1]]$wis_out, in2 = blast_hit[[2]]$wis_out, index.tables = T, col.indices = list(qseqid=12,sseqid=1), header = T,n_threads = params_list$numWorkers) #,weight.col=c(22,8) ), unique.hit.weights = T, process.weights.func = max) #RBH(in1 = in1, in2 = in2, index.tables = T, col.indices = list(qseqid=12,sseqid=1), header = T,n_threads = params_list$numWorkers) #,weight.col=c(22,8) ), unique.hit.weights = T, process.weights.func = max)

        #data.table::fwrite(x = list(RBH_out$in1), file = out1,quote = F,col.names = T,row.names = F, sep = "\t", sep2 = c("\n","\t","\n"), nThread = params_list$numWorkers,compress = "auto")
        #data.table::fwrite(x = list(RBH_out$in2), file = out2,quote = F,col.names = T,row.names = F, sep = "\t", sep2 = c("\n","\t","\n"), nThread = params_list$numWorkers,compress = "auto")
        if(all(!is.null(unlist(sapply(RBH_out,is.null))))){
          write.table(x = as.data.frame(RBH_out$in1), file = gzfile(out1, open = "w"),quote = F,col.names = T,row.names = F, sep = "\t")
          write.table(x = as.data.frame(RBH_out$in2), file = gzfile(out2, open="w"),quote = F,col.names = T,row.names = F, sep = "\t")
          #data.table::fwrite(x = list(RBH_out$in1), file = out1,quote = F,col.names = T,row.names = F, sep = "\t", nThread = params_list$numWorkers,compress = "auto", verbose = verbose)
          #data.table::fwrite(x = list(RBH_out$in2), file = out2,quote = F,col.names = T,row.names = F, sep = "\t", nThread = params_list$numWorkers,compress = "auto", verbose = verbose)
        }
        #save(RBH_out, tx_CDS_lengths, file = "tmp.RData")
        #calculate_HSP_coverage(RBH_out$in1,transcript_region_lengths = tx_CDS_lengths, col.indices=list(qseqid=12,sseqid=1,qstart=10,qend=11,sstart=2,send=3), COMPLETE.format.ids = T,params_list = params_list)
        #calculate_HSP_coverage(RBH_out$in2,transcript_region_lengths = tx_CDS_lengths, col.indices=list(qseqid=12,sseqid=1,qstart=10,qend=11,sstart=2,send=3), COMPLETE.format.ids = T,params_list = params_list)
      } #else{
      #   stop(paste(out1,"exists! (OR)",in1," & ",in2,"does not exist"))
      # }
    )
    # }, mc.cores = params_list$numWorkers )
    #})
  }, .options = furrr::furrr_options(seed = seed, scheduling=F)) #params_list$numWorkers

  # all2all_rbh_out <- parallel::mclapply(all2all_wis_out, function(blast_hit){
  #
  #   RBH_out <- RBH(in1 = blast_hit[[1]]$wis_out, in2 = blast_hit[[2]]$wis_out, index.tables = T, col.indices = list(qseqid=12,sseqid=1), header = T,n_threads = params_list$numWorkers) #,weight.col=c(22,8) ), unique.hit.weights = T, process.weights.func = max)
  #   #data.table::fwrite(x = list(RBH_out$in1),file = out1,quote = F,col.names = T,row.names = F, sep = "\t", nThread = params_list$numWorkers)
  #   #data.table::fwrite(x = list(RBH_out$in2),file = out2,quote = F,col.names = T,row.names = F, sep = "\t", nThread = params_list$numWorkers)
  #   r_name <- unique(unlist(lapply(blast_hit, function(x){
  #     y <- stringi::stri_split(str = x$run_name,fixed = ".",simplify=T)
  #     return(paste(y[,c(1,2)], collapse = "."))
  #   })))
  #
  #   return(list(rbh_out=RBH_out,run_name=paste(r_name,collapse = "-"),BLAST_files=c(blast_hit[[1]]$BLAST_file,blast_hit[[2]]$BLAST_file)))
  # }, mc.cores = params_list$numWorkers, mc.silent = !verbose, mc.set.seed = seed)
  # #rm(all2all_wis_out)
  # save(all2all_rbh_out, file="all2all_rbh_out.RData")

  furrr::future_map2(.x=group_combinations$clusters_left,.y=group_combinations$clusters_right ,.f=function(query,subject){
    #lapply(clusters_left, function(query){
    #  parallel::mclapply(clusters_right,function(subject){ #grep(clusters_left,clusters_right,ignore.case = T,invert = T,value=T)
    in1 <- paste(output_dir,query,".",subject,".fw.all2all.rbh_out.gz",sep="")
    in2 <- paste(output_dir,query,".",subject,".bk.all2all.rbh_out.gz",sep="")
    out1 <- paste(output_dir,query,".",subject,".all2all.final_out.gz",sep="")
    #out2 <- paste(output_dir,query,".",subject,".bk.final_out.gz",sep="")
    out_cov_data <- paste(output_dir,query,".",subject,".all2all.coverage.gz",sep="")
    run_name <-  paste(query,subject,sep=".")
    try(
      if (!file.exists(out_cov_data) && !file.exists(out1) && file.exists(in1) && file.exists(in2) && file.info(in1)$size > 0 && file.info(in2)$size > 0 || clean_extract ) {

        final_blast_table <- NULL

        #all2all_final_out <- parallel::mclapply(all2all_rbh_out, function(rbh_tup){
        try(final_blast_table <- calculate_HSP_coverage(fw_blast_table =in1,bk_blast_table = in2,col.indices=list(qseqid=12,sseqid=1,query_len=13,subject_len=14,align_len=23), group=run_name,COMPLETE.format.ids = T,params_list = params_list, header = T,verbose = verbose, min_coverage_filter = params_list$MIN_COVERAGE_THRESHOLD, run.mode = run.mode,n_threads=ceiling(params_list$numWorkers/2))) #transcript_region_lengths = tx_CDS_lengths
        if(!is.null(final_blast_table) && all(unlist(sapply(final_blast_table, nrow)) > 0)){
          write.table(x = final_blast_table$blast_table$BLAST_hits,file = gzfile(out1, open = "w"),quote = F,col.names = T,row.names = F, sep = "\t") #,nThread = params_list$numWorkers,compress = "auto", verbose = verbose)
          #data.table::fwrite(x = list(final_blast_table$blast_table[[2]]),file = out2_data,quote = F,col.names = T,row.names = F, sep = "\t",nThread = params_list$numWorkers,compress = "auto", verbose = verbose)
          write.table(x = final_blast_table$coverage, file = gzfile(out_cov_data, open = "w"),quote = F,col.names = T,row.names = F, sep = "\t") #,nThread = params_list$numWorkers,compress = "auto", verbose = verbose)
        }
        #return(list(coverage=final_blast_table, run_name=run_name, BLAST_files=c(in1,in2)))
      } ) #, mc.cores = params_list$numWorkers, mc.silent = !verbose, mc.set.seed = seed)
    #rm(all2all_rbh_out)
    #save(all2all_final_out, file="all2all_final_out.RData")
    #return(all2all_final_out)
  }, .options = furrr::furrr_options(seed = seed, scheduling=F)) #params_list$numWorkers
  #gc()

  # final_blast_tables <- furrr::future_map2(.x=group_combinations$clusters_left,.y=group_combinations$clusters_right ,.f=function(query,subject){
  #   in1_data <- paste(output_dir,query,".",subject,".fw.rbh_out",sep="")
  #   in2_data <- paste(output_dir,query,".",subject,".bk.rbh_out",sep="")
  #   out1_data <- paste(output_dir,query,".",subject,".fw.final_out",sep="")
  #   out2_data <- paste(output_dir,query,".",subject,".bk.final_out",sep="")
  #   out_cov_data <- paste(output_dir,query,".",subject,".coverage",sep="")
  #   try(if(file.exists(in1_data) && file.exists(in2_data)){ #&& !file.exists(out1_data) && !file.exists(out2_data) && !file.exists(out_cov_data)
  #     final_blast_table <- calculate_HSP_coverage(fw_blast_table = in1_data,bk_blast_table = in2_data,col.indices=list(qseqid=12,sseqid=1,query_len=13,subject_len=14,align_len=23), group=paste(query,subject,sep="."),COMPLETE.format.ids = T,params_list = params_list, header = T,verbose = T, min_coverage_filter = params_list$MIN_COVERAGE_THRESHOLD) #transcript_region_lengths = tx_CDS_lengths
  #     if(!is.null(unlist(sapply(final_blast_table, nrow))) && unlist(sapply(final_blast_table, nrow)) > 0){
  #       data.table::fwrite(x = list(final_blast_table$blast_table[[1]]),file = out1_data,quote = F,col.names = T,row.names = F, sep = "\t",nThread = params_list$numWorkers)
  #       data.table::fwrite(x = list(final_blast_table$blast_table[[2]]),file = out2_data,quote = F,col.names = T,row.names = F, sep = "\t",nThread = params_list$numWorkers)
  #       data.table::fwrite(x = list(final_blast_table$coverage,file), file = out_cov_data,quote = F,col.names = T,row.names = F, sep = "\t",nThread = params_list$numWorkers)
  #     }
  #     return(final_blast_table)
  #   } #else{
  #   #   stop(paste(out1_data,",",out2_data,"&",out_cov_data,"exists! (OR)",in1_data,",",in2_data, "does not exist!"))
  #   # }
  #   )
  # }, .options = furrr::furrr_options(seed = TRUE, scheduling=params_list$numWorkers))

  cat(print_toc(tictoc::toc(quiet = T, log = T)))

  #save(final_blast_tables, file="final_blast_tables.RData") #DEBUG

  # ##PLOT_CODE
  # tmp2 <- c()
  # tmp3 <- c()
  # tmp4 <- c()
  # tmp5 <- c()
  # tmp6 <- c()
  # plot_filename <- paste(ref_org,org,sep = "-")
  # tmp3 <- sapply(HSP_fw, USE.NAMES = T ,function(x){
  #   return(data.frame(query=x$query,subject=x$subject,min_cov=x$min_cov,max_cov=x$max_cov,cov_q=x$cov_q,cov_s=x$cov_s,same_CDS_count=as.logical(x$same_CDS_count),hit_from=x$hit_from,hit_to=x$hit_to,query_len=x$query_len,subject_len=x$subject_len,q_align_len=x$q_align_len,s_align_len=x$s_align_len,query_from=x$query_from,query_to=x$query_to, pident=x$pident))
  # })
  # tmp4 <- sapply(HSP_bk, USE.NAMES = T ,function(x){
  #   return(data.frame(query=x$query,subject=x$subject,min_cov=x$min_cov,max_cov=x$max_cov,cov_q=x$cov_q,cov_s=x$cov_s,same_CDS_count=as.logical(x$same_CDS_count),hit_from=x$hit_from,hit_to=x$hit_to,query_len=x$query_len,subject_len=x$subject_len,q_align_len=x$q_align_len,s_align_len=x$s_align_len,query_from=x$query_from,query_to=x$query_to, pident=x$pident))
  # })
  # #print(str(tmp3))
  # if(length(tmp3)>0){ ##IF length(tmp3) == 0 then there are no genes matching between organisms
  #   for(i in 1:ncol(tmp3)){
  #     #print(data.frame(tmp3[,i],colnames(tmp3)[i]))
  #     #print(colnames(tmp3)[i])
  #     #print(data.frame(t(tmp3[[i]]),colnames(tmp3)[i]))
  #     tmp2 <- rbind(data.frame(tmp3[,i],gene=colnames(tmp3)[i]),tmp2)
  #   }
  #   tmp2 <- tmp2 %>% mutate(ref_org=rep(ref_org)) %>% mutate(org=rep(org)) %>% mutate(direction=rep("forward"))
  #   if(length(tmp4)>0){
  #     for(i in 1:ncol(tmp4)){
  #       #print(data.frame(tmp3[,i],colnames(tmp3)[i]))
  #       #print(colnames(tmp3)[i])
  #       #print(data.frame(t(tmp3[[i]]),colnames(tmp3)[i]))
  #       tmp5 <- rbind(data.frame(tmp4[,i],gene=colnames(tmp4)[i]),tmp5)
  #     }
  #     tmp5 <- tmp5 %>% mutate(ref_org=rep(org)) %>% mutate(org=rep(ref_org)) %>% mutate(direction=rep("backward"))
  #     #tmp6 <- rbind(tmp2,tmp5)
  #     if(mean(tmp2$min_cov)>=mean(tmp5$min_cov)){
  #       tmp6 <- tmp2
  #     }else{
  #       tmp6 <- tmp5
  #     }}else{
  #       tmp6 <- tmp2
  #     }
  #   #print(tmp6[which(tmp6$min_cov > 1),])
  #   ##Remove outliers
  #   tmp6 <- tmp6[!is.na(tmp6$min_cov),]
  #   tmp6 <- tmp6[!is.na(tmp6$same_CDS_count),]
  #   #tmp6 <- tmp6[which(tmp6$min_cov <= 1),] #tmp6[!which(tmp6$min_cov > 1),]
  #   #tmp6 <- tmp6[which(tmp6$max_cov <= 1),]
  #   tmp6 <- unique(tmp6)
  #
  #   #print(head(tmp6))
  #   #ggplot(tmp6, aes(x=min_cov*100, fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion") + facet_wrap(~direction)
  #   if(nrow(tmp6)>0){
  #     #print(head(tmp6))
  #
  #     #print(tmp2[which(is.na(tmp2)),])
  #     #DENSITY PLOT FOR ALL GENES COLOURED BY SAME_CDS_COUNT
  #     density_plot <- ggplot(tmp6, aes(x=min_cov*100, fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion") + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))
  #     density_plot_max <- ggplot(tmp6, aes(x=max_cov*100, fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Max.Coverage of all transcripts(pairwise)") + ylab("Proportion") + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))
  #
  #     #GROUPED BY GENE
  #     genewise_plot <- ggplot(tmp6, aes(x=min_cov*100, group=gene ,fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion") + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))
  #     genewise_plot_max <- ggplot(tmp6, aes(x=max_cov*100, group=gene ,fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Max.Coverage of all transcripts(pairwise)") + ylab("Proportion") + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))
  #
  #     # facet_wrap(~gene, drop=FALSE, scales=c("free"))
  #     total_pages <- n_pages(genewise_plot + facet_wrap_paginate(facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y")))
  #     total_pages_max <- n_pages(genewise_plot_max + facet_wrap_paginate(facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y")))
  #     #GROUPED BY SAME_CDS_COUNT
  #     cds_count_plot <- ggplot(tmp6, aes(x=min_cov*100 ,group=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity') + xlab("Same_CDS_Count") + ylab("Proportion") + facet_wrap(~same_CDS_count, drop=FALSE, scales=c("free")) + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))
  #     cds_count_plot_max <- ggplot(tmp6, aes(x=max_cov*100 ,group=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity') + xlab("Same_CDS_Count") + ylab("Proportion") + facet_wrap(~same_CDS_count, drop=FALSE, scales=c("free")) + ggtitle(paste(gsub(pattern = delimiter,x = ref_org,replacement = " "),gsub(pattern = delimiter,x =org,replacement = " "),sep = "-"))
  #
  #     ##Do a boxplot with facet_wrap(~ref organisms)
  #     #print(head(tmp6))
  #     if(file.exists(paste(plot_out_path,paste(plot_filename,".pdf",sep=""),sep="/"))){
  #       file.remove(paste(plot_out_path,paste(plot_filename,".pdf",sep=""),sep="/"))
  #     }
  #
  #     pdf(file =paste(plot_out_path,paste(plot_filename,".pdf",sep=""),sep="/"),title = plot_filename)
  #     plot.new()
  #     text(.5, .5, "MINIMUM_COVERAGE")
  #     print(density_plot)
  #     #print(total_pages)
  #     if(total_pages > 0){
  #       for(i in 1:total_pages){
  #         try(print(genewise_plot + facet_wrap_paginate(facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y"), page=i))) #facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,scales=c("free_y"),
  #       }
  #     }
  #     #print(genewise_plot)
  #     print(cds_count_plot)
  #     plot.new()
  #     text(.5, .5, "MAXIMUM_COVERAGE")
  #     print(density_plot_max)
  #     #print(total_pages)
  #     if(total_pages_max > 0){
  #       for(i in 1:total_pages_max){
  #         try(print(genewise_plot_max + facet_wrap_paginate(facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y"), page=i))) #facets=c("gene","same_CDS_count"),nrow=2,ncol=2,shrink=F,scales=c("free_y"),
  #       }
  #     }
  #     #print(genewise_plot)
  #     print(cds_count_plot_max)
  #     plot.new()
  #     text(.5, .5, "BLAST Coverage Plots")
  #     for(row_q in unique(tmp6$query)){
  #       for(row_s in unique(tmp6$subject)){
  #         subset_data <- tmp6[which(tmp6$query==row_q & tmp6$subject==row_s),]
  #         #print(subset_data)
  #         if(nrow(subset_data)>0){
  #           subset_data <- subset_data %>% mutate(group=1:nrow(subset_data))
  #           line_coords <- pivot_longer(subset_data[,c("query","subject","hit_from","hit_to","query_from","query_to","pident","min_cov","group","max_cov")], cols = c("hit_from","hit_to","query_from","query_to") ,
  #                                       names_to = "direction", values_to = "coords")
  #           line_coords$direction[which(line_coords$direction=="hit_from" | line_coords$direction=="hit_to")] <- "hit"
  #           line_coords$direction[which(line_coords$direction=="query_from" | line_coords$direction=="query_to")] <- "query"
  #           data_rows=nrow(subset_data)
  #           #print(data_rows)
  #           line_coords <- data.frame(query=rep(unique(line_coords$query),2*data_rows),subject=rep(unique(line_coords$subject),2*data_rows),from=line_coords$coords[which(line_coords$direction=="query")],to=line_coords$coords[which(line_coords$direction=="hit")], groups=rep(1:data_rows,each=2),pident=rep(subset_data$pident,each=2),min_cov=rep(subset_data$min_cov,each=2),max_cov=rep(subset_data$max_cov,each=2))
  #           line_coords$groups <- factor(line_coords$groups)
  #           line_coords$min_cov <- as.numeric(line_coords$min_cov)
  #           line_coords$to <- as.numeric(line_coords$to)
  #           line_coords$from <- as.numeric(line_coords$from)
  #           plot_labels <- as.vector(t(apply(subset_data[,c("min_cov","max_cov")], MARGIN = c(1,2),FUN=function(x){
  #             return(round(x*100,2))
  #           })))
  #           try(print(ggplot(line_coords,aes(x=from,y=to,group=groups,color=pident)) + ylim(0,unique(subset_data$query_len)) + xlim(0,unique(subset_data$subject_len)) +
  #                       geom_line(na.rm = T) + geom_point(na.rm = T) + geom_text(aes(label=plot_labels)) + ylab(paste(unique(line_coords$query),"(",unique(subset_data$query_len),")")) + xlab(paste(unique(line_coords$subject),"(",unique(subset_data$subject_len),")"))))
  #           HSP <- rbind(HSP,subset_data)
  #         }
  #       }
  #     }
  #
  #     dev.off()
  #   }
  # }
  # ###PLOT-CODE OVER

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
  run_status <- processx::run( command = COMPLETE_vars$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"select_ref_org_groups",params_list$REF_ORGS_FILE, params_list$GROUPS_PATH,COMPLETE_vars$parallel, params_list$SELECT_REF_ORG_GROUPS_METHOD, params_list$numWorkers) ,spinner = T,stdout = "",stderr = "")
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

#@param fix.invalid.labels Fix labels of files which do not have COMPLETE.format.ids?. Default - !remove.invalid.files (TRUE).  (Refer ?COMPLETE_PIPELINE_DESIGN to know about COMPLETE.format.ids)
#' Group FASTA Sequences Based on a column index of COMPLETE.format.ids (?COMPLETE_PIPELINE_DESIGN)
#'
#' group_FASTA_genes() for run.more="gene" & group_FASTA_clusters() for run.mode="cluster" wrap this function. run.mode is used in FIND_TRANSCRIPT_ORTHOLOGS(). This function write groupings (clusters/genes) of each organism (org_name) into paste(params_list$OUT_PATH,"/genes/",org_name,"/ORG_CLUSTERS.", id.col.index,sep="")
#'
#' @param gene_list Vector/Filename of Gene names (FASTA Filenames where the sequences are stored)
#' @param params_list Output of load_params()
#' @param id.col.index The index of Column of COMPLETE.format.ids (?COMPLETE_PIPELINE_DESIGN) to groups sequences into. Use id.col.index=1 for grouping sequences based on Transcript IDs, id.col.index=2 to group sequences into Organisms, id.col.index=3 for grouping sequences based on Gene Names and id.col.index=4 to groups sequences into Ortholog Clusters. Check COMPLETE_vars$FORMAT_ID_INDEX for indices
#' @param remove.invalid.files Remove files which do not have COMPLETE.format.ids? Default - FALSE . (Refer ?COMPLETE_PIPELINE_DESIGN to know about COMPLETE.format.ids)
#' @param verbose Print DEBUG Messages?
#' @export
group_FASTA <- function(gene_list ,params_list, id.col.index, remove.invalid.files=F, verbose=T){ #fix.invalid.labels=!remove.invalid.files
  if(!any(grepl(x = class(params_list), pattern = "COMPLETE-options"))){
    stop("Error: params_list is not a COMPLETE-options class. Use load_params()")
  }

  if (length(gene_list) == 1 && file.exists(gene_list)) {
    genes <- tolower(gsub('[[:punct:]]+','_', factor(scan(gene_list, character(),quiet = T))))
    genes <- genes[grep("gene",genes, invert = T, fixed = T)]
  }else{
    genes <- tolower(as.vector(gene_list))
  }
  # if(remove.invalid.files==T && fix.invalid.labels==T){
  #   stop("Both remove.invalid.files and fix.invalid.labels cannot be TRUE")
  # }
  dir.create(params_list$GROUPS_PATH,showWarnings = F, recursive = T)
  id.col.index <- as.numeric(id.col.index)
  grouping_by <- names(COMPLETE_vars$FORMAT_ID_INDEX[id.col.index])
  tictoc::tic(msg = paste("Grouping FASTA into", grouping_by,"..."))
  unlink(x = params_list$GROUPS_PATH,recursive = T,force = T,expand = T)
  try(unlink(x = paste(params_list$GROUPS_PATH,"/../groups_noncds",sep=""),recursive = T,force = T,expand = T))
  dir.create(path = params_list$GROUPS_PATH,showWarnings = F,recursive = T)
  #fasta_files <- list.files(path = params_list$FASTA_OUT_PATH,pattern = genes,all.files = F,full.names = T,recursive = T,include.dirs = F)

  tmp_p1 <- processx::process$new(command = Sys.which("find"),args = c(params_list$FASTA_OUT_PATH),stdout = "|", supervise = T, cleanup = T)
  odb_file <- paste(params_list$OUT_PATH,"/genes/*/odb.final_map",sep="") #list.files(path =paste(params_list$OUT_PATH,"/genes/",sep=""), pattern = "odb.final_map", recursive = T,include.dirs = F,full.names = T) #processx::process$new(command = Sys.which("find"),args = c(paste(params_list$OUT_PATH,"/genes/*/odb.final_map",sep="")),stdout = "")
  tmp_p3 <- processx::process$new(command = COMPLETE_vars$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"get_all_odb_genes", gene_list,odb_file), stdout = "|",stderr=NULL)
  tmp_p3$wait()
  genes <- unique(c(genes,tolower(gsub('[[:punct:]]+','_', tmp_p3$read_all_output_lines()))))
  safe_genes_tmp <- tempfile(pattern="safe_names")
  data.table::fwrite(file = safe_genes_tmp,x = list(genes),col.names = F,row.names = F,nThread = params_list$numWorkers)
  #tmp_p1$wait()
  tmp_p2 <- processx::process$new(command = Sys.which("grep"), args=c("-i","-f", safe_genes_tmp), stdin= tmp_p1$get_output_connection(), stdout = "|")
  #tmp_p2$wait()
  fasta_files <- tmp_p2$read_all_output_lines()
  fasta_files <- grep(pattern = "*.fai", x = fasta_files,ignore.case = T,value = T,invert = T)
  #fasta_files <- c(fasta_files,tmp_p3$read_all_output())
  #print(fasta_files)

  # furrr::future_map(fasta_files, function(x){ #furrr::future_map #parallel::mclapply #lapply
  #   if(!file.exists(x) && file.info(x)$size ==0){
  #     return()
  #   }
  #   fasta_recs <- Biostrings::readDNAStringSet(filepath = x,use.names = T, format = "fasta")
  #   #print(names(fasta_recs)) #DEBUG
  #   split_recs <- stringi::stri_split(str = names(fasta_recs), fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)
  #   tryCatch({
  #     org_name <- unique( split_recs[, COMPLETE_vars$FORMAT_ID_INDEX$ORG] )
  #     #gene_name <- unique( split_recs[, COMPLETE_vars$FORMAT_ID_INDEX$GENE] )
  #     #odb_clusters <- unique( split_recs[, COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS] )
  #   }, error=function(cond){
  #     org_name <- basename(dirname(x))
  #     if(verbose){
  #       message(paste("Storing output in ",org_name))
  #     }
  #   })
  #   #if(is.null(org_name)){
  #   #  org_name <- basename(dirname(x))
  #   #}
  #
  #   if(ncol(split_recs)!=length(COMPLETE_vars$FORMAT_ID_INDEX)){ ##CHECKING IF FASTA IDs are COMPLETE.format.ids
  #     #print(x) #DEBUG
  #     if(verbose){
  #       message(paste(x," : does not have COMPLETE.format.ids"))
  #     }
  #     if(remove.invalid.files){
  #       unlink(x = x,force = T, expand = T)
  #       return(NULL)
  #     }
  #
  #     # if(fix.invalid.labels){
  #     #   tryCatch({
  #     #     label_sequenceIDs(fasta_path = x, org = org_name, gene = gene_name, odb_clusters = odb_clusters, params_list = params_list)
  #     #     fasta_recs <- Biostrings::readDNAStringSet(filepath = x,use.names = T, format = "fasta")
  #     #   }, error=function(cond){
  #     #     if(verbose){
  #     #       message(paste("Cannot label",x,". Missing information in the FASTA IDs!"))
  #     #     }
  #     #   })
  #     # }
  #
  #   }
  #
  #   #print(org_name) #DEBUG
  #   all_clusters <- unique(unlist(furrr::future_map2(seq_along(fasta_recs),names(fasta_recs), function(rec_num, rec_name){
  #     split_rec <- stringi::stri_split(str = rec_name, fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)
  #     rec_clusters <- unique(stringi::stri_split(str = split_rec[, id.col.index], fixed = ",", simplify=T)) #ncol(split_rec)
  #     #print(rec_clusters) #DEBUG
  #     furrr::future_map(rec_clusters, function(each_cluster){ #furrr::future_map #parallel::mclapply
  #       #print(fasta_recs[rec_num]) #DEBUG
  #       #print(each_cluster) #DEBUG
  #       #print( paste(params_list$GROUPS_PATH,"/",each_cluster,".",tools::file_ext(x),sep = ""),)  #DEBUG
  #       if(!stringi::stri_isempty(each_cluster)){
  #         Biostrings::writeXStringSet(x = fasta_recs[rec_num],filepath = paste(params_list$GROUPS_PATH,"/",each_cluster,".",tools::file_ext(x),sep = ""), append = T,format = "fasta")
  #       }
  #     }, .options = furrr::furrr_options(seed = TRUE, scheduling=params_list$numWorkers)) #, mc.cores = params_list$numWorkers,mc.preschedule = T,mc.silent = !verbose)
  #     return(rec_clusters)
  #   }, .options = furrr::furrr_options(seed = TRUE, scheduling=params_list$numWorkers))))
  #   data.table::fwrite(x = list(all_clusters), file = paste(params_list$OUT_PATH,"/genes/",org_name,"/ORG_CLUSTERS.", grouping_by,sep=""),quote = F,row.names = F,col.names = F, nThread = params_list$numWorkers)
  #
  # }, .options = furrr::furrr_options(seed = TRUE, scheduling=params_list$numWorkers)) #, mc.core = params_list$numWorkers,mc.preschedule = T,mc.silent = !verbose)
  # cat(print_toc(tictoc::toc(quiet = T)))

  furrr::future_map(fasta_files, function(x){
    processx::run( command = COMPLETE_vars$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"group_FASTA_seqs",COMPLETE_vars$parallel,x,params_list$GROUPS_PATH,params_list$SEQUENCE_ID_DELIM,params_list$numWorkers,params_list$OUT_PATH,id.col.index,params_list$FASTA_OUT_PATH ) ,spinner = T,stdout = NULL,stderr = NULL)
  }, .options = furrr::furrr_options(seed = TRUE, scheduling=params_list$numWorkers))
  # all_groups_list <- parallel::mclapply(list.files(path = paste(params_list$OUT_PATH,"/genes/",sep=""),include.dirs=TRUE, full.names=TRUE),function(x){
  #   if(file.exists(paste(x,"/ORG_CLUSTERS.",grouping_by,sep="")) && file.info(paste(x,"/ORG_CLUSTERS.",grouping_by,sep=""))$size > 0 ){
  #     return(scan(paste(x,"/ORG_CLUSTERS.",grouping_by,sep=""), character(), quiet = T))
  #   }
  # }, mc.cores =  params_list$numWorkers, mc.preschedule = T, mc.silent = !verbose)
  all_groups <- unique(basename(tools::file_path_sans_ext(list.files(path=params_list$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F)))) #unique(unlist(all_groups_list,recursive = T)) #purrr::reduce(all_groups_list, union)
  data.table::fwrite(x = list(all_groups), file = paste(params_list$OUT_PATH,"/ALL_GROUPS.txt",sep=""), quote = F, row.names = F,col.names = F,na = "-", append = F, nThread = params_list$numWorkers)
  return(all_groups)
}

#@param gene_list Vector or File with a list of genes to extract data for(check the github repo for an example).
#' (2) - Find Transcript Orthologs
#'
#' This function can be executed after COMPLETE::EXTRACT_DATA() and is the continuation of R-COMPLETE pipeline
#'
#' This is the main function which calls all the other functions and performs and end-end execution of finding transcript level orthologs. It runs the iterative Transcript Ortholog Extraction pipeline which is used to reduce the pool of genes, reduce the pool of organisms, find transcript level orthologs (check ?extract_transcript_orthologs)
#'
#' @note ONLY USE THIS FUNCTION WHEN RUNNING THE PIPELINE OF R-COMPLETE. Use other helper function to work with custom BLAST files not generated by this R package. run.mode="gene" is NOT RECOMMENDED because the sequences are grouped based on gene names
#'
#' @param gene_list Vector/Filename of gene names. (essentially Filenames of FASTA sequences)
#' @param params_list Filename of a formatted parameter file (check the github repo for an example) or Output of load_params().
#' @param blast_program Give the name of the BLAST program to use (if in $PATH) or give the absolute path to the BLAST program. BLAST options are taken from params_list. Default is Sys.which("tblastx")
#' @param run.mode A value from COMPLETE_vars$FORMAT_ID_INDEX. Default - COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS. Find transcript orthologs in the level of Orgs, Genes or Ortholog Clusters. Genes have more tight orthology and fewer transcript orthologs which may be very similar. Ortholog Clusters are a level higher than Genes (Because an Ortholog Cluster can have more than one gene) and have highest number of transcript orthologs with a lot of dissimilarity. run.mode=COMPLETE_vars$FORMAT_ID_INDEX$GENE is NOT RECOMMENDED because the sequences are grouped based on gene names, while run.mode=COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS groups sequences based on protein identity.
#' @param verbose Print DEBUG Messages?
#' @param seed Seed value for reproducibility. Default - 123
#' @export
FIND_TRANSCRIPT_ORTHOLOGS <- function(gene_list, params_list, blast_program=Sys.which("tblastx"), run.mode=COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS, verbose=F, seed=123){ #gene_list
  set.seed(seed)

  if(stringi::stri_isempty(blast_program)){
    stop("BLAST+ not found in $PATH. Provide blast_program")
  }
  else{
    BLAST_BIN <- dirname(blast_program)
  }

  if(is.na(match(run.mode,COMPLETE_vars$FORMAT_ID_INDEX)) || is.null(run.mode)) {
    stop(paste("run.mode must be one of COMPLETE_vars$FORMAT_ID_INDEX"))
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

  tictoc::tic(msg="Finding transcript orthologs :")

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

    all_gtf_stats <- NULL
    tryCatch({
       all_gtf_stats <- data.table::fread(file = paste(loaded_PARAMS$OUT_PATH,"/all_gtf_stats.csv",sep=""),sep = ",",header = T,quote = "",fill = T,na.strings = "-", nThread = loaded_PARAMS$numWorkers)
     }, error= function(cond){
       message(paste("Coerced GTF stats file", paste(loaded_PARAMS$OUT_PATH,"/all_gtf_stats.csv",sep=""),"not found/invalid. Plots will not be generated. Please rerun EXTRACT_DATA()"))
       warning(paste("Coerced GTF stats file", paste(loaded_PARAMS$OUT_PATH,"/all_gtf_stats.csv",sep=""),"not found/invalid. Plots will not be generated. Please rerun EXTRACT_DATA()"))
     })

    # available_genes_list <- parallel::mclapply(paste(loaded_PARAMS$OUT_PATH,"/genes/",REF_ORGS,sep=""),function(x){
    #   if(file.exists(paste(x,"/AVAILABLE_GENES",sep="")) && file.info(paste(x,"/AVAILABLE_GENES",sep=""))$size > 0 ){
    #     return(scan(paste(x,"/AVAILABLE_GENES",sep=""), character(),quiet = T))
    #   }
    # }, mc.cores =  loaded_PARAMS$numWorkers, mc.preschedule = T)
    # available_genes <- unique(purrr::reduce(available_genes_list, union))

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
      all_clusters <- unique(all_clusters)
      files_in_dir <- list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F)
      clusters_in_dir <- unique(basename(tools::file_path_sans_ext(files_in_dir))) #unique(stringi::stri_split(str = list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = F,recursive = F,include.dirs = F), simplify=T, fixed = ".")[,1])
      #print(head(all_clusters)) #DEBUG
      #print(head(files_in_dir)) #DEBUG
      #print(head(clusters_in_dir)) #DEBUG
      #print(match(all_clusters,clusters_in_dir)) #DEBUG
      if( any(is.na(match(all_clusters,clusters_in_dir))) || length(clusters_in_dir) == 0 || length(all_clusters) == 0 ){ #all_clusters[which(!is.na(match(all_clusters,clusters_in_dir)))] #|| length(which(!is.na(match(all_clusters,clusters_in_dir)))) != length(which(!is.na(match(clusters_in_dir,all_clusters)))) #any(is.na(match(all_clusters[which(!is.na(match(all_clusters,clusters_in_dir)))] ,clusters_in_dir[which(!is.na(match(clusters_in_dir,all_clusters)))])))
        stop("Regrouping sequences...\n")
      }
    },error=function(cond){
      message(cond) #paste(cond,":",loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""))
      tictoc::tic(msg = paste("Storing sequences into groups (",names(COMPLETE_vars$FORMAT_ID_INDEX)[as.numeric(run.mode)],") ..."))
      #print("here1")
      unlink(paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""),force = T,expand = T)
      all_clusters <- group_FASTA(gene_list=gene_list, params_list = loaded_PARAMS, id.col.index = as.numeric(run.mode),remove.invalid.files=T,verbose = T) # fix.invalid.labels = T
      #print(all_clusters) #DEBUG
      cat(print_toc(tictoc::toc(quiet = T, log=T)))
      #write.table(x = all_clusters,file = paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), quote = F, row.names = F,col.names = F,na = "-", append = F)
    }, finally = {
      #print(all_clusters) #DEBUG
      #unlink(x = grep(x = list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F) , pattern = "cds", ignore.case = T,invert = T, value = T), recursive = F,force = T,expand = T)
      files_in_dir <- list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F)
      non_cds_file_list <- grep(x = files_in_dir , pattern = "cds", ignore.case = T,invert = T, value = T)
      dir.create(path = file.path(loaded_PARAMS$GROUPS_PATH,"/../groups_noncds"), showWarnings = F,recursive = T)
      if(length(non_cds_file_list) > 0){
        file.rename(from = non_cds_file_list,to = paste(loaded_PARAMS$GROUPS_PATH,"/../groups_noncds/",basename(non_cds_file_list),sep = ""))
        files_in_dir <- list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F)
        data.table::fwrite(x = list(unique(tools::file_path_sans_ext(basename(files_in_dir)))),file = paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), quote = F, row.names = F,col.names = F,na = "-", append = F, nThread = loaded_PARAMS$numWorkers)
      }
      available_clusters <- unique(tools::file_path_sans_ext(basename(files_in_dir)))#all_clusters
    })
    #}

    #available_clusters <- scan(paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), character(),quiet = T) #unique((tools::file_path_sans_ext(list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = F,recursive = F,include.dirs = F))))
    if(length(available_clusters)==0){
      #available_clusters <- "ungrouped"
      stop("No clusters were found!. Try other values for run.mode")
    }

    #STEP 1 - ONLY for run.mode="cluster" - Place ungrouped sequences into groups (all2allblast BLAST ungrouped cluster againts all clusters) ##EG - ungrouped.1013114at2759.fw.all2all, 1013114at2759.ungrouped.bk.all2all
    #all2allblast and then wisard and then RBH for grouping ungrouped clusters
    if(any(grepl(pattern = "ungrouped",x = available_clusters,ignore.case = T))){
      message("STEP 1 - Placing ungrouped sequences into groups\n")
      tictoc::tic(msg = "Placing ungrouped sequences into groups...")

      parallel::mclapply(available_clusters, function(x){
        extract_transcript_orthologs(blast_program = blast_program, blast_options =  loaded_PARAMS$BLAST_OPTIONS, params_list = loaded_PARAMS, transcript_region = ".cds",clusters_left = "ungrouped",clusters_right = x, input_dir = loaded_PARAMS$GROUPS_PATH,output_dir = all2all_out,clean_extract = loaded_PARAMS$CLEAN_EXTRACT,verbose = verbose, seed=seed, run.mode="both")
      }, mc.cores = ceiling(loaded_PARAMS$numWorkers/2), mc.silent = !verbose, mc.set.seed = seed) #loaded_PARAMS$numWorkers

      if(length(available_clusters)>1){
    parallel::mclapply(list.files(path = all2all_out, pattern = ".final_out.gz", ignore.case = T,include.dirs = F,full.names = T), function(x){
      final_hits <- LoadBLASTHits(infile = x,header = T)
      #final_hits <- split(final_hits,f=factor(final_hits$direction))
      #lapply(final_hits, function(y){
      #  id_clusters_all <- unique(stringi::stri_split(str = unique(c(y$query_id,y$subject_id)), fixed = loaded_PARAMS$SEQUENCE_ID_DELIM, simplify=T))
      id_clusters_all <- unique(stringi::stri_split(str = unique(c(final_hits$query_id,final_hits$subject_id)), fixed = loaded_PARAMS$SEQUENCE_ID_DELIM, simplify=T))
        id_clusters_sub <- NULL
         id_clusters_sub <- unique(purrr::reduce(stringi::stri_split(str = grep(pattern = "ungrouped", ignore.case = T, x=id_clusters_all, value = T, invert = T), fixed=",", simplify=F),intersect))
         if(is.null(id_clusters_sub)){
           id_clusters_sub <- unique(purrr::reduce(stringi::stri_split(str = grep(pattern = "ungrouped", ignore.case = T, x=id_clusters_all, value = T, invert = T), fixed=",", simplify=F),union))
         }
         if(is.null(id_clusters_sub)){
           return(NULL)
         }
        ids_nogrp <- grep(pattern = "ungrouped", x = unique(c(y$query_id,y$subject_id)),ignore.case = T,value = T)
        #ids_grp <- grep(pattern = paste(id_clusters_sub,collapse ="|"), x = unique(c(y$query_id,y$subject_id)),ignore.case = T,value = T)

          split_id <- stringi::stri_split(str = ids_nogrp, fixed = loaded_PARAMS$SEQUENCE_ID_DELIM, simplify=T)
          old_id <- split_id
          in_file <- paste(loaded_PARAMS$FASTA_OUT_PATH,"/",split_id[,COMPLETE_vars$FORMAT_ID_INDEX$ORG],"/",gsub('[[:punct:]]+','_', split_id[,COMPLETE_vars$FORMAT_ID_INDEX$GENE]), sep="")
          split_id[,COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS] <- paste(paste("-",id_clusters_sub,sep=""),collapse = ",") #paste("(,",paste(id_clusters_sub,collapse = ","),",)",sep="")
          #new_id <- paste(split_id,sep=loaded_PARAMS$SEQUENCE_ID_DELIM)
          invisible( lapply(seq_along(1:nrow(split_id)), function(idx){
            parallel::mclapply(loaded_PARAMS$TRANSCRIPT_REGIONS, function(reg){
              if(stringi::stri_cmp_eq(paste(old_id[idx,],collapse=loaded_PARAMS$SEQUENCE_ID_DELIM), paste(split_id[idx,],collapse=loaded_PARAMS$SEQUENCE_ID_DELIM))){
                return(NULL)
              }
              processx::run( command = COMPLETE_vars$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"sed_replace", paste(in_file[idx],".",reg,sep=""), paste(old_id[idx,],collapse=loaded_PARAMS$SEQUENCE_ID_DELIM), paste(split_id[idx,],collapse=loaded_PARAMS$SEQUENCE_ID_DELIM) ) ,spinner = T,stdout = NULL,stderr = NULL)
            }, mc.cores = 3)
          })) #, mc.cores = loaded_PARAMS$numWorkers)

      #})
    })
      }else{#create custom clusters for the ungrouped sequences
      #  parallel::mclapply(list.files(path = all2all_out, pattern = ".final_out.gz", ignore.case = T,include.dirs = F,full.names = T), function(x){
          final_hits <- LoadBLASTHits(infile = paste(all2all_out,"/","ungrouped.ungrouped.all2all.final_out.gz",sep=""),header = T)
          final_coverage <- read.table(file = paste(all2all_out,"/","ungrouped.ungrouped.all2all.coverage.gz",sep=""),header = T)
          final_coverage <- final_coverage[final_coverage$min>loaded_PARAMS$MIN_COVERAGE_THRESHOLD,]
          final_coverage <- final_coverage %>% mutate(cluster=kmeans(final_coverage$min,centers = ceiling(sd(final_coverage$min)))$cluster)
          #final_coverage <- split(final_coverage,f=factor(final_coverage$cluster))
            id_clusters_q <- stringi::stri_split(str = final_coverage$query_id, fixed = loaded_PARAMS$SEQUENCE_ID_DELIM, simplify=T)
            id_clusters_s <- stringi::stri_split(str = final_coverage$subject_id, fixed = loaded_PARAMS$SEQUENCE_ID_DELIM, simplify=T)
            id_clusters_q[,COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS] <- final_coverage$cluster
            id_clusters_q <- unique(id_clusters_q)
            id_clusters_s[,COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS] <- final_coverage$cluster
            id_clusters_s <- unique(id_clusters_s)
            id_clusters_q <- data.frame(seq_id=paste(id_clusters_q[,COMPLETE_vars$FORMAT_ID_INDEX$TRANSCRIPT_ID],id_clusters_q[,COMPLETE_vars$FORMAT_ID_INDEX$ORG],id_clusters_q[,COMPLETE_vars$FORMAT_ID_INDEX$GENE],sep=loaded_PARAMS$SEQUENCE_ID_DELIM),cluster=id_clusters_q[,COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS])
            id_clusters_s <- data.frame(seq_id=paste(id_clusters_s[,COMPLETE_vars$FORMAT_ID_INDEX$TRANSCRIPT_ID],id_clusters_s[,COMPLETE_vars$FORMAT_ID_INDEX$ORG],id_clusters_s[,COMPLETE_vars$FORMAT_ID_INDEX$GENE],sep=loaded_PARAMS$SEQUENCE_ID_DELIM),cluster=id_clusters_s[,COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS])
            id_clusters_q <- unique(id_clusters_q %>% group_by(seq_id) %>% mutate(cluster=paste("-cluster",cluster,sep="",collapse=",")))
            id_clusters_s <- unique(id_clusters_s %>% group_by(seq_id) %>% mutate(cluster=paste("-cluster",cluster,sep="",collapse=",")))
            id_clusters_q <- unlist(apply(id_clusters_q, MARGIN = 1, FUN=function(x){
              return(paste(x,collapse = loaded_PARAMS$SEQUENCE_ID_DELIM,sep=loaded_PARAMS$SEQUENCE_ID_DELIM))
            }))
            id_clusters_s <- unlist(apply(id_clusters_s, MARGIN = 1, FUN=function(x){
              return(paste(x,collapse = loaded_PARAMS$SEQUENCE_ID_DELIM,sep=loaded_PARAMS$SEQUENCE_ID_DELIM))
            }))
            id_clusters_all <- unique(stringi::stri_split(str = unique(c(id_clusters_q,id_clusters_s)), fixed = loaded_PARAMS$SEQUENCE_ID_DELIM, simplify=T))
            split_id <- id_clusters_all
            old_id <- id_clusters_all
            old_id[,COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS] <- "ungrouped"
            in_file <- paste(loaded_PARAMS$FASTA_OUT_PATH,"/",split_id[,COMPLETE_vars$FORMAT_ID_INDEX$ORG],"/",gsub('[[:punct:]]+','_', split_id[,COMPLETE_vars$FORMAT_ID_INDEX$GENE]), sep="")
            split_id[,COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS] <- paste(paste("-",id_clusters_sub,sep=""),collapse = ",") #paste("(,",paste(id_clusters_sub,collapse = ","),",)",sep="")
            #new_id <- paste(split_id,sep=loaded_PARAMS$SEQUENCE_ID_DELIM)
            invisible( lapply(seq_along(1:nrow(split_id)), function(idx){
              parallel::mclapply(loaded_PARAMS$TRANSCRIPT_REGIONS, function(reg){
                if(stringi::stri_cmp_eq(paste(old_id[idx,],collapse=loaded_PARAMS$SEQUENCE_ID_DELIM), paste(split_id[idx,],collapse=loaded_PARAMS$SEQUENCE_ID_DELIM))){
                  return(NULL)
                }
                processx::run( command = COMPLETE_vars$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"sed_replace", paste(in_file[idx],".",reg,sep=""), paste(old_id[idx,],collapse=loaded_PARAMS$SEQUENCE_ID_DELIM), paste(split_id[idx,],collapse=loaded_PARAMS$SEQUENCE_ID_DELIM) ) ,spinner = T,stdout = NULL,stderr = NULL)
              }, mc.cores = 3)
            })) #, mc.cores = loaded_PARAMS$numWorkers)

          #})
         #})
      }
    all_clusters <- group_FASTA(gene_list=gene_list, params_list = loaded_PARAMS, id.col.index = as.numeric(run.mode),remove.invalid.files=T,verbose = T)

    files_in_dir <- list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F)
    non_cds_file_list <- grep(x = files_in_dir , pattern = "cds", ignore.case = T,invert = T, value = T)
    dir.create(path = file.path(loaded_PARAMS$GROUPS_PATH,"/../groups_noncds"), showWarnings = F,recursive = T)
    if(length(non_cds_file_list) > 0){
      file.rename(from = non_cds_file_list,to = paste(loaded_PARAMS$GROUPS_PATH,"/../groups_noncds/",basename(non_cds_file_list),sep = ""))
      files_in_dir <- list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F)
      data.table::fwrite(x = list(unique(tools::file_path_sans_ext(basename(files_in_dir)))),file = paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), quote = F, row.names = F,col.names = F,na = "-", append = F, nThread = loaded_PARAMS$numWorkers)
    }
    available_clusters <- unique(tools::file_path_sans_ext(basename(files_in_dir)))#all_clusters

      cat(print_toc(tictoc::toc(quiet = T, log=T)))
    }else{
      message("STEP 1 - Skipped because all sequences are grouped\n") # or run.mode != COMPLETE_vars$FORMAT_ID_INDEX$CLUSTERS
    }
    # }else{
    #   message("STEP 1 - Skipped because run.mode='gene'\n")
    #   all_genes <- group_FASTA(params_list = loaded_PARAMS, id.col.index = COMPLETE_vars$FORMAT_ID_INDEX$GENE)
    #   available_clusters <- list.files(path=loaded_PARAMS$GROUPS_PATH,full.names = T,recursive = F,include.dirs = F)
    # }

    #     if(loaded_PARAMS$SELECT_REF_ORG_GROUPS){
    #       all_clusters <- select_ref_org_groups(loaded_PARAMS)
    #       write.table(x = all_clusters,file = paste(loaded_PARAMS$OUT_PATH,"/ALL_GROUPS.txt",sep=""), quote = F, row.names = F,col.names = F,na = "-")
    #     }
    #     available_clusters <- all_clusters

    #available_cluster_combinations <- unique(tidyr::crossing(available_clusters, available_clusters)) #combination is not necessary because we only want orthologous transcript sequences from within the groups and not across
    #STEP 2 - Select transcript level orthologs with minimum coverage between clusters/genes
    message(paste("STEP 2 - Select transcript level orthologs between orgs/genes/clusters, based on minimum coverage\n"))
    tictoc::tic(msg = "Extracting Transcript Orthologs...")
    parallel::mclapply(available_clusters, function(x){
      extract_transcript_orthologs(blast_program = blast_program, blast_options =  loaded_PARAMS$BLAST_OPTIONS, params_list = loaded_PARAMS, transcript_region = "cds", clusters_left = x,clusters_right = x, input_dir = loaded_PARAMS$GROUPS_PATH,output_dir = all2allfinal_out, clean_extract = loaded_PARAMS$CLEAN_EXTRACT, verbose = verbose, seed=seed, run.mode="both")
    }, mc.cores = ceiling(loaded_PARAMS$numWorkers), mc.silent = !verbose, mc.set.seed = seed) #loaded_PARAMS$numWorkers
    cat(print_toc(tictoc::toc(quiet = T,log=T)))

    cat(paste("Run-Time Log:\n"))
    cat(paste(tictoc::tic.log(),collapse = "\n"))

    if(!is.null(all_gtf_stats) && nrow(all_gtf_stats) > 0){
     # parallel::mclapply(available_clusters, function(x){
        #load *.coverage and *.fw.final_out, *.bk.final_out from all2allfinal_out
        tmp_cove$coverage <- tmp_cove$coverage %>% mutate(query_org=unlist(purrr::map(tmp_cove$coverage[,1],function(x){
          #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
          org <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
          #print(org)
          return(org[,COMPLETE_vars$FORMAT_ID_INDEX$ORG])
        })) )
        tmp_cove$coverage <- tmp_cove$coverage %>% mutate(subject_org=unlist(purrr::map(tmp_cove$coverage[,2],function(x){
          #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
          org <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
          return(org[,COMPLETE_vars$FORMAT_ID_INDEX$ORG])
        })) )
        tmp_cove$coverage <- tmp_cove$coverage %>% mutate(query_transcript_id=unlist(purrr::map(tmp_cove$coverage[,1],function(x){
          #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
          tx_id <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
          return(stringi::stri_split_fixed(tx_id[,COMPLETE_vars$FORMAT_ID_INDEX$TRANSCRIPT_ID],pattern = params_list$TRANSCRIPT_ID_DELIM, simplify=T)[,1])
        })) )
        tmp_cove$coverage <- tmp_cove$coverage %>% mutate(subject_transcript_id=unlist(purrr::map(tmp_cove$coverage[,2],function(x){
          #genes <- unlist(stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM))
          tx_id <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
          return(stringi::stri_split_fixed(tx_id[,COMPLETE_vars$FORMAT_ID_INDEX$TRANSCRIPT_ID],pattern = params_list$TRANSCRIPT_ID_DELIM, simplify=T)[,1])
        })) )
        query_gene=unlist(purrr::map(tmp_cove$coverage[,1],function(x){
          gene <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
          return(gene[,COMPLETE_vars$FORMAT_ID_INDEX$GENE])
        }))
        subject_gene=unlist(purrr::map(tmp_cove$coverage[,2],function(x){
          gene <- stringi::stri_split_fixed(x,pattern = params_list$SEQUENCE_ID_DELIM, simplify=T)
          return(gene[,COMPLETE_vars$FORMAT_ID_INDEX$GENE])
        }))
        tmp_cove$coverage$gene <- NA
        tmp_cove$coverage$gene[query_gene==subject_gene] <- factor(query_gene[query_gene==subject_gene])
        tmp_cove$blast_table <- split(tmp_cove$blast_table$BLAST_hits,f=tmp_cove$blast_table$BLAST_hits$direction)
        tmp_cove$blast_table[[1]] <- dplyr::full_join(x = tmp_cove$coverage,y = tmp_cove$blast_table[[1]]) #, by=c("query_id"="query_id","subject_id"="subject_id","query_transcript_id"="query_transcript_id","subject_transcript_id"="subject_transcript_id", "query_len"="query_len", "subject_len"="subject_len", "query_org"="query_org", "subject_org"="subject_org","direction"="direction"))
        tmp_cove$blast_table[[2]] <- dplyr::full_join(x = tmp_cove$coverage,y = tmp_cove$blast_table[[2]]) #, by=c("query_id"="query_id","subject_id"="subject_id", "query_len"="query_len", "subject_len"="subject_len", "query_org"="query_org", "subject_org"="subject_org"))
        tmp_cove$coverage <- tmp_cove$coverage[tmp_cove$coverage$query_transcript_id!=tmp_cove$coverage$subject_transcript_id,]
        tmp_cove$coverage <- tmp_cove$coverage[tmp_cove$coverage$query_org!=tmp_cove$coverage$subject_org,]
        tmp_cove$coverage <- tmp_cove$coverage %>% mutate(query_cds_count=all_gtf_stats[match(tmp_cove$coverage$query_transcript_id,all_gtf_stats$transcript_id),c("cds_count")])
        tmp_cove$coverage <- tmp_cove$coverage %>% mutate(subject_cds_count=all_gtf_stats[match(tmp_cove$coverage$subject_transcript_id,all_gtf_stats$transcript_id),c("cds_count")])
        tmp_cove$coverage <- tmp_cove$coverage %>% mutate(same_cds_count=query_cds_count==subject_cds_count)
        #tmp_cove$coverage <- tmp_cove$coverage[complete.cases(tmp_cove$coverage),] #IF you do not want any NAs

        #PLOTTING
        density_plot <- ggplot(tmp_cove$coverage, aes(x=tmp_cove$coverage$min_raw_cov*100, fill=tmp_cove$coverage$same_cds_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion") #+ ggtitle(paste("Min. Coverage between ",x,"-",x,sep = ""))
        density_plot_max <- ggplot(tmp_cove$coverage, aes(x=tmp_cove$coverage$max_raw_cov*100, fill=tmp_cove$coverage$same_cds_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Max.Coverage of all transcripts(pairwise)") + ylab("Proportion") #+ ggtitle(paste("Max. Coverage between ",x,"-",x,sep = ""))
        if(length(unique(tmp_cove$coverage$gene)) > 1){
          genewise_plot <- ggplot(tmp_cove$coverage, aes(x=tmp_cove$coverage$min_raw_cov*100, group= tmp_cove$coverage$gene ,fill= tmp_cove$coverage$same_cds_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion")# + ggtitle(paste("Min. Coverage between ",x,"-",x," (grouped by gene)",sep = ""))
          genewise_plot_max <- ggplot(tmp_cove$coverage, aes(x=tmp_cove$coverage$max*100, group= tmp_cove$coverage$gene ,fill= tmp_cove$coverage$same_cds_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Max.Coverage of all transcripts(pairwise)") + ylab("Proportion") #+ ggtitle(paste("Max. Coverage between ",x,"-",x," (grouped by gene)",sep = ""))
          total_pages <- n_pages(genewise_plot + facet_wrap_paginate(facets=c("gene","same_cds_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y")))
          total_pages_max <- n_pages(genewise_plot_max + facet_wrap_paginate(facets=c("gene","same_cds_count"),nrow=2,ncol=2,shrink=F,drop=F,scales=c("free_y")))
        }

      #  hit_combinations <- unique(tidyr::crossing(unique(tmp_cove$coverage$query_id),unique(tmp_cove$coverage$subject_id)))
      #  furrr::future_map2( .x=hit_combinations[,1],.y=hit_combinations[,2],.f=function(row_q,row_s){
      #    print(row_q,row_s)   #DEBUG
      #    if(!stringi::stri_cmp_eq(row_q,row_s)){

            subset_data1 <- tmp_cove$blast_table[[1]] #tmp_cove$blast_table[[1]][which(!is.na(match(tmp_cove$blast_table[[1]]$query_id,row_q)) & !is.na(match(tmp_cove$blast_table[[1]]$subject_id,row_s))),] %>% mutate(blast_dir="forward") #dplyr::bind_rows( , tmp_cove$blast_table[[1]][which(!is.na(match(tmp_cove$blast_table[[1]]$subject,row_q)) & !is.na(match(tmp_cove$blast_table[[1]]$query,row_s))),]

            subset_data2 <- tmp_cove$blast_table[[2]] #tmp_cove$blast_table[[2]][which(!is.na(match(tmp_cove$blast_table[[2]]$query_id,row_q)) & !is.na(match(tmp_cove$blast_table[[2]]$subject_id,row_s))),] %>% mutate(blast_dir="backward") #dplyr::bind_rows( , tmp_cove$blast_table[[2]][which(!is.na(match(tmp_cove$blast_table[[2]]$subject,row_q)) & !is.na(match(tmp_cove$blast_table[[2]]$query,row_s))),]

            subset_data <- dplyr::bind_rows(subset_data1,subset_data2)
            b_table <- subset_data
            #subset_data_plots <- parallel::mclapply(list(subset_data1,subset_data2),function(b_table){
              if(nrow(b_table) > 0){

                b_table <- b_table %>% mutate(group=1:nrow(b_table))

                b_table <- b_table[order(b_table$start),]
                b_table1_GO <- GenomicRanges::makeGRangesFromDataFrame(b_table,keep.extra.columns = T,seqnames.field = c("subject_id"), start.field = "start",end.field = "end",na.rm=TRUE)
                b_table$start <- NULL
                b_table$end <- NULL
                b_table2_GO <- GenomicRanges::makeGRangesFromDataFrame(b_table,keep.extra.columns = T,seqnames.field = c("query_id"), start.field = "query_HSP_from",end.field = "query_HSP_to",na.rm=TRUE)

                line_coords <- tidyr::pivot_longer(subset_data[,c("query_id","subject_id","start","end","query_HSP_from","query_HSP_to","pidentity","min","group","max")], cols = c("start","end","query_HSP_from","query_HSP_to") ,
                                                   names_to = "position", values_to = "coords")
                line_coords$position[which(line_coords$position=="start" | line_coords$position=="end")] <- "subject"
                line_coords$position[which(line_coords$position=="query_HSP_from" | line_coords$position=="query_HSP_to")] <- "query"
                data_rows=nrow(subset_data)
                #print(data_rows)
                line_coords <- unique(data.frame(query=rep(line_coords$query_id,each=2),subject=rep(line_coords$subject_id,each=2),from=line_coords$coords[which(line_coords$position=="query")],to=line_coords$coords[which(line_coords$position=="subject")], raw_cov_q=rep(subset_data$raw_cov_q,each=2), raw_cov_s=rep(subset_data$raw_cov_s,each=2),pident=rep(subset_data$pidentity,each=2),query_len=rep(subset_data$query_len,each=2),subject_len=rep(subset_data$subject_len,each=2),hit_src=line_coords$position,min_cov=rep(subset_data$min_raw_cov,each=2),max_cov=rep(subset_data$max_raw_cov,each=2))) #groups=rep(1:data_rows,each=2),
                line_coords_grp <- tidyr::crossing(unique(line_coords$query),unique(line_coords$subject) )
                line_coords_grp <- line_coords_grp[line_coords_grp[,1]!=line_coords_grp[,2],]
                line_coords_grp <- line_coords_grp %>% mutate(groups=1:nrow(line_coords_grp))
                colnames(line_coords_grp) <- c("query","subject","groups")
                line_coords <- dplyr::inner_join(line_coords,line_coords_grp)
                line_coords$groups <- factor(line_coords$groups)
                #line_coords$min <- as.numeric(line_coords$min)
                line_coords$to <- as.numeric(line_coords$to)
                line_coords$from <- as.numeric(line_coords$from)
                # plot_labels <- as.vector(t(apply(line_coords[,c("min_cov","max_cov")], MARGIN = c(1,2),FUN=function(x){
                #   return(round(x*100,2))
                # })))

                line_coords_plots <- lapply(split(line_coords,line_coords$groups),function(line_coords_split){
                  #line_coords_split <- line_coords_split %>% mutate(plot_labels=NA)
                  #print(which.min(line_coords_split$from)) #DEBUG
                  #print(which.max(line_coords_split$to)) #DEBUG
                  #line_coords_split$plot_labels[which.min(line_coords_split$from)] <- round(line_coords_split$min_cov[which.min(line_coords_split$from)] * 100)
                  #line_coords_split$plot_labels[which.max(line_coords_split$to)] <- round(line_coords_split$max_cov[which.max(line_coords_split$to)] * 100)
                 # print(line_coords_split)
                  line_coords_split <- unique(line_coords_split[complete.cases(line_coords_split),])
                  if(nrow(line_coords_split)>0){
                  line_coords_split$groups <- rep(1:(nrow(line_coords_split)/2), each=2)
                  x <-  ggplot(line_coords_split,aes(x=from,y=to,group=groups,color=pident)) +geom_line(na.rm = T) + geom_point(aes(shape=factor(hit_src)),na.rm = T) + geom_abline(color=c("red")) + ylim(0,max(line_coords_split$subject_len)) + xlim(0,max(line_coords_split$query_len)) + xlab(unique(line_coords_split$query)) + ylab(unique(line_coords_split$subject))
                  return(x)
                  }
                 # return(line_coords_split)
                })

                #return(ggplot(line_coords,aes(x=from,y=to,group=groups,color=pident)) +geom_line(na.rm = T) + geom_point(na.rm = T)) #+ geom_label(aes(label=plot_labels, group=factor(groups)), size=3.6, angle=45) ) #=as.numeric(plot_labels), fill# + ylim(0,unique(subset_data$query_len)) + xlim(0,unique(subset_data$subject_len)) +
                #+ ylab(paste(unique(line_coords$query),"(",unique(subset_data[[1]]$query_len),")")) + xlab(paste(unique(line_coords$subject),"(",unique(subset_data[[1]]$subject_len),")"))))
              }
            #}, mc.cores = 2)

            #print(subset_data)

          }
        #}, .options = furrr::furrr_options(seed = TRUE, scheduling=loaded_PARAMS$numWorkers))

    #  }, mc.cores = loaded_PARAMS$numWorkers)

    #calculate gene conservation - calculate_gene_conservation.R - probably not needed
    ##Maybe write one for cluster conservation/coverage across organisms

  }else{
    stop(paste(blast_program," NOT found (blast_program). Give the right path to blast_program."))
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

#@param sep2 Delimiter 2 of the BLAST File columns. Default - c("","|",""). Check ?data.table::fwrite or ?data.table::fread
#' Find Reciprocal Blast Hits (RBH)
#'
#' Find RBH between BLAST results of different organisms/genes/transcripts (FASTA/FASTQ). The BLAST results must be of the format 6 and can be converted from BLAST format 11 with convert_BLAST_format().
#' The command with the required column names are given below.
#'
#' convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
#'
#' Optional : You can provide the file/table with indexed Transcsript IDs. The format must be "file"[tab]"long_id"[tab]"index" (without a header). It can be generated with the  index_FASTA_IDs() (check ?index_FASTA_IDs or index_fastaIDs() in fs::path_package("COMPLETE","exec","functions.sh"))
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
RBH <- function(in1,in2,sep="\t",header=F, transcript_ID_metadata=NULL, col.names=NULL, index.tables=T,col.indices, unique.hit.weights=F, process.weights.func=max, n_threads=tryCatch(parallel::detectCores(all.tests = T, logical = T), error=function(cond){return(2)})){ #sep2=c("","|","")

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

  #in_data <- parallel::mclapply(list(in1,in2), function(in_file){
  in_data <- furrr::future_map(list(in1,in2), function(in_file){
    #in_file_data <- NULL
    if(!is.null(in_file) && grepl(pattern = "tbl|df|data.frame", x = class(in_file),ignore.case = T)){
      in_file_data <- in_file
      if(!is.null(col.names)){
        colnames(in_file_data) <- col.names
      }
    }else if(!is.null(in_file) && length(in_file)==1 && grepl(pattern = "character", x = class(in_file),ignore.case = T)){
      if(file.exists(in_file)){
        in_file_data <-  try(LoadBLASTHits(infile = in_file, transcript_ID_metadata = transcript_ID_metadata, col.names = col.names, sep=sep, header=header))
      }else{
        message(paste(in_file,"does not exist!"))
      }
    }
    return(in_file_data)
  }, .options = furrr::furrr_options(seed = TRUE, scheduling=F)) #n_threads#,mc.cores = 2,mc.silent = F)
  if(any(unlist(lapply(in_data, is.null)))){
    stop("RBH() - Either in1 or in2 does not exist or has invalid format")
  }
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

  in1_g <- data.frame(from=in_data[[1]][, col.indices[["qseqid"]] ], to=in_data[[1]][,col.indices[["sseqid"]]], stringsAsFactors = T)
  in1_g <- in1_g[apply(in1_g, MARGIN=1,FUN=function(x){return(!any(is.na(x)))}),]

  in2_g <- data.frame(from=in_data[[2]][,col.indices[["qseqid"]]], to=in_data[[2]][,col.indices[["sseqid"]]], stringsAsFactors = T)
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
        furrr::future_map2(.x=hit_combinations$in1_valid_hits,.y=hit_combinations$in2_valid_hits,.f= function(x,y){
          #parallel::mclapply(in1_valid_hits, function(x){
          #return( parallel::mclapply(in2_valid_hits, function(y){
          if(!stringi::stri_cmp_eq(x,y)){
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
        } , .options = furrr::furrr_options(seed = TRUE, scheduling=F)) ) ) #n_threads
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
    #in_data[[1]] <- in_data[[1]][in1_RBH_rows,]
    #in_data[[2]] <- in_data[[2]][in2_RBH_rows,]
  }
  if(index.tables){
    in_data[[1]] <- deindex_BLAST_table(in_data[[1]], col.indices[["qseqid"]])
    in_data[[1]] <- deindex_BLAST_table(in_data[[1]], col.indices[["sseqid"]])
    in_data[[2]] <- deindex_BLAST_table(in_data[[2]], col.indices[["qseqid"]])
    in_data[[2]] <- deindex_BLAST_table(in_data[[2]], col.indices[["sseqid"]])
  }

  return(list(in1=in_data[[1]][in1_RBH_rows,],in2=in_data[[2]][in2_RBH_rows,]))
}

