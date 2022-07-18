#.libPaths("/data/meyer/viz/tools/miniconda3/envs/local_root2/lib/R/library")

library(stringi)
library(GenomicRanges)
library(purrr)
library(dplyr)
library(checkmate)
library(Rgb)
library(ggplot2)
require(parallel)
library(ggforce)
library(tidyverse)

###FUNCTIONS
# function to combine two data frames with unequal set of columns. missing
# columns in one frame are filled with NA. empty data frames allowed
get_param <- function(file, param){
  file_in <- read.delim(file,header = F,sep = c("="))
  return(file_in[grep(pattern = param,x = file_in$V1, ignore.case = T),][2])
}

rbind.all.columns <- function(x, y) {
  if ((ncol(x) > 0) & (ncol(y) > 0)) {
    x.diff <- setdiff(colnames(x), colnames(y))
    y.diff <- setdiff(colnames(y), colnames(x))
    x[, c(as.character(y.diff))] <- NA
    y[, c(as.character(x.diff))] <- NA
    return(rbind(x, y))
  } else {
    return(rbind(x, y))
  }
}

swap_columns <- function(x, col_x, col_y, swap_names=T){
  if(!is.null(x[col_x]) & !is.null(x[col_y])){
    tmp <- x[col_x]
    x[col_x] <- x[col_y]
    x[col_y] <- tmp
    if(swap_names){
      x_index <- which(names(x)==col_x)
      y_index <- which(names(x)==col_y)
      colnames(x)[x_index] <- col_y 
      colnames(x)[y_index] <- col_x
    }
    return(x)
  }else{
    print("invalid dataframes/col_names")
  }
}

merge_wis <- function(x, y){
  new_list <- list()
  new_list$alignments <- c(x$alignments, y$alignments)
  #new_list$alignments$max_score <- c(rep(x$max_score, length(x$alignments)),rep(y$max_score, length(y$alignments)))
  new_list$max_score <- sum(x$max_score, y$max_score)
  #print(new_list)
  return(unlist(new_list))
}

create_GRObject <- function(blast_output){
  tmp_blast <- blast_output
  col_names <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand", "qlen", "slen", "qseq","sseq","nident","positive")
  names(tmp_blast) <- col_names
  
  tmp_blast$subject_gene <- unlist(map(tmp_blast$sseqid,function(x){
    genes <- unlist(stri_split_fixed(x,pattern = seqID_delimiter))
    return(paste(genes[2],genes[4],sep=seqID_delimiter))
  }))
  
  tmp_blast$query_gene <- unlist(map(tmp_blast$qseqid,function(x){
    genes <- unlist(stri_split_fixed(x,pattern = seqID_delimiter))
    return(paste(genes[2],genes[4],sep=seqID_delimiter))
  }))
  
  #gene_check <- map2(tmp_blast$subject_gene,tmp_blast$query_gene, function(x,y){
  # print(x)
  #  print(unlist(stri_split_fixed(x,pattern = delimiter))[2])
  #  print(unlist(stri_split_fixed(y,pattern = delimiter)[2]))
  #  return(unlist(stri_split_fixed(x,pattern = delimiter)[2])==unlist(stri_split_fixed(y,pattern = delimiter)[2]))
  #})
  
  #tmp_blast <- tmp_blast[which(gene_check == T),]
  tmp_blast <- tmp_blast[which(tmp_blast$evalue < e_cutoff),] 
  ##ADDING FILTERING OPTIONS HERE
  #filter_param <- get_param(param_file,"filter")
  
  #gene_filter <- map2(tmp_blast$qseqid,tmp_blast$sseqid, function(x,y){
  # g_s <- unlist(stri_split_fixed(x,pattern = "::"))[4]
  #  g_a <- unlist(stri_split_fixed(y,pattern = "::"))[4]
  # switch (tolower(filter_param),
  #  "strict" = return(g_s==g_a),
  #  "moderate" = grepl(g_s,g_a,ignore.case = T),
  #  "none" = return(g_a==g_a),
  #  )
  #})
  
  #tmp_blast <- tmp_blast[which(gene_filter == T),]
  
  ##CHANGE strands
  change_sstrand <- which(tmp_blast$sstart > tmp_blast$send)
  if (length(change_sstrand) > 0) {
    # print("changing sstart")
    tmp <- tmp_blast[change_sstrand,]$sstart
    tmp_blast[change_sstrand,]$sstart <- tmp_blast[change_sstrand,]$send
    tmp_blast[change_sstrand,]$send <- tmp
    tmp <- NULL
  }
  change_qstrand <- which(tmp_blast$qstart > tmp_blast$qend)
  if (length(change_qstrand) > 0) {
    # print("changing qstart")
    tmp <- tmp_blast[change_qstrand,]$qstart
    tmp_blast[change_qstrand,]$qstart <- tmp_blast[change_qstrand,]$qend
    tmp_blast[change_qstrand,]$qend <- tmp
    tmp <- NULL
  }
  
  
  frame_check <- unlist(map(tmp_blast$frames,function(x){
    frames <- as.integer(unlist(stri_split_fixed(x,pattern = "/")))
    ##ONLY FORWARD STRANDS are taken
    #return(frames[1]==frames[2] && (frames[1]==1 || (frames[1]%%3==0 && frames[1]>=0 )) && (frames[2]==1 || (frames[2]%%3==0 && frames[2]>=0)))
    return(frames[1]==frames[2] && (frames[1]==1 && frames[2]==1))
  }))
  
  tmp_blast <- tmp_blast[which(frame_check == T),]
  
  tmp_gr <- GRanges(Rle(tmp_blast$sseqid), ranges =  IRanges(tmp_blast$sstart, end = tmp_blast$send)) #, names = orths$sseqid))
  
  tmp_gr$Hsp_num <- c(1:nrow(tmp_blast)) #seq(1,nrow(tmp_blast),1)
  tmp_gr$Hsp_bit.score <- tmp_blast$bitscore
  tmp_gr$Hsp_score <- tmp_blast$qcovhsp
  tmp_gr$Hsp_evalue <- tmp_blast$evalue
  tmp_gr$Hsp_query.from <- tmp_blast$qstart
  tmp_gr$Hsp_query.to <- tmp_blast$qend
  tmp_gr$query_id <- tmp_blast$qseqid
  tmp_gr$query_len <- tmp_blast$qlen
  tmp_gr$subject_len <- tmp_blast$slen
  tmp_gr$Hsp_hit.from <- tmp_blast$sstart
  tmp_gr$Hsp_hit.to <- tmp_blast$send
  tmp_gr$Hsp_query.frame <-  unlist(map(tmp_blast$frames,function(x){
    frames <- as.integer(unlist(stri_split_fixed(x,pattern = "/")))
    return(frames[1])
  }))
  tmp_gr$Hsp_hit.frame <-  unlist(map(tmp_blast$frames,function(x){
    frames <- as.integer(unlist(stri_split_fixed(x,pattern = "/")))
    return(frames[2])
  }))
  tmp_gr$query_org <- unlist(map(tmp_blast$qseqid,function(x){
    org <- unlist(stri_split_fixed(x,pattern = seqID_delimiter))
    return(org[3])
  }))
  tmp_gr$subject_org <- unlist(map(tmp_blast$sseqid,function(x){
    org <- unlist(stri_split_fixed(x,pattern = seqID_delimiter))
    return(org[3])
  }))
  tmp_gr$Hsp_pidentity <-tmp_blast$pident
  #gr$Hsp_positive <- 
  tmp_gr$Hsp_gaps <- tmp_blast$gaps
  tmp_gr$Hsp_align.len <-  tmp_blast$length
  tmp_gr$subject_gene <-  tmp_blast$subject_gene
  tmp_gr$query_gene <-  tmp_blast$query_gene
  return(tmp_gr)
}

run_WISARD <- function(query_id,subject_ids, gr,score_col){
  result_list <- c()
  for(sub_id in subject_ids){
  result <- get_wis(gr[which(gr$query_id == query_id & seqnames(gr) == sub_id),],max_score = score_col,overlap = 0) 
  if(result$max_score > 0 && length(result$alignments) > 0){
    #print(result)
    #g_name <- unique(result$alignments$subject_gene)
    #all_results[[g_name]] <- merge_wis(all_results[[g_name]], result)
    #return(merge_wis(all_results[[g_name]], result))
    #print(result)
    #g_name <- unique(unlist(stri_split_fixed(result$alignments$subject_gene,pattern=seqID_delimiter,n = 1,tokens_only = T)))
    #print(g_name)
    #if(!is.null(g_name)){
      #if(is.null(result_list)){
      #  result_list[[g_name]] <- result
      #}else{
      #  result_list[[g_name]] <- merge_wis(result_list[[g_name]], result)
      #}
    #}
    result_list <- c(result_list, result)
    #print(result_list)
  }
  #print(result_list)
  
  }
  return(result_list)
}

run_WISARD_gene <- function(subject_gene, gr,score_col){
  result <- get_wis(gr[which(gr$subject_gene==subject_gene),],max_score = score_col, overlap = 0) 
  if(result$max_score > 0 && length(result$alignments)>0){
    #print(result)
    #g_name <- unique(result$alignments$subject_gene)
    #all_results[[g_name]] <- merge_wis(all_results[[g_name]], result)
    #return(merge_wis(all_results[[g_name]], result))
    return(result)
  }
}

exec_WISARD <- function(gr, score_col){
  all_results <- list()
  #for (sseq in unique(gr$query_id)){   
  #  result <- get_wis(gr[which(sseq == gr$query_id),],max_score = "Hsp_score", overlap = 0) 
  #  if(result$max_score > 0){
  #    #print(result)
  #    g_name <- unique(result$alignments$subject_gene)
  #    all_results[[g_name]] <- merge_wis(all_results[[g_name]], result)
  #  }
  #}
  n_cores <- detectCores(all.tests = TRUE, logical = TRUE)
  if(is.na(n_cores)){
    n_cores=2
  }
  
  query_vector <- unique(gr$query_id)
  #print(runValue(seqnames(gr)))
  subject_vector <- unique(runValue(seqnames(gr)))
  child_results <- mcmapply(FUN=run_WISARD, query_vector, MoreArgs=list(subject_ids=subject_vector,gr=gr,score_col=score_col) ,mc.cores = n_cores-1, SIMPLIFY = FALSE,USE.NAMES = F)
  for (child in child_results) {
    #print(child)
  #  #g_name <- unique(child$alignments$subject_gene)
    g_name <- unique(unlist(stri_split_fixed(child$alignments$subject_gene,pattern=seqID_delimiter,n = 1,tokens_only = T)))
    #g_name <- names(child)
    print(g_name)
    
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
  #return(all_results)
  return(all_results)
  #return(tmp_cr)
}

exec_WISARD_gene <- function(gr,score_col){
  all_results <- list()
  
  n_cores <- try(detectCores(all.tests = TRUE, logical = TRUE))
  if(is.na(n_cores)){
    n_cores=2
  }
  subject_gene <- unique(gr$subject_gene)
  child_results <- mcmapply(FUN=run_WISARD_gene, subject_gene, MoreArgs=list(gr=gr,score_col=score_col) ,mc.cores = n_cores-1, SIMPLIFY = FALSE)
  for (child in child_results) {
    #g_name <- unique(child$alignments$subject_gene)
    g_name <- unique(unlist(stri_split_fixed(child$alignments$subject_gene,pattern=seqID_delimiter,n = 1,tokens_only = T)))
    if(!is.null(g_name)){
      if(is.null(all_results)){
        all_results[[g_name]] <- child
      }else{
        all_results[[g_name]] <- merge_wis(all_results[[g_name]], child)
      }
    }
  }
  return(all_results)
}

###CODE BASED ON IM's INPUT
feature_length <- function(transcript_id,gtf_data,feature="CDS"){
  ###gtf_data SHOULD BE READ USING Rgb::read.gtf
  ###With the followiing options, 
  ###attr = c("split"),features = c("CDS")
  #gtf_sel <- gtf_data[which(gtf_data$transcript_id == transcript_id & gtf_data$feature == feature),]
  #print(transcript_id)
  gtf_sel <- gtf_data[which(stri_detect(fixed = transcript_id,str = gtf_data$attributes)),]
  gtf_sel <- gtf_sel[which(gtf_sel$feature == feature & gtf_sel$frame==0),]
  #print(gtf_sel[order(gtf_sel$start),])
  
  return(sum(gtf_sel[order(gtf_sel$start),]$end - gtf_sel[order(gtf_sel$start),]$start))
}

feature_count <- function(transcript_id,gtf_data,feature="CDS"){
  ###gtf_data SHOULD BE READ USING Rgb::read.gtf
  ###With the followiing options, 
  ###attr = c("split"),features = c("CDS")
  
  #gtf_sel <- gtf_data[which(gtf_data$transcript_id == transcript_id & gtf_data$feature == feature),]
  gtf_sel <- gtf_data[which(stri_detect(fixed = transcript_id,str = gtf_data$attributes)),]
  gtf_sel <- gtf_sel[which(gtf_sel$feature == feature),]
  #print(transcript_id)
  #print(gtf_sel)
  return(nrow(gtf_sel))
  #return(max(gtf_sel$exon_number))
}

dissolve_GR_Overlaps <- function(GR_object){
  ##Essentially checks each GR element within the GR object with others to find and dissolve overlaps(only the longest overlapping item is kept)
  overlaps <- findOverlapPairs(GR_object, type="within",select = "all")
  all_overlaps <- findOverlapPairs(GR_object, type="any",select = "all")
  all_overlaps <- all_overlaps[-which(all_overlaps@first==overlaps@first & all_overlaps@second==overlaps@second)]
  all_overlaps <- all_overlaps[-which(all_overlaps@first==all_overlaps@second)]
  #overlaps <- overlaps[-which(overlaps@first==overlaps@second)]
  #overlaps <- unlist(zipup(overlaps))
  
  map2(ranges(tmp_overlaps),ranges(tmp_overlaps),function(x,y){
    if(x$start < y$start & x$end > y$end){
      return(x$end-x$start)
    }
    if(x$start > y$start & x$end < y$end){
      return(y$end-y$start)
    }
    if(x$start > y$start & x$end > y$end){
      overhang
    } 
  })
  
  discard_list <- c()
  overhang_list <- c()
  for (i in 1:length(GR_object)) {
    ##Find type of overlap (fully overlapped, (left/right) overhang)
    left_range <- ranges(GR_object[i]@first)
    right_range <- ranges(GR_object[i]@second)
    
    if(left_range$start < right_range$start & left_range$end > right_range$end){
      discard_list <- c(discard_list, GR_object[i]@second) #right engulfed, discarding it
    }
    if(left_range$start > right_range$start & left_range$end < right_range$end){
      discard_list <- c(discard_list, GR_object[i]@first) #left engulfed
    }
    if(left_range$start > right_range$start & left_range$end > right_range$end){
      overhang
    } 
    
    #print(GR_object[which(GR_object$start <= element$start & GR_object$start <= element$end),])
  }
}

#HSP_Coverage(query_mRNA = query_mRNA, org_fw = q_fw,org_bk=s_fw, s_gtf=ref_gtf, q_gtf=org_gtf)
HSP_Coverage <- function(query, org_fw, org_bk, gtf_stats, feature="CDS"){
  # org_fw <- q_fw
  # org_bk <- s_fw
  tmp_df <- c()
  #for(query in unique(query_mRNA)){
    #print(query)
    q_result <- org_fw[which(org_fw$query_id == query),]
    #print(q_result)
    ##q_result <- q_result[which(q_result$Hsp_query.frame == 1),]
    if(length(q_result) > 0){
      id_query <- unlist(stri_split_fixed(query,pattern = delimiter))[1]
      q_CDS_length <- unique(q_result$query_len) #feature_length(id_query, q_gtf, feature=feature)
      subject_mRNA <- unique(unique(seqnames(org_fw)@values),unique(seqnames(org_bk)@values))
      #subject_mRNA <- unique(seqnames(org_fw)@values)
      #print(subject_mRNA)
      #print(length(q_result))
      
      for(subject in subject_mRNA){
        result <- org_bk[which(org_bk$query_id == subject),]
        subject_result <- result[which(seqnames(result) == query)]
        #subject_result <- subject_result[which(subject_result$Hsp_hit.frame==1),]
        query_result <- q_result[intersect(which(q_result$query_id == query),which(seqnames(q_result) == subject)),]
        
        if(length(query_result) > 0 && length(subject_result) > 0){
          id_subject <- unlist(stri_split_fixed(subject,pattern = delimiter))[1]
          if(length(subject_result) > 0){
            s_CDS_length <- unique(subject_result$query_len) #feature_length(id_subject, s_gtf, feature=feature)
            ##OLD code in the comment to the right
            s_align_length <- sum(subject_result$Hsp_hit.to - subject_result$Hsp_hit.from) #sum(subject_result$Hsp_align.len)
            q_align_length <- sum(query_result$Hsp_hit.to - query_result$Hsp_hit.from) #sum(query_result$Hsp_align.len)
            #seqlengths(query_result)[which(names(seqlengths(query_result))==query)]
            #print(q_align_length)
            #print(s_align_length)
            #print(q_CDS_length)
            #print(s_CDS_length)
            #s_overlaps <- dissolve_GR_Overlaps(subject_result)
            cov_q <- q_align_length/q_CDS_length
            cov_s <- s_align_length/s_CDS_length
            print(paste(paste(query,"(",cov_q,")",sep = ""),paste(subject,"(",cov_s,")",sep = ""),sep="->"),)
            print(paste(query,"[","q_align_len/q_CDS_length:",q_align_length,"/",q_CDS_length,"]",sep=""))
            print(paste(subject,"[","s_align_len/s_CDS_length:",s_align_length,"/",s_CDS_length,"]",sep=""))
            #print(subject)
            #print(cov_q)
            #print(cov_s)
            feature_comp <- gtf_stats[gtf_stats$rna_id==id_subject,]$cds_count == gtf_stats[gtf_stats$rna_id==id_query,]$cds_count #(feature_count(id_query,q_gtf)==feature_count(id_subject,s_gtf))
            print(feature_comp)
            tmp_df <- rbind(tmp_df, data.frame(query=query, subject=subject,min_cov=min(cov_q, cov_s),max_cov=max(cov_q, cov_s), same_CDS_count=feature_comp,cov_q=cov_q, cov_s=cov_s,hit_from=subject_result$Hsp_hit.from,hit_to=subject_result$Hsp_hit.to,query_len=q_CDS_length,subject_len=s_CDS_length,q_align_len=q_align_length,s_align_len=s_align_length,query_from=subject_result$Hsp_query.from,query_to=subject_result$Hsp_query.to, pident=subject_result$Hsp_pidentity))
          }
        }
      }
    
    }
  #}
  return(tmp_df)
}

calc_delta_cov <- function(query_mRNA,org_bk_result,org_fw_result){
  for(query in unique(query_mRNA)){
    #print(query)
    query_result <- org_fw_result[which(org_fw_result$query_id == query),]
    #query_result <- org_fw_result[which(org_fw_result$query_id == query_mRNA),]
    subject_mRNA <- unique(seqnames(org_fw_result)@values)
    for(subject in subject_mRNA){
      ##FOR NOW OVERLAPs are not considered, eg
      ##[1] rna36069_cds(+)::tpm3::xenopus_tropicalis   337-591
      ##[2] rna36069_cds(+)::tpm3::xenopus_tropicalis    14-349
      ##[3] rna36069_cds(+)::tpm3::xenopus_tropicalis     4-141
      ## in the above sequences, results 2 & 3 overlap, which would technically create a bias when calculating delta_cov
      ## for now I am not doing anything to remove them. I should ask Rob if he can do anythin about overlap within results
      #print(subject)
      subject_result <- org_bk_result[which(org_bk_result$query_id == subject),]
      result_fw <- query_result[intersect(which(query_result$query_id == query),which(seqnames(query_result)@values == subject)),]
      result_bk <- subject_result[intersect(which(subject_result$query_id == subject),which(seqnames(subject_result)@values == query)),]
      
      id_query <- unlist(stri_split_fixed(query,pattern = delimiter))[1]
      id_subject <- unlist(stri_split_fixed(subject,pattern = delimiter))[1]
      
      feature_comp <- (feature_count(id_query,ref_gtf)==feature_count(id_subject,org_gtf))
      cov_fw <- sum(result_fw$Hsp_score) #sum(result_fw$Hsp_score)/sum(result_fw$Hsp_align.len)
      cov_bk <- sum(result_bk$Hsp_score) #sum(result_bk$Hsp_score)/sum(result_bk$Hsp_align.len)
      del_cov <-abs(cov_fw-cov_bk) #^2 #sum(result_fw$Hsp_align.len)/CDS_length(subject,org_gtf) #sum(result_fw$Hsp_score)/sum(result_fw$Hsp_align.len)
      #print(result_fw)
      #print(result_bk)
      #print(paste(cov_fw,cov_bk,del_cov))
      if(!is.na(del_cov) && cov_bk!=0 && cov_fw!=0){
        #print("here3")
        tmp_delcov <- data.frame(query=query,subject=subject,gene=gene,cov_fw=cov_fw,cov_bk=cov_bk,delta_cov=del_cov,same_CDS_count=feature_comp, query_org=ref_org, subject_org=org)
        return(tmp_delcov)
      }
    }
  }
}
reverse_BLAST <- function(gr_results){
    #Reversing BLAST direction
    tmp_org_fw <- c()
    tmp_org_fw <- as.data.frame(gr_results)
    tmp_org_fw <- swap_columns(tmp_org_fw,"seqnames", "query_id", swap_names = F)
    tmp_org_fw <- swap_columns(tmp_org_fw,"start", "Hsp_query.from", swap_names = F)
    tmp_org_fw <- swap_columns(tmp_org_fw,"end", "Hsp_query.to", swap_names = F)
    tmp_org_fw <- swap_columns(tmp_org_fw,"query_len", "subject_len", swap_names = F)
    tmp_org_fw <- swap_columns(tmp_org_fw,"query_org", "subject_org", swap_names = F)
    tmp_org_fw <- swap_columns(tmp_org_fw,"query_gene", "subject_gene", swap_names = F)
    tmp_org_fw <- swap_columns(tmp_org_fw,"Hsp_query.frame", "Hsp_hit.frame", swap_names = F)
    return(GRanges(tmp_org_fw))
}

###ENTRYPOINT
set.seed(123)
#source("/run/user/1000/gvfs/sftp:host=max-login.mdc-berlin.net/data/meyer/rob/wisard/wisard/R/WIS_functions.R")
#source("/run/user/1000/gvfs/sftp:host=max-login.mdc-berlin.net/data/meyer/rob/wisard/wisard/R/greedy.R")
#source("/run/user/1000/gvfs/sftp:host=max-login.mdc-berlin.net/data/meyer/rob/wisard/wisard/R/read_blast.R")

#source("/run/user/25654/gvfs/sftp:host=max-login/data/meyer/rob/wisard/wisard/R/WIS_functions.R")
#source("/run/user/25654/gvfs/sftp:host=max-login/data/meyer/rob/wisard/wisard/R/greedy.R")
#source("/run/user/25654/gvfs/sftp:host=max-login/data/meyer/rob/wisard/wisard/R/read_blast.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) Org list path (2) Ref org list path (3) Path to *.orths files and (3) (Combined) gtf_stats.csv", call.=FALSE)
}

#orgs_list_path <- "files/oneway/set.tmp"
#orgs_ref_path <- "files/reference_ORGS.txt"
#orths_path <- "files/all2all_final"
#e_cutoff <- as.numeric("1e-05")
#gtf_path <- "files/annos"
param_file <- "parameters.txt"
param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^") #Convert multibyte seperator to one byte sep #read.table(param_file,sep="==")
delimiter <- as.character(param_table[which(param_table=="transcript_delimiter"),c(2)])
seqID_delimiter <- as.character(param_table[which(param_table=="seqID_delimiter"),c(2)])
orgs_ref_path <- as.character(param_table[which(param_table=="ref_orgs"),c(2)])
e_cutoff <- as.numeric(param_table[which(param_table=="e_value"),c(2)])
clean_extract <- as.logical(param_table[which(param_table=="clean_extract"),c(2)])
#gtf_path <- as.character(param_table[which(param_table=="annos_path"),c(2)])
wisard_path <- as.character(param_table[which(param_table=="wisard_path"),c(2)])
plot_out_path <- as.character(param_table[which(param_table=="plot_path"),c(2)])
mincov_threshold <- as.numeric(param_table[which(param_table=="mincov_threshold"),c(2)])

tryCatch(require(wisard), error = function(){
  source(paste(wisard_path,"/WIS_functions.R",sep = ""))
  source(paste(wisard_path,"/greedy.R",sep = ""))
  source(paste(wisard_path,"/read_blast.R",sep = ""))
  
})

#delimiter<- "_"
#setwd("/run/user/1000/gvfs/sftp:host=max-login.mdc-berlin.net/data/meyer/viz/mrna_loc")

orgs_list_path <- args[1]
#orgs_ref_path <- args[2]
orths_path <- args[2]
gtf_stats_file <- args[3]
#e_cutoff <- as.numeric(args[4])
#gtf_path <- args[5]
#param_file <- args[6]

dir.create(plot_out_path)
orgs.list <- factor(scan(orgs_list_path, character()))
orgs.ref <-  factor(scan(orgs_ref_path, character()))

all_fw_results <- c()
all_bk_results <- c()

if (clean_extract==T) {
  try(file.remove("files/all_bk_results.RData"))
  try(file.remove("files/all_fw_results.RData"))
}

if(file.exists("files/all_bk_results.RData") & file.exists("files/all_fw_results.RData")  & file.info("files/all_bk_results.RData")$size > 0 & file.info("files/all_fw_results.RData")$size > 0){
  try(load("files/all_bk_results.RData"))
  try(load("files/all_fw_results.RData"))
  #print(str(all_bk_results))
  #print(str(all_fw_results))
}
if(is.null(all_fw_results) || is.null(all_bk_results)){
  for (ref in orgs.ref) {
    for (org in orgs.list) {
      if(ref != org){
        file0 <- paste(orths_path,"/",ref,"-",org,".orths", sep="")
        file1 <- paste(orths_path,"/",org,"-",ref,".orths", sep="")
        if(file.exists(file0) && file.exists(file1) && file.info(file0)$size > 0 && file.info(file1)$size > 0 ){
          orths0 <- read.table(paste(orths_path,"/",ref,"-",org,".orths", sep="")) #"danio_rerio-xenopus_laevis.orths")
          orths1 <- read.table(paste(orths_path,"/",org,"-",ref,".orths", sep="")) #"xenopus_laevis-danio_rerio.orths")
          #orths0 <- read.table("danio_rerio-xenopus_tropicalis.orths")
          #orths1 <- read.table("xenopus_tropicalis-danio_rerio.orths")
          
          ###OLD_CODE
          ###OMITTING THIS BECAUSE WE ARE PRETTY CONFIDENT ABOUT our orthologs    
          ###merge the two orths files
          ###this ensures bi-directionality and high confidence orthology
          ##high_q_orthologs <- c(intersect(unique(orths0[,c("sseqid")]), unique(orths1[,c("qseqid")])),intersect(unique(orths0[,c("qseqid")]), unique(orths1[,c("sseqid")])))
          
          ##orths0 <- orths0[duplicated(na.omit(match(orths0$qseqid,high_q_orthologs)), na.omit(match(orths0$sseqid,high_q_orthologs))),]
          ##orths1 <- orths1[duplicated(na.omit(match(orths1$qseqid,high_q_orthologs)), na.omit(match(orths1$sseqid,high_q_orthologs))),]
          
          ##swap qseqid and sseqid columns in either of the dataframes before merging them
          #orths1<- swap_columns(orths1, "qseqid", "sseqid", swap_names = F)
          #orths <- full_join(orths0,orths1)
          
          #orths <- orths[which(orths$evalue < e_cutoff),] 
          
          gr_fw <- create_GRObject(orths0)
          gr_bk <- create_GRObject(orths1)
          
          all_fw_results[[ref]][[org]] <- exec_WISARD(gr_fw,"Hsp_score")
          all_bk_results[[ref]][[org]] <- exec_WISARD(gr_bk,"Hsp_score")
          
          #for (sseq in unique(gr$subject_gene)){   
          #    result <- get_wis(gr[which(sseq == gr$subject_gene),],max_score = "Hsp_score", overlap = 0) 
          #    if(result$max_score > 0){
          #    #print(result)
          #    g_name <- unique(result$alignments$subject_gene)
          #    all_fw_results[[g_name]] <- merge_wis(all_fw_results[[g_name]], result)
          #    }
          #  }
          
        }else{
          print(paste("WARNING:",file0," OR ", file1," doesn't exist"))
        }
      }
    }
  }
  
  save(all_fw_results,file="files/all_fw_results.RData")
  save(all_bk_results,file="files/all_bk_results.RData")
}
#map(all_fw_results, as.data.frame)

##transcript_grah code.txt goes here###

#gtf_dir <- dir(as.character(gtf_path),pattern = "*.gtf.gz", full.names = T)
if(file.exists(gtf_stats_file) && file.info(gtf_stats_file)$size > 0){
  full_gtf_stats <- read.csv(gtf_stats_file, header=F)
  names(full_gtf_stats) <- c("org","search_group","gene_name","gene_id","rna_id","total_exon_len","total_cds_len","five_len","three_len","exon_count","cds_count")
}else{
  stop("(Combined) GTF stats file (gtf_stats.csv) doesn't exist!")
}
gene_list <- unique(union(unique(unlist(map(all_fw_results,function(x){map(x,names)}))),unique(unlist(map(all_bk_results,function(x){map(x,names)})))))
print(gene_list)
#delta_cov_fw <- c()

full_orgs <- union(orgs.list,orgs.ref)
n_cores <- detectCores(all.tests = TRUE, logical = TRUE)
if(is.na(n_cores)){
  n_cores=2
}

HSP <- c()

for(ref_org in orgs.ref){
  ##Read ref-org GTF with "CDS" features
#  ref_gtf_name <- gtf_dir[grep(pattern = paste("\\b",ref_org,"\\b",sep=""),x = gtf_dir, ignore.case = T)]
#  if(length(ref_gtf_name)>1){ref_gtf_name <- ref_gtf_name[length(ref_gtf_name)]}
#  if(!identical(ref_gtf_name,character(0)) && file.exists(ref_gtf_name)){
#    ref_gtf_con <- gzfile(ref_gtf_name,"r")
#    ref_gtf <- read.gtf(file = ref_gtf_con,attr = c("intact"))#,features = c("CDS"))
  # ref_gtf_file <- paste("files/genes/",ref_org,"/gtf_stats.csv")  
  
    HSP_fw <- c()
    HSP_bk <- c()
    tmp_cov <- c()
    for(org in full_orgs){
      if(org!=ref_org){
        ##Read org GTF with "CDS" features
#        org_gtf_name <- gtf_dir[grep(pattern = paste("\\b",org,"\\b",sep=""),x = gtf_dir, ignore.case = T)]
#        if(length(org_gtf_name)>1){org_gtf_name <- org_gtf_name[length(org_gtf_name)]}
#        print(paste(ref_gtf_name,org_gtf_name,sep="-"))
#        if(!identical(org_gtf_name,character(0)) && file.exists(org_gtf_name)){
#          org_gtf_con <- gzfile(org_gtf_name,"r")
#          org_gtf <- read.gtf(file = org_gtf_name ,attr = c("intact"))#,features = c("CDS"))
          
          #orths0 <- read.table("danio_rerio-xenopus_tropicalis.orths")
          #orths1 <- read.table("xenopus_tropicalis-danio_rerio.orths")
          for(gene in gene_list){
            print(gene)
            
            #fw_org_fw_indices <- which(names(all_fw_results[[ref_org]][[org]]) == gene) # org(A)::transcript(i) -> org(B)::transcript(j)
            #bk_org_fw_indices <- which(names(all_fw_results[[org]][[ref_org]]) == gene) # org(A)::transcript(i) -> org(B)::transcript(j)
            #fw_org_bk_indices <- which(names(all_bk_results[[ref_org]][[org]]) == gene) # org(A)::transcript(i) <- org(B)::transcript(j)
            #bk_org_bk_indices <- which(names(all_bk_results[[org]][[ref_org]]) == gene) # org(A)::transcript(i) <- org(B)::transcript(j)
            
            #if(length(fw_org_fw_indices) > 0 && length(fw_org_bk_indices) > 0){ # && length(bk_org_fw_indices) > 0 && length(bk_org_bk_indices) > 0){
              #gene_fw_result <- all_fw_results[[ref_org]][[org]][[gene]]$alignments
              #s_gene_fw <- all_fw_results[[org]][[ref_org]][[gene]]$alignments
              #gene_bk_result <- all_bk_results[[ref_org]][[org]][[gene]]$alignments
              #s_gene_bk <- all_bk_results[[org]][[ref_org]][[gene]]$alignments
              #query_mRNA <- unique(org_fw_result$query_id)
              
              #technically, q_fw == s_bk, & s_fw == q_bk
              #fw_org_fw_indices point to q_fw
              #fw_org_bk_indices point to q_bk
              #bk_org_fw_indices point to s_fw
              #bk_org_bk_indices point to s_bk
              
              gtf_stats <- full_gtf_stats[ ( full_gtf_stats$org==ref_org | full_gtf_stats$org==org ) & full_gtf_stats$search_group==gene,]
            
              q_fw <- all_fw_results[[ref_org]][[org]][[gene]]$alignments #gene_fw_result[which(gene_fw_result$subject_org == ref_org & gene_fw_result$query_org == org),]
              q_bk <- all_bk_results[[ref_org]][[org]][[gene]]$alignments #gene_bk_result[which(gene_bk_result$query_org == ref_org & gene_bk_result$subject_org == org),]
              s_fw <- all_fw_results[[org]][[ref_org]][[gene]]$alignments #s_gene_fw[which(s_gene_fw$subject_org == org & s_gene_fw$query_org == ref_org),]
              s_bk <- all_bk_results[[org]][[ref_org]][[gene]]$alignments #s_gene_bk[which(s_gene_bk$subject_org == ref_org & s_gene_bk$query_org == org),]

              if(!is.null(q_bk) && is.null(q_fw)){
                q_fw <- reverse_BLAST(q_bk)
              }
              if(is.null(q_bk) && !is.null(q_fw)){
                q_bk <- reverse_BLAST(q_fw)
              }
              if(is.null(s_fw) && !is.null(s_bk)){
                s_fw <- reverse_BLAST(s_bk)
              }
              if(!is.null(s_fw) && is.null(s_bk)){
                s_bk <- reverse_BLAST(s_fw)
              }
              queries <- unique(unique(q_fw$query_id),unique(s_bk$query_id))
              subjects <- unique(unique(s_fw$query_id),unique(q_bk$query_id))
              #subjects <- unique(unique(seqnames(q_fw)@values),unique(seqnames(q_bk)@values))
              
              #print(queries)
              #print(subjects)
              
              ##DELTA_COV code
              # child_results <- mcmapply(FUN=calc_delta_cov, query_mRNA, MoreArgs=list(org_fw_result=q_fw,org_bk_result=s_fw) ,mc.cores = n_cores-1, SIMPLIFY = FALSE)
              # for (child in child_results) {
              #   #print(child)           
              #   if(is.null(delta_cov_fw)){
              #     #|| delta_cov_fw[which(delta_cov_fw$query == query & delta_cov_fw$subject == subject),] == -1
              #     delta_cov_fw <- child
              #   }else{
              #     #delta_cov_fw[which(delta_cov_fw$query == query & delta_cov_fw$subject == subject),c("delta_cov")] <- min(del_cov, delta_cov_fw[which(delta_cov_fw$query == query & delta_cov_fw$subject == subject),c("delta_cov")])
              #     existing_del <- delta_cov_fw[c(delta_cov_fw$query == child$query & delta_cov_fw$subject == child$subject),]
              #     if(nrow(existing_del) == 1){
              #       #if(existing_del$delta_cov != child$delta_cov)
              #       delta_cov_fw[c(delta_cov_fw$query == child$query & delta_cov_fw$subject == child$subject),c("delta_cov")] <- min(existing_del$delta_cov,child$delta_cov)
              #     }else{
              #       delta_cov_fw <- rbind(delta_cov_fw,child)
              #       #print(head(delta_cov_fw))
              #     }
              #   }
              # }
              
              ## fw's go together and bk's go together
              
              HSP_fw_results <- c()
              HSP_bk_results <- c()
              
              #sum of length of HSPs / length of CDS
              if(!is.null(queries)){
                HSP_fw_results <- mcmapply(FUN=HSP_Coverage, queries, MoreArgs=list(org_bk=q_bk,org_fw=q_fw,gtf_stats=gtf_stats) ,mc.cores = n_cores-1, SIMPLIFY = FALSE, USE.NAMES=FALSE)
                for(fw in HSP_fw_results){
                  HSP_fw[[gene]] <- rbind(HSP_fw[[gene]],fw)
                }
              }
              if(!is.null(subjects)){
                HSP_bk_results <- mcmapply(FUN=HSP_Coverage, subjects, MoreArgs=list(org_bk=s_bk,org_fw=s_fw,gtf_stats=gtf_stats) ,mc.cores = n_cores-1, SIMPLIFY = FALSE, USE.NAMES=FALSE)
                for(bk in HSP_bk_results){
                  HSP_bk[[gene]] <- rbind(HSP_bk[[gene]],bk)
                }
                }
              #print(head(HSP_fw))
              #print(str(HSP_fw))
            #}
          }
          
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
              ###PLOT-CODE OVER
              
            }
          }
          HSP_fw <- c()
          HSP_bk <- c()
          
         # close(org_gtf_con)
        }
        
      } 
      
    }
    #close(ref_gtf_con)  
  #}
#}

#save(delta_cov_fw,file="delta_cov_fw.RData")
#save(HSP_fw,file="HSP_fw.RData")
#save(HSP_bk,file="HSP_bk.RData")
save(HSP,file="files/HSP.RData")
HSP <- HSP[HSP$min_cov>=mincov_threshold,]
HSP <- HSP[HSP$same_CDS_count==TRUE,]
HSP <- HSP[order(HSP$min_cov,decreasing = T),]
write.csv(HSP,file="files/HSP.csv",quote = F, sep = "\t", row.names = F)

##PLOT
if(file.exists(paste(plot_out_path,"ORGWISE_PLOTS.pdf",sep="/"))){
  file.remove(paste(plot_out_path,"ORGWISE_PLOTS.pdf",sep="/"))
}
##REF_ORGwise plot

pdf(file =paste(plot_out_path,"ORGWISE_PLOTS.pdf",sep="/"),title = "Organism-wise plots")
print(ggplot(HSP, aes(x=min_cov*100, group=ref_org ,fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion") + facet_wrap_paginate(~same_CDS_count + ref_org))
##ORGwise plot
orgwise_plot <- ggplot(HSP, aes(x=min_cov*100, group=same_CDS_count ,fill=same_CDS_count)) + geom_density(aes(y=..density..),position = 'identity', alpha=0.65) + xlab("Min.Coverage of all transcripts(pairwise)") + ylab("Proportion") 
total_pages <- n_pages(orgwise_plot + facet_wrap_paginate(~org + ref_org,nrow=2,ncol=2))
for(i in 1:total_pages){
  try(print(orgwise_plot + facet_wrap_paginate(~org + ref_org,nrow=2,ncol=2,page=i)))
}
dev.off()



sessionInfo()
#delta_cov %>% ggplot(aes(delta_cov$delta_cov, fill = delta_cov$same_CDS_count)) + theme_bw(base_size = 20) + geom_histogram( aes(y = ..density..), position = 'identity') + scale_fill_manual(values=c("black", "orange"))

