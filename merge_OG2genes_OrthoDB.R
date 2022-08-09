#/run/user/25654/gvfs/sftp:host=max-login.mdc-berlin.net/
###Takes optional input of genelist to reduce ODB filesize
require(dplyr)
require(tidyverse)
require(parallel)

args = commandArgs(trailingOnly=TRUE)

if (length(args)>1) {
  stop("Give either the (1) Genelist (to map only the OG clusters containing the genes) (OR) no arguments at all (to map all OG clusters)", call.=FALSE)
}

###FUNCTIONS
reduce_OG_to_genelist <- function(gene_list, OG_table, n_cores, mode="hard") {
    genes <- gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- genes[grep("gene",tolower(genes),invert = T,fixed=T)]
    ###Accounting for gene names having punctuations or version numbers
    genes <- c(genes, unlist(lapply(strsplit(genes,"_"), function(x){return(x[[1]])})))
    if (tolower(SEARCH_MODE)!="hard") {
      matched <- unique(unlist(mclapply(FUN=function(x){grep(x,OG_table$gene_name,ignore.case = T,value = F)}, genes,mc.cores = n_cores-1, mc.cleanup=T)))
    }else{
      matched <- unique(which(!is.na(match(OG_table$gene_name,genes))))
    }
    #matched <- which(!is.na(match(OG_table$gene_name,genes)))
    ##unmatched <- OG_table[is.na(matched),]
    #unmatched_genes <- genes[is.na(match(genes,OG_table$gene_name))]
    #matched_genes <- genes[!is.na(match(genes,OG_table$gene_name))]
    #print(paste("Unmatched Genes:",length(unmatched_genes)))
    #print(unmatched_genes)
    #print(paste("Matched Genes:",length(matched_genes)))
    #print(matched_genes)
    ###Accounting for gene names having punctuations or version numbers
    #stripped_unmatched <- unlist(lapply(strsplit(unmatched_genes,"_"), function(x){return(x[[1]])})) 
    #if(!is.null(stripped_unmatched) && length(stripped_unmatched) > 0){
      #if (tolower(SEARCH_MODE)!="hard") {
        #matched <- unique(c(matched, unlist(lapply(stripped_unmatched, function(x){grep(x,OG_table$gene_name,ignore.case = T,value = F)}))))
     #   matched <- unique(c(matched, unlist(mclapply(FUN=function(x){grep(x,OG_table$gene_name,ignore.case = T,value = F)}, stripped_unmatched,mc.cores = n_cores-1, mc.cleanup=F))))
      #}else{
      #  matched <- unique(c(matched, which(!is.na(match(OG_table$gene_name,stripped_unmatched)))))
      #}
    #}
    
    if(!is.null(matched) && length(matched) > 0){
      return(OG_table[matched,])  
    }else{
      stop("Error: No genes were found, check gene list...Exiting")
    }
    
    ###
    
    #OG_table <- OG_table[matched,]
    # OG_groups <- OG_table  %>% group_rows()
    # out.OG_table <- bind_rows(mclapply(FUN=function(x){
    #  curr_group <- OG_table[x,]
    #  g_names <- paste(unique(curr_group %>% pull("gene_name")),collapse = ",")
    #  g_ids <- paste(unique(curr_group %>% pull("gene_id")),collapse = ",")
    #  return(data.frame(OG= unique(curr_group %>% pull("OG")),gene_ids=g_ids,gene_names=g_names))
    #}, OG_groups, mc.cores = n_cores-1, mc.cleanup=F)) 
    #write.table(out.OG_table, file =  gzfile(odb.user_output_file, open = "w") ,quote = F,sep = "\t", row.names = F, col.names = F)
}
###

###ENTRYPOINT

param_file <- "parameters.txt"
if(!file.exists(param_file) || file.info(param_file)$size < 0){
  stop("ERROR: parameters.txt is missing and is required")
}
param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^", header = T) #Convert multibyte seperator to one byte sep #read.table(param_file,sep="==")

max_concurrent_jobs <- as.numeric(param_table[which(param_table=="max_concurrent_jobs"),c(2)])
ODB_PREFIX <- param_table[which(param_table=="orthodb_path_prefix"),c(2)]
SEARCH_MODE <- param_table[which(param_table=="gene_search_mode"),c(2)]
CLEAN_EXTRACT <- as.logical(param_table[which(param_table=="clean_extract"),c(2)])
n_cores <- detectCores(all.tests = TRUE, logical = TRUE)
if(is.na(n_cores) || n_cores > max_concurrent_jobs){
  n_cores=max_concurrent_jobs
}

odb.OG2genes_file <- paste(file.path(ODB_PREFIX),"_OG2genes.tab.gz",sep = "")
odb.genes_file <- paste(file.path(ODB_PREFIX),"_genes.tab.gz",sep = "")
odb.output_file <-paste(file.path(ODB_PREFIX),"_OGgenes_fixed.tab.gz",sep = "")
odb.user_output_file <-paste(file.path(ODB_PREFIX),"_OGgenes_fixed_user.tab.gz",sep = "")

if(CLEAN_EXTRACT){
  unlink(c(odb.output_file, odb.user_output_file),force = T,expand = T)
}

if(!file.exists(odb.output_file)  || file.info(odb.output_file)$size==0){
  OG_table <- tibble(read.table(gzfile(odb.OG2genes_file,open = "r"), sep="\t", header = F, fill = T, col.names = c("OG","gene_id"))) #read_table("/fast/AG_Meyer/viz/OrthoDB/odb10v1_OG2genes.tab", col_names = c("OG","gene_id"), col_types= c("c","c"))
  #genes_table <- tibble(read.table(gzfile(odb.genes_file, open = "r"), sep="\t", header = F, fill = T)) #read_table("/fast/AG_Meyer/viz/OrthoDB/odb10v1_genes.tab", col_names = c("gene_id","org_id","V3","gene_name","V5","V6","V7","gene_description"), col_types= c("c","c","c","c","c","c","c"))
  #genes_table <- tibble(read.table(textConnection(gsub("\t", "^", readLines(gzfile(odb.genes_file,open = "r"))),open = "r"),sep="^", header = F, fill = T, col.names =  c("gene_id","org_id","V3","gene_name","V5","V6","V7","gene_description")))
  genes_table_lines <- readLines(gzfile(odb.genes_file,open = "r"))
  #genes_table <- tibble(read.table(textConnection(unlist(mcmapply(function(x){
  #  return(gsub("\t", "^",x))
  #},readLines(gzfile(odb.genes_file,open = "r")),USE.NAMES = F,  mc.cores = n_cores-1))),sep="^", header = F, fill = T))
  genes_table <- mcmapply(function(x){
    return(gsub("\t", "^",x))
  },genes_table_lines,USE.NAMES = F,mc.cores = n_cores-1,mc.cleanup=T)
  #save("genes_table", file="genes_table.RData")
  rm("genes_table_lines")
  genes_table <- tibble(read.table(textConnection(genes_table ,open = "r"),sep="^", header = F, fill = T, col.names =  c("gene_id","org_id","V3","gene_name","V5","V6","V7","gene_description")))
  #save("genes_table", file="genes_table.RData")
  #stop()
  #head(OG_table) 
  #head(genes_table)
  #colnames(OG_table)  <- c("OG","gene_id")
  #colnames(genes_table) <- c("gene_id","org_id","V3","gene_name","V5","V6","V7","gene_description")
  
  genes_table <- genes_table[,c("gene_id","gene_name")]
  save("genes_table", file="genes_table.RData")
  #OG_table_rows <- OG_table %>% group_by(gene_id) %>% group_rows()
  #genes_table_rows <- genes_table %>% group_by(gene_id) %>% group_rows()
  
  OG2genes_names <- full_join(OG_table,genes_table, by= "gene_id")
  #OG2genes_table <- OG2genes_table[!any(is.na(OG2genes_table[,c("OG","gene_id","gene_description")])),]
  #OG2genes_table <- OG2genes_table[!any(is.na(OG2genes_table[,c(1,2,5)])),]
  #OG2genes_table <- unique(OG2genes_table)
  #OG2genes_names <- OG2genes_table[,c("OG","gene_id","gene_name")] #OG2genes_table[,c(1,2,5)]
  rm("OG_table","genes_table")
  gc(full = T)
  #OG2genes_names <- OG2genes_names[5000:10000,]
  OG2genes_names <- OG2genes_names[!is.na(OG2genes_names$gene_name),]
  OG2genes_names <- OG2genes_names %>% group_by(OG)
  save("OG2genes_names", file="OG2genes_names.RData")
  #out.OG_table <- bind_rows(lapply(OG_groups, function(x){
  #  curr_group <- OG2genes_names[x,]
  #  g_names <- paste(unique(curr_group %>% pull("gene_name")),collapse = ",")
  #  g_ids <- paste(unique(curr_group %>% pull("gene_id")),collapse = ",")
  #  return(data.frame(OG= unique(curr_group %>% pull("OG")),gene_ids=g_ids,gene_names=g_names))
  #}))
  OG_groups <- OG2genes_names  %>% group_rows() #OG2genes_names  %>% attr("group")
  out.OG_table <- c()
  out.OG_table <- mclapply(FUN=function(x){
    curr_group <- OG2genes_names[x,]
    g_names <- paste(unique(curr_group %>% pull("gene_name")),collapse = ",")
    g_ids <- paste(unique(curr_group %>% pull("gene_id")),collapse = ",")
    return(data.frame(OG= unique(curr_group %>% pull("OG")),gene_ids=g_ids,gene_names=g_names))
  }, OG_groups, mc.cores = n_cores-1, mc.cleanup=T)
  out.OG_table <- bind_rows(out.OG_table)
  save("out.OG_table", file="out.OG_table.RData")
  write.table(out.OG_table, file =  gzfile(odb.output_file, open = "w") ,quote = F,sep = "\t", row.names = F, col.names = F)
}

if (length(args)==1 && (!file.exists(odb.user_output_file)  || file.info(odb.user_output_file)$size==0)) {
  if (!exists("out.OG_table")) {
    out.OG_table <- tibble(read.table(gzfile(odb.output_file,open = "r"), sep="\t", header = F, fill = T))
    names(out.OG_table) <- c("OG","gene_id","gene_name") 
  }
  print(args)
  gene_list <- args[1]
  print(paste("Gene list is provided : ",gene_list,sep = ""))
  #out.OG_table <- c()
  out.OG_table <- reduce_OG_to_genelist(gene_list = gene_list,OG_table = out.OG_table,n_cores = n_cores, mode=SEARCH_MODE)
    #OG_groups <- OG2genes_names  %>% group_rows() #OG2genes_names  %>% attr("group")
    #out.OG_table <- bind_rows(mclapply(FUN=function(x){
    #  curr_group <- OG2genes_names[x,]
    #  g_names <- paste(unique(curr_group %>% pull("gene_name")),collapse = ",")
    #  g_ids <- paste(unique(curr_group %>% pull("gene_id")),collapse = ",")
    #  return(data.frame(OG= unique(curr_group %>% pull("OG")),gene_ids=g_ids,gene_names=g_names))
    #}, OG_groups, mc.cores = n_cores-1, mc.cleanup=F)) 
    #names(OG_groups) <- OG2genes_names %>% group_keys()
    #save("OG2genes_table", file="OG2genes_table.RData")
    #save("out.OG_table", file="out.OG_table.RData")
  write.table(out.OG_table, file =  gzfile(odb.user_output_file, open = "w") ,quote = F,sep = "\t", row.names = F, col.names = F)
}

#sessionInfo()