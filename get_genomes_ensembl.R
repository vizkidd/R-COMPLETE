require(biomartr)
require(dplyr)
require(biomaRt)
require(parallel)
require(tools) ##FOR removing file extensions (gz)
require(RCurl)
require(curl)
require(purrr)
require(fs)
require(processx)
require(ps)
require(stringi)
require(data.table)
require(stringr)

##Credits to https://nbisweden.github.io/workshop-RNAseq/2011/lab_download.html for GTF download through biomart with biomaRt

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) Gene List", call.=FALSE)
}
set.seed(123)

options(RCurlOptions = list(ssl.verifyhost=0, ssl.verifypeer=0))

#FUNCTIONS
add_to_process <- function(p_cmd,p_args=list(),p_wait){
  process_list[sapply(process_list, function(x){
    return(!x$is_alive())
  })] <<- NULL
  tryCatch(expr = function(){
    process_list[sapply(process_list, is.null)] <<- NULL
  }, error = function(){
    if(is.null(process_list)){
      process_list <<- c()
    }
  })
  print(paste("Adding process to list...(",length(process_list),")",sep=""))
  if(length(process_list)>=numWorkers){
    #for (p_id in seq_along(process_list)) {
    #  if(process_list[[p_id]]$is_alive()){
    #    print(paste("Process Q Full...Waiting for a process to end(",length(process_list),")",sep=""))
    #    save(process_list, files="proces_list.RData")
    #    process_list[[p_id]]$wait(timeout=-1)
    #    process_list[[p_idx]] <<- NULL
    #    break;
    #  }
    print(paste("Process Q Full...Waiting for a process to end(",length(process_list),")",sep=""))
    save(process_list, file="proces_list.RData")
    dead_procs <- c()
    for (x in seq_along(process_list)) {
      if(process_list[[x]]$is_alive()){
        process_list[[x]]$wait(timeout=-1)
        break;
      }else{
        dead_procs <- c(dead_procs,x)
      }
    }
    #lapply(seq_along(process_list), function(x){
    #  if(process_list[[x]]$is_alive()){
    #    process_list[[x]]$wait(timeout=-1)
    #    break;
    #  }else{
    #    dead_procs <- c(dead_procs,x)
    #  }
    #})
    
    if (length(dead_procs)>0) {
      dead_procs <- unique(dead_procs)
      process_list[[dead_procs]] <<- NULL 
    }
  }
  
  process_list <<- append(process_list,process$new(command=p_cmd,args = p_args, supervise = TRUE,stdout = "",stderr = "2>&1")) #stderr = T, stdout =  T
  if(p_wait){
    process_list[[1]]$wait(timeout=-1)
  }
}

check_files <-function(fasta_path,org,genes){
  if(dir.exists(fasta_path) && length(dir(fasta_path)) > 0){
    missing_genes <- c()
    available_genes <- c()
    if(file.exists(paste("files/genes/",org,"/","MISSING_GENES",sep=""))){
      missing_genes <- gsub('[[:punct:] ]+','_', factor(scan(paste("files/genes/",org,"/","MISSING_GENES",sep=""), character())))
    }
    if(file.exists(paste("files/genes/",org,"/","AVAILABLE_GENES",sep=""))){
      available_genes <- gsub('[[:punct:] ]+','_', factor(scan(paste("files/genes/",org,"/","AVAILABLE_GENES",sep=""), character())))
    }
    files_in_dir <- unique(sapply(list.files(fasta_path,no.. = T,recursive = F), FUN=function(x){stri_split_fixed(string = x, pattern = fixed('.'),n=2,simplify = T)[1]})) 
    missing_genes <- missing_genes[is.na(match(missing_genes, available_genes))] 
    print(paste("Genes in Dir:",length(files_in_dir),", Missing:",length(missing_genes),", Available:",length(available_genes),", User Genes:",length(genes)))
    if(!is.na(all(match(files_in_dir,available_genes) && all(match(available_genes,files_in_dir)))) && is.na(all(match(missing_genes, files_in_dir)))){
      print(paste(fasta_path,": Check PASSED!"))
      return(TRUE)
    }
  }
  print(paste(fasta_path,": Check FAILED!"))
  return(FALSE)
}

get_genome_mart <- function(org) {
  mart.dataset <- grep(x = org.meta.list$dataset, pattern=regex(stri_split_fixed(org,"_",2,simplify = T)[2],ignore_case = T),fixed=F, value = T)
  if(length(mart.dataset)==1){
    using.mart.data <- useMart(ENSEMBL_MART, mart.dataset)
    
    #coords <- biomaRt::select(x=using.mart.data, keys=keys(using.mart.data, "chromosome_name"), columns=c("start_position","end_position","chromosome_name"), keytype = "chromosome_name" )
    #seq <- data.frame(sequences=biomaRt::select(x=using.mart.data, keys=keys(using.mart.data, "chromosome_name"), columns=c("cdna","chromosome_name"), keytype = "chromosome_name" ))
    curl_handle <- getCurlMultiHandle()
    dir.create(path = paste(GENOMES_PATH,"/",org,"/",sep=""),showWarnings = F,recursive = T)
    chr_names <- keys(using.mart.data, "chromosome_name")
    mclapply(chr_names, function(x){
      seq <- getBM(mart=using.mart.data, values=x, curl=curl_handle, attributes=c("cdna"), filters = "chromosome_name", uniqueRows=T, useCache=F) 
      fwrite(x = list(seq[,1]),file = paste(GENOMES_PATH,"/",org,"/",x,".fa.gz",sep=""), compress = "gzip", quote = F,nThread=numWorkers, col.names=F)
    }, mc.cleanup = T, mc.cores = numWorkers,mc.silent = T)
    
    
    lapply(chr_names, function(x){
      if (file.exists(paste(GENOMES_PATH,"/",org,"/",x,".fa.gz",sep="")) && file.info(paste(GENOMES_PATH,"/",org,"/",x,".fa.gz",sep=""))$size > 0) {
        fwrite(list(paste(">",x,sep="")),file = paste(GENOMES_PATH,"/",org,".fa.gz",sep=""),append = T,quote = F,compress = "gzip", col.names = F, nThread = numWorkers)
        fwrite(list(gsub(pattern = "\n", replacement = "",fixed=T,x = readLines(con = gzfile(paste(GENOMES_PATH,"/",org,"/",x,".fa.gz",sep=""))))),file = paste(GENOMES_PATH,"/",org,".fa.gz",sep=""),append = T,quote = F,compress = "gzip", col.names = F, nThread = numWorkers) 
      }
    })
    
    unlink(x = paste(GENOMES_PATH,"/",org,"/"),force = T,expand = T)
    
    #martDisconnect(using.mart.data)
    
    return(paste(GENOMES_PATH,"/",org,".fa.gz",sep=""))
  }else{
    print(paste("Error : Genome of",org,"not found/cannot be downloaded. Manually download the files or give the links in user list!"))
    return(FALSE)
  }
}

get_gtf_mart <- function(org){
  ##Credits to https://nbisweden.github.io/workshop-RNAseq/2011/lab_download.html for GTF download through biomart with biomaRt
  gtf_attributes <- c("strand",
                    "ensembl_gene_id",
                    "ensembl_transcript_id",
                    "external_gene_name",
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
  mart.dataset <- grep(x = org.meta.list$dataset, pattern=regex(stri_split_fixed(org,"_",2,simplify = T)[2],ignore_case = T),fixed=F, value = T)
  if(length(mart.dataset)==1){
    using.mart.data <- useMart(ENSEMBL_MART, mart.dataset)
    
    bm_gtf <- getBM(mart=using.mart.data,attributes=gtf_attributes,uniqueRows=T, useCache=F) #bm_gtf <- getBM(mart=using.mart.data,attributes=gtf_attributes,uniqueRows=T, useCache=F, filters = c("external_gene_name"), values = genes)
    bm_gtf <- dplyr::arrange(bm_gtf,chromosome_name,start_position)
    
    bm_gtf$strand[bm_gtf$strand==1] <- "+"
    bm_gtf$strand[bm_gtf$strand== -1] <- "-"
    bm_gtf$source_name <- "Ensembl"
    
    write.table(x=bm_gtf,file = gzfile(paste(ANNOS_PATH,"/",org,".gtf.gz",sep="")),sep="\t",row.names=F,quote=F)
    
    #martDisconnect(using.mart.data)
    return(paste(ANNOS_PATH,"/",org,".gtf.gz",sep=""))
  }else{
    print(paste("Error : GTF of",org,"not found/cannot be downloaded. Manually download the files or give the links in user list!"))
    return(FALSE)
  }
}

try_genome_download <- function(org, release, assembly_type) {
  g_name <- paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")
  
  if(!file.exists(g_name) || !file.info(g_name)$size > 0){
    g_path <- try(getGenome(db = tolower(GENOMES_SOURCE), organism=org, reference = F, gunzip = F, path = GENOMES_PATH, release = as.numeric(release),assembly_type = assembly_type))
  
    if(is.logical(g_path)){
      g_path <- get_genome_mart(org)
    }
  }else{
    if(system2("gzip",args = c("-t", g_name), wait = T,stdout = NULL, stderr = NULL) == 0){
      g_path<-g_name
    }else{
      lapply(list.files(GENOMES_PATH, full.names = T, ignore.case = T, no.. = T, pattern = regex(stri_split_fixed(org,"_",2,simplify = T),ignore_case = T)),function(x){
        if(grepl(x,pattern = "fa|gz")){
          try(file.remove(x, showWarnings=F))
        }
      })
      g_path <- try(getGenome(db = tolower(GENOMES_SOURCE), organism=org, reference = F, gunzip = F, path = GENOMES_PATH, release = as.numeric(release),assembly_type = assembly_type))
      
      if(is.logical(g_path)){
        g_path <- get_genome_mart(org)
      }
    }
  }
  
  if(!is.logical(g_path) && g_path!=g_name){
    try(file_move(g_path, g_name))
    g_path <- g_name
  }
  
  return(g_path)
}

try_gtf_download <- function(org, release, assembly_type) {
    a_name <- paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")
    if(!file.exists(a_name) || !file.info(a_name)$size > 0){   
      a_path <- try(getGTF(db = tolower(GENOMES_SOURCE), organism=org, path = file.path(ANNOS_PATH),release = as.numeric(release),assembly_type = assembly_type))

      if (is.logical(a_path) || !all(grepl("gene_name",readLines(con = gzfile(a_path,open = "r"))))) { 
        a_path <- get_gtf_mart(org)
      }
    }else{
      if(system2("gzip",args = c("-t", a_name), wait = T,stdout = NULL, stderr = NULL) == 0){
        a_path<-a_name
      }else{
        lapply(list.files(ANNOS_PATH, full.names = T, ignore.case = T, no.. = T, pattern = regex(stri_split_fixed(org,"_",2,simplify = T),ignore_case = T)),function(x){
          if(grepl(x,pattern = "gtf|gz")){
            try(file.remove(x, showWarnings=F))
          }
        })
        a_path <- try(getGTF(db = tolower(GENOMES_SOURCE), organism=org, path = file.path(ANNOS_PATH),release = as.numeric(release),assembly_type = assembly_type))
        
        if (is.logical(a_path) || !any(grepl("gene_name",readLines(con = gzfile(a_path,open = "r"))))) { 
          a_path <- get_gtf_mart(org)
        }
      }
    }
    
    if(!is.logical(a_path) && a_path!=a_name){
      try(file_move(a_path, a_name))
      # a_path <- a_name
    }
    
  return(a_path)
} 

fetch_genome_db <- function(org_row, genes) {
  #print(org_row)
  org <- org_row["name"]
  release <- as.numeric(org_row["release"])
  assembly_type = "primary_assembly" #c("primary_assembly","toplevel")
  #print(paste(org,release,assembly_type))

  fasta_path <- file.path(OUT_PATH ,org) 
  #print(fasta_path)
  details <- file.info(list.files(path = fasta_path, full.names = TRUE))
  if(CLEAN_DOWNLOAD==TRUE){
    try(file.remove(paste(GENOMES_PATH, "/",org,".fa.gz",sep = ""),showWarnings=F))
    try(file.remove(paste(ANNOS_PATH, "/",org,".gtf.gz",sep = ""),showWarnings=F))
  }
  
  genome_path <- try_genome_download(org, release, assembly_type)
  
  gtf_path <- try_gtf_download(org, release, assembly_type)
  
  print(paste(genome_path,gtf_path,org))
  
  if(CLEAN_EXTRACT || !check_files(fasta_path,org,genes) && !is.logical(gtf_path) && !is.logical(genome_path)){
    do.call(add_to_process,list(p_cmd = c("./jobhold.sh"), p_args = c(paste("extract",org,sep="_"), "./extract_genomic_regions.sh", genome_path, gtf_path, gene_list, org), p_wait = SUBPROCESS_WAIT))
  }
  Sys.sleep(5) 
}
fetch_genome_user <- function(data, genes){
  genome_path <- c()
  gtf_path <- c()
  #print(data)
  genome <- data["genome"]
  gtf <- data["gtf"]
  org <- data["org"]
  
  if(genome == "-" || gtf == "-"){
    fetch_genome_db(org, genes)
    break;
  }
  genome_path<-paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")
  gtf_path<-paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")  
  fasta_path <- file.path(OUT_PATH ,org)
  
  if(grepl("http",genome)){
    if(!file.exists(genome_path) || !file.info(genome_path)$size > 0){
      curl_fetch_disk(genome, paste(GENOMES_PATH, "/",basename(URLdecode(genome)),sep=""))
      genome_path <- paste(GENOMES_PATH, "/",basename(URLdecode(genome)),sep="")
    }else{
      if(system2("gzip",args = c("-t", genome_path), wait = T,stdout = NULL, stderr = NULL) == 0){
        genome_path<-genome_path 
      }else{
        lapply(list.files(GENOMES_PATH, full.names = T, ignore.case = T, no.. = T, pattern = regex(stri_split_fixed(org,"_",2,simplify = T),ignore_case = T)),function(x){
          if(grepl(x,pattern = "fa|gz")){
            try(file.remove(x, showWarnings=F))
          }
        })
        curl_fetch_disk(genome, basename(URLdecode(genome)))
        genome_path <- paste(GENOMES_PATH, "/",basename(URLdecode(genome)),sep="") #basename(URLdecode(genome))
      }
    }
    
  }else{
    genome_path <- genome
  }
  
  if(grepl("http",gtf)){
    if(!file.exists(gtf_path) || !file.info(gtf_path)$size > 0){
      curl_fetch_disk(gtf, paste(ANNOS_PATH, "/",basename(URLdecode(gtf)),sep=""))
      gtf_path <- paste(ANNOS_PATH, "/",basename(URLdecode(gtf)),sep="")
    }else{
      if(system2("gzip",args = c("-t", gtf_path), wait = T,stdout = NULL, stderr = NULL) == 0){
        gtf_path<-gtf_path        
      }else{
        lapply(list.files(ANNOS_PATH, full.names = T, ignore.case = T, no.. = T, pattern = regex(stri_split_fixed(org,"_",2,simplify = T),ignore_case = T)),function(x){
          if(grepl(x,pattern = "gtf|gz")){
            try(file.remove(x, showWarnings=F))
          }
        })
        curl_fetch_disk(gtf, paste(ANNOS_PATH, "/",basename(URLdecode(gtf)),sep=""))
        gtf_path <- paste(ANNOS_PATH, "/",basename(URLdecode(gtf)),sep="")
      }
    }
    
  }else{
    gtf_path <- gtf
  }
  
  if(!is.logical(genome_path) && file_ext(genome_path) == "gz"){
    try(file_move(genome_path, paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")))
    genome_path <- paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")
  }else if(!is.logical(genome_path)){
    g_ext <- file_ext(genome_path)
    try(file_move(genome_path, paste(GENOMES_PATH, "/",org,".",g_ext,sep = "")))
    genome_path <- paste(GENOMES_PATH, "/",org,".",g_ext,sep = "")
  }
  if(!is.logical(gtf_path) && file_ext(gtf_path) == "gz"){
    try(file_move(gtf_path, paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")))
    gtf_path <- paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")
  }else if(!is.logical(gtf_path)){
    a_ext <- file_ext(gtf_path)
    try(file_move(gtf_path, paste(ANNOS_PATH, "/",org,".",a_ext,sep = "")))
    gtf_path <- paste(ANNOS_PATH, "/",org,".",a_ext,sep = "")
  }
  
  print(paste(genome_path,gtf_path,org))
  
  if(CLEAN_EXTRACT || !check_files(fasta_path,org,genes) && !is.logical(gtf_path) && !is.logical(genome_path)){
    do.call(add_to_process,list(p_cmd = c("./jobhold.sh"), p_args = c(paste("extract",org,sep="_"), "./extract_genomic_regions.sh",genome_path, gtf_path, gene_list, org), p_wait = SUBPROCESS_WAIT))
  }
}

####ENTRY POINT
process_list <<- c()

param_file <- "parameters.txt"

if(!file.exists(param_file) || file.info(param_file)$size < 0){
  stop("ERROR: parameters.txt is missing and is required")
}

param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^", header = T) #Convert multibyte seperator to one byte sep

GENOMES_PATH <<- param_table[which(param_table=="genomes_path"),c(2)]
ANNOS_PATH <<- param_table[which(param_table=="annos_path"),c(2)]
TEMP_PATH <<- param_table[which(param_table=="temp_path"),c(2)]
GROUPS_PATH <<- param_table[which(param_table=="groups_path"),c(2)]
OUT_PATH <<- param_table[which(param_table=="fasta_path"),c(2)]
BED_PATH <<- as.character(param_table[which(param_table=="bed_path"),c(2)])
CLEAN_DOWNLOAD <<- as.logical(param_table[which(param_table=="clean_download"),c(2)])
CLEAN_EXTRACT <<- as.logical(param_table[which(param_table=="clean_extract"),c(2)])
SUBPROCESS_WAIT <<- as.logical(param_table[which(param_table=="subprocess_wait"),c(2)])
BLAST_REGION <<- tolower(as.character(param_table[which(param_table=="blast_region"),c(2)]))
GENOMES_SOURCE <<- as.character(param_table[which(param_table=="genomes_source"),c(2)])
REMOVE_DOWNLOADS <<- as.character(param_table[which(param_table=="remove_downloads"),c(2)])
USER_GENOMES <<- as.character(param_table[which(param_table=="user_genomes"),c(2)])
max_concurrent_jobs <<- as.numeric(param_table[which(param_table=="max_concurrent_jobs"),c(2)])

gene_list <- args[1]

print(paste("CLEAN_DOWNLOAD:",CLEAN_DOWNLOAD))
print(paste("CLEAN_EXTRACT:",CLEAN_EXTRACT))
print(paste("REMOVE DOWNLOADS:",REMOVE_DOWNLOADS))
print(paste("GENOMES_SOURCE:",GENOMES_SOURCE))
print(paste("MAX PROCESSES:",max_concurrent_jobs))
print(paste("SUBPROCESS_WAIT:",SUBPROCESS_WAIT))

tryCatch(numWorkers <<- detectCores(all.tests = T, logical = T), error=function(){numWorkers <<- 2})
if(is.na(numWorkers) || numWorkers > max_concurrent_jobs){
  numWorkers <<- max_concurrent_jobs
}

dir.create("files",showWarnings = F, recursive = T)
dir.create(BED_PATH,showWarnings = F, recursive = T)
unlink(TEMP_PATH, recursive = T,force = T,expand = T)
dir.create(TEMP_PATH,showWarnings = F, recursive = T)
unlink(GROUPS_PATH, recursive = T,force = T,expand = T)
dir.create(GROUPS_PATH,showWarnings = F, recursive = T)
dir.create(OUT_PATH,showWarnings = F, recursive = T)

genes <- gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
genes <- genes[grep("gene",tolower(genes), invert = T, fixed = T)]

ENSEMBL_MART <<- "ENSEMBL_MART_ENSEMBL"
using.mart <- useMart(ENSEMBL_MART) #For biomaRt
org.meta.list <<- listDatasets(mart=using.mart) #For biomaRt
org.meta <- listGenomes(db = tolower(GENOMES_SOURCE), type = "all", details = T) #For biomartr
#for(org in getENSEMBLGENOMESInfo()$name){
#  org.meta <- rbind(org.meta, is.genome.available(db = "tolower(GENOMES_SOURCE)", org, details = T))
#}

check_list <- list()

#Read user data
#if(tolower(GENOMES_SOURCE)=="both" || tolower(GENOMES_SOURCE)=="user"){
if(!stri_isempty(USER_GENOMES) && !is.null(USER_GENOMES)){
  user_data <- read.csv(USER_GENOMES,header = F)
  names(user_data) <- c("org","genome","gtf")
}

print("Checking and downloading organisms...")

#if(tolower(GENOMES_SOURCE)=="both" || tolower(GENOMES_SOURCE)=="user"){
if(nrow(user_data)!=0 && !is.null(user_data)){
  #for (row in 1:nrow(user_data)) {
  #  #print(user_data[row,])
  #  fetch_genome_user(user_data[row,], genes)
  #}
  apply(user_data, MARGIN = 1, function(x){
    fetch_genome_user(x, genes)
  })
}
#}

#if(tolower(GENOMES_SOURCE)=="both" || tolower(GENOMES_SOURCE)!="user"){
#mclapply(org.meta$name, fetch_genome_db, mc.cores = 1)
apply(org.meta, MARGIN = 1, FUN = function(x){
  fetch_genome_db(x, genes)
})
#}

print(process_list)

p_pid <- c()
for(proc in process_list){
  if(proc$is_alive()){
    print(proc$get_status())
    if(ps_is_running(proc$as_ps_handle())){
      ## if SIGCHLD is overwritten the process is lost, so we try to get pid and check if the process is running
      proc$wait(timeout = -1) 
    }else{
      print(proc$print())
      print(proc$get_status())
      proc$interrupt()
    }
  }
  #print(proc$get_exit_status())
}


write.table(x = org.meta,file = "files/org_meta.txt", quote = F,sep=",", row.names = F)
try(write.table(x = user_data,file = "files/org_meta.txt", quote = F,sep=",", row.names = F, append = T))

sessionInfo()