require(biomartr)
require(parallel)
require(tools) ##FOR removing file extensions (gz)
require(RCurl)
require(curl)
require(purrr)
require(fs)
require(processx)
require(ps)
require(stringr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) folder path for storing genomes (2) folder path for storing annotations (3) the filename with gene symbols(without header) (4) FASTA OUTPUT PATH(absolute path)", call.=FALSE)
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
    save(process_list, files="proces_list.RData")
    dead_procs <- c()
    lapply(seq_along(process_list), function(x){
      if(process_list[[x]]$is_alive()){
        process_list[[x]]$wait(timeout=-1)
      }
      dead_procs <- c(dead_procs,x)
    })
    dead_procs <- unique(dead_procs)
    process_list[[dead_procs]] <<- NULL
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
    files_in_dir <- unique(sapply(list.files(fasta_path,no.. = T,recursive = F), FUN=function(x){str_split_fixed(string = x, pattern = fixed('.'),n=2)[1]})) 
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
fetch_genome_ensembl <- function(org, genes) {
  
  if(tolower(GENOMES_SOURCE)=="both"){
    if(!is.na(match(org,user_data$org))){
      return()
    }
  }
  fasta_path <- file.path(OUT_PATH ,org) 
  #print(fasta_path)
  details <- file.info(list.files(path = fasta_path, full.names = TRUE))
  if(CLEAN_DOWNLOAD==TRUE){
    try(file.remove(paste(GENOMES_PATH, "/",org,".fa.gz",sep = ""),showWarnings=F))
    try(file.remove(paste(ANNOS_PATH, "/",org,".gtf.gz",sep = ""),showWarnings=F))
  }
  g_name <- c()
  a_name <- c()
  g_name <- paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")
  a_name <- paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")
  if(!file.exists(g_name) || !file.info(g_name)$size > 0){
    genome_path <-  try(getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH)) 
    if(system2("gzip",args = c("-t", genome_path), wait = T,stdout = NULL, stderr = NULL)!=0){ #IF file didnt download properly
      lapply(list.files(GENOMES_PATH, full.names = T, ignore.case = T, no.. = T, pattern = str_split_fixed(org,"_",2)[1]),function(x){
        if(grepl(x,pattern = "fa|gz")){
          try(file.remove(x, showWarnings=F))
        }
      })
      genome_path <-  try(getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH))
    }
  }else{
    if(system2("gzip",args = c("-t", g_name), wait = T,stdout = NULL, stderr = NULL) == 0){
      genome_path<-g_name
    }else{
      lapply(list.files(GENOMES_PATH, full.names = T, ignore.case = T, no.. = T, pattern = str_split_fixed(org,"_",2)[1]),function(x){
        if(grepl(x,pattern = "fa|gz")){
          try(file.remove(x, showWarnings=F))
        }
      })
      genome_path <-  try(getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH)) 
    }
  }
  
  if(!file.exists(a_name) || !file.info(a_name)$size > 0){
    gtf_path <- try(getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH))) 
    if(system2("gzip",args = c("-t", gtf_path), wait = T,stdout = NULL, stderr = NULL)!=0){ #IF file didnt download properly
      lapply(list.files(ANNOS_PATH, full.names = T, ignore.case = T, no.. = T, pattern = str_split_fixed(org,"_",2)[1]),function(x){
        if(grepl(x,pattern = "gtf|gz")){
          try(file.remove(x, showWarnings=F))
        }
      })
      gtf_path <- try(getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH)))
      }
  }else{
    if(system2("gzip",args = c("-t", a_name), wait = T,stdout = NULL, stderr = NULL) == 0){
      gtf_path<-a_name
    }else{
      lapply(list.files(ANNOS_PATH, full.names = T, ignore.case = T, no.. = T, pattern = str_split_fixed(org,"_",2)[1]),function(x){
        if(grepl(x,pattern = "gtf|gz")){
          try(file.remove(x, showWarnings=F))
        }
      })
      gtf_path <- try(getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH))) 
    }
    
  }
  if(CLEAN_EXTRACT || !check_files(fasta_path,org,genes)){
    do.call(add_to_process,list(p_cmd = c("./jobhold.sh"), p_args = c(paste("extract",org,sep="_"), "./extract_genomic_regions.sh", genome_path, gtf_path, gene_list, org), p_wait = SUBPROCESS_WAIT))
  }else{
    if(file.exists(genome_path) && file_ext(genome_path) == "gz"){
      try(file_move(genome_path, paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")))
    }
    if(file.exists(gtf_path) && file_ext(gtf_path) == "gz"){
      try(file_move(gtf_path, paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")))
    }}
  Sys.sleep(5) 
}
fetch_genome_user <- function(data, genes){
  genome_path <- c()
  gtf_path <- c()
  print(data)
  genome <- data["genome"]
  gtf <- data["gtf"]
  org <- data["org"]
  
  if(genome == "-" || gtf == "-"){
    fetch_genome_ensembl(org, genes)
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
        lapply(list.files(GENOMES_PATH, full.names = T, ignore.case = T, no.. = T, pattern = str_split_fixed(org,"_",2)[1]),function(x){
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
        lapply(list.files(ANNOS_PATH, full.names = T, ignore.case = T, no.. = T, pattern = str_split_fixed(org,"_",2)[1]),function(x){
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
  
  print(genome_path)
  print(gtf_path)
  
  if(CLEAN_EXTRACT || !check_files(fasta_path,org,genes)){
    do.call(add_to_process,list(p_cmd = c("./jobhold.sh"), p_args = c(paste("extract",org,sep="_"), "./extract_genomic_regions.sh",genome_path, gtf_path, gene_list, org), p_wait = SUBPROCESS_WAIT))
  }else{
    if(file.exists(genome_path) && file_ext(genome_path) == "gz"){
      try(file_move(genome_path, paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")))
    }
    if(file.exists(gtf_path) && file_ext(gtf_path) == "gz"){
      try(file_move(gtf_path, paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")))
    }}
}

####ENTRY POINT
process_list <<- c()

param_file <- "parameters.txt"

if(!file.exists(param_file) || file.info(param_file)$size < 0){
  stop("ERROR: parameters.txt is missing and is required")
}

param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^", header = T) #Convert multibyte seperator to one byte sep

GENOMES_PATH <- param_table[which(param_table=="genomes_path"),c(2)]
ANNOS_PATH <- param_table[which(param_table=="annos_path"),c(2)]
TEMP_PATH <- param_table[which(param_table=="temp_path"),c(2)]
GROUPS_PATH <- param_table[which(param_table=="groups_path"),c(2)]
OUT_PATH <- param_table[which(param_table=="fasta_path"),c(2)]
BED_PATH <- as.character(param_table[which(param_table=="bed_path"),c(2)])
CLEAN_DOWNLOAD <- as.logical(param_table[which(param_table=="clean_download"),c(2)])
CLEAN_EXTRACT <- as.logical(param_table[which(param_table=="clean_extract"),c(2)])
SUBPROCESS_WAIT <- as.logical(param_table[which(param_table=="subprocess_wait"),c(2)])
BLAST_REGION <- tolower(as.character(param_table[which(param_table=="blast_region"),c(2)]))
GENOMES_SOURCE <- as.character(param_table[which(param_table=="genomes_source"),c(2)])
REMOVE_DOWNLOADS <- as.character(param_table[which(param_table=="remove_downloads"),c(2)])
USER_GENOMES <- as.character(param_table[which(param_table=="user_genomes"),c(2)])
max_concurrent_jobs <- as.numeric(param_table[which(param_table=="max_concurrent_jobs"),c(2)])

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

org.meta <- c()
for(org in getENSEMBLGENOMESInfo()$name){
  org.meta <- rbind(org.meta, is.genome.available(db = "ensembl", org, details = T))
}

check_list <- list()

#Read user data
if(tolower(GENOMES_SOURCE)=="both" || tolower(GENOMES_SOURCE)=="user"){
  user_data <- read.csv(USER_GENOMES,header = F)
  names(user_data) <- c("org","genome","gtf")
}

print("Checking and downloading organisms...")

if(tolower(GENOMES_SOURCE)=="both" || tolower(GENOMES_SOURCE)=="user"){
  if(nrow(user_data)!=0){
    #for (row in 1:nrow(user_data)) {
    #  #print(user_data[row,])
    #  fetch_genome_user(user_data[row,], genes)
    #}
    apply(user_data, MARGIN = 1, function(x){
      fetch_genome_user(x, genes)
    })
  }
}

if(tolower(GENOMES_SOURCE)=="both" || tolower(GENOMES_SOURCE)=="ensembl"){
  #mclapply(org.meta$name, fetch_genome_ensembl, mc.cores = 1)
  lapply(org.meta$name, FUN = function(x){
    fetch_genome_ensembl(x, genes)
  })
}

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