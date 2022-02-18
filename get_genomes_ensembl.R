require(biomartr)
require(parallel)
require(tools) ##FOR removing file extensions (gz)
require(RCurl)
require(curl)
require(purrr)
require(fs)
require(processx)
require(stringr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  #stop("Give the (1) folder path for storing genomes (2) folder path for storing annotations (3) the filename with gene symbols(without header)", call.=FALSE)
  stop("Give the (1) folder path for storing genomes (2) folder path for storing annotations (3) the filename with gene symbols(without header) (4) FASTA OUTPUT PATH(absolute path)", call.=FALSE)
}
set.seed(123)

options(RCurlOptions = list(ssl.verifyhost=0, ssl.verifypeer=0)) #qsub -l data -V ./run_pipeline_ensembl.sh files/trimlist.txt 

#DECLARATIONS
#GENOMES_PATH <- "files/genomes"
#ANNOS_PATH <- "files/annos"
#gene_list <- "genelist.txt"

#FUNCTIONS
add_to_process <- function(p_cmd,p_args=list(),p_wait){
  #system2("./extract_genomic_regions.sh",args = c(genome_path, gtf_path, file.path("files","bed",paste(org,"_genes",sep = "")), gene_list, org), stderr = T, stdout =  T,wait = SUBPROCESS_WAIT)
  #print(p_cmd)
  #print(p_args)
  #print(p_wait)
  print(paste("Adding process to list...(",length(process_list),")",sep=""))
  #process_list <<- append(process_list,list(tryCatch(system2(p_cmd,args = p_args, wait = p_wait, stderr = F, stdout = F),finally = function(){process_count <<- process_count+1}))) #stderr = T, stdout =  T
  #process_list <<- append(process_list,tryCatch(system2(p_cmd,args = p_args, wait = p_wait),finally = function(){process_count <<- process_count+1})) #stderr = T, stdout =  T
  if(length(process_list)==numWorkers){
    p_idx <- 0
    for (p_id in process_list) {
      p_idx <- p_idx + 1
      if(p_id$is_alive()){
        print(paste("Process Q Full...Waiting for a process to end(",length(process_list),")",sep=""))
        p_id$wait(timeout=-1)
        break;
      }
    }
    if(p_idx > 0){
      process_list[[p_idx]] <- NULL
    }
  }
  
  process_list <<- append(process_list,process$new(command=p_cmd,args = p_args, supervise = TRUE,stdout = "",stderr = "2>&1",post_process = function(){process_count <<- process_count+1})) #stderr = T, stdout =  T
  if(p_wait){
    process_list[[length(process_list)]]$wait(timeout=-1)
  }
}

check_files <-function(fasta_path){
  if(dir.exists(fasta_path)){
    missing_genes <- c()
    available_genes <- c()
    if(file.exists(paste(fasta_path,"/","MISSING_GENES",sep=""))){
      missing_genes <- gsub('[[:punct:] ]+','_', factor(scan(paste(fasta_path,"/","MISSING_GENES",sep=""), character())))
    }
    if(file.exists(paste(fasta_path,"/","AVAILABLE_GENES",sep=""))){
      available_genes <- gsub('[[:punct:] ]+','_', factor(scan(paste(fasta_path,"/","AVAILABLE_GENES",sep=""), character())))
    }
    files_in_dir <- list.files(fasta_path,pattern=paste(".",BLAST_REGION,sep = ""))
    #print(files_in_dir)
    #sum(unlist(map2(genes, files_in_dir,function(x,y){
    #print(paste(y,x))
    #grepl(y,pattern = x,ignore.case = T)
    #}))) == length(genes)
    
    if(length(missing_genes)==length(genes) || length(available_genes)==length(files_in_dir) || length(missing_genes)==0){
      print(paste(fasta_path,": Check PASSED!"))
      return(TRUE)
      #}
    }
  }
  print(paste(fasta_path,": Check FAILED!"))
  return(FALSE)
}
fetch_genome_ensembl <- function(org) {
  
  if(tolower(GENOMES_SOURCE)=="both"){
    if(!is.na(match(org,user_data$org))){
      return()
    }
  }
  fasta_path <- file.path(OUT_PATH ,org) ## file_path_as_absolute(file.path(OUT_PATH ,org))
  print(fasta_path)
  details <- file.info(list.files(path = fasta_path, full.names = TRUE))
  #print(details)
  #print(rownames(details[details$size == 0, ]))
  if(CLEAN_DOWNLOAD==TRUE){
    try(file.remove(paste(GENOMES_PATH, "/",org,".fa.gz",sep = ""),showWarnings=F))
    try(file.remove(paste(ANNOS_PATH, "/",org,".gtf.gz",sep = ""),showWarnings=F))
    #unlink(paste(OUT_PATH,"/",org,sep=""), recursive=TRUE)
  }
  g_name <- c()
  a_name <- c()
  g_name <- paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")
  a_name <- paste(ANNOS_PATH, "/",org,".gtf.gz",sep = "")
  #if(!dir.exists(paths = fasta_path) || length(dir(path = fasta_path)) <= 1 || length(rownames(details[details$size == 0, ])) != 0  )
  #{
  if(!file.exists(g_name) || !file.info(g_name)$size > 0){
    genome_path <-  try(getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH)) # tryCatch(getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH),error = function(e){print("Sleeping...");Sys.sleep(360);},finally=getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH))
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
      genome_path <-  try(getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH)) # tryCatch(getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH),error = function(e){print("Sleeping...");Sys.sleep(360);},finally=getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH))
    }
  }
  
  if(!file.exists(a_name) || !file.info(a_name)$size > 0){
    gtf_path <- try(getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH))) # tryCatch(getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH)),error = function(e){print("Sleeping...");Sys.sleep(360)},finally=getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH)))
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
      gtf_path <- try(getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH))) # tryCatch(getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH)),error = function(e){print("Sleeping...");Sys.sleep(360)},finally=getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH)))
    }
    
  }
  if(CLEAN_EXTRACT || !check_files(fasta_path)){
    ##g_unzip.str <- paste("gunzip","-d", genome_path, sep = " ")
    #system2(command = "gunzip",args = c("-d","-f", genome_path),wait=T)
    ##unzip.str <- paste("gunzip","-d", gtf_path, sep = " ")
    #system2(command = "gunzip",args = c("-d","-f", gtf_path),wait=T)
    #new_gtf_path <- paste(ANNOS_PATH,paste(org,".gtf",sep=""), sep="/")
    #system2(command = "mv",args = c(tools::file_path_sans_ext(gtf_path),new_gtf_path),wait=T)
    ##faidx.str <- paste("samtools","faidx", genome_path, sep = " ")
    #system2("samtools", args = c("faidx", tools::file_path_sans_ext(genome_path)),wait=T)
    
    #script.str <- paste("./extract_genomic_regions.sh", genome_path, gtf_path, file.path("files","bed",paste(org,"_genes",sep = "")), gene_list, org, sep=" ")
    #system2("./extract_genomic_regions.sh",args = c(tools::file_path_sans_ext(genome_path), new_gtf_path, file.path("files","bed",paste(org,"_genes",sep = "")), gene_list, org), wait = SUBPROCESS_WAIT)
    #process_list <- append(process_list,
    #system2("./extract_genomic_regions.sh",args = c(genome_path, gtf_path, file.path("files","bed",paste(org,"_genes",sep = "")), gene_list, org), stderr = T, stdout =  T,wait = SUBPROCESS_WAIT)
    do.call(add_to_process,list(p_cmd = c("./extract_genomic_regions.sh"), p_args = c(genome_path, gtf_path, file.path("files","bed",paste(org,"_genes",sep = "")), gene_list, org), p_wait = SUBPROCESS_WAIT))
    #process_count = process_count + 1
  }else{
    if(file.exists(genome_path) && file_ext(genome_path) == "gz"){
      file_move(genome_path, paste(GENOMES_PATH, "/",org,".fa.gz",sep = ""))
    }
    if(file.exists(gtf_path) && file_ext(gtf_path) == "gz"){
      file_move(gtf_path, paste(ANNOS_PATH, "/",org,".gtf.gz",sep = ""))
    }}
  #}
  Sys.sleep(5) 
}
fetch_genome_user <- function(data){
  genome_path <- c()
  gtf_path <- c()
  print(data)
  if(data$genome == "-" || data$gtf == "-"){
    fetch_genome_ensembl(data$org)
    break;
  }
  genome_path<-paste(GENOMES_PATH, "/",data$org,".fa.gz",sep = "")
  gtf_path<-paste(ANNOS_PATH, "/",data$org,".gtf.gz",sep = "")  
  fasta_path <- file.path(OUT_PATH ,data$org)
  #print(gtf_path)

  if(grepl("http",data$genome)){
    if(!file.exists(genome_path) || !file.info(genome_path)$size > 0){
      curl_fetch_disk(data$genome, paste(GENOMES_PATH, "/",basename(URLdecode(data$genome)),sep=""))
      genome_path <- paste(GENOMES_PATH, "/",basename(URLdecode(data$genome)),sep="")
    }else{
      if(system2("gzip",args = c("-t", genome_path), wait = T,stdout = NULL, stderr = NULL) == 0){
        genome_path<-genome_path #paste(GENOMES_PATH, "/",org,".fa.gz",sep = "")        
      }else{
        lapply(list.files(GENOMES_PATH, full.names = T, ignore.case = T, no.. = T, pattern = str_split_fixed(data$org,"_",2)[1]),function(x){
          if(grepl(x,pattern = "fa|gz")){
            try(file.remove(x, showWarnings=F))
          }
        })
        curl_fetch_disk(data$genome, basename(URLdecode(data$genome)))
        genome_path <- paste(GENOMES_PATH, "/",basename(URLdecode(data$genome)),sep="") #basename(URLdecode(data$genome))
      }
    }
    
  }else{
    genome_path <- data$genome
  }
  
  if(grepl("http",data$gtf)){
    #print(paste(ANNOS_PATH, "/",basename(URLdecode(data$gtf)),sep=""))
    if(!file.exists(gtf_path) || !file.info(gtf_path)$size > 0){
      curl_fetch_disk(data$gtf, paste(ANNOS_PATH, "/",basename(URLdecode(data$gtf)),sep=""))
      gtf_path <- paste(ANNOS_PATH, "/",basename(URLdecode(data$gtf)),sep="")
    }else{
      if(system2("gzip",args = c("-t", gtf_path), wait = T,stdout = NULL, stderr = NULL) == 0){
        gtf_path<-gtf_path        
      }else{
        lapply(list.files(ANNOS_PATH, full.names = T, ignore.case = T, no.. = T, pattern = str_split_fixed(data$org,"_",2)[1]),function(x){
          if(grepl(x,pattern = "gtf|gz")){
            try(file.remove(x, showWarnings=F))
          }
        })
        curl_fetch_disk(data$gtf, paste(ANNOS_PATH, "/",basename(URLdecode(data$gtf)),sep=""))
        gtf_path <- paste(ANNOS_PATH, "/",basename(URLdecode(data$gtf)),sep="")
      }
    }
    
  }else{
    gtf_path <- data$gtf
  }
  ##CONVERT GFFS to GTFS
  #gffread files/annos/XENTR_9.1_Xenbase.v10-lift.gff3 -T -O -E -o files/annos/XENTR_9.1_Xenbase.v10-lift.gtf
  print(genome_path)
  print(gtf_path)
  
  #print(paste(file_path_sans_ext(gtf_path),".gtf",sep = ""))
  #print(tolower(file_ext(gtf_path)))
  #if(tolower(file_ext(gtf_path))!="gtf"){
  #  system2("gffread",args = c(gtf_path, "-T", "-O", "-E", "-o", paste(file_path_sans_ext(gtf_path),".gtf",sep = "")),wait = T)
  #  gtf_path <- paste(file_path_sans_ext(gtf_path),".gtf",sep = "")
  #}
  if(CLEAN_EXTRACT || !check_files(fasta_path)){
    #process_list <- append(process_list,
    #system2("./extract_genomic_regions.sh",args = c(genome_path, gtf_path, file.path("files","bed",paste(data$org,"_genes",sep = "")), gene_list, data$org),stderr = T, stdout =  T,wait = SUBPROCESS_WAIT)
    do.call(add_to_process,list(p_cmd = c("./extract_genomic_regions.sh"), p_args = c(genome_path, gtf_path, file.path("files","bed",paste(data$org,"_genes",sep = "")), gene_list, data$org), p_wait = SUBPROCESS_WAIT))
    #process_count = process_count + 1
  }else{
    if(file.exists(genome_path) && file_ext(genome_path) == "gz"){
      file_move(genome_path, paste(GENOMES_PATH, "/",org,".fa.gz",sep = ""))
    }
    if(file.exists(gtf_path) && file_ext(gtf_path) == "gz"){
      file_move(gtf_path, paste(ANNOS_PATH, "/",org,".gtf.gz",sep = ""))
    }}
}
#ENTRY POINT
process_list <<- c()
process_count <<- 0

#try(file.remove("files/extraction_logs.e"))
#try(file.remove("files/extraction_logs.o"))

param_file <- "parameters.txt"
param_table <- read.table(param_file,sep="=")

GENOMES_PATH <- param_table[which(param_table=="genomes_path"),c(2)]
ANNOS_PATH <- param_table[which(param_table=="annos_path"),c(2)]
gene_list <- args[1]
OUT_PATH <- param_table[which(param_table=="fasta_path"),c(2)]
CLEAN_DOWNLOAD <- as.logical(param_table[which(param_table=="clean_download"),c(2)])
CLEAN_EXTRACT <- as.logical(param_table[which(param_table=="clean_extract"),c(2)])
SUBPROCESS_WAIT <- as.logical(param_table[which(param_table=="subprocess_wait"),c(2)])
BLAST_REGION <- as.character(param_table[which(param_table=="blast_region"),c(2)])
GENOMES_SOURCE <- as.character(param_table[which(param_table=="genomes_source"),c(2)])
REMOVE_DOWNLOADS <- as.character(param_table[which(param_table=="remove_downloads"),c(2)])
USER_GENOMES <- as.character(param_table[which(param_table=="user_genomes"),c(2)])

print(paste("CLEAN_DOWNLOAD:",CLEAN_DOWNLOAD))
print(paste("CLEAN_EXTRACT:",CLEAN_EXTRACT))
print(paste("REMOVE DOWNLOADS:",REMOVE_DOWNLOADS))
print(paste("GENOMES_SOURCE:",GENOMES_SOURCE))
print(paste("SUBPROCESS_WAIT:",SUBPROCESS_WAIT))

tryCatch(numWorkers <- detectCores(all.tests = T, logical = T), error=function(){numWorkers=2})

genes <- gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))

org.meta <- c()
for(org in getENSEMBLGENOMESInfo()$name){
  org.meta <- rbind(org.meta, is.genome.available(db = "ensembl", org, details = T))
}
#getGenomeSet(db = "ensembl", org.meta$name, reference = F, clean_retrieval = T, gunzip = T, update = T, path = "files/genomes")
check_list <- list()
user_data <- read.csv(USER_GENOMES,header = F)
names(user_data) <- c("org","genome","gtf")

if(!CLEAN_DOWNLOAD && !CLEAN_EXTRACT){
  #print("SKIPPING download and extraction of genomes. Check clean_extract and clean_download parameters in parameters.txt")
  print("Checking files...")
  for (fasta_dir in list.dirs(OUT_PATH)) {
    check_result <- check_files(fasta_dir)
    if(!check_result){
      check_list <- append(check_list,list(c(basename(fasta_dir))))
    }
    print(paste(fasta_dir," : ", check_result))
  }
  print("Checked files..Organisms which failed")
  print(unlist(check_list)) #
  #check_ens <- na.omitmatch(unlist(check_list),org.meta$name)
  #if(!is.null(unlist(check_list))){
  org.meta <- org.meta[na.omit(match(unlist(check_list),org.meta$name)),]
  user_data <- user_data[na.omit(match(unlist(check_list),user_data$org)),]
  #}
}

if(tolower(GENOMES_SOURCE)=="both" || tolower(GENOMES_SOURCE)=="user"){
  if(nrow(user_data)!=0){
    for (row in 1:nrow(user_data)) {
      #print(user_data[row,])
      fetch_genome_user(user_data[row,])
    }
  }
}

if(tolower(GENOMES_SOURCE)=="both" || tolower(GENOMES_SOURCE)=="ensembl"){
  mclapply(org.meta$name, fetch_genome_ensembl, mc.cores = 1)
}

for(proc in process_list){
  if(proc$is_alive()){
    proc$wait(timeout = -1)
  }
  print(proc$get_exit_status())
}


#print(process_list)
#print(paste("Processes compeleted successfully : ",is.null(process_list[which(process_list!=0)])))
#if(SUBPROCESS_WAIT==F){
#  Sys.sleep(300)
#}
sessionInfo()