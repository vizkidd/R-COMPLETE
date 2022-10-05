suppressMessages(require(parallel))
suppressMessages(require(Biostrings))
suppressMessages(require(stringi))

##FUNCTIONS
label_sequenceIDs <- function(fasta_path,org,gene_list,odb_gene_map=NULL,params_list){
  if(!is.null(odb_gene_map)){
    if(file.exists(odb_gene_map) && file.info(odb_gene_map)$size > 0){
      odb_gene_map <- read.table(file = odb_gene_map,header = F,quote = "",sep = "\t")
      #local ortho_cluster=$(grep -w $gene_name $odb_clusters | awk -F'\t' '{if (length(c) == 0){c=$1;}else{c=c","$1;}}END{print c}')
    }
  }
  if (is.character(gene_list)) {
    genes <- factor(scan(gene_list, character(),quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- genes[grep("gene",tolower(genes), invert = T, fixed = T)]
  }else{
    genes <- as.vector(gene_list)
  }

  invisible(mclapply(genes, function(x){
    safe_gene <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(x))
    file_path <- paste(fasta_path,safe_gene,sep="")
    if(!is.null(odb_gene_map)){
      odb_clusters <- paste(odb_gene_map[grep(pattern = x,x = odb_gene_map[,2],ignore.case = T,value = F),1],collapse = ",")
    }else{
      odb_clusters <- "ungrouped"
    }
    lapply(list.files(path = fasta_path,pattern = safe_gene,full.names = T), function(y){
      seq_set <- readDNAStringSet(filepath = y,format = "fasta",use.names = T)
      names(seq_set) <- paste(names(seq_set),org,x,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)
      writeXStringSet(x = seq_set,filepath = y,append = F,format = "fasta")
    })

  }, mc.cores = params_list$numWorkers,mc.silent = T))

}

check_param <- function(param_table,param_id,optional=F,CAST_FUN=as.character,create_dir=F){
  if(is.character(param_table)){
    param_table <- read.table(textConnection(gsub("==", "^", readLines(param_table))),sep="^", header = T) #Convert multibyte seperator to one byte sep
  }
  param_index <- which(param_table==param_id)
  if(is.null(param_index) || is.na(param_index)){
    stop(paste("Parameter :",param_id,"is missing!"))
  }
  param_value <- param_table[param_index,c(2)]
  #print(CAST_FUN(param_value))
  if(create_dir){
    dir.create(CAST_FUN(param_value),showWarnings = F,recursive = T)
  }
  if(!stri_isempty(param_value) || optional){
    return(CAST_FUN(param_value))
  }else{
    stop(paste("Parameter :",param_id,"is empty and is not optional!"))
  }
}
##ENTRYPOINT
set.seed(123)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Give the (1) Organism FASTA Path (2) Name of the Organism (3) Gene List (4) Parameter List (5) Filename of OrthoDB Ortholog Clusters mapped to gene names (Optional). If not provided, sequences cluster ID/name 'ungrouped' is used", call.=FALSE)
}

org_fasta_path <- args[1]
org_name <- args[2]
gene_list <- args[3]
param_file <- args[4] #"parameters.txt"
if(length(args)==5){
  odb_gene_map <- args[5]
}else{
  odb_gene_map <- NULL
}

if(!file.exists(param_file) || file.info(param_file)$size <= 4){
  stop("ERROR: parameters.txt is missing and is required")
}

param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^", header = T) #Convert multibyte seperator to one byte sep #read.table(param_file,sep="==")

params_list <- c()

params_list$max_concurrent_jobs <- check_param(param_table,"max_concurrent_jobs",optional=T,CAST_FUN=as.numeric)
if(is.na(params_list$max_concurrent_jobs) || length(params_list$max_concurrent_jobs) == 0){
  params_list$max_concurrent_jobs <- detectCores(all.tests = T, logical = T)
}
params_list$TRANSCRIPT_ID_DELIM <- check_param(param_table,"transcript_delimiter",optional=F,CAST_FUN=as.character) #param_table[which(param_table=="transcript_delimiter"),c(2)]
params_list$SEQUENCE_ID_DELIM <- check_param(param_table,"seqID_delimiter",optional=F,CAST_FUN=as.character)
params_list$numWorkers <- detectCores(all.tests = TRUE, logical = TRUE)
if(is.na(params_list$numWorkers) || params_list$numWorkers > params_list$max_concurrent_jobs){
  params_list$numWorkers <- params_list$max_concurrent_jobs
}

label_sequenceIDs(fasta_path = org_fasta_path,org = org_name,gene_list = gene_list,odb_gene_map = odb_gene_map,params_list = params_list)
