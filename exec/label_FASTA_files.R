suppressMessages(require(parallel))
suppressMessages(require(Biostrings))
suppressMessages(require(stringi))
suppressMessages(require(furrr))
suppressMessages(require(purrr))
#suppressMessages(require(COMPLETE))

##FUNCTIONS
deduplicate_FASTA <- function(fasta_file, duplicates.method, n_threads=4) {
  seq_set <- Biostrings::readDNAStringSet(filepath = fasta_file,format = "fasta",use.names = T)

  seq_set <- dplyr::bind_rows(x = parallel::mclapply(split(seq_set,factor(names(seq_set))),function(x){
    if(stringi::stri_cmp_eq(duplicates.method,"merge")){

      merged_seq <- Biostrings::DNAString(gsub("[[:space:]]", "", paste(x,collapse="")))
      merged_seq_name <- unique(names(x))
      return(data.frame(seq_name=merged_seq_name, seq=paste(merged_seq)))
    }else if(stringi::stri_cmp_eq(duplicates.method,"make_unique")){
      if(length(x)==1){
        return(data.frame(seq_name=names(x), seq=paste(x)))
      }
      uniq_seqs <- unique(paste(x))
      seq_match_val <- match(uniq_seqs,x) #match(x,uniq_seqs)

      if(any(!is.na(seq_match_val)) && length(unique(seq_match_val)) == 1){ #length(uniq_seqs) == 1
        uniq_seqs <- gsub("[[:space:]]", "", paste(uniq_seqs)) # Biostrings::DNAString(
        uniq_seq_name <- unique(names(uniq_seqs))
        names(uniq_seqs) <- uniq_seq_name
        return(data.frame(seq_name=uniq_seq_name, seq=paste(uniq_seqs)))
      }else{
        uniq_seq_name <- names(x)[seq_match_val]
        if(any(duplicated(uniq_seq_name))){
          uniq_seq_name <- paste(paste("block",1:length(uniq_seq_name),sep=""),uniq_seq_name,sep=".")
        }
        #uniq_seqs <- list(uniq_seqs)
        names(uniq_seqs) <- uniq_seq_name
        return(data.frame(seq_name=uniq_seq_name, seq=paste(uniq_seqs)))
      }
    }else if(stringi::stri_cmp_eq(duplicates.method,"delete")){
      x <- x[!which(duplicated(paste(x)))]
      x <- x[!which(duplicated(names(x)))]
      if(nrow(x) > 0){
        return(data.frame(seq_name=names(x), seq=paste(x)))
      }
    }
  }, mc.cores = n_threads,mc.silent = T))

  seq_set_tmp <- Biostrings::DNAStringSet(seq_set[,2], use.names = F)
  names(seq_set_tmp) <- seq_set[,1]
  seq_set <- seq_set_tmp
  return(seq_set)
}

label_sequenceIDs <- function(fasta_path,org,gene_list,odb_gene_map=NULL,params_list, duplicates.method="merge"){

  if(is.null(duplicates.method) || !grepl(pattern = c("merge|delete|make_unique"), x = duplicates.method,ignore.case = T)){
    stop(paste("duplicates.method must be one of merge|delete|make_unique"))
  }

  if(!is.null(odb_gene_map)){
    if(file.exists(odb_gene_map) && file.info(odb_gene_map)$size > 0){
      odb_gene_map <- read.table(file = odb_gene_map,header = F,quote = "",sep = "\t")
      #local ortho_cluster=$(grep -w $gene_name $odb_clusters | awk -F'\t' '{if (length(c) == 0){c=$1;}else{c=c","$1;}}END{print c}')
    }else{
      warning(paste(odb_gene_map,"does not exist!"))
      odb_gene_map <- NULL
    }
  }
  if (is.character(gene_list)) {
    genes <- factor(scan(gene_list, character(),quiet = T)) #gsub('[[:punct:] ]+','_', factor(scan(gene_list, character())))
    genes <- tolower(genes[grep("gene",tolower(genes), invert = T, fixed = T)])
  }else{
    genes <- tolower(as.vector(gene_list))
  }

  invisible( future_map(genes, function(x){ #invisible(mclapply
    #print(x) #DEBUG
    safe_gene <- gsub(pattern="[[:punct:]]", replacement = "_",tolower(x))
    #file_path <- paste(fasta_path,safe_gene,sep="")
    odb_clusters <- "ungrouped"
    if(!is.null(odb_gene_map)){
      odb_clusters <- paste(odb_gene_map[grep(pattern = x,x = odb_gene_map[,2],ignore.case = T,value = F),1],collapse = ",")
      if(stringi::stri_isempty(odb_clusters)){
        odb_clusters <- "ungrouped"
      }
    }
    #print(paste(x, safe_gene,odb_clusters)) #DEBUG
    future_map(list.files(path = fasta_path,pattern = safe_gene,full.names = T,ignore.case = T), function(y){ #mclapply
      #if(file.exists(y)){
      seq_set <- deduplicate_FASTA(fasta_file=y, duplicates.method=duplicates.method, n_threads=params_list$numWorkers)

        split_seq_names <- stringi::stri_split(str = names(seq_set), fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)
        if( ncol(split_seq_names) == 1 ){ #length(COMPLETE$FORMAT_ID_INDEX)
          #print(paste(names(seq_set),org,x,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)) #DEBUG
          #print(names(seq_set)) #DEBUG
          names(seq_set) <- paste(names(seq_set),org,x,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)
          writeXStringSet(x = seq_set,filepath = y,append = F,format = "fasta")
          #return(paste(names(seq_set),org,x,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM)) #DEBUG
        }else if( ncol(split_seq_names) < 4 && ncol(split_seq_names) >= 1){
          names(seq_set) <- paste(split_seq_names[,1] ,org,x,odb_clusters,sep = params_list$SEQUENCE_ID_DELIM) #stringi::stri_split(str = split_seq_names, fixed = params_list$SEQUENCE_ID_DELIM, simplify=T)[,1]
          writeXStringSet(x = seq_set,filepath = y,append = F,format = "fasta")
        }else if(ncol(split_seq_names) != 4){
          message(paste("Error : Check sequence ID format of", y))
        }
      #}
    },.options = furrr::furrr_options(seed=T, scheduling=params_list$numWorkers))#, mc.cores = params_list$numWorkers,mc.silent = T)

  },.options = furrr::furrr_options(seed=T, scheduling=params_list$numWorkers)) ) #, mc.cores = params_list$numWorkers,mc.silent = T))

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

invisible(label_sequenceIDs(fasta_path = org_fasta_path,org = org_name,gene_list = gene_list,odb_gene_map = odb_gene_map,params_list = params_list,duplicates.method="merge"))
