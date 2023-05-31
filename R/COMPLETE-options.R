
#' Internal Function - Check and Install GNU parallel
#'
#' Install GNU parallel if not available
#'
#' @return Path to GNU parallel
install_parallel <- function(){
  if (stringi::stri_isempty(Sys.which("parallel")) && !file.exists(paste(fs::path_home(),"/bin/parallel",sep=""))) {
    install_status <- processx::run( command = COMPLETE_env$SHELL ,args=c(fs::path_package("COMPLETE","exec","functions.sh"),"install_parallel") ,spinner = T,stdout = "",stderr = "")
    if(install_status$status>0){
      stop("Problem installing GNU parallel. Is $SHELL set? check permissions in $HOME directory and check https://www.gnu.org/software/parallel/parallel_tutorial.html")
    }else{
      return(paste(fs::path_home(),"/bin/parallel",sep=""))
    }
  }else if(!stringi::stri_isempty(Sys.which("parallel"))){
    return(Sys.which("parallel"))
  }else if(file.exists(paste(fs::path_home(),"/bin/parallel",sep=""))){
    return(paste(fs::path_home(),"/bin/parallel",sep=""))
  }else{
    stop("Problem with GNU parallel installation. Is $SHELL set? check permissions in $HOME/bin directory and check https://www.gnu.org/software/parallel/parallel_tutorial.html")
  }
}

#' Connect to a biomaRt repo/dataset
#'
#' This function tries to connect to BIOMART services (using biomaRt) and retries if the connection fails with a stepping timeout of 5secs.
#' The MART_FUN is called with do.call()
#'
#' @param MART_FUN Fuction to be called
#' @param args Named list of arguments to be used by the function
#' @return Name of the dataset in biomaRt
mart_connect <- function(MART_FUN=NULL,args=c(),verbose=F){
  if (!is.null(MART_FUN)) {
    #print(args)
    time_out <- 0
    while (time_out < 600) { #max time out is 10mins, code will fibonacci to it stepping up 5 seconds
      time_out <- time_out+5
      catch_value <- tryCatch(do.call(MART_FUN,args),
                              error=function(cond){
                                if(verbose){
                                  message(cond)
                                  #print(class(cond))
                                  message(paste("\nWaiting for",time_out,"s & trying again...\n"))
                                }
                                Sys.sleep(time_out)
                              })
      if(!any(grepl(pattern = "error|exception|try-error|try|fail|timeout",ignore.case = T,x = class(catch_value))) && !is.null(catch_value)){
        #print("Done!")
        return(catch_value)
      }
    }
    stop("\nEnsembl is unresponsive, please check connection\n")
  }else{
    return(0)
  }
}

#' Checks for parameter and return the value
#'
#' If a parameter ID does not exist or if the value is empty and if the parameter is NOT optional then the execution is stopped
#'
#' @param param_table parameter table or parameter file name
#' @param param_id ID of the paramter to be extracted from parameter file
#' @param gene_list File with a list of genes to extract data for(check the github repo for an example)
check_param <- function(param_table,param_id,optional=F,CAST_FUN=as.character,create_dir=F){
  if(length(param_table) == 1 && file.exists(param_table)){
    param_table <- read.table(textConnection(gsub("==", "^", readLines(param_table))),sep="^", header = T) #Convert multibyte seperator to one byte sep
  }
  param_index <- which(param_table==param_id)
  if(is.null(param_index) || is.na(param_index)){
    stop(paste("Parameter :",param_id,"is missing!"))
  }

  if(nrow(param_table[param_index,]) > 1){
    stop(paste("Parameter :",param_id,"has multiple entries!"))
  }
  param_value <- param_table[param_index,c(2)]

  if(create_dir){
    dir.create(CAST_FUN(param_value),showWarnings = F,recursive = T)
  }
  if(!stringi::stri_isempty(param_value) || optional){
    return(CAST_FUN(param_value))
  }else{
    stop(paste("Parameter :",param_id,"is empty and is not optional!"))
  }
}

#' Loads the parameters from param file
#'
#' The function loads parameters from the param file into Global R variables and returns a named parameter table. The parameters are checked and an error is output if there are any missing parameters or if the folders are not accesible
#'
#' @examples
#'     params_list <- load_params(param_file)
#'
#' @note Parameters are case sensitive and are de-limited with '=='.
#'
#' @param param_file parameter file name
#' @return Named list of parameters
#' @export
load_params <- function(param_file){
  if(!is.character(param_file)){
    message("Give the file name of the parameters")
    return(NULL)
  }

  # if (grepl(pattern = "bash",ignore.case = T,x = Sys.getenv("SHELL"))) {
  #   SHELL <<- Sys.getenv("SHELL")
  #   print(paste("SHELL :",SHELL))
  # }else if(file.exists("/bin/bash")){
  #   SHELL <<- "/bin/bash"
  #   print(paste("SHELL :",SHELL))
  # }else{
  #   stop(paste("SHELL (",SHELL,") : bash not available, or not in $PATH or SHELL=/bin/bash not set"))
  # }

  #param_file <<- param_file
  param_table <- read.table(textConnection(gsub("==", "^", readLines(param_file))),sep="^", header = T) #Convert multibyte seperator to one byte sep

  param_dups <- duplicated(param_table$param_id)
  if(any(param_dups)){
    stop(paste("Duplicate entries in parameter file :", paste(unique(param_table$param_id[param_dups]),collapse=",") ))
  }

  GENOMES_PATH <- tools::file_path_as_absolute(check_param(param_table,"genomes_path",optional=F,CAST_FUN=as.character,create_dir=T))  #param_table[which(param_table=="genomes_path"),c(2)]
  ANNOS_PATH <- tools::file_path_as_absolute(check_param(param_table,"annos_path",optional=F,CAST_FUN=as.character,create_dir=T))  #param_table[which(param_table=="annos_path"),c(2)]
  TEMP_PATH <- tools::file_path_as_absolute(check_param(param_table,"temp_path",optional=T,CAST_FUN=as.character,create_dir=T)) #param_table[which(param_table=="temp_path"),c(2)]
  if(stringi::stri_isempty(TEMP_PATH)){
    TEMP_PATH <- tempdir()
  }
  GROUPS_PATH <- tools::file_path_as_absolute(check_param(param_table,"groups_path",optional=F,CAST_FUN=as.character,create_dir=T)) #param_table[which(param_table=="groups_path"),c(2)]
  FASTA_OUT_PATH <- tools::file_path_as_absolute(check_param(param_table,"fasta_path",optional=F,CAST_FUN=as.character,create_dir=T))  #param_table[which(param_table=="fasta_path"),c(2)]
  OUT_PATH <- tools::file_path_as_absolute(check_param(param_table,"out_path",optional=F,CAST_FUN=as.character,create_dir=T))
  # if(stringi::stri_isempty(OUT_PATH)){
  #   OUT_PATH <<- "."
  # }
  CLEAN_EXTRACT <- check_param(param_table,"clean_extract",optional=F,CAST_FUN=as.logical) #as.logical(param_table[which(param_table=="clean_extract"),c(2)])
  TRANSCRIPT_ID_DELIM <- check_param(param_table,"transcript_delimiter",optional=F,CAST_FUN=as.character) #param_table[which(param_table=="transcript_delimiter"),c(2)]
  SEQUENCE_ID_DELIM <- check_param(param_table,"seqID_delimiter",optional=F,CAST_FUN=as.character) #param_table[which(param_table=="seqID_delimiter"),c(2)]
  # DATA_SOURCE <- check_param(param_table,"data_source",optional=F,CAST_FUN=as.character) #tolower(param_table[which(param_table=="data_source"),c(2)])
  # if(grepl(pattern="both|user",DATA_SOURCE,ignore.case = T)){
  #   USER_GENOMES <-  tools::file_path_as_absolute(check_param(param_table,"user_genomes",optional=F,CAST_FUN=as.character)) #as.character(param_table[which(param_table=="user_genomes"),c(2)])
  # }else{
  #   USER_GENOMES <- ""
  # }
  TRANSCRIPT_REGIONS <- tolower(gsub("[[:space:]]","",x = unlist(stringi::stri_split( check_param(param_table,"transcript_regions",optional=F,CAST_FUN=as.character) ,fixed = ",")))) #param_table[which(param_table=="transcript_regions"),c(2)]
  STRAND <-  check_param(param_table,"strand",optional=F,CAST_FUN=as.character) #param_table[which(param_table=="strand"),c(2)]
  max_concurrent_jobs <- check_param(param_table,"max_concurrent_jobs",optional=T,CAST_FUN=as.numeric) #as.numeric(param_table[which(param_table=="max_concurrent_jobs"),c(2)])
  if(is.na(max_concurrent_jobs) || length(max_concurrent_jobs) == 0){
    max_concurrent_jobs <- parallel::detectCores(all.tests = T, logical = T)
  }
  tryCatch(numWorkers <- parallel::detectCores(all.tests = T, logical = T), error=function(){numWorkers <- 2})
  if(is.na(numWorkers) || numWorkers > max_concurrent_jobs){
    numWorkers <- max_concurrent_jobs
  }

  GENE_SEARCH_MODE <- check_param(param_table,"gene_search_mode",optional=T,CAST_FUN=as.character)
  if(is.na(GENE_SEARCH_MODE) || length(GENE_SEARCH_MODE) == 0){
    GENE_SEARCH_MODE <- "HARD"
  }
  SELECT_ALL_GENES <- check_param(param_table,"select_all_genes_from_cluster",optional=T,CAST_FUN=as.logical)
  if(is.na(SELECT_ALL_GENES) || length(SELECT_ALL_GENES) == 0){
    SELECT_ALL_GENES <- FALSE
  }

  SELECT_REF_ORG_GROUPS <- check_param(param_table,"select_groups_with_ref_orgs",optional=T,CAST_FUN=as.logical)
  if(is.na(SELECT_REF_ORG_GROUPS) || length(SELECT_REF_ORG_GROUPS) == 0){
    SELECT_REF_ORG_GROUPS <- FALSE
  }

  SELECT_REF_ORG_GROUPS_METHOD <- toupper(check_param(param_table,"select_groups_with_ref_orgs_method",optional=T,CAST_FUN=as.character))
  if(is.na(SELECT_REF_ORG_GROUPS_METHOD) || length(SELECT_REF_ORG_GROUPS_METHOD) == 0){
    SELECT_REF_ORG_GROUPS_METHOD <- "ANY"
  }

  E_VALUE_THRESH <- check_param(param_table,"e_value_thresh",optional=F,CAST_FUN=as.double)
  MIN_IDENT_THRESH  <- check_param(param_table,"minIdent_thresh",optional=F,CAST_FUN=as.double)
  BLAST_OPTIONS  <- check_param(param_table,"blast_options",optional=F,CAST_FUN=as.character)
  TRANSCRIPT_REGIONS_DELIMITER <- check_param(param_table,"transcript_regions_delimiter",optional=F,CAST_FUN=as.character)
  #BLAST_DB_PATH <- tools::file_path_as_absolute(check_param(param_table,"blastdb_path",optional=F,CAST_FUN=as.character,create_dir=T))
  REF_ORGS_FILE <- tools::file_path_as_absolute(check_param(param_table,"ref_orgs",optional=F,CAST_FUN=as.character))
  #CLUSTER_OPTIONS  <- check_param(param_table,"cluster_options",optional=T,CAST_FUN=as.character)
  BED_PATH <- tools::file_path_as_absolute(check_param(param_table,"bed_path",optional=F,CAST_FUN=as.character,create_dir=T))
  PLOT_PATH <- tools::file_path_as_absolute(check_param(param_table,"plot_path",optional=F,CAST_FUN=as.character,create_dir=T))
  GENE_DROP_THRESHOLD <- check_param(param_table,"gene_drop_thresh",optional=F,CAST_FUN=as.double)
  ORTHODB_PREFIX <- check_param(param_table,"orthodb_path_prefix",optional=T,CAST_FUN=as.character)
  MACSE_PATH <- tools::file_path_as_absolute(check_param(param_table,"macse_path",optional=F,CAST_FUN=as.character))
  MAFFT_PATH <- tools::file_path_as_absolute(check_param(param_table,"mafft_path",optional=F,CAST_FUN=as.character))
  TRANSAT_PATH <- tools::file_path_as_absolute(check_param(param_table,"transat_path",optional=F,CAST_FUN=as.character))
  RNADECODER_PATH <- tools::file_path_as_absolute(check_param(param_table,"rnadecoder_path",optional=F,CAST_FUN=as.character))
  #python3_path  <- check_param(param_table,"python3_path",optional=F,CAST_FUN=as.character)
  FASTTREE_PATH <- tools::file_path_as_absolute(check_param(param_table,"fasttree_path",optional=F,CAST_FUN=as.character))
  ALIGNMENTS_PATH <- tools::file_path_as_absolute(check_param(param_table,"alignments_path",optional=F,CAST_FUN=as.character,create_dir=T))
  ALIGNMENT_GAP_THRESHOLD  <- check_param(param_table,"aln_gap_thres",optional=F,CAST_FUN=as.double)
  MIN_COVERAGE_THRESHOLD <- check_param(param_table,"min_coverage_thres",optional=F,CAST_FUN=as.double)

  COMPLETE_env$numWorkers <- max_concurrent_jobs
COMPLETE_env$USE_ORTHODB <- !stringi::stri_isempty(ORTHODB_PREFIX)

  #COMPLETE_env <- new.env(parent=emptyenv())
  #COMPLETE_env$PARAMS <-
  param_list <- list(param_file=param_file,
                     GENOMES_PATH=GENOMES_PATH,
                     ANNOS_PATH=ANNOS_PATH,
                     TEMP_PATH=TEMP_PATH,
                     GROUPS_PATH=GROUPS_PATH,
                     FASTA_OUT_PATH=FASTA_OUT_PATH,
                     OUT_PATH=OUT_PATH,
                     #USER_GENOMES=USER_GENOMES,
                     CLEAN_EXTRACT=CLEAN_EXTRACT,
                     TRANSCRIPT_ID_DELIM=TRANSCRIPT_ID_DELIM,
                     SEQUENCE_ID_DELIM=SEQUENCE_ID_DELIM,
                     #DATA_SOURCE=DATA_SOURCE,
                     TRANSCRIPT_REGIONS=TRANSCRIPT_REGIONS,
                     STRAND=STRAND,
                     max_concurrent_jobs=max_concurrent_jobs,
                     numWorkers=numWorkers,
                     GENE_SEARCH_MODE=GENE_SEARCH_MODE,
                     SELECT_ALL_GENES=SELECT_ALL_GENES,
                     SELECT_REF_ORG_GROUPS=SELECT_REF_ORG_GROUPS,
                     SELECT_REF_ORG_GROUPS_METHOD=SELECT_REF_ORG_GROUPS_METHOD,
                     E_VALUE_THRESH=E_VALUE_THRESH,
                     MIN_IDENT_THRESH=MIN_IDENT_THRESH,
                     BLAST_OPTIONS=BLAST_OPTIONS,
                     TRANSCRIPT_REGIONS_DELIMITER=TRANSCRIPT_REGIONS_DELIMITER,
                     #BLAST_DB_PATH=BLAST_DB_PATH,
                     REF_ORGS_FILE=REF_ORGS_FILE,
                     BED_PATH=BED_PATH,
                     PLOT_PATH=PLOT_PATH,
                     GENE_DROP_THRESHOLD=GENE_DROP_THRESHOLD,
                     ORTHODB_PREFIX=ORTHODB_PREFIX,
                     MACSE_PATH=MACSE_PATH,
                     MAFFT_PATH=MAFFT_PATH,
                     TRANSAT_PATH=TRANSAT_PATH,RNADECODER_PATH=RNADECODER_PATH,
                     FASTTREE_PATH=FASTTREE_PATH,
                     ALIGNMENTS_PATH=ALIGNMENTS_PATH,
                     ALIGNMENT_GAP_THRESHOLD=ALIGNMENT_GAP_THRESHOLD,
                     MIN_COVERAGE_THRESHOLD=MIN_COVERAGE_THRESHOLD,
                     PARAMETERS_LOADED=TRUE)

  #return(COMPLETE_env$PARAMS)
  if(stringi::stri_isempty(COMPLETE_env$BLAST_BIN)){
    message("Warning: NCBI-BLAST path is empty")
  }
  
  class(param_list) <- c(class(param_list), "COMPLETE-options")
  return(param_list)
}

#' INITIALIZES THE R-COMPLETE LIBRARY
#'
#' First function to run during package load
#' 
INITIALIZE <- function() {
  options(RCurlOptions = list(ssl.verifyhost=0, ssl.verifypeer=0, timeout=200,maxconnects=200,connecttimeout=200)) #ssl.verifyhost=0, ssl.verifypeer=0,
  #options(download.file.method="curl")
  
  COMPLETE_env <<- new.env(parent=baseenv())
  
  COMPLETE_env$ENSEMBL_MART <- "ENSEMBL_MART_ENSEMBL"
  
  if (curl::has_internet()) {
    COMPLETE_env$using.mart <- mart_connect(biomaRt::useMart,args=list(COMPLETE_env$ENSEMBL_MART)) #For biomaRt
    COMPLETE_env$org.meta.list <- mart_connect(biomaRt::listDatasets,args=list(mart=COMPLETE_env$using.mart)) #For biomaRt
    COMPLETE_env$org.meta <- mart_connect(biomartr::listGenomes,args=list(db = "ensembl", type = "all", details = T)) #For biomartr #db = tolower(GENOMES_SOURCE)
  }else{
    warning("Warning : Check if there is an internet connection")
  }
  
  COMPLETE_env$curl_handle <- RCurl::getCurlMultiHandle()
  COMPLETE_env$process_list <- c()
  COMPLETE_env$SKIP_USER_DATA <- FALSE
  COMPLETE_env$USE_ORTHODB <- FALSE
  
  #if(!grepl(x=Sys.info()["sysname"],pattern="linux",ignore.case = T)){
  #  stop("R-COMPLETE Pipeline only supports Linux (and bash) :(")
  #}
  
  if (grepl(pattern = "bash",ignore.case = T,x = Sys.getenv("SHELL"))) {
    COMPLETE_env$SHELL <- Sys.getenv("SHELL")
    cat(paste("SHELL :",COMPLETE_env$SHELL,"\n"))
  }else if(file.exists("/bin/bash")){
    COMPLETE_env$SHELL <- "/bin/bash"
    cat(paste("SHELL :",COMPLETE_env$SHELL,"\n"))
  }else{
    stop(paste("SHELL : bash not available, or not in $PATH or SHELL=/bin/bash not set"))
  }
  
  COMPLETE_env$user_home <- fs::path_home()
  COMPLETE_env$parallel <- install_parallel()
  COMPLETE_env$numWorkers <- tryCatch(parallel::detectCores(all.tests = T, logical = T), error=function(cond){return(2)})
  COMPLETE_env$max_file_handles <- as.numeric(processx::run(command = COMPLETE_env$SHELL, args = c("-c","ulimit -n"))$stdout)
  COMPLETE_env$BLAST_BIN <- dirname(Sys.which("tblastx"))
  COMPLETE_env$FORMAT_ID_INDEX <- list(TRANSCRIPT_ID=1,ORG=2,GENE=3,CLUSTERS=4)
  Sys.chmod(fs::path_package("COMPLETE","exec","functions.sh"), "777", use_umask = FALSE)
  
  if(stringi::stri_isempty(COMPLETE_env$BLAST_BIN)){
    message("Warning: NCBI-BLAST path is empty")
  }
}

#.onAttach <- function(libname, pkgname) {
#  # to show a startup message
#  packageStartupMessage("Initializing R-COMPLETE")
#  tryCatch({INITIALIZE()}, error=function(cond) {stop(cond)})
#  packageStartupMessage("R-COMPLETE loaded successfully!")
#}

.onLoad <- function(libname, pkgname) {
  # something to run
  packageStartupMessage("Initializing R-COMPLETE")
  tryCatch({
    INITIALIZE()
    }, error=function(cond) {stop(cond)})
  packageStartupMessage("R-COMPLETE loaded successfully!")
}

.onUnload <- function(libname, pkgname){
  ##rm(COMPLETE_env)
  if(exists("COMPLETE_env", mode="environment")){
    COMPLETE_env <<- new.env(parent=emptyenv())
  }
}

#' @importFrom S4Vectors mcols
NULL
