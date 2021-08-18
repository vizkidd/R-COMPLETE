require(biomartr)
require(parallel)
require(tools) ##FOR removing file extensions (gz)
require(RCurl)
require(curl)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  #stop("Give the (1) folder path for storing genomes (2) folder path for storing annotations (3) the filename with gene symbols(without header)", call.=FALSE)
  stop("Give the (1) folder path for storing genomes (2) folder path for storing annotations (3) the filename with gene symbols(without header) (4) FASTA OUTPUT PATH(absolute path)", call.=FALSE)
}
set.seed(123)

options(RCurlOptions = list(ssl.verifyhost=0, ssl.verifypeer=0)) #qsub -l data -V ./run_pipeline_ensembl.sh files/trimlist.txt 
sessionInfo()

#DECLARATIONS
#GENOMES_PATH <- "files/genomes"
#ANNOS_PATH <- "files/annos"
#gene_list <- "genelist.txt"

GENOMES_PATH <- args[1]
ANNOS_PATH <- args[2]
gene_list <- args[3]
OUT_PATH <- args[4]

numWorkers <- detectCores(all.tests = T, logical = T) 

#FUNCTIONS
fetch_genome_ensembl <- function(org) {
  fasta_path <- file.path(OUT_PATH ,org) ## file_path_as_absolute(file.path(OUT_PATH ,org))
  print(fasta_path)
  details <- file.info(list.files(path = ANNOS_PATH, full.names = TRUE))
  #print(details)
  #print(rownames(details[details$size == 0, ]))
  if(!file.exists(paths = paste(ANNOS_PATH,org,".gtf", sep="/")))
  {
    #genome_path <- tryCatch(getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH),error = function(e){print("Sleeping...");Sys.sleep(360);},finally=getGenome(db = "ensembl", org, reference = F, gunzip = F, path = GENOMES_PATH))
    gtf_path <- tryCatch(getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH)),error = function(e){print("Sleeping...");Sys.sleep(360)},finally=getGTF(db = "ensembl", org, path = file.path(ANNOS_PATH)))
    #g_unzip.str <- paste("gunzip","-d", genome_path, sep = " ")
    #system2(command = "gunzip",args = c("-d", genome_path))
    #unzip.str <- paste("gunzip","-d", gtf_path, sep = " ")
    system2(command = "gunzip",args = c("-d","-f", gtf_path),wait=T)
    new_gtf_path <- paste(ANNOS_PATH,paste(org,".gtf",sep=""), sep="/")
    system2(command = "mv",args = c(tools::file_path_sans_ext(gtf_path),new_gtf_path),wait=T)
    #faidx.str <- paste("samtools","faidx", genome_path, sep = " ")
    #system2("samtools", args = c("faidx", tools::file_path_sans_ext(genome_path)))
    #script.str <- paste("./extract_genomic_regions.sh", genome_path, gtf_path, file.path("files","bed",paste(org,"_genes",sep = "")), gene_list, org, sep=" ")
    #system2("./extract_genomic_regions.sh",args = c(tools::file_path_sans_ext(genome_path), tools::file_path_sans_ext(gtf_path), file.path("files","bed",paste(org,"_genes",sep = "")), gene_list, org, OUT_PATH), wait = T)
    Sys.sleep(5)
  }
}

#ENTRY POINT
org.meta <- c()
for(org in getENSEMBLGENOMESInfo()$name){
  org.meta <- rbind(org.meta, is.genome.available(db = "ensembl", org, details = T))
}
#getGenomeSet(db = "ensembl", org.meta$name, reference = F, clean_retrieval = T, gunzip = T, update = T, path = "files/genomes")
mclapply(org.meta$name, fetch_genome_ensembl, mc.cores = 1)
