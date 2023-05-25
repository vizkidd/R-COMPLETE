# R - COMPLETE
Pipeline for extracting localization elements/motifs using a comparitive approach. Library can be installed and tested on R(=4.1). For a list of genes the pipeline downloads full transcript sequences from all the organisms (or selected organisms) from ENSEMBL or NCBI, Formats the headers and stores them according to the OrthoDB clusters. Sequences without clusters are clustered based on sequence identity/sequence coverage of Reciprocal BLAST Hits of CDS regions. Next step would be to stitch the transcript regions into full length transcripts and align them. Pipeline is format agnostic(to my knowledge - TMK), packaged with a multi-threaded BLAST framework, can load and convert BLAST formats (internally - to GRanges) and sports functions for performing Reciprocal Bidirectional Hits. The multithreaded BLAST framework can handle many-many BLAST Hits across organisms. The package is interfaced with bash, R transforms the data & calls the shots and bash handles files & BLASTing. Check requirements and installation instuctions before proceeding.

*Ironically this repo is incomplete but the functionality in it works. Under Construction Indefinitely. Documentation can be found within the package, Play around with the functions for the rest*

## Installation (on R - Linux or Docker with WSL in Windows):
    sudo apt-get update && sudo apt-get install curl bzip2 parallel 
    BiocManager::install(c("Rhtslib", "devtools", "BiocManager", "Biostrings", "biomaRt", "S4Vectors", "IRanges", "rtracklayer", "GenomicRanges", "BiocGenerics"))
    devtools::install_github("https://github.com/vizkidd/R-COMPLETE/")

## REQUIRES:
   * Linux with BASH ($SHELL must be set or /bin/bash must exist) (export SHELL="/bin/bash")
   * Parameters File (fs::path_package("COMPLETE","inst","parameters.txt"))
   * GNU parallel (in $PATH - BASH functions)
   * [Samtools](http://www.htslib.org/download/) (in $PATH - BASH functions)
   * [Bedtools](https://github.com/arq5x/bedtools2/releases) (in $PATH - BASH functions)
   * [ncbi-blast+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Compile from .src.tar.gz with ./configure && make all && sudo make install) (or sudo alien -i ncbi-blast-X.XX.X+-2.src.rpm) (or download binaries)
        * Check if you have the binaries for #*blastdb_path*# and #*makeblastdb*#
   * [Zlib](https://zlib.net/) - (Compile from sources) (or) (sudo apt install libz-dev or yum install zlib-devel)
   * LZMA SDK - (sudo apt-get install liblzma-dev or yum install xz-devel)
   * BZLIB - (sudo apt-get install libbz2-dev libclang-dev or yum install bzip2-devel.x86_64)
   * [OrthoDB (ODB) Flat Files (>= v10.1)](https://www.orthodb.org/?page=filelist) (Pipeline is tested with ODB v10.1) 
        * [odb10v1_species.tab.gz](https://v101.orthodb.org/download/odb10v1_species.tab.gz) - Ortho DB organism ids based on NCBI taxonomy ids (mostly species level) 
        * [odb10v1_genes.tab.gz](https://v101.orthodb.org/download/odb10v1_genes.tab.gz)  -Ortho DB genes with some info 
        * [odb10v1_OG2genes.tab.gz](https://v101.orthodb.org/download/odb10v1_OG2genes.tab.gz) - OGs to genes correspondence 
    (OR)
        * odb10v1_OGgenes_fixed.tab.gz - Merged & Transformed ODB file (Done within pipeline)
        * odb10v1_OGgenes_fixed_user.tab.gz - Merged & Transformed ODB file BASED on user gene list (Done within pipeline)

## Documentation
* ?COMPLETE_PIPELINE_DESIGN (in R docs)

### PARAMETERS :
   The pipeline takes a single parameter file. This design was chosen,
   + To expose as many options as possible to the end-user
   + The pipeline uses BASH to BLAST and handle files (which is significantly faster than R) and the parameter file is shared between R and BASH.    

    The file is of the format [param_id==value==comment] where param_id and value columns are CASE-SENSITIVE (because its unnecessarily hard to check and convert param types in BASH). 
    A default/example file is in fs::path_package("COMPLETE","inst","parameters.txt")

### USER DATA :
     Columns Org, genome, gtf

### COMPLETE.format.ids :
  + The Ordering of ID labels can be referred from COMPLETE_env$FORMAT_ID_INDEX
  + Sequences are labelled with the following long ID format of R-COMPLETE (specific to this pipeline and referred to as COMPLETE.format.ids) (seqID_delimiter & transcripID_delimiter set in parameters, "::" & "||" respectively in this context )
  
        >$transcript_id $transcripID_delimiter $transcript_region ($strand) $seqID_delimiter $seqID_delimiter $org_name $gene_name $seqID_delimiter $ortho_cluster
        >SOME_TRANSCRIPT||cds(+)::SOMEORG::RANDOMGENE::ORTHOLOG_CLUSTERS
        >ENSDART00000193157||cds(+)::danio_rerio::sulf1::18335at7898,51668at7742,360590at33208

### FLOW :

   1) EXTRACT_DATA() - Extracts the transcript regions for Protein Coding Transcripts (provided in parameters, pipeline requires cds,5utr,3utr)
     from BIOMART and/or User provided genomes & GTFs. This functions uses biomaRt/biomartr for extracting data from BIOMART
     and BASH function extract_genomic_regions() for user provided data. Extraction priority/flow : User Data > biomaRt > biomartr
          - ODB Files are merged and transformed with BASH function merge_OG2genes_OrthoDB()
          - Orthologous genes are found for genes which are not present in the organism with BASH function check_OrthoDB()
          - Flank lengths are calculated from GTF data for missing UTRs (with variance correction, check ?calculate_gtf_stats)
          - FASTA Nucleotide Sequences for given TRANSCRIPT_REGIONS are fetched from BIOMART/Genome

   2) FIND_ORTHOLOGS() -
