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

