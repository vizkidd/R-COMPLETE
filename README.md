# mrna_zipcodes
Pipeline for extracting localization elements/motifs using a comparitive approach. Library can be installed and tested on R(=4.1) (ironically this repo is incomplete but the functionality in it works)..Under Construction Indefinitely. *Documentation can be found within the package, Play around with the functions for the rest*

## Installation (on R - Linux or Docker with WSL in Windows):
* sudo apt-get update && sudo apt-get install curl bzip2
* Zlib - (Compile from sources - https://zlib.net/) (or) (sudo apt install libz-dev or yum install zlib-devel)
* LZMA SDK - (sudo apt-get install liblzma-dev or yum install xz-devel)
* BZLIB - (sudo apt-get install libbz2-dev libclang-dev or yum install bzip2-devel.x86_64)
* R-Packages : Rhtslib, devtools, BiocManager, Biostrings, biomaRt, S4Vectors, IRanges, Rgb, rtracklayer, GenomicRanges, BiocGenerics
~~* Rgb (Archived R Package) - https://cran.r-project.org/src/contrib/Archive/Rgb/
~  * Download the archive and extract it in the R library
* devtools::install_github("https://github.com/vizkidd/R-COMPLETE/")

## Documentation
* ?COMPLETE_PIPELINE_DESIGN (in R docs)
