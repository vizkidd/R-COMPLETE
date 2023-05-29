## ToC
+ [About](#about)
+ [Installation](#install)
+ [Requirements](#requires)
+ [Run Examples](#examples)
+ [Documentation](#docs)
    + [Parameters](#params)
    + [User Data](#user_data)
    + [COMPLETE.format.ids](#ids)
    + [Flow](#flow)
        + [EXTRACT_DATA()](#fun1)      

# R - COMPLETE

<a name="about"/>

Pipeline for extracting localization elements/motifs using a comparitive approach. Library can be installed and tested on R(=4.1). For a list of genes the pipeline downloads full transcript sequences from all the organisms (or selected organisms) from ENSEMBL or NCBI, Formats the headers and stores them according to the OrthoDB clusters. Sequences without clusters are clustered based on sequence identity/sequence coverage of Reciprocal BLAST Hits of CDS regions. Next step would be to stitch the transcript regions into full length transcripts and align them. Pipeline is format agnostic(to my knowledge - TMK), packaged with a multi-threaded BLAST framework, can load and convert BLAST formats (internally - to GRanges) and sports functions for performing Reciprocal Bidirectional Hits. The multithreaded BLAST framework can handle many-many BLAST Hits across organisms. The package is interfaced with bash, R transforms the data & calls the shots and bash handles files & BLASTing. Check requirements and installation instuctions before proceeding.

*Ironically this repo is incomplete but the functionality in it works. Under Construction Indefinitely. Documentation can be found within the package, Play around with the functions for the rest*

## Installation (on R - Linux or RStudio Docker with WSL in Windows) :

<a name="install"/>

```diff
sudo apt-get update && sudo apt-get install curl bzip2 parallel liblmdb-dev ncbi-blast+ samtools bedtools libz-dev liblzma-dev libbz2-dev libclang-dev gffread curl lsof
BiocManager::install(c("Rhtslib", "devtools", "BiocManager", "Biostrings", "biomaRt", "S4Vectors", "IRanges", "rtracklayer", "GenomicRanges", "BiocGenerics"))
devtools::install_github("https://github.com/vizkidd/R-COMPLETE/")
```

<a name="requires"/>

## REQUIRES :
+ Linux with BASH ($SHELL must be set or /bin/bash must exist) (export SHELL="/bin/bash")
+ **[Config Files](#files)**
+ Lot of space in `genomes_path`, `fasta_path` and `annos_path` path locations (in [parameters file](#params))
+ GNU parallel (in $PATH - BASH functions)
+ [GffRead](https://github.com/gpertea/gffread)
+ [Samtools](http://www.htslib.org/download/) (in $PATH - BASH functions)
+ [Bedtools](https://github.com/arq5x/bedtools2/releases) (in $PATH - BASH functions)
+ [ncbi-blast+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Compile from .src.tar.gz with ./configure && make all_r && sudo make install) (or sudo alien -i ncbi-blast-X.XX.X+-2.src.rpm) (or download binaries) ([Docs](https://www.ncbi.nlm.nih.gov/books/NBK52640/) & [Compilation](http://www.ncbi.nlm.nih.gov/books/NBK279671/#_introduction_Source_tarball))
     + Check if you have the binaries for <b>*blastdb_path*</b> and <b>*makeblastdb*</b>
+ [Zlib](https://zlib.net/) - (Compile from sources) (or) (sudo apt install libz-dev or yum install zlib-devel)
+ LZMA SDK - (sudo apt-get install liblzma-dev or yum install xz-devel)
+ BZLIB - (sudo apt-get install libbz2-dev libclang-dev or yum install bzip2-devel.x86_64)

<a name="odb"/>

### OrthoDB : (Optional)
+ [OrthoDB (ODB) Flat Files (>= v10.1)](https://www.orthodb.org/?page=filelist) (Pipeline is tested with ODB v10.1) 
     + [odb10v1_species.tab.gz](https://v101.orthodb.org/download/odb10v1_species.tab.gz) - Ortho DB organism ids based on NCBI taxonomy ids (mostly species level) 
     + [odb10v1_genes.tab.gz](https://v101.orthodb.org/download/odb10v1_genes.tab.gz)  -Ortho DB genes with some info 
     + [odb10v1_OG2genes.tab.gz](https://v101.orthodb.org/download/odb10v1_OG2genes.tab.gz) - OGs to genes correspondence 
    (OR)
     + odb10v1_OGgenes_fixed.tab.gz - Merged & Transformed ODB file (Done within pipeline - Only once)
     + odb10v1_OGgenes_fixed_user.tab.gz - Merged & Transformed ODB file BASED on user gene list (Done within pipeline - For different gene sets)

**NOTE : Set ``orthodb_path_prefix`` in parameters file if you use OrthoDB files**

<a name="tools"/>

### Tools - (Paths for parameters file) :
+ [MACSE](https://bioweb.supagro.inra.fr/macse/) (Path to the .jar)
+ [MAFFT](https://mafft.cbrc.jp/alignment/software/installation_without_root.html) (Compile from sources <b>with extensions</b> because <b>*mafft-qinsi*</b> is required)
+ [TRANSAT](https://e-rna.org/transat/help.cgi#data) (Download preferred tarball and check INSTALL file)
+ [RNADECODER](https://github.com/jujubix/rnadecoder) (Compiled program is in the bin/ of the repo) <b>(Give path to the folder containing the binary)</b>
+ [Python3](https://www.python.org/download/releases/3.0/)
+ [FastTree](www.microbesonline.org/fasttree)

<a name="files"/>

### Files (Config) :
+ [Parameters](inst/pkg_data/parameters.txt)
+ [User Data](inst/pkg_data/user_data.txt) (Optional)
+ [Reference Organisms](inst/pkg_data/reference_ORGS.txt) (`COMPLETE_env$org.meta` has the list of organisms available)

<a name="examples"/>

## Run Example :
To run the example, from the context of your current working directory, 
+ Download [OrthoDB(ODB) files](#odb) (optional) and [Tools](#tools)
+ Check [config files](#files)
+ Provide paths and options in the parameters file
    + **NOTE : Default parameters file ([parameters.txt](inst/pkg_data/parameters.txt)) is at** ``fs::path_package("COMPLETE","pkg_data","parameters.txt")``

```{R}
params_list <- COMPLETE::load_params(fs::path_package("COMPLETE","pkg_data","parameters.txt"))
COMPLETE::EXTRACT_DATA(params_list = params_list, gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt"), user_data = fs::path_package("COMPLETE","pkg_data", "user_data.txt"), only.user.data = F )
```

<span style="color: #ff0000">**NOTE : First run will take some time due to conversion of ODB file structure (if OrthoDB is used)**</span>

## Documentation

<a name="docs"/>

```{R}
?COMPLETE_PIPELINE_DESIGN (in R docs)
```

### PARAMETERS :

<a name="params"/>

   The pipeline takes a single parameter file. This design was chosen,
   + To expose as many options as possible to the end-user
   + The pipeline uses BASH to BLAST and handle files (significantly faster than R) and the parameter file is shared between R and BASH.    
   + Parameters tagged as output in comment column are outputs from COMPLETE
   
```diff
* Delimited by '=='
* Inputs and Ouputs are specified in the comments
* The file is of the format [param_id==value==comment] where param_id and value columns are CASE-SENSITIVE (hard to check and convert param types in BASH). 
* A default/example file is in fs::path_package("COMPLETE","pkg_data","parameters.txt")
```
[EXTRACT_DATA()](#fun1) - `genomes_path`, `annos_path`, 

### USER DATA : (Oprtional)

<a name="user_data"/>

```diff
* Columns Org, genome, gtf
* Can accept empty or '-' in genome and/or gtf column. If empty or '-', the genome/gtf is looked up in ENSEMBL or NCBI DBs 
* A default/example file is in fs::path_package("COMPLETE","pkg_data","user_data.txt")
```

### COMPLETE.format.ids :

<a name="ids"/>

  + Order of FASTA ID labels are stored in ``COMPLETE_env$FORMAT_ID_INDEX``
  + Sequences are labelled with the following long ID format of R-COMPLETE (specific to this pipeline and referred to as COMPLETE.format.ids) (seqID_delimiter & transcripID_delimiter set in parameters, "::" & "||" respectively in this context )
  + COMPLETE.format.ids are indexed(internally) with `COMPLETE::index_BLAST_tables()` for compatibility

```diff
>$transcript_id $transcripID_delimiter $transcript_region ($strand) $seqID_delimiter $seqID_delimiter $org_name $gene_name $seqID_delimiter $ortho_cluster
>SOME_TRANSCRIPT||cds(+)::SOMEORG::RANDOMGENE::ORTHOLOG_CLUSTERS
>ENSDART00000193157||cds(+)::danio_rerio::sulf1::18335at7898,51668at7742,360590at33208
```

### FLOW :

<a name="flow"/>
<a name="fun1"/>

   1) <b>EXTRACT_DATA()</b> - Extracts the transcript regions for Protein Coding Transcripts (provided in parameters, pipeline requires cds,5utr,3utr)
     from BIOMART and/or User provided genomes & GTFs. This functions uses biomaRt/biomartr for extracting data from BIOMART
     and BASH function *extract_transcript_regions()* for user provided data.
     Extraction priority/flow : User Data > biomaRt > biomartr
        + ODB Files are merged and transformed with BASH function *merge_OG2genes_OrthoDB()*
        + Orthologous genes are found for genes which are not present in the organism with BASH function *check_OrthoDB()*
        + Flank lengths are calculated from GTF data for missing UTRs (with variance correction, check ?*calculate_gtf_stats*)
        + FASTA Nucleotide Sequences for given TRANSCRIPT_REGIONS are fetched from BIOMART/Genome

   2) FIND_ORTHOLOGS() -
