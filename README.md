## ToC
+ [About](#about)
+ [Installation](#install)
+ [Requirements](#requires)
+ [Run Examples](#examples)
+ [Documentation](#docs)
    + [Parameters](#params)
    + [User Data](#user_data)
    + [COMPLETE.format.ids](#ids)
    + [BLAST functions](#blast)
    + [Flow](#flow)
        + [EXTRACT_DATA()](#fun1)
        + [FIND_TRANSCRIPT_ORTHOLOGS()](#fun2)      

# R - COMPLETE

<a name="about"/>

Pipeline for extracting localization elements/motifs using a comparitive approach. Library can be installed and tested on R(>=4.1). For a list of genes, the pipeline downloads full transcript sequences for all organisms (or selected organisms) from ENSEMBL or NCBI, Formats the headers and stores them according to the clusters (optionally taken from [OrthoDB](https://www.orthodb.org/) or clustered otherwise). Sequences without cluster information are clustered based on sequence identity/sequence coverage of Reciprocal Bidirectional BLAST Hits (RBH) of CDS regions. Next step would be to stitch the transcript regions into full length transcripts and align them. Pipeline is format agnostic(to my knowledge - TMK), packaged with a multi-threaded BLAST framework, can load and convert BLAST formats (internally - to GRanges) and sports functions for performing Reciprocal Bidirectional Hits. The multithreaded BLAST framework can handle many-many BLAST Hits across organisms. The package is interfaced with bash, R transforms the data & calls the shots and bash handles files & BLASTing. Check requirements and installation instuctions before proceeding.

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
     + [odb10v1_OG2genes.tab.gz](https://v101.orthodb.org/download/odb10v1_OG2genes.tab.gz) - OGs to genes correspondence <br> **(OR)**
     + odb10v1_OGgenes_fixed.tab.gz - Merged & Transformed ODB file (Done within pipeline - Only once)
     + odb10v1_OGgenes_fixed_user.tab.gz - Merged & Transformed ODB file BASED on user gene list (Done within pipeline - For different gene sets)

**NOTE : Set ``orthodb_path_prefix`` in parameters file if you use OrthoDB files**

<a name="tools"/>

### Tools - (Paths in parameters file) :
+ [MACSE](https://bioweb.supagro.inra.fr/macse/) (Path to the .jar)
+ [MAFFT](https://mafft.cbrc.jp/alignment/software/installation_without_root.html) (Compile from sources <b>with extensions</b> because <b>*mafft-qinsi*</b> is required)
+ [TRANSAT](https://e-rna.org/transat/help.cgi#data) (Download preferred tarball and check INSTALL file)
+ [RNADECODER](https://github.com/jujubix/rnadecoder) (Compiled program is in the bin/ of the repo) <b>(Give path to the folder containing the binary)</b>
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
gene_list = fs::path_package("COMPLETE","pkg_data","genelist.txt")
user_data = fs::path_package("COMPLETE","pkg_data", "user_data.txt")
COMPLETE::EXTRACT_DATA(params_list = params_list, gene_list = gene_list, user_data = user_data, only.user.data = F )
COMPLETE::FIND_TRANSCRIPT_ORTHOLOGS(params_list = params_list, gene_list = gene_list, blast_program = Sys.which("tblastx"), group.mode=COMPLETE_env$FORMAT_ID_INDEX$CLUSTERS, run.mode="both", verbose=F, seed=123)
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

### USER DATA : (Optional)

<a name="user_data"/>

```diff
* Columns Org, genome, gtf
* Can accept empty or '-' in genome and/or gtf column. If empty or '-', the genome/gtf is looked up in ENSEMBL or NCBI DBs 
* A default/example file is in fs::path_package("COMPLETE","pkg_data","user_data.txt")
```

### COMPLETE.format.ids :

<a name="ids"/>

  + Order of FASTA ID labels are stored in ``COMPLETE_env$FORMAT_ID_INDEX``
  + Sequences are labelled with the following long ID format of R-COMPLETE (specific to this pipeline and referred to as COMPLETE.format.ids) <i>(seqID_delimiter & transcripID_delimiter set in parameters, `::` & `||` respectively in this context)</i>
  + COMPLETE.format.ids are indexed(internally) with `COMPLETE::index_BLAST_tables()` for compatibility

```diff
>$transcript_id $transcripID_delimiter $transcript_region ($strand) $seqID_delimiter $seqID_delimiter $org_name $gene_name $seqID_delimiter $ortho_cluster
>SOME_TRANSCRIPT||cds(+)::SOMEORG::RANDOMGENE::ORTHOLOG_CLUSTERS
>ENSDART00000193157||cds(+)::danio_rerio::sulf1::18335at7898,51668at7742,360590at33208
```
<a name="blast"/>

### BLAST Functions :
+ one2one
+ all2all
    
<a name="flow"/>

### FLOW :

<a name="fun1"/>

   1) <b>EXTRACT_DATA()</b> - Extracts the transcript regions for Protein Coding Transcripts (provided in parameters, pipeline requires cds,5utr,3utr)
     from BIOMART and/or User provided genomes & GTFs. This functions uses biomaRt/biomartr for extracting data from BIOMART
     and BASH function *extract_transcript_regions()* for user provided data.
     Extraction priority/flow : User Data > biomaRt > biomartr
        + ODB Files are merged and transformed with BASH function *merge_OG2genes_OrthoDB()*
        + Orthologous genes are found for genes which are not present in the organism with BASH function *check_OrthoDB()*
        + Flank lengths are calculated from GTF data for missing UTRs (with variance correction, check ?*calculate_gtf_stats*)
        + FASTA Nucleotide Sequences for given TRANSCRIPT_REGIONS are fetched from BIOMART/Genome

<a name="fun2"/>

   2) <b>FIND_TRANSCRIPT_ORTHOLOGS()</b> - Finds transcript-level orthologs based on minimum coverage and/or maximum sequence identity (check `?extract_transcript_orthologs`). This function reduces pool of organisms ~~based on the common genes between reference species~~ and reduces the . **Has different grouping modes(`group.mode`) and run modes(`run.mode`),** to group transcript orthologs at the level of <i>organisms, genes or Ortholog Clusters</i>, sequences are grouped into any level of `COMPLETE_env$FORMAT_ID_INDEX` (Default - `COMPLETE_env$FORMAT_ID_INDEX$CLUSTERS`) and select transcript orthologs. Gene level grouping has more tight orthology and fewer transcript orthologs. Ortholog Cluster level grouping is a level higher than Genes (Because an Ortholog Cluster can have more than one gene) and have highest number of transcript orthologs with a lot of dissimilarity. Different run modes determine how transcript-level orthologs are selected by their HSP coverages after grouping. After grouping, non-overlapping BLAST hits which maximize coverage for each transcript are chosen with [WISARD](https://github.com/robbueck/wisard/). Transcripts which do not have bi-directional hits are discarded with RBH(). Finally, HSP coverage is calculated with `COMPLETE::calculate_HSP_coverage()` and transcripts are processed according to `run.mode` option which can be one of, 
        + "coverage_distance" - Hits are filtered based on distance between bi-directional minimum HSP coverages (coverage_distance <= min_coverage_filter). This option selects more BLAST hits and should be used when the coverage values are very low (and the BLAST Hits/sequences are distant). `coverage_distance = 1 - (2 * aligned_length) / (query_length + subject_length)`. ``(coverage_distance >= min_coverage_filter)`
        + "coverage_filter" - Filters Hits based on minimum coverage of HSPs from either direction. Use this option when the coverage values are high (and the BLAST Hits/sequences are closely related). `(cov_q >= min_coverage_filter && cov_s >= min_coverage_filter)`
        + "both" - Uses both "coverage_distance" and "coverage_filter" and is very strict. (Default) . `(coverage_distance >= min_coverage_filter && cov_q >= min_coverage_filter && cov_s >= min_coverage_filter)`
        + "no_filter" - Only calculates HSP coverages and does not filter any Hits . 
Execution Flow (for each group): all2all BLAST -> [wisard](https://github.com/robbueck/wisard/) -> RBH -> calculate_HSP_coverage()
        + (ONLY for run.mode="cluster") - Place ungrouped sequences into groups 
        + Select transcript level orthologs with minimum coverage between clusters/genes
        + Creates compatible sets of organisms based on compatible sets of genes

**NOTE : ONLY USE THIS FUNCTION WHEN RUNNING THE PIPELINE OF R-COMPLETE. Use other helper function to work with custom BLAST files not generated by this R package.**
