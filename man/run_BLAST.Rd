% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMPLETE-BLAST_wrapper.R
\name{run_BLAST}
\alias{run_BLAST}
\title{Perform Nucleotide BLAST}
\usage{
run_BLAST(
  query_path,
  subject_path,
  blast_DB_dir = file.path(tempdir(), run_name),
  blast_out = NULL,
  blast_program = Sys.which("tblastx"),
  run_name = "BLAST",
  blast_options = "",
  COMPLETE.format.ids = F,
  params_list = NULL,
  blast.sequence.limit = 0,
  n_threads = 8,
  gzip.output = F,
  clean_extract = F,
  verbose = F,
  seed = 123,
  return_data = TRUE,
  use.feather = T,
  header = F
)
}
\arguments{
\item{query_path}{Path to Query FASTA or Biostrings::DNAString|AAString|RNAString}

\item{subject_path}{Path to Subject FASTA or or Biostrings::DNAString|AAString|RNAString}

\item{blast_DB_dir}{Path to BLAST DBs, if provided, The Query and Subject FASTA are copied into this directory and then BLASTed. Default is NULL}

\item{blast_out}{Path to BLAST output file, Default BLAST FORMAT is 11. It is converted internally to BLAST Format 6 and returned as a GRanges Object. Default is tempdir(). Can also be NULL}

\item{blast_program}{Give path to the BLAST program. eg, Default - Sys.which("tblastx") if tblastx is in SHELL $PATH.}

\item{run_name}{Name of the BLAST run. Only for logging (Optional)}

\item{blast_options}{Extra Options to be passed to the BLAST program}

\item{COMPLETE.format.ids}{Do BLAST Hit IDs of BLAST Hits (query and subject) have R-COMPLETE's long format IDs? (TRUE if using BLAST results from this package, FALSE (Default) otherwise) (Refer ?COMPLETE_PIPELINE_DESIGN) (Optional)}

\item{params_list}{Output from load_params() (Optional)}

\item{blast.sequence.limit}{Maximum number of sequences allowed in each BLAST file. Default - 0. If the query FASTA sequences > blast.sequence.limit, the sequences are split into multiple files and BLASTed. Use 0 to not split & copy the files into temporary files}

\item{n_threads}{Number of threads. Default - 8}

\item{gzip.output}{Should the output files be Gzipped? Default - FALSE}

\item{clean_extract}{Delete Output file if exists? (Optional). Default - F}

\item{verbose}{Print DEBUG Messages?}

\item{seed}{Seed value for reproducibility. Default - 123}

\item{return_data}{Logical. Should GRanges Object of BLAST Hits be returned? Default - TRUE}

\item{use.feather}{Should feather API be used to read BLAST Hits? - Default - TRUE}

\item{header}{Does the BLAST format output have header? (Some do, check it, by default fmt 6 is used which does not have a header) - FALSE}
}
\value{
BLAST Hits as GRanges Object
}
\description{
BLAST between two FASTA files. BLAST output is saved temporarily in BLAST format 11, internally converted to format 6 with convert_BLAST_format() and BLAST Hits are returned as a GRanges Object. If blast_out is provided then the intermediate BLAST Format 11 output can be found in the same directory with the same name and *.blast11 extension. If blast_DB_dir is provided then the Subject FASTA is copied into this directory and then BLASTed.
}
\note{
ASSUMES Nucleotide FASTA sequences. Files are copied to tempdir() and split by blast.sequence.limit before blasting. Not checking/Not working for Protein/Peptide sequences and FASTQ files. Columns of GRanges BLAST 6 output are c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive")
}
\examples{
    run_BLAST(query_path = "query.fasta",subject_path = "subject.fasta",blast_DB_dir = "files/blastdb", blast_program="tblastx", run_name = "blast_positive",blast_options = "-strand plus", gzip.output=F, use.feather=T)
    run_BLAST(query_path = Biostrings::DNAString(gsub("[\r\n]", "", "ATTGTCGAAGTTGTCGCTCGAGAGGCGGGAGTTTACCGACACTTTTCCTCAGAAGTTTAC
    CGTGAAGCTGACCGGAGAACGGCGGAGTCTGTGCTGAATCTGCCATCATGTCCAGGCGGA
    ")), subject_path=Biostrings::DNAString(gsub("[\r\n]", "", "TTACCGACACTTTTCCTCAGAAGTTTACCGTGAAGCTGACCGGAGAACGGCGGAGTCTGT
    GCTGAATCTGCCATCATGTCCAGGCGGAGCTCACTCATAGTGCCGATGAAGAGTATTGAG
    ")))

}
