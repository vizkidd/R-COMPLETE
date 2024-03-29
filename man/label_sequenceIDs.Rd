% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMPLETE-data_extraction.R
\name{label_sequenceIDs}
\alias{label_sequenceIDs}
\title{Label Sequence IDs of A FASTA File (Refer ?COMPLETE_PIPELINE_DESIGN about COMPLETE.format.IDs)}
\usage{
label_sequenceIDs(
  fasta_path,
  org,
  gene,
  odb_clusters,
  duplicates.method,
  params_list
)
}
\arguments{
\item{fasta_path}{Path to FASTA File}

\item{org}{Name of the organism}

\item{gene}{Filename of gene list or Vector of gene names}

\item{odb_clusters}{ODB Clusters}

\item{duplicates.method}{merge/delete/make_unique. How to handle sequences with duplicate names?. For CDS/UTR blocks - merge (Concatenate sequences with duplicate names), for Exon blocks - make_unique (Sequence names/IDs are made unique), delete - for deleting all duplicate sequences (first seq is kept).}

\item{params_list}{Output of load_params()}
}
\description{
Convert short sequence IDs into longer COMPLETE format IDs
}
