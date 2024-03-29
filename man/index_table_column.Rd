% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMPLETE-find_orthologs.R
\name{index_table_column}
\alias{index_table_column}
\title{Index Table Column}
\usage{
index_table_column(in_table, index_col, offset = 0)
}
\arguments{
\item{in_table}{Input table}

\item{index_col}{Column index of BLAST table to index (shorten IDs)}

\item{offset}{Offset value to add to the indices}
}
\value{
BLAST table with indexed IDs from index_col (original index_col is attached to the table)
}
\description{
Convert longer FASTA IDs into short indices. BLAST table IDs need to be short for downstream processing (eg, graph package does not accept long names). The IDs are stored as an extra column and returned along with the BLAST table.
}
\note{
WARNING : Do not index it more than once (although you may do index it as much as you like)
}
\examples{
 #Indexing columns 1 & 2
in1_data <- index_BLAST_table(in1_data,2,offset = 0)
in1_data <- index_BLAST_table(in1_data,1,offset = length(levels(factor( in1_data[,1] ))))
#TO INDEX TWO TABLES
in1_data <- index_BLAST_table(in1_data,2,offset = 0)
in1_data <- index_BLAST_table(in1_data,1,offset = length(levels(factor( in2_data[,1] ))) + length(levels(factor( in1_data[,2] ))))
in2_data <- index_BLAST_table(in2_data,2,offset = length(levels(factor( in1_data[,1] ))) + length(levels(factor( in2_data[,2] ))))
in2_data <- index_BLAST_table(in2_data,1,offset = length(levels(factor( in2_data[,1] ))) + length(levels(factor( in1_data[,2] ))))

}
