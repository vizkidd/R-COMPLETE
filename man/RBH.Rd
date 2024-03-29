% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMPLETE-find_orthologs.R
\name{RBH}
\alias{RBH}
\title{Find Reciprocal Blast Hits (RBH)}
\usage{
RBH(
  in1,
  in2,
  sep = "\\t",
  header = F,
  transcript_ID_metadata = NULL,
  col.names = NULL,
  index.tables = T,
  col.indices,
  unique.hit.weights = F,
  process.weights.func = max,
  n_threads = tryCatch(parallel::detectCores(all.tests = T, logical = T), error =
    function(cond) {
     return(2)
 }),
  use.feather = FALSE
)
}
\arguments{
\item{in1}{Input (query->subject) BLAST/WISARD hits table/filename}

\item{in2}{Input (query<-subject) BLAST/WISARD hits table/filename}

\item{sep}{Delimiter for the input Files. Only valid if in1 and in2 are files. Default - '\\t'}

\item{header}{Does the input files have header?. Only valid if in1 and in2 are files. Default - FALSE}

\item{transcript_ID_metadata}{Tab-delimited File with the filenames, indexed transcript IDs and the long transcript IDs.}

\item{col.names}{Columns names for the BLAST/WISARD tables/files}

\item{index.tables}{Should the IDs in the tables be indexed? (TRUE (Default) if COMPLETE.format.ids/Long BLAST Sequence IDs are used)}

\item{col.indices}{A Named List with indices of columns Query sequence ID (qseqid), Subject sequence ID (sseqid), and columns to be used as edge weights (eg, "Hsp_score","max_score" etc). Eg col.indices=list(qseqid=1,sseqid=2,weight.col=c(8,22)) OR col.indices=list(qseqid=1,sseqid=2)}

\item{unique.hit.weights}{Should only the unique Weights be taken for all Query->Subject Hits? (TRUE/FALSE (Default)). Only valid if weight.col is given in col.indices}

\item{process.weights.func}{Pass a function name to process the weights (eg, sum/max/min etc) (Default - max). Only valid if weight.col is given in col.indices}

\item{n_threads}{Number of Threads (Optional)}

\item{use.feather}{Should feather API be used to read BLAST Hits? - Default - FALSE}
}
\value{
Named List list(in1,in2) with the selected Hits which are RBHs between in1 and in2 data
}
\description{
Find RBH between BLAST results of different organisms/genes/transcripts (FASTA/FASTQ). The BLAST results must be of the format 6 and can be converted from BLAST format 11 with convert_BLAST_format().
The command with the required column names are given below.
}
\details{
convert_BLAST_format(in_file,outfile = out_file,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))

Optional : You can provide the file/table with indexed Transcsript IDs. The format must be "file"[tab]"long_id"[tab]"index" (without a header). It can be generated with the  index_FASTA_IDs() (check ?index_FASTA_IDs or index_fastaIDs() in fs::path_package("COMPLETE","exec","functions.sh"))
}
\note{
Both the input files are expected to have a header, the same number of columns with matching column order (and names). Give only indices for col.indices. Order of execution is unique.hit.weights followed by process.weights.func (if any/all these options are set), i.e Unique Weights are chosen for each hit (if unique.hit.weights=T) and then weights are processed using process.weights.func (if process.weights.func is set). If the function fails, try indexing the tables with index.tables=T. Weight Columns in col.indices are optional, when not given all the reciprocal hits are returned
}
\examples{
   convert_BLAST_format(in_file1,outfile = out_file1,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
   convert_BLAST_format(in_file2,outfile = out_file2,outformat=6,cols=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","score","gaps","frames","qcovhsp","sstrand","qlen","slen","qseq","sseq","nident","positive"))
   blast_GO1 <- GRObject_from_BLASTHits(blast_input = out_file1, COMPLETE.format.ids = T, col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
   blast_GO2 <- GRObject_from_BLASTHits(blast_input = out_file2, COMPLETE.format.ids = T, col.indices=list(qseqid=1,sseqid=2,evalue=11,qstart=7,qend=8,sstart=9,send=10,bitscore=12,qcovhsp=16,qlen=18,slen=19,frames=15,pident=3,gaps=14,length=4,sstrand=17))
   wis1 <- run_WISARD(blast_hits = blast_GO1,score_col = "Hsp_score",COMPLETE.format.ids = T) #score_col=16
   wis2 <- run_WISARD(blast_hits = blast_GO2,score_col = "Hsp_score",COMPLETE.format.ids = T) #score_col=16
   wis1 <- melt_wisard_list(wis1)
   wis2 <- melt_wisard_list(wis2)
   RBH(in1 = wis1, in2 = wis2, index.tables = T, col.indices = list(qseqid=12,sseqid=1,weight.col=c(8,22)),col.names = c("subject_id","start","end","width","strand","Hsp_num","Hsp_bit.score","Hsp_score","Hsp_evalue","Hsp_query.from","Hsp_query.to","query_id","query_len","subject_len","Hsp_hit.from","Hsp_hit.to","Hsp_query.frame","Hsp_hit.frame","Hsp_pidentity","Hsp_gaps","Hsp_align.len","max_score"),unique.hit.weights = T, process.weights.func = max)

}
