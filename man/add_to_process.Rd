% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMPLETE-data_extraction.R
\name{add_to_process}
\alias{add_to_process}
\title{Internal Function - Add process to a list of process}
\usage{
add_to_process(
  p_cmd,
  p_args = list(),
  verbose = F,
  logfile = NULL,
  params_list
)
}
\arguments{
\item{p_cmd}{Command to be executed}

\item{p_args}{Arguments to be passed to the command}

\item{verbose}{Print status messages?}

\item{logfile}{Redirect output (stdout & stderr) to this file}

\item{params_list}{Output of load_params()}
}
\value{
Process ID from processx::new() which can be used for further monitoring of the process
}
\description{
This function adds process to a list and makes sure hte number of processes do not
exceed a user defined limit (< max_concurrent_jobs/cores/threads). Uses processx library to call
a child process. Process is appended to COMPLETE_env$process_list
}
