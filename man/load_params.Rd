% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMPLETE-options.R
\name{load_params}
\alias{load_params}
\title{Loads the parameters from param file}
\usage{
load_params(param_file)
}
\arguments{
\item{param_file}{parameter file name}
}
\value{
Named list of parameters
}
\description{
The function loads parameters from the param file into Global R variables and returns a named parameter table. The parameters are checked and an error is output if there are any missing parameters or if the folders are not accesible
}
\note{
Parameters are case sensitive and are de-limited with '=='.
}
\examples{
    params_list <- load_params(param_file)

}
