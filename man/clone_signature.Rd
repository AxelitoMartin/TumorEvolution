% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clone_signature.R
\name{clone_signature}
\alias{clone_signature}
\title{clone_signature}
\usage{
clone_signature(
  clone_mut,
  maf,
  sample.names,
  num_sign = 5,
  projectname = "",
  path = ".",
  num_parallelCores = 1
)
}
\arguments{
\item{clone_mut}{A list containing the variants observed in each variant of interest}

\item{maf}{The corresponding maf file to the samples used in canopy.}

\item{sample.names}{Sample names to consider in the maf file.}

\item{num_sign}{Number of signatures to look for.}

\item{projectname}{Name of the project to save results.}

\item{path}{Relative path to where results should be saved. Default is current directory.}

\item{num_parallelCores}{number of cores to use}
}
\value{
Results pdf
}
\description{
Get the mutational signature of each clone in a canopy run.
}
