% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{read.varietal}
\alias{read.varietal}
\title{Read a .csv file with varietal data}
\usage{
read.varietal(
  fname,
  snames = NULL,
  sep = ",",
  dec = ".",
  method = "KS",
  check.names = FALSE,
  row.names = 1,
  ...
)
}
\arguments{
\item{fname}{file name (character string)}

\item{snames}{either a vector of sample names, an integer marking
the length of the sample name prefix, or
\code{NULL}. \code{read.varietal} assumes that the row names of
the \code{.csv} file consist of character strings marking the
sample names, followed by a number.}

\item{sep}{the field separator character.  Values on each line of
the file are separated by this character.}

\item{dec}{the character used in the file for decimal points.}

\item{method}{an optional string specifying the dissimilarity
measure which should be used for comparing this with other
datasets. Should be one of either \code{"KS"} (for
Kolmogorov-Smirnov) or \code{"Kuiper"} (for Kuiper)}

\item{check.names}{logical.  If \code{TRUE} then the names of the
variables in the frame are checked to ensure that they are
syntactically variable names.}

\item{row.names}{logical. See the documentation for the
\code{read.table} function.}

\item{...}{optional arguments to the built-in \code{read.csv}
function}
}
\value{
an object of class \code{varietal}, i.e. a list with the
    following items:

\code{x}: a compositional data table

\code{snames}: a vector of strings corresponding to the sample names

\code{name}: the name of the dataset, extracted from the file path
}
\description{
Reads a data table containing compositional data (e.g. chemical
concentrations) for multiple grains and multiple samples
}
\examples{
fn <- system.file("SNSM/Ttn_chem.csv",package="provenance")
Ttn <- read.varietal(fname=fn,snames=3)
plot(MDS(Ttn))
}
