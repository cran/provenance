% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolbox.R
\name{varietal2distributional}
\alias{varietal2distributional}
\title{Convert varietal to distributional data}
\usage{
varietal2distributional(x, bycol = FALSE, plot = FALSE)
}
\arguments{
\item{x}{an object of class \code{varietal}.}

\item{bycol}{logical. If \code{TRUE}, returns a list of
distributional objects (one for each element). If \code{FALSE},
returns a single distributional object (containing the PC1
scores for each sample).}

\item{plot}{logical. If \code{TRUE}, shows the PCA biplot that is
used when \code{bycol} is \code{FALSE}.}
}
\description{
Convert an object of class \code{varietal} either to a list of
distributional objects by breaking it up into separate elements, or
to a single distributional object corresponding to the first
principal component.
}
\examples{
Ttn_file <- system.file("SNSM/Ttn_chem.csv",package="provenance")
Ttn <- read.varietal(fn=Ttn_file,snames=3)
varietal2distributional(Ttn,bycol=FALSE,plot=TRUE)
}
