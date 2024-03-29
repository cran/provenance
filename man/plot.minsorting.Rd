% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.minsorting}
\alias{plot.minsorting}
\title{Plot inferred grain size distributions}
\usage{
\method{plot}{minsorting}(x, cumulative = FALSE, components = NULL, ...)
}
\arguments{
\item{x}{an object of class \code{minsorting}}

\item{cumulative}{boolean flag indicating whether the grain size
distribution should be plotted as a density or cumulative
probability curve.}

\item{components}{string or list of strings with the names of a
subcomposition that needs plotting}

\item{...}{optional parameters to be passed on to graphics::matplot
(see \code{?par} for details)}
}
\description{
Plot the grain size distributions of the different minerals under
consideration
}
\examples{
data(endmembers,densities)
OPH <- subset(endmembers,select="ophiolite")
distribution <- minsorting(OPH,densities,phi=2,sigmaphi=1,
                           medium="air",by=0.05)
plot(distribution,components=c('F','px','opaques'))
}
\seealso{
minsorting
}
