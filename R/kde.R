# Abramson
# get geometric mean pilot density
getG <- function(pdens) {
    fpos <- pdens[pdens>0]
    N <- length(fpos)
    out <- exp(sum(log(fpos))/N)
    return(out)
}

# Abramson
# get fixed bandwidth pilot density
pilotdensity <- function(dat,bw){
    n <- length(dat)
    dens <- rep(0,n)
    for (i in 1:n){
        dens[i] <- mean(stats::density(dat,bw,from=(dat[i]-1e-10),
                                to=(dat[i]+1e-10),n=2)$y)
    }
    return(dens)
}

# adaptive KDE algorithm of Abramson (1982) as summarised by Jahn (2007)
Abramson <- function(dat,from,to,bw,...){
    n <- length(dat)
    lambda <- rep(0,n,length.out=n)
    pdens <- pilotdensity(dat,bw)
    G <- getG(pdens)
    dens <- rep(0,512)
    for (i in 1:n){
        lambda = sqrt(G/pdens[i])
        dens <- dens + stats::density(dat[i],bw*lambda,
                       from=from,to=to,...)$y
    }
    return(dens)
}

#' Create a kernel density estimate
#'
#' Turns a vector of numbers into an object of class \code{KDE} using
#' a combination of the Botev (2010) bandwidth selector and the
#' Abramson (1982) adaptive kernel bandwidth modifier. 
#' @param x a vector of numbers
#' @param from minimum age of the time axis. If NULL, this is set
#' automatically
#' @param to maximum age of the time axis. If NULL, this is set
#' automatically
#' @param bw the bandwidth of the KDE. If NULL, bw will be calculated
#' automatically using \code{botev()}
#' @param adaptive boolean flag controlling if the adaptive KDE
#' modifier of Abramson (1982) is used
#' @param log transform the ages to a log scale if TRUE
#' @param ... optional arguments to be passed on to \code{density}
#' @return an object of class \code{KDE}, i.e. a list
#' containing the following items:
#'
#' x: horizontal plot coordinates
#'
#' y: vertical plot coordinates
#'
#' bw: the base bandwidth of the density estimate
#'
#' ages: the data values from the input to the \code{KDE} function
#' @examples
#' data(Namib)
#' samp <- Namib$DZ$x[['N1']]
#' dens <- KDE(samp,0,3000,kernel="epanechnikov")
#' plot(dens)
#' @seealso KDEs
#' @export
KDE <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,...){
    out <- list()
    class(out) <- "KDE"
    out$log <- log
    if (is.na(from) | is.na(to)) {
        mM <- setmM(x,from,to,log)
        from <- mM$m
        to <- mM$M
    }
    if (log) {
        d <- log10(x)
        from <- log10(from)
        to <- log10(to)
        bw <- bw/(stats::median(x)*log(10))
    } else {
        d <- x
    }
    out$x <- seq(from=from,to=to,length.out=512)
    if (is.na(bw)){ bw <- botev(d) }
    if (adaptive){
        out$y <- Abramson(d,from,to,bw,...)
    } else {
        out$y <- stats::density(d,bw,from=from,to=to,...)$y
    }
    if (log) out$x <- 10^(out$x)
    out$y <- out$y/(sum(out$y)*(to-from)/512)
    out$bw <- bw
    out$ages <- x
    return(out)
}

#' Generate an object of class \code{KDEs}
#'
#' Convert a dataset of class \code{distributional} into an object of
#' class \code{KDEs} for further processing by the
#' \code{summaryplot} function.
#' @param x an object of class \code{distributional}
#' @param from minimum limit of the x-axis.
#' @param to maximum limit of the x-axis.
#' #' @param bw the bandwidth of the kernel density estimates. If bw =
#' NULL, the bandwidth will be set automatically using \code{botev()}
#' @param samebandwidth boolean flag indicating whether the same
#' bandwidth should be used for all samples. If samebandwidth = TRUE
#' and bw = NULL, then the function will use the median bandwidth of
#' all the samples.
#' @param adaptive boolean flag switching on the adaptive bandwidth
#' modifier of Abramson (1982)
#' @param pch (optional) symbol to be used to mark the sample points along the x-axis
#' @param normalise boolean flag indicating whether or not the KDEs
#' should all integrate to the same value.
#' @param log boolean flag indicating whether the data should by
#' plotted on a logarithmic scale.
#' @param ... optional parameters to be passed on to \code{density}
#' @return an object of class \code{KDEs}, i.e. a list containing the
#' following items:
#'
#' kdes: a named list with objects of class \code{KDE}
#'
#' from: the beginning of the common time scale
#'
#' to: the end of the common time scale
#' 
#' themax: the maximum probability density of all the KDEs
#'
#' pch: the plot symbol to be used by \code{plot.KDEs}
#'
#' xlabel: the x-axis label to be used by \code{plot.KDEs}
#' @examples
#' data(Namib)
#' KDEs <- KDEs(Namib$DZ,0,3000,pch=NA)
#' summaryplot(KDEs,ncol=3)
#' @seealso KDE
#' @export
KDEs <- function(x,from=NA,to=NA,bw=NA,samebandwidth=TRUE,
                 adaptive=TRUE,pch=NA,normalise=FALSE,log=FALSE,...){
    if (is.na(from) | is.na(to)) {
        mM <- setmM(unlist(x$x),from,to,log)
        from <- mM$m
        to <- mM$M
    }
    snames <- names(x$x)
    thekdes <- list()
    n <- list()
    themax <- -1
    if (is.na(bw) & samebandwidth) bw <- commonbandwidth(x)
    for (name in snames){
        thekdes[[name]] <- KDE(x$x[[name]],from,to,bw,adaptive,log,...)
        if (normalise){
            maxval <- max(thekdes[[name]]$y)
            if (themax < maxval) {themax <- maxval}
        }
    }
    out <- list()
    class(out) <- "KDEs"
    out$kdes <- thekdes
    out$from <- from
    out$to <- to
    out$themax <- themax
    out$pch <- pch
    out$xlabel <- x$xlabel
    return(out)
}