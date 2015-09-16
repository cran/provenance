#' Read a .csv file with continuous (detrital zircon) data
#'
#' Reads a data table containing continuous data (e.g. detrital zircon
#' ages)
#' @param datafile the path of a .csv file with the input data,
#' arranged in columns.
#' @param errorfile the (optional) path of a .csv file with the
#' standard errors of the input data, arranged by column in the same
#' order as \code{datafile}. Must be specified if the data are to be
#' compared with the Sircombe-Hazelton dissimilarity.
#' @param method an optional string specifying the dissimilarity
#' measure which should be used for comparing this with other
#' datasets. Should be one of either \code{"KS"} (for
#' Kolmogorov-Smirnov) or \code{"SH"} (for Sircombe and Hazelton). If
#' \code{method = "SH"}, then \code{errorfile} should be specified. If
#' \code{method = "SH"} and \code{errorfile} is unspecified, then the
#' program will default back to the Kolmogorov-Smirnov dissimilarity.
#' @param xlabel an optional string specifying the nature and units of
#' the data.  This string is used to label kernel density estimates.
#' @param colmap an optional string with the name of one of R's
#' built-in colour palettes (e.g., heat.colors, terrain.colors,
#' topo.colors, cm.colors), which are to be used for plotting the data.
#' @return an object of class \code{distributional}, i.e. a list with the
#' following items:
#' 
#' \code{x}: a named list of vectors containing the numerical data for each sample
#' 
#' \code{err}: an (optional) named list of vectors containing the standard errors of \code{x}
#'
#' \code{method}: either "KS" (for Kolmogorov-Smirnov) or "SH" (for Sircombe Hazelton)
#' 
#' \code{breaks}: a vector with the locations of the histogram bin edges
#' 
#' \code{xlabel}: a string containing the label to be given to the x-axis on all plots
#' @examples
#' agefile <- system.file("DZ.csv",package="provenance")
#' errfile <- system.file("DZerr.csv",package="provenance")
#' DZ <- read.distributional(agefile,errfile)
#' plot(KDE(DZ$x$N1))
#' @export
read.distributional <- function(datafile,errorfile=NULL,method="KS",xlabel="age [Ma]",colmap='rainbow') {
    out <- list()
    if (method=="SH" & is.null(errorfile)) method <- "KS"
    class(out) <- "distributional"
    out$method <- method
    out$x <- list()
    out$err <- list()
    out$colmap <- colmap
    dat <- utils::read.csv(datafile,header=TRUE)
    ns = length(dat)
    for (i in 1:ns){
        out$x[[names(dat)[i]]] = dat[!is.na(dat[,i]),i]
    }
    if (!is.null(errorfile)){
        err <- utils::read.csv(errorfile,header=TRUE)
        for (i in 1:ns) {
            out$err[[names(dat)[i]]] = dat[!is.na(err[,i]),i]
        }
    }
    d <- unlist(out$x)
    ng <- length(d) # number of grains
    nb <- log(ng/ns,base=2)+1
    out$breaks <- seq(min(d),max(d),length.out=nb+1)
    out$xlabel <- xlabel
    return(out)
}

#' Read a .csv file with categorical data
#'
#' Reads a data table containing categorical data (e.g. petrographic,
#' heavy mineral or geochemical data)
#'
#' @param fname a string with the path to the .csv file
#' @param method either "bray" (for the Bray-Curtis distance) or
#' "aitchison" (for Aitchison's central logratio distance). If
#' omitted, the function defaults to 'aitchison', unless there are
#' zeros present in the data.
#' @param colmap an optional string with the name of one of R's
#' built-in colour palettes (e.g., heat.colors, terrain.colors,
#' topo.colors, cm.colors), which are to be used for plotting the data.
#' @return an object of class \code{compositional}, i.e. a list with the
#' following items:
#' 
#' \code{x}: a data frame with the samples as rows and the categories as columns
#'
#' \code{method}: either "aitchison" (for Aitchison's centred logratio
#' distance) or "bray" (for the Bray-Curtis distance)
#' @examples
#' fname <- system.file("Major.csv",package="provenance")
#' Major <- read.compositional(fname)
#' plot(PCA(Major))
#' @export
read.compositional <- function(fname,method=NULL,colmap='rainbow') {
    out <- list()
    class(out) <- "compositional"
    out$x <- utils::read.csv(fname,header=TRUE,row.names=1)
    if (is.null(method)){
        if (any(out$x==0)) { method <- "bray" }
        else { method <- "aitchison" }
    }
    out$method <- method
    out$colmap <- colmap
    if (any(out$x==0) & method=="aitchison"){
        stop(paste("This dataset contains zeros and is",
                   "incompatible with the 'aitchison' distance"))
    }
    return(out)
}
#' Read a .csv file with mineral and rock densities
#'
#' Reads a data table containing densities to be used for
#' hydraulic sorting corrections (minsorting and srd functions)
#'
#' @param fname a string with the path to the .csv file
#' @return a vector with mineral and rock densities
#' @examples
#' data(Namib,densities)
#' N8 <- subset(Namib$HM,select="N8")
#' distribution <- minsorting(N8,densities,phi=2,sigmaphi=1,medium="air",by=0.05)
#' plot(distribution)
#' @export
read.densities <- function(fname){
    return(utils::read.csv(fname,header=TRUE))
}

#' Kolmogorov-Smirnov dissimilarity
#'
#' Returns the Kolmogorov-Smirnov dissimilarity between two samples
#'
#' @param x the first sample as a vector
#' @param y the second sample as a vector
#' @return a scalar value representing the maximum vertical distance
#' between the two cumulative distributions
#' @examples
#' data(Namib)
#' print(KS.diss(Namib$DZ$x[['N1']],Namib$DZ$x[['T8']]))
#' @export
KS.diss <- function(x,y) {
    xx = sort(x)
    cdftmp = stats::ecdf(xx)
    cdf1 = cdftmp(xx)
    xy = sort(y)
    cdftmp = stats::ecdf(xy)
    cdfEstim = cdftmp(xy)
    cdfRef = stats::approx(xx, cdf1, xy, yleft = 0, yright = 1, ties = "mean")
    dif = cdfRef$y - cdfEstim
    dif = abs(dif)
    out = max(dif)
    return(out)
}

#' Calculate the dissimilarity matrix between two \code{distributional} or
#' \code{compositional} datasets
#'
#' Calculate the dissimilarity matrix between two datasets of class
#' \code{distributional} or \code{compositional} using the Kolmogorov-Smirnov,
#' Sircombe-Hazelton, Aitchison or Bray Curtis distance
#' 
#' @param x an object of class \code{distributional} or \code{compositional}
#' @param method (optional) either "KS", "SH", "aitchison" or "bray"
#' @examples
#' data(Namib)
#' print(round(100*diss(Namib$DZ)))
#' @return an object of class \code{diss}
#' @rdname diss
#' @export
diss <- function(x,method){ UseMethod("diss",x) }
#' @rdname diss
#' @export
diss.distributional <- function(x,method=NULL) {
    if (!is.null(method)) x$method <- method
    n <- length(x$x)
    d <- mat.or.vec(n,n)
    rownames(d) <- names(x$x)
    colnames(d) <- names(x$x)
    if (x$method=="SH") c2 <- getc2(x)
    for (i in 1:n){
        for (j in 1:n){
            if (x$method=="SH"){
                d[i,j] <- SH.diss(x,i,j,c.con=c2)
            }
            if (x$method=="KS"){
                d[i,j] <- KS.diss(x$x[[i]],x$x[[j]])
            }
        }
    }
    out <- stats::as.dist(d)
    class(out) <- append("diss",class(out))
    return(out)
}
#' @rdname diss
#' @export
diss.compositional <- function(x,method=NULL){
    if (!is.null(method)) x$method <- method
    if (x$method=="aitchison"){
        out <- stats::dist(clr(x)$x)
    } else {
        snames <- names(x)
        ns <- length(snames)
        d <- mat.or.vec(ns,ns)
        rownames(d) <- snames
        colnames(d) <- snames
        for (i in 1:ns){
            for (j in 1:ns){
                d[i,j] <- bray.diss(x$x[i,],x$x[j,])
            }
        }
        out <- stats::as.dist(d)
    }
    class(out) <- append("diss",class(out))
    return(out)
}

#' Bray-Curtis dissimilarity
#'
#' Calculates the Bray-Curtis dissimilarity between two samples
#' @param x a vector containing the first compositional sample
#' @param y a vector of length(x) containing the second compositional sample
#' @return a scalar value
#' @examples
#' data(Namib)
#' print(bray.diss(Namib$HM$x["N1",],Namib$HM$x["N2",]))
#' @export
bray.diss <- function(x,y){
    return(as.numeric(sum(abs(x-y))/sum(x+y)))
}

#' Multidimensional Scaling
#'
#' Performs classical or nonmetric Multidimensional Scaling analysis
#' @param x an object of class \code{distributional}, \code{compositional} or \code{diss}
#' @param classical boolean flag indicating whether classical (TRUE)
#' or nonmetric (FALSE) MDS should be used
#' @param ... optional arguments to be passed onto \code{diss} (if
#' \code{x} is of class \code{compositional} or \code{distributional})
#' or onto \code{cmdscale} or \code{isoMDS} (if \code{x} is of class
#' \code{dist}).
#' @return an object of class \code{MDS}, i.e. a list containing the
#' following items:
#'
#' \code{points}: a two column vector of the fitted configuration
#'
#' \code{classical}: a boolean flag indicating whether the MDS
#' configuration was obtained by classical (\code{TRUE}) or nonmetric
#' (\code{FALSE}) MDS.
#'
#' \code{diss}: the dissimilarity matrix used for the MDS analysis
#' 
#' \code{stress}: (only if \code{classical=TRUE}) the final stress
#' achieved (in percent)
#' @examples
#' data(Namib)
#' plot(MDS(Namib$Major,classical=TRUE))
#' @rdname MDS
#' @export
MDS <- function(x,...){ UseMethod("MDS",x) }
#' @rdname MDS
#' @export
MDS.compositional <- function(x,classical=FALSE,...){
    d <- diss.compositional(x,...)
    return(MDS.diss(d,classical=classical))
}
#' @rdname MDS
#' @export
MDS.distributional <- function(x,classical=FALSE,...){
    d <- diss.distributional(x,...)
    return(MDS.diss(d,classical=classical))
}
#' @rdname MDS
#' @export
MDS.diss <- function(x,classical=FALSE,...){
    out <- list() 
    if (classical){
        out$points <- stats::cmdscale(x)
    } else {
        out <- MASS::isoMDS(d=x,...)
    }
    out$classical <- classical
    out$diss <- x
    class(out) <- "MDS"
    return(out)
}

#' Centred logratio transformation
#'
#' Calculates Aitchison's centered logratio transformation for a
#' dataset of class \code{compositional}
#' @param x an object of class \code{compositional}
#' @param ... optional arguments of the generic function
#' @return a matrix of clr coordinates
#' @examples
#' # The following code shows that applying provenance's PCA function
#' # to compositional data is equivalent to applying R's built-in
#' # princomp function to the clr transformed data.
#' data(Namib)
#' plot(PCA(Namib$Major))
#' dev.new()
#' clrdat <- clr(Namib$Major)$x
#' biplot(princomp(clrdat))
#' @export
clr <- function(x,...){ UseMethod("clr",x) }
clr.default <- function(x,...){stop()}
#' @rdname clr
#' @export
clr.compositional <- function(x,...){
    out <- x
    g <- apply(log(x$x),1,mean)
    nc <- ncol(x$x)
    gg <- matrix(rep(g,nc),ncol=nc,byrow=FALSE)
    out$x <- log(x$x) - gg
    return(out)
}

#' Principal Component Analysis
#'
#' Performs PCA of compositional data using a centred logratio distance
#' @param x an object of class \code{compositional}
#' @param ... optional arguments to R's \code{princomp function}
#' @return an object of classes \code{PCA}, which is synonymous to
#' the stats packages' \code{princomp} class.
#' @examples
#' data(Namib)
#' plot(MDS(Namib$Major,classical=TRUE))
#' dev.new()
#' plot(PCA(Namib$Major),asp=1)
#' print("This example demonstrates the equivalence of classical MDS and PCA")
#' @export
PCA <- function(x,...){
    clrdat <- clr(x)
    pc <- stats::princomp(clrdat$x,...)
    class(pc) <- append("PCA",class(pc))
    return(pc)
}

#' Get a subset of a distributional dataset
#'
#' Return a subset of provenance data according to some specified indices
#' @param x an object of class \code{distributional}
#' @param subset logical expression indicating elements or rows to keep:
#' missing values are taken as false.
#' @param select a vector of sample names
#' @param ... optional arguments for the generic subset function
#' @return an object of class \code{distributional}
#' @seealso read.distributional
#' @examples
#' data(Namib)
#' coast <- subset(Namib$HM,select=c("N1","N2","T8","T13","N12","N13"))
#' summaryplot(coast,ncol=2)
#' @export
subset.distributional <- function(x,subset=NULL,select=NULL,...){
    out <- x
    if (!is.null(subset)){
        i <- which(subset,arr.ind=TRUE)
    } else if (!is.null(select)){
        i <- which(names(x) %in% select)
    } else {
        return(out)
    }
    if (length(x$err)==length(x$x)) out$err <- x$err[i]
    out$x <- x$x[i]
    return(out)
}
#' Get a subset of a compositional dataset
#'
#' Return a subset of provenance data according to some specified indices
#' @param x an object of class \code{compositional}
#' @param subset logical expression indicating elements or rows to keep:
#' missing values are taken as false.
#' @param select a vector of sample names.
#' @param components a vector specifying a subcomposition
#' @param ... optional arguments for the generic subset function
#' @return an object of class \code{compositional}
#' @seealso read.compositional
#' @export
subset.compositional <- function(x,subset=NULL,select=NULL,components=NULL,...){
    out <- x
    if (!is.null(subset)){
        i <- which(subset,arr.ind=TRUE)
    } else if (!is.null(select)){
        i <- which(names(x) %in% select)
    } else {
        i <- 1:length(names(x))
    }
    if (!is.null(components)){
        j <- which(colnames(x$x) %in% components,arr.ind=TRUE)
    } else {
        j <- 1:length(colnames(x$x))
    }
    out$x <- x$x[i,j]
    if (methods::is(x,"SRDcorrected")){
        out$restoration <- x$restoration[i]
        for (sname in rownames(out$x)){
            out$restoration[[sname]] <- x$restoration[[sname]][,j]
        }
    }
    return(out)
}

# returns list of dissimilarities between common items
getdisslist <- function(slist){
    dnames <- names(slist)
    lablist <- lapply(slist,function(x) names(x))
    commonlabels <- Reduce(intersect,lablist)
    for (name in dnames){
        slist[[name]] <- subset(slist[[name]],select=commonlabels)
    }
    disslist <- slist
    for (name in dnames){
        disslist[[name]] <- diss(slist[[name]])
    }
    return(disslist)
}

#' Generalised Procrustes Analysis (GPA)
#'
#' Given a number of input datasets, this function performs an MDS
#' analysis on each of these and the feeds the resulting
#' configurations into a GPA algorithm, which uses a combination of
#' transformations (reflections, rotations, translations and scaling)
#' to find a 'consensus' configuration which best matches all the
#' component configurations in a least-squares sense.
#'
#' @param ... a sequence of datasets of classes \code{distributional}
#' and \code{compositional}
#' @return an object of class \code{GPA}, i.e. a list containing the
#' following items:
#' 
#' \code{points}: a two column vector with the coordinates of the group configuration
#'
#' \code{labels}: a list with the sample names
#' @author Ian L. Dryden
#' @references Dryden, Ian, and Maintainer Ian Dryden. "Shapes
#' package." Vienna, Austria: R Foundation for Statistical Computing
#' (2012).
#' @examples
#' data(Namib)
#' GPA <- procrustes(Namib$DZ,Namib$HM)
#' plot(GPA)
#' @export
procrustes <- function(...) {
    dnames <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
    slist <- list(...)
    names(slist) <- dnames
    disslist <- getdisslist(slist)
    n <- length(labels(disslist[[1]]))
    m <- length(disslist)
    X <- array(dim=c(n,2,m))
    for (i in 1:m){
        md <- MASS::isoMDS(disslist[[i]],k=2)
        X[,,i] <- md$points
    }
    result <- shapes::procGPA(X, scale=TRUE, reflect=TRUE)
    out <- list()
    out$points <- result$mshape
    out$labels <- labels(disslist[[1]])
    class(out) <- "GPA"
    return(out)
}

#' Individual Differences Scaling of provenance data
#'
#' Performs 3-way Multidimensional Scaling analysis using Carroll and
#' Chang (1970)'s INdividual Differences SCALing method as implemented
#' using De Leeuw and Mair (2011)'s stress majorization algorithm.
#' @param ... a sequence of datasets of class \code{distributional} or
#' \code{compositional}
#' @param type is either "ratio", "interval", or "ordinal"
#' @return an object of class \code{INDSCAL}, i.e. a list containing
#' the following items:
#' 
#' delta: Observed dissimilarities
#'
#' obsdiss: List of observed dissimilarities, normalized
#'
#' confdiss: List of configuration dissimilarities
#'
#' conf: List of matrices of final configurations
#' 
#' gspace: Joint configurations aka group stimulus space
#' 
#' cweights: Configuration weights
#'
#' stress: Stress-1 value
#'
#' spp: Stress per point
#' 
#' sps: Stress per subject (matrix)
#' 
#' ndim: Number of dimensions
#'
#' model: Type of smacof model
#' 
#' niter: Number of iterations
#'
#' nobj: Number of objects
#' @author
#' Jan de Leeuw and Patrick Mair
#' @references
#' de Leeuw, J., & Mair, P. (2009). Multidimensional scaling using
#' majorization: The R package smacof. Journal of Statistical
#' Software, 31(3), 1-30, < http://www.jstatsoft.org/v31/i03/>
#' @examples
#' data(Namib)
#' plot(indscal(Namib$DZ,Namib$HM))
#' @export
indscal <- function(...,type='ordinal'){
    dnames <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
    slist <- list(...)
    names(slist) <- dnames
    disslist <- getdisslist(slist)
    out <- smacof::smacofIndDiff(disslist, constraint = "indscal",type=type)
    class(out) <- "INDSCAL"
    return(out)
}

# set minimum and maximum values of a dataset
setmM <- function(x,from=NA,to=NA,log=FALSE){
    if (is.na(from)) { from <- min(x); setm <- TRUE }
    else { setm <- FALSE }
    if (is.na(to)) { to <- max(x); setM <- TRUE }
    else { setM <- FALSE }
    if (setm) {
        if (log) { from <- from/2 }
        else {
            if (2*from-to<0) {from <- 0}
            else {from <- from-(to-from)/10}
        }
    }
    if (setM) {
        if (log) { to <- 2*to }
        else { to <- to+(to-from)/10 }
    }
    return(list(m=from,M=to))
}

#' @export
names.distributional <- function(x){
    return(names(x$x))
}
#' @export
names.compositional <- function(x){
    return(rownames(x$x))
}
#' @export
names.KDEs <- function(x){
    return(names(x$kdes))
}
#' @export
names.ternary <- function(x){
    return(rownames(x$xyz))
}

#' Calculate the number of grains required to achieve a desired level of sampling resolution
#'
#' Returns the number of grains that need to be analysed to decrease
#' the likelihood of missing any fraction greater than a given size
#' below a given level.
#' @param f the size of the smallest resolvable fraction (0<f<1)
#' @param p the probability that all n grains in the sample have missed
#' at least one fraction of size f
#' @param n, the number of grains in the sample
#' @return the number of grains needed to reduce the chance of missing
#' at least one fraction f of the total population to less than p
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' # number of grains required to be 99% that no fraction greater than 5% was missed:
#' print(get.n(0.01))
#' # number of grains required to be 90% that no fraction greater than 10% was missed:
#' print(get.n(p=0.1,f=0.1))
#' @export
get.n <- function(p=0.05,f=0.05){
    n <- 1
    while(T){
        pp <- get.p(n,f)
        if (pp<p){ break }
        else {n <- n+1}
    }
    return(n)
}

#' Calculate the probability of missing a given population fraction
#'
#' For a given sample size, returns the likelihood of missing any
#' fraction greater than a given size
#' @param n the number of grains in the detrital sample
#' @param f the size of the smallest resolvable fraction (0<f<1)
#' @return the probability that all n grains in the sample have missed
#' at least one fraction of size f
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' print(get.p(60))
#' print(get.p(117))
#' @export
get.p <- function(n,f=0.05){
    p <- 0
    M <- 1/f
    for (i in 1:M){
        p <- p + (-1)^(i-1) * choose(M,i)*(1-i*f)^n
    }
    return(p)
}

#' Calculate the largest fraction that is likely to be missed
#'
#' For a given sample size, returns the largest fraction which has
#' been sampled with p x 100 % likelihood.
#' @param n the number of grains in the detrital sample
#' @param p the required level of confidence
#' @return the largest fraction that is sampled with at least 100 x p%
#' certainty
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' print(get.f(60))
#' print(get.f(117))
#' @export
get.f <- function(n,p=0.05){
    fmin <- 0
    fmax <- 1
    for (i in 1:100){
        f <- (fmax+fmin)/2
        if (get.p(n,f)<p) { fmax <- f }
        else { fmin <- f }
    }
    return((fmin+fmax)/2)
}

get.densities <- function(X,dtable){
    if (!methods::is(X,"compositional")) stop("input is not of class compositional")
    minerals <- colnames(X$x)
    i <- which(colnames(dtable) %in% colnames(X$x), arr.ind=TRUE)
    return(dtable[i])
}

ndim <- function(X){
    return(length(dim(X)))
}

sumcols <- function(X,x){
    if (length(x)>1 & ndim(X[,x])>0) # >1 class, >1 sample
        out <- apply(X[,x],1,sum)
    if (length(x)>1 & ndim(X[,x])==0) # >1 class, 1 sample
        out <- sum(X[,x])
    if (length(x)==1 & ndim(X[,x])>0) # 1 class, >1 sample
        out <- sum(X[,x])
    if (length(x)==1 & ndim(X[,x])==0) # 1 class, 1 sample
        out <- X[,x]
    names(out) <- rownames(X)
    return(out)
}

# X is a vector of strings
sumlabels <- function(X){
    out <- X[1]
    n <- length(X)
    if (n==1) return(out)
    for (i in 2:length(X)){
        out <- paste(out,X[i],sep='+')
    }
    return(out)
}

#' Group components of a composition
#'
#' Adds several components of a composition together into a single component
#' @param X a compositional dataset
#' @param ... a series of new labels assigned to strings or vectors of strings
#' denoting the components that need amalgamating
#' @return an object of the same class as X with fewer components
#' @examples
#' data(Namib)
#' HMcomponents <- c("zr","tm","rt","TiOx","sph","ap","ep",
#'                   "gt","st","amp","cpx","opx")
#' am <- amalgamate(Namib$PTHM,feldspars=c("KF","P"),
#'                  lithics=c("Lm","Lv","Ls"),heavies=HMcomponents)
#' plot(ternary(am))
#' @rdname amalgamate
#' @export
amalgamate <- function(X,...){ UseMethod("amalgamate",X) }
#' @rdname amalgamate
#' @export
amalgamate.default <- function(X,...){
    groups <- list(...)
    ng <- length(groups)
    labels <- names(groups)
    out <- NULL
    for (i in 1:ng){
        colsum <- sumcols(X,groups[[i]])
        out <- cbind(out,colsum)
    }
    colnames(out) <- labels
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.compositional <- function(X,...){
    out <- X
    out$x <- amalgamate(X$x,...)
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.SRDcorrected <- function(X,...){
    out <- X
    out$x <- amalgamate.default(X$x,...)
    for (sname in names(X$restoration)){
        out$restoration[[sname]] <-
            amalgamate.default(X$restoration[[sname]],...)
    }
    return(out)
}

ternaryclosure <- function(X,x,y,z){ 
    xlab <- sumlabels(x)
    ylab <- sumlabels(y)
    zlab <- sumlabels(z)
    out <- cbind(sumcols(X,x),sumcols(X,y),sumcols(X,z))
    den <- rowSums(out)
    out <- apply(out,2,'/',den)
    if (methods::is(out,"matrix")) {
        colnames(out) <- c(xlab,ylab,zlab)
    } else {
        names(out) <- c(xlab,ylab,zlab)
    }
    return(out)
}

#' Define a ternary composition
#'
#' Create an object of class \code{ternary}
#' @param X an object of class \code{compositional}
#' @param x string or a vector of strings indicating the variables making up
#' the first subcomposition of the ternary system. If omitted, the first
#' component of X is used instead.
#' @param y second (set of) variables
#' @param z third (set of) variables
#' @return an object of class \code{ternary}, i.e. a
#' list containing:
#'
#' xyz: a three column matrix (or vector) of ternary compositions.
#'
#' and (if X is of class \code{SRDcorrected})
#'
#' restoration: a list of intermediate ternary compositions inherited
#' from the SRD correction
#'
#' @seealso restore
#' @examples
#' data(Namib)
#' tern <- ternary(Namib$PT,c('Q'),c('KF','P'),c('Lm','Lv','Ls'))
#' plot(tern,type="QFL")
#' @export
ternary <- function(X,x=NULL,y=NULL,z=NULL){
    if (!methods::is(X,"compositional")) stop("Input does not have class compositional")
    if (is.null(x)) x <- colnames(X$x)[1]
    if (is.null(y)) y <- colnames(X$x)[2]
    if (is.null(z)) z <- colnames(X$x)[3]
    arg <- deparse(substitute(x))
    out <- list()
    out$xyz <- ternaryclosure(X$x,x,y,z)
    class(out) <- "ternary"
    if (methods::is(X,"SRDcorrected")){
        out$restoration <- list()
        snames <- names(X$restoration)
        for (sname in snames){
            out$restoration[[sname]] <-
                ternaryclosure(X$restoration[[sname]],x,y,z)
        }
        class(out) <- append("SRDcorrected",class(out))
    }
    return(out)
}
