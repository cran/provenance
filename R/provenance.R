# save plot d as a pdf file with name f
saveplot <- function(f, d){
    dev.set(d)
    dev.copy2pdf(file=f)
}

# Botev:
# computes the discrete cosine transform of the column vector data
dct1d <- function(dat){
    n <- length(dat);
    # Compute weights to multiply DFT coefficients
    weight <- c(1,2*exp(-1i*(1:(n-1))*pi/(2*n)));
    # Re-order the elements of the columns of x
    dat <- c(dat[seq(1,n-1,2)], dat[seq(n,2,-2)]);
    # Multiply FFT by weights:
    out <- Re(weight* fft(dat));
    return(out)
}

# Botev:
# this implements the function t-zeta*gamma^[l](t)
fixedpoint <-  function(tt,N,II,a2){
    l <- 7
    f <- 2*(pi^(2*l))*sum((II^l)*a2*exp(-II*(pi^2)*tt))
    for (s in (l-1):2){
        K0 <- prod(seq(1,2*s-1,2))/sqrt(2*pi)
        const <- (1+(1/2)^(s+1/2))/3
        TT <- (2*const*K0/N/f)^(2/(3+2*s))
        f <- 2*pi^(2*s)*sum(II^s*a2*exp(-II*pi^2*TT))
    }
    out <- tt-(2*N*sqrt(pi)*f)^(-2/5)
    return(out)
}

#' Compute the optimal kernel bandwidth
#'
#' Uses the diffusion algorithm of Zdravko Botev (2011)
#' to calculate the bandwidth for kernel density estimation
#'
#' @param x a vector of ordinal data
#' @return a scalar value with the optimal bandwidth
#' @author Dzdravko Botev
#' @references Botev, Z. I., J. F. Grotowski, and
#' D. P. Kroese. "Kernel density estimation via diffusion." The Annals
#' of Statistics 38.5 (2010): 2916-2957.
#' @examples
#' fname <- system.file("DZ.csv",package="provenance")
#' bw <- botev(read.DZdata(fname)$x$N1)
#' print(bw)
#' @export
botev <- function(x){
    n <- 512
    minimum <- min(x)
    maximum <- max(x)
    Range <- maximum - minimum
    MIN <- minimum-Range/10
    MAX <- maximum+Range/10
    R <- MAX-MIN
    dx <- R/n
    xmesh <- MIN+seq(0,R,dx)
    if (anyDuplicated(x)){
        N <- length(as.numeric(names(table(x))))
    } else {
        N <- length(x)
    }
    w <- hist(x,xmesh,plot=FALSE)
    initialdata <- (w$counts)/N
    initialdata <- initialdata/sum(initialdata)
    a <- dct1d(initialdata)
    II <- (1:(n-1))^2
    a2 <- (a[2:n]/2)^2
    tstar <- tryCatch(uniroot(fixedpoint,c(0,.1),
                      N=N,II=II,a2=a2,tol=10^(-14))$root,
                      error=function(e).28*N^(-2/5))
    bandwidth <- sqrt(tstar)*R
    return(bandwidth)
}

# returns the median bandwidth for plotting several KDEs together
getcommonbandwidth <- function(x){
    dnames <- names(x$x)
    n <- length(dnames)
    bw <- rep(0,n)
    for (i in 1:n){
        bw[i] <- botev(x$x[[i]])
    }
    return(median(bw))
}

# Sircombe and Hazelton:
# Density estimator with sample-point adaptive bandwidth
adapt.f <- function(x,h,eval){
    f <- eval*0
    for (i in 1:length(eval)){f[i] <- mean(dnorm(x-eval[i],0,sd=h))}
    return(f)
}

# Sircombe and Hazelton:
sig2.con <- function(x,sigx,UCV=TRUE){
    sig2max <- max(sigx^2)
    h <- 1.06*min(sd(x),IQR(x)/1.34)/length(x)^0.2
    if (UCV) h <- bw.ucv(x,lower=h/20,upper=h)
    sig2.con <- sig2max+h^2
    return(sig2.con)
}

# Sircombe and Hazelton:
# Calculate integrated squared KDEs
Rf <- function(x,sig2,c.con){
    h1 <- sqrt(c(outer(c.con-sig2,c.con-sig2,"+")))
    xdiff <- c(outer(x,x,"-"))
    rf <- mean(dnorm(xdiff,sd=h1))
    return(rf)
}

#' Sircombe and Hazelton distance
#'
#' Calculates Sircombe and Hazelton's L2 distance between the Kernel
#' Functional Estimates (KFEs) of two samples with specified
#' analytical uncertainties
#'
#' @param x an object of class \code{DZdata}
#' @param i index of the first sample
#' @param j index of the second sample
#' @param c.con smoothing bandwidth of the kernel functional estimate
#' @author Keith Sircombe and Martin Hazelton
#' @references Sircombe, K. N., and M. L. Hazelton. "Comparison of
#' detrital zircon age distributions by kernel functional estimation."
#' Sedimentary Geology 171.1 (2004): 91-111.
#' @examples
#' datfile <- system.file("DZ.csv",package="provenance")
#' errfile <- system.file("DZerr.csv",package="provenance")
#' DZ <- read.DZdata(datfile,errfile)
#' d <- SH.dist(DZ,1,2)
#' print(d)
#' @export
SH.dist <- function(x,i,j,c.con=0){
    X <- x$x[[i]]
    Y <- x$x[[j]]
    sigX <- x$err[[i]]
    sigY <- x$err[[j]]
    sig2X <- sigX^2
    sig2Y <- sigY^2
    if (c.con <= 0) c.con <- max(sig2.con(X,sigX),sig2.con(Y,sigY))
    rfX <- Rf(X,sig2X,c.con)
    rfY <- Rf(Y,sig2Y,c.con)
    h1 <- sqrt(c(outer(c.con-sig2X,c.con-sig2Y,"+")))
    XYdiff <- c(outer(X,Y,"-"))
    Itmp <- mean(dnorm(XYdiff,sd=h1))
    dXY <- rfX+rfY-2*Itmp
    dXY <- sqrt(dXY)
    return(dXY)
}

# Sircombe and Hazelton:
# calculate 'c' value for data (this returns c^2)
getc2 <- function(x){
    n <- length(x$x)
    c2 <- 0
    for (i in 1:n){
        foo <- sig2.con(x$x[[i]],x$err[[i]])
        c2 <- max(c2,foo)
    }
    return(c2)
}

#' Read a .csv file with continuous (detrital zircon) data
#'
#' Reads a data table containing continuous data (e.g. detrital zircon
#' ages)
#' @param datafile the path of a .csv file with the input data,
#' arranged in columns.
#' @param errorfile the (optional) path of a .csv file with the
#' standard errors of the input data, arranged by column in the same
#' order as \code{datafile}. Must be specified if the data are to be
#' compared with the Sircombe-Hazelton distance.
#' @param metric an optional string specifying the dissimilarity
#' measure which should be used for comparing this with other
#' datasets. Should be one of either \code{"KS"} (for
#' Kolmogorov-Smirnov) or \code{"SH"} (for Sircombe and Hazelton). If
#' \code{metric = "SH"}, then \code{errorfile} should be specified. If
#' \code{metric = "SH"} and \code{errorfile} is unspecified, then the
#' program will default back to the Kolmogorov-Smirnov dissimilarity.
#' @param xlabel an optional string specifying the nature and units of
#' the data.  This string is used to label kernel density estimates.
#' @return an object of class \code{DZdata}, i.e. a list with the
#' following items:
#' 
#' \code{x}: a named list of vectors containing the numerical data for each sample
#' 
#' \code{err}: an (optional) named list of vectors containing the standard errors of \code{x}
#'
#' \code{metric}: either "KS" (for Kolmogorov-Smirnov) or "SH" (for Sircombe Hazelton)
#' 
#' \code{breaks}: a vector with the locations of the histogram bin edges
#' 
#' \code{xlabel}: a string containing the label to be given to the x-axis on all plots
#' @examples
#' fname <- system.file("DZ.csv",package="provenance")
#' DZ <- read.DZdata(fname)
#' plot(getKDE(DZ$x$N1))
#' @export
read.DZdata <- function(datafile,errorfile=NULL,metric="KS",xlabel="age [Ma]") {
    out <- list()
    if (metric=="SH" & is.null(errorfile)){
        stop(paste0("The input file must include the analytical errors",
                    "of the input data to use the Sircombe-Hazelton dissimilarity"))
    }
    class(out) <- "DZdata"
    out$metric <- metric
    out$x <- list()
    out$err <- list()
    dat <- read.csv(datafile,header=TRUE)
    ns = length(dat)
    for (i in 1:ns){
        out$x[[names(dat)[i]]] = dat[!is.na(dat[,i]),i]
    }
    if (!is.null(errorfile)){
        err <- read.csv(errorfile,header=TRUE)
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
#' @param metric either "bray" (for the Bray-Curtis distance) or
#' "aitchison" (for Aitchison's central logratio distance). If
#' omitted, the function defaults to 'aitchison', unless there are
#' zeros present in the data.
#' @return an object of class \code{HMdata}, i.e. a list with the
#' following items:
#' 
#' \code{x}: a data frame with the samples as rows and the categories as columns
#'
#' \code{metric}: either "aitchison" (for Aitchison's centred logratio
#' distance) or "bray" (for the Bray-Curtis distance)
#' @examples
#' fname <- system.file("Major.csv",package="provenance")
#' Major <- read.HMdata(fname)
#' plot(getPCA(Major))
#' @export
read.HMdata <- function(fname,metric=NULL) {
    out <- list()
    class(out) <- "HMdata"
    out$x <- read.csv(fname,header=TRUE,row.names=1)
    if (is.null(metric)){
        if (any(out$x==0)) { metric <- "bray" }
        else { metric <- "aitchison" }
    }
    out$metric <- metric
    if (any(out$x==0) & metric=="aitchison"){
        stop(paste("This dataset contains zeros and is",
                   "incompatible with the 'aitchison' distance"))
    }
    return(out)
}

#' Kolmogorov-Smirnov distance
#'
#' Returns the Kolmogorov-Smirnov distance between two samples
#'
#' @param x the first sample as a vector
#' @param y the second sample as a vector
#' @return a scalar value representing the maximum vertical distance
#' between the two cumulative distributions
#' @examples
#' fname <- system.file("DZ.csv",package="provenance")
#' DZ <- read.DZdata(fname)
#' print(KS.dist(DZ$x[['N1']],DZ$x[['T8']]))
#' @export
KS.dist <- function(x,y) {
    xx = sort(x)
    cdftmp = ecdf(xx)
    cdf1 = cdftmp(xx)
    xy = sort(y)
    cdftmp = ecdf(xy)
    cdfEstim = cdftmp(xy)
    cdfRef = approx(xx, cdf1, xy, yleft = 0, yright = 1, ties = "mean")
    dif = cdfRef$y - cdfEstim
    dif = abs(dif)
    out = max(dif)
    return(out)
}

#' Calculate the dissimilarity matrix between two \code{DZdata} or
#' \code{HMdata} datasets
#'
#' Calculate the dissimilarity matrix between two datasets of class
#' \code{DZdata} or \code{HMdata} using the Kolmogorov-Smirnov,
#' Sircombe-Hazelton, Aitchison or Bray Curtis distance
#' 
#' @param x an object of class \code{DZdata} or \code{HMdata}
#' @param metric (optional) either "KS", "SH", "aitchison" or "bray"
#' @param ... optional arguments to the generic dist function
#' @examples
#' fname <- system.file("DZ.csv",package="provenance")
#' DZ <- read.DZdata(fname)
#' print(round(100*dist(DZ)))
#' @return a matrix of pairwise distances
#' @rdname dist
#' @export
dist <- function(x,...){ UseMethod("dist",x) }
#' @rdname dist
#' @export
dist.default <- function(x,...){stats::dist(x,...)}
#' @rdname dist DZdata
#' @export
dist.DZdata <- function(x,metric=NULL,...) {
    if (!is.null(metric)) x$metric <- metric
    n = length(x$x)
    diss = mat.or.vec(n,n)
    rownames(diss) = names(x$x)
    colnames(diss) = names(x$x)
    if (x$metric=="SH") c2 <- getc2(x)
    for (i in 1:n){
        for (j in 1:n){
            if (x$metric=="SH"){
                diss[i,j] = SH.dist(x,i,j,c.con=c2)
            }
            if (x$metric=="KS"){
                diss[i,j] = KS.dist(x$x[[i]],x$x[[j]])
            }
        }
    }
    return (as.dist(diss))
}
#' @rdname dist
#' @export
dist.HMdata <- function(x,metric=NULL,...){
    if (!is.null(metric)) x$metric <- metric
    if (x$metric=="aitchison"){
        return(dist(clr(x)$x))
    } else {
        return(vegan::vegdist(x$x,x$metric))
    }
}

# a function to plot the nearest neighbour lines
plotlines <- function(conf,diss) {
    # rank the samples according to their pairwise proximity
    i = t(apply(as.matrix(diss),1,function(x) order(x))[2:3,])
    # coordinates for the lines
    x1 = as.vector(conf[i[,1],1]) # calculate (x,y)-coordinates ...
    y1 = as.vector(conf[i[,1],2]) # ... of nearest neighbours
    x2 = as.vector(conf[i[,2],1]) # calculate (x,y)-coordinates ...
    y2 = as.vector(conf[i[,2],2]) # ... of second nearest neighbours
    for (j in 1:nrow(conf)) {
        lines(c(conf[j,1],x1[j]),c(conf[j,2],y1[j]),lty=1) # solid line
        lines(c(conf[j,1],x2[j]),c(conf[j,2],y2[j]),lty=2) # dashed line
    }
}

#' Plot an MDS configuration
#'
#' Plots the coordinates of a multidimensional scaling analysis as an
#' X-Y scatter plot or 'map' and, if x$classical = FALSE, a Shepard
#' plot.
#' 
#' @param x an object of class \code{MDS}
#' @param nnlines if TRUE, draws nearest neighbour lines
#' @param pch plot character (see ?plot for details)
#' @param cex magnification of the plot character (see ?par for details)
#' @param xlab a string with the label of the x axis
#' @param ylab a string with the label of the y axis
#' @param xaxt if = 'y', adds ticks to the x axis
#' @param yaxt if = 'y', adds ticks to the y axis
#' @param ... optional arguments to the generic \code{plot} function
#' @method plot MDS
#' @export
plot.MDS <- function(x,nnlines=FALSE,pch=NA,cex=NA,xlab="",ylab="",xaxt='n',yaxt='n',...){
    plot(x$points, type="n", asp=1, xlab=xlab, ylab=ylab,xaxt=xaxt,yaxt=yaxt,...)
    # draw lines between closest neighbours
    if (nnlines) {
        if (is.na(pch)) pch=21
        if (is.na(cex)) cex=2.5
        plotlines(x$points,x$diss)
    }
    points(x$points, pch=pch, cex=cex, col='red', bg='white')
    text(x$points, labels = labels(x$diss))
    if (!x$classical){
        dev.new()
        shep <- MASS::Shepard(x$diss, x$points)
        plot(shep, pch=".")
        lines(shep$x, shep$yf, type="S")
        title(paste0("Stress = ",x$stress))
    }
}

#' Multidimensional Scaling
#'
#' Performs classical or nonmetric Multidimensional Scaling analysis
#' @param x an object of class \code{DZdata}, \code{HMdata} or \code{dist}
#' @param classical boolean flag indicating whether classical (TRUE)
#' or nonmetric (FALSE) MDS should be used
#' @param ... optional arguments to be passed onto \code{cmdscale} or
#' \code{isoMDS}
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
#' fname <- system.file("Major.csv",package="provenance")
#' Major <- read.HMdata(fname)
#' plot(getMDS(Major,classical=TRUE))
#' @rdname getMDS
#' @export
getMDS <- function(x,...){ UseMethod("getMDS",x) }
#' @rdname getMDS
#' @export
getMDS.HMdata <- function(x,classical=FALSE,...){
    diss <- dist.HMdata(x)
    return(getMDS.dist(diss,classical=classical,...))
}
#' @rdname getMDS
#' @export
getMDS.DZdata <- function(x,classical=FALSE,...){
    diss <- dist.DZdata(x)
    return(getMDS.dist(diss,classical=classical,...))
}
#' @rdname getMDS
#' @export
getMDS.dist <- function(x,classical=FALSE,...){
    out <- list() 
    if (classical){
        out$points <- cmdscale(x)
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
#' dataset of class \code{HMdata}
#' @param x an object of class \code{HMdata}
#' @param ... optional arguments of the generic function
#' @return a matrix of clr coordinates
#' @examples
#' fname <- system.file("Major.csv",package="provenance")
#' Major <- read.HMdata(fname)
#' clrdat <- clr(Major)$x
#' diss <- dist(clrdat)
#' plot(getMDS(diss,classical=TRUE))
#' @export
clr <- function(x,...){ UseMethod("clr",x) }
clr.default <- function(x,...){stop()}
#' @rdname clr
#' @export
clr.HMdata <- function(x,...){
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
#' @param x an object of class \code{HMdata}
#' @param ... optional arguments to R's \code{princomp function}
#' @return an object of class \code{PCA}
#' @examples
#' fname <- system.file("Major.csv",package="provenance")
#' Major <- read.HMdata(fname)
#' plot(getMDS(Major,classical=TRUE))
#' dev.new()
#' plot(getPCA(Major),asp=1)
#' print("This example demonstrates the equivalence of classical MDS and PCA")
#' @export
getPCA <- function(x,...){
    clrdat <- clr(x)
    pc <- princomp(clrdat$x,...)
    class(pc) <- append("PCA",class(pc))
    return(pc)
}

#' Compositional biplot
#'
#' Plot the results of a principal components analysis as a biplot
#' @param x an object of class \code{PCA}
#' @param ... optional arguments of the \code{biplot} function
#' @examples
#' fname <- system.file("Major.csv",package="provenance")
#' Major <- read.HMdata(fname)
#' plot(getPCA(Major))
#' @method plot PCA
#' @export
plot.PCA <- function(x,...){
    biplot(x,...)
}

#' Get a subset of a provenance dataset
#'
#' Return a subset of provenance data according to some specified
#' indices
#' @param x an object of class \code{DZdata} or \code{HMdata}
#' @param i the indices of the samples to be returned
#' @param ... optional arguments for the generic subset function
#' @return an object of class \code{DZdata} or \code{HMdata}
#' @examples
#' fname <- system.file("HM.csv",package="provenance")
#' HM <- read.HMdata(fname)
#' i <- match(c("N1","N2","T8","T13","N12","N13"),names(HM))
#' foo <- subset(HM,i)
#' summaryplot(foo,ncol=2)
#' @export
subset <- function(x,...){ UseMethod("subset",x) }
#' @rdname subset
#' @export
subset.default <- function(x,...){base::subset(x,...)}
#' @rdname subset
#' @export
subset.DZdata <- function(x,i,...){
    out <- x
    if (length(x$err)==length(x$x)) out$err <- x$err[i]
    out$x <- x$x[i]
    return(out)
}
#' @rdname subset
#' @export
subset.HMdata <- function(x,i,...){
    out <- x
    out$x <- x$x[i,]
    return(out)
}

# returns list of distances between common items
getdistlist <- function(slist){
    dnames <- names(slist)
    lablist <- lapply(slist,function(x) names(x))
    commonlabels <- Reduce(intersect,lablist)
    for (name in dnames){
        i <- match(commonlabels,names(slist[[name]]))
        slist[[name]] <- subset(slist[[name]],i)
    }
    distlist <- slist
    for (name in dnames){
        distlist[[name]] <- dist(slist[[name]])
    }
    return(distlist)
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
#' @param ... a sequence of datasets of classes \code{DZdata}
#' and \code{HMdata}
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
#' DZ <- read.DZdata(system.file("DZ.csv",package="provenance"))
#' HM <- read.HMdata(system.file("HM.csv",package="provenance"))
#' GPA <- procrustes(DZ,HM)
#' plot(GPA)
#' @export
procrustes <- function(...) {
    dnames <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
    slist <- list(...)
    names(slist) <- dnames
    distlist <- getdistlist(slist)
    n <- length(labels(distlist[[1]]))
    m <- length(distlist)
    X <- array(dim=c(n,2,m))
    for (i in 1:m){
        md <- MASS::isoMDS(distlist[[i]],k=2)
        X[,,i] <- md$points
    }
    result <- shapes::procGPA(X, scale=TRUE, reflect=TRUE)
    out <- list()
    out$points <- result$mshape
    out$labels <- labels(distlist[[1]])
    class(out) <- "GPA"
    return(out)
}

#' Plot a Procrustes configuration
#'
#' Plots the group configuration of a Generalised Procrustes Analysis
#'
#' @param x an object of class \code{GPA}
#' @param ... optional arguments to the generic \code{plot} function
#' @examples
#' DZ <- read.DZdata(system.file("DZ.csv",package="provenance"))
#' HM <- read.HMdata(system.file("HM.csv",package="provenance"))
#' GPA <- procrustes(DZ,HM)
#' plot(GPA)
#' @method plot GPA
#' @export
plot.GPA <- function(x,...){
    plot(x$points[,1],x$points[,2],type="n",asp=1,...)
    text(x$points[,1],x$points[,2],x$labels)
}

#' Individual Differences Scaling of provenance data
#'
#' Performs 3-way Multidimensional Scaling analysis using Carroll and
#' Chang (1970)'s INdividual Differences SCALing method as implemented
#' using De Leeuw and Mair (2011)'s stress majorization algorithm.
#' @param ... a sequence of datasets of class \code{DZdata} or
#' \code{HMdata}
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
#' DZ <- read.DZdata(system.file("DZ.csv",package="provenance"))
#' HM <- read.HMdata(system.file("HM.csv",package="provenance"))
#' plot(indscal(DZ,HM))
#' @export
indscal <- function(...,type='ordinal'){
    dnames <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
    slist <- list(...)
    names(slist) <- dnames
    distlist <- getdistlist(slist)
    out <- smacof::smacofIndDiff(distlist, constraint = "indscal",type=type)
    class(out) <- "INDSCAL"
    return(out)
}

#' Plot an INDSCAL group configuration and source weights
#'
#' Given an object of class \code{INDSCAL}, generates two plots: the
#' group configuration and the subject weights. Together, these
#' describe a 3-way MDS model.
#' @param x an object of class \code{INDSCAL}
#' @param asp the aspect ratio of the plot
#' @param xlab a string with the label of the x axis
#' @param ylab a string with the label of the y axis
#' @param xaxt if = 'y', adds ticks to the x axis
#' @param yaxt if = 'y', adds ticks to the y axis
#' @param ... optional arguments to the generic plot function
#' @examples
#' DZ <- read.DZdata(system.file("DZ.csv",package="provenance"))
#' HM <- read.HMdata(system.file("HM.csv",package="provenance"))
#' plot(indscal(DZ,HM))
#' @method plot INDSCAL
#' @export
plot.INDSCAL <- function(x,asp=1,xlab="",ylab="",
                         xaxt='n',yaxt='n',...){
    par(mar=c(1,1,1,1))
    graphics::plot(x$gspace,type="n",asp=asp,xlab=xlab,
         ylab=ylab,xaxt=xaxt,yaxt=yaxt,...)
    text(x$gspace,labels=labels(x$obsdiss[[1]]))
    title('Group Configuration')
    X <- unlist(lapply(x$cweights,function(foo) foo[1,1]))
    Y <- unlist(lapply(x$cweights,function(foo) foo[2,2]))
    dev.new()
    graphics::plot(X,Y,type="n",asp=1,...)
    text(X,Y,names(x$cweights))
    title('Source Weights')
}

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
        dens[i] <- mean(density(dat,bw,from=(dat[i]-1e-10),
                                to=(dat[i]+1e-10),n=2)$y)
    }
    return(dens)
}

# adaptive KDE algorithm of Abramson (1982) as summarised by Jahn (2007)
Abramson <- function(dat,from,to,bw){
    n <- length(dat)
    lambda <- rep(0,n,length.out=n)
    pdens <- pilotdensity(dat,bw)
    G <- getG(pdens)
    dens <- rep(0,512)
    for (i in 1:n){
        lambda = sqrt(G/pdens[i])
        dens <- dens + density(dat[i],bw*lambda,from=from,to=to)$y
    }
    return(dens)
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
#' #' @param adaptive boolean flag controlling if the adaptive KDE
#' modifier of Abramson (1982) is used
#' @param log transform the ages to a log scale if TRUE
#' @examples
#' fname <- system.file("DZ.csv",package="provenance")
#' DZ <- read.DZdata(fname)
#' samp <- DZ$x[['N1']]
#' dens <- getKDE(samp,0,3000)
#' plot(dens)
#' @export
getKDE <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE){
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
        bw <- bw/(median(x)*log(10))
    } else {
        d <- x
    }
    out$x <- seq(from=from,to=to,length.out=512)
    if (is.na(bw)){ bw <- botev(d) }
    if (adaptive){
        out$y <- Abramson(d,from,to,bw)
    } else {
        out$y <- density(d,bw,from=from,to=to)$y
    }
    if (log) out$x <- 10^(out$x)
    out$y <- out$y/(sum(out$y)*(to-from)/512)
    out$bw <- bw
    out$ages <- x
    return(out)
}

#' Generate an object of class \code{KDEs}
#'
#' Convert a dataset of class \code{DZdata} into an object of
#' class \code{KDEs} for further processing by the
#' \code{summaryplot} function.
#' @param x an object of class \code{DZdata}
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
#' @examples
#' fname <- system.file("DZ.csv",package="provenance")
#' DZ <- read.DZdata(fname)
#' KDEs <- getKDEs(DZ,0,3000,pch=NA)
#' summaryplot(KDEs,ncol=3)
#' @export
getKDEs <- function(x,from=NA,to=NA,bw=NA,samebandwidth=TRUE,
                 adaptive=TRUE,pch=NA,normalise=FALSE,log=FALSE){
    if (is.na(from) | is.na(to)) {
        mM <- setmM(unlist(x$x),from,to,log)
        from <- mM$m
        to <- mM$M
    }
    snames <- names(x$x)
    thekdes <- list()
    n <- list()
    themax <- -1
    if (is.na(bw) & samebandwidth) bw <- getcommonbandwidth(x)
    for (name in snames){
        thekdes[[name]] <- getKDE(x$x[[name]],from,to,bw,adaptive,log)
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

emptyplot <- function(){
    plot(c(0,1),c(0,1),type='n',axes=FALSE,xlab="",ylab="")
}

#' @export
names.DZdata <- function(x){
    return(names(x$x))
}
#' @export
names.HMdata <- function(x){
    return(rownames(x$x))
}
#' @export
names.KDEs <- function(x){
    return(names(x$kdes))
}

#' Plot a kernel density estimate
#'
#' Plots an object of class \code{KDE}
#' @param x an object of class \code{KDE}
#' @param pch the symbol used to show the samples. Set \code{pch = NA}
#' to turn them off
#' @param xlab the label of the x-axis
#' @param ylab the label of the y-axis
#' @param ... optional parameters to be passed on to the graphics object
#' @examples
#' fname <- system.file("DZ.csv",package="provenance")
#' DZ <- read.DZdata(fname)
#' samp <- DZ$x[['N1']]
#' dens <- getKDE(samp,0,3000)
#' plot(dens)
#' @method plot KDE
#' @export
plot.KDE <- function(x,pch='|',xlab="age [Ma]",ylab="",...){
    if (x$log) {
        plot(x$x,x$y,type='l',log="x",xlab=xlab,ylab=ylab,...)
    } else {
        plot(x$x,x$y,type='l',xlab=xlab,ylab=ylab,...)
    }
    points(x$ages,rep(par("usr")[3]/2,length(x$ages)),pch=pch)
    text(tail(x$x,n=1),.9*max(x$y),paste0("n=",length(x$ages)),pos=2)
}

# Plot multiple KDEs. Used by summaryplot function
plot.KDEs <- function(x,sname,annotate=TRUE,...){
    if (x$themax>0){ # normalise
        M <- x$themax
    } else {
        M <- max(x$kdes[[sname]]$y)
    }
    if (annotate){
        plot(x$kdes[[sname]],pch=NA,...)
    } else {
        plot(x$kdes[[sname]],pch=x$pch,
             axes=FALSE,xlab="",ylab="",...)
    }
}

#' Plot continuous data as histograms or cumulative age distributions
#'
#' Plot one or several samples from a \code{DZdata} dataset as a
#' histogram or Cumulative Age Distributions (CAD).
#' @param x an object of class \code{DZdata}
#' @param snames a string or a vector of string with the names of the
#' samples that need plotting if \code{snames} is a vector, then the
#' function will default to a CAD.
#' @param annotate boolean flag indicating whether the x- and y-axis
#' should be labeled
#' @param CAD boolean flag indicating whether the data should be
#' plotted as a cumulative age distribution or a histogram. For
#' multi-sample plots, the function will override this value with
#' \code{TRUE}.
#' @param pch an optional the symbol to mark the sample points along
#' the CAD
#' @param verticals boolean flag indicating if the horizontal lines of
#' the CAD should be connected by vertical lines
#' @param col a colour map
#' @param ... optional arguments to the generic \code{plot} function
#' @examples
#' DZ <- read.DZdata(system.file("DZ.csv",package="provenance"))
#' plot(DZ,c('N1','N2'),CAD=TRUE)
#' @method plot DZdata
#' @export
plot.DZdata <- function(x,snames=NULL,annotate=TRUE,CAD=FALSE,
                        pch=NA,verticals=TRUE,col=par("col"),...){
    if (is.null(snames)) snames <- names(x)
    if (length(snames)>1){
        n <- length(snames)
        if (col==par("col")) col <- heat.colors(n)
        plot(ecdf(x$x[[snames[1]]]),pch=pch,verticals=verticals,
             col=col[1],xlab=x$xlabel,main="",...)
        for (i in 2:n){
            lines(ecdf(x$x[[snames[i]]]),pch=pch,verticals=verticals,
                  col=col[i],xlab=x$xlabel,main="",...)
        }
        legend("bottomright",legend=snames,lwd=1,col=col)
    } else {
        if (annotate){
            if (CAD) { plot(ecdf(x$x[[snames]]),pch=pch,
                            verticals=verticals,col=col,main=snames,...) }
            else { hist(x$x[[snames]],x$breaks,col=col) }
        } else {
            if (CAD) { plot(ecdf(x$x[[snames]]),pch=pch,verticals=verticals,
                        axes=FALSE,xlab="",ylab="",main="",col=col,...)
            } else {hist(x$x[[snames]],x$breaks, axes=FALSE,xlab="",ylab="",main="",col=col) }
        }
    }
}

#' Plot a pie chart
#'
#' Plots an object of class \code{HMdata} as a pie chart
#' @param x an object of class \code{HMdata}
#' @param sname the sample name
#' @param annotate a boolean flag controlling if the pies of the
#' pie-chart should be labeled
#' @param col colour scheme to fill the pie charts
#' @param ... optional parameters to be passed on to the graphics object
#' @examples
#' fname <- system.file("HM.csv",package="provenance")
#' HM <- read.HMdata(fname)
#' plot(HM,'N1')
#' @method plot HMdata
#' @export
plot.HMdata <- function(x,sname,annotate=TRUE,col='default',...){
    i <- which(names(x) %in% sname)
    if (col=='default') col <- heat.colors(length(colnames(x$x)))
    pos <- x$x[i,]>0
    if (annotate){
        pie(unlist(x$x[i,]),col=col,...)
    } else {
        pie(unlist(x$x[i,]),labels=NA,col=col,...)
    }
}

# annotation of various plots used in summaryplot function
annotation <- function(x,...){ UseMethod("annotation",x) }
annotation.default <- function(x,...){stop()}
annotation.KDEs <- function(x,height=NA,...){
    if (is.na(height)){ par(mar=c(2,0,0,0)) }
    else { par(mai=c(height/2,0,0,0)) }
    plot(c(x$from,x$to),c(0,1),type='n',axes=FALSE,xlab="",ylab="")
    Axis(side=1)
    text(x=.5*(x$from+x$to),.5,label=x$xlabel)
    par(mar=c(0,0,0,0))
}
annotation.HMdata <- function(x,...){
    labels <- colnames(x$x)
    comp <- rep(1,length(labels))
    pie(comp,labels=labels,col=heat.colors(length(labels)))
}
annotation.DZdata <- function(x,height=NA,...){
    if (is.na(height)){ par(mar=c(2,0,0,0)) }
    else { par(mai=c(height/2,0,0,0)) }
    m <- min(x$breaks)
    M <- max(x$breaks)
    plot(c(m,M),c(0,1),
         type='n',axes=FALSE,xlab="",ylab="")
    Axis(side=1)
    text(x=.5*(m+M),.5,label=x$xlabel)
    par(mar=c(0,0,0,0))
}

#' Joint plot of several provenance datasets
#'
#' Arranges kernel density estimates and pie charts in a grid format
#' 
#' @param ... a sequence of datasets of class \code{HMdata},
#' \code{KDEs}, or \code{DZdata}
#' @param ncol the number of columns
#' @return a summary plot of all the data comprised of KDEs for the
#' datasets of class \code{KDEs}, pie charts for those of class
#' \code{HMdata} and histograms for those of class \code{DZdata}.
#' @examples
#' DZ <- read.DZdata(system.file("DZ.csv",package="provenance"))
#' KDEs <- getKDEs(DZ,0,3000)
#' HM <- read.HMdata(system.file("HM.csv",package="provenance"))
#' PT <- read.HMdata(system.file("PT.csv",package="provenance"))
#' summaryplot(KDEs,HM,PT,ncol=2)
#' @export
summaryplot <- function(...,ncol=1){
    dnames <- sapply(match.call(expand.dots=TRUE)[-1:-2], deparse)
    dlist <- list(...)
    names(dlist) <- dnames
    classes <- unlist(lapply(dlist,class))
    nd <- length(dlist) # number of datasets
    snames <- unique(unlist(lapply(dlist,names)))
    ns <- length(snames)
    w <- rep(1,nd) # column widths
    w[which( classes %in% c("KDEs","DZdata") )] <- 2
    w <- rep(c(1,w),ncol)
    nppc <- ceiling(ns/ncol)
    np <- (nppc+1)*ncol*(nd+1) # number of subpanels
    layout(matrix(1:np,nppc+1,length(w)),w,rep(1,nppc+1))
    si <- ceiling(seq(from=0,to=ns,length.out=ncol+1)) # sample index
    par(xpd=TRUE,mar=c(0,0,0,0))
    for (i in 1:ncol){ # loop through columns
        for (j in (si[i]+1):(si[i+1])){
            sname <- snames[j]
            emptyplot()
            text(0,0.5,labels=sname,pos=4)
        }
        if (si[i+1]-si[i]<nppc) emptyplot()
        emptyplot()
        for (d in dlist){ # loop through datasets
            for (j in (si[i]+1):si[i+1]){
                sname <- snames[j]
                if (sname %in% names(d)){
                    plot(d,sname,annotate=FALSE)
                } else {
                    emptyplot()
                }
            }
            if (si[i+1]-si[i]<nppc) emptyplot()
            if (i==ncol){
                ds <- dev.size()[2]/(nppc+1)
                annotation(d,height=ds)
            } else {
                emptyplot()
            }
        }
    }
}

lines.ternary <- function(type='empty'){
    if (type=='QFL'){
        xy1 <- xyz2xy(c(90,90),c(10,0),c(0,10))
        xy2 <- xyz2xy(c(75,75),c(25,0),c(0,25))
        xy3 <- xyz2xy(c(0,75),c(25,25/4),c(75,25*3/4))
        xy4 <- xyz2xy(c(0,75),c(75,25*3/4),c(25,25/4))
        xy5 <- xyz2xy(c(0,90),c(50,5),c(50,5))
        lines(xy1); lines(xy2); lines(xy3); lines(xy4); lines(xy5)
        text(xyz2xy(20,-1,2),labels='quartarenite',adj=0)
        text(xyz2xy(40,-2,10),labels='sublitharenite',adj=0)
        text(xyz2xy(40,10,-2),labels='subarkose',adj=1)
        text(xyz2xy(30,70,10),labels='arkose',srt=68)
        text(xyz2xy(30,50,30),labels='lithic arkose',srt=85)
        text(xyz2xy(30,30,50),labels='feldspathic litharenite',srt=-85)
        text(xyz2xy(30,10,70),labels='litharenite',srt=-68)
        return(c('Q','F','L'))
    }
    if (type=='QFL.dickinson') {
        xy1 <- xyz2xy(c(97,0),c(0,85),c(3,15))
        xy2 <- xyz2xy(c(25,51.6),c(0,40),c(75,8.4))
        lines(xy1);lines(xy2)
        text(xyz2xy(20,40,40),labels='magmatic arc')
        text(xyz2xy(55,15,30),labels='recycled orogen')
        text(xyz2xy(50,45,5),labels='continental block',srt=65)
        return(c('Q','F','L'))
    }
}

#' Plot a ternary diagram
#'
#' Plots triplets of compositional data on a ternary diagram
#' @param x vector with the first variable
#' @param y vector with the second variable
#' @param z vector with the third variable
#' @param type adds annotations to the ternary diagram, one of either
#' \code{empty}, \code{QFL}, \code{QFL.dickinson} or
#' \code{QmFLt.dickinson}
#' @param pch plot symbol
#' @param labels vector of strings to be added to the plot symbols
#' @param ... optional arguments to the generic \code{points} function
#' @examples
#' PT <- read.HMdata(system.file("PT.csv",package="provenance"))
#' q <- PT$x[,'Q']
#' f <- PT$x[,'KF'] + PT$x[,'P']
#' l <- PT$x[,'Lm'] + PT$x[,'Lv'] + PT$x[,'Ls']
#' ternaryplot(q,f,l,type='QFL.dickinson',labels=names(PT))
#' @export
ternaryplot <- function(x,y,z,type='empty',pch=NA,labels=names(x),...){
    par(mar=c(1,1,1,1))
    plot(c(0,1),c(0,1),type='n',xaxt='n',yaxt='n',
         xlab='',ylab='',asp=1,bty='n')
    xyzlabels <- lines.ternary(type=type)
    corners <- xyz2xy(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0))
    lines(corners)
    if (any(is.na(xyzlabels))){
        xlab <- deparse(substitute(x))
        ylab <- deparse(substitute(y))
        zlab <- deparse(substitute(z))
    } else {
        xlab <- xyzlabels[1]
        ylab <- xyzlabels[2]
        zlab <- xyzlabels[3]
    }
    text(corners[1:3,],labels=c(xlab,ylab,zlab),pos=c(3,1,1))
    xy <- xyz2xy(x,y,z)
    if (is.na(pch) & is.null(labels)){ pch <- 1 }
    if (!is.na(pch)) points(xy,pch=pch,...)
    if (!is.null(labels)){ text(xy,labels=labels,pos=1) }
}

xyz2xy <- function(x,y,z){
    n <- length(x)
    bl <- (1-sin(pi/3))/2 # baseline
    xy <- matrix(0,nrow=n,ncol=2)
    xy[,1] <- 0.5*(x+2*z)/(x+y+z)
    xy[,2] <- bl + sin(pi/3)*x/(x+y+z)
    return(xy)
}

#' Calculate the probability of missing a given population fraction
#'
#' For a given sample size, returns the likelihood of missing any
#' fraction greater than a given size
#' @param n the number of grains in the detrital sample
#' @param f the size of the smallest resolvable fraction (0<f<1)
#' @return the probability that all N grains in the sample have missed
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

test <- function(){
    setwd("/home/pvermees/Dropbox/provenance/provenance")
    graphics.off()

    DZ <- read.DZdata("./inst/DZ.csv","./inst/DZerr.csv")
    DZerr <- read.DZdata("./inst/DZerr.csv","./inst/DZerr.csv")
    HM <- read.HMdata("./inst/HM.csv")
    PT <- read.HMdata("./inst/PT.csv")
    Major <- read.HMdata("./inst/Major.csv")
    Trace <- read.HMdata("./inst/Trace.csv")

#    plot(getKDE(DZ$x[['N1']],log=TRUE))
    
#    Q <- PT$x[,'Q']
#    F <- PT$x[,'KF'] + PT$x[,'P']
#    L <- PT$x[,'Lm'] + PT$x[,'Lv'] + PT$x[,'Ls']
#    ternaryplot(60,25,15,type='QFL')
    
#    plot(DZ,'N1',CAD=TRUE)

    KDEs <- getKDEs(DZ)
    summaryplot(KDEs,ncol=4)
    
#    KDEs <- getKDEs(DZ,0,3000)
#    summaryplot(KDEs,HM,PT,Major,Trace,ncol=2)

#    ind <- indscal(DZ,HM)
#    plot(ind)
#    graphics.off()
#    plot(getMDS(Major,classical=TRUE),nnlines=TRUE)
    
#    i <- match(c("N1","N2","T8","T13","N12","N13"),names(HM))
#    foo <- subset(HM,i)
#    summaryplot(foo,ncol=2)

#    plot(indscal(DZ,HM))
#    saveplot("TarimPCA.pdf",2)

#    print(get.f(60))
}

#test()
