#' Adaptive Neyman test in DEGraph package
#'
#' @param X1 nxT (number of reps x width of interval). A matrix represents coverage for condition 1.
#' @param X2 nxT (number of reps x width of interval). A matrix represents coverage for condition 2.
#' @param candK a vector, candidate values for the true number of Fourier components.
#' @param na.rm A logical value indicating whether variables with NA in at least one of the n1
#   + n2 observations should be discarder before the test is performed.
#' @param pool
#' @param poolVarX
#' @param poolVarY
#'
#' @return The list of adaptive Neyman statistic \code{T}, p.value \code{p}, \code{kstar}, \code{varX}, and  \code{varY}.
#' @export
#'
#' @examples
AN.test <- function (X1, X2, candK = 1:ncol(X1), na.rm = FALSE, pool = FALSE, poolVarX, poolVarY)
{
    varX <- varY <- T_stat <- c()
    if (na.rm) {
        na1 <- apply(X1, 2, FUN = function(x) sum(is.na(x)))
        na2 <- apply(X2, 2, FUN = function(x) sum(is.na(x)))
        idxs <- which((na1 == 0) & (na2 == 0))
        X1 <- X1[, idxs]
        X2 <- X2[, idxs]
    }
    # Avoid var = 0 when two counts at a position t are the same
    id1_same <- c()
    id2_same <- c()
    id_same <- c()
    for (t in 1:ncol(X1)) {
        if (max(X1[,t]) == min(X1[,t]) ) {
            id1_same <- c(id1_same,t)
        }

        if (max(X2[,t]) == min(X2[,t]) ) {
            id2_same <- c(id2_same,t)
        }
    }
    id_same <- union(id1_same, id2_same)
    if (length(id_same) > 0 ) {
        X1 <- X1[, -id_same]
        X2 <- X2[, -id_same]
    }
    p <- ncol(X1)
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    Tstar <- -Inf
    kstar <- NA

    if (pool == FALSE) {
        for (kk in 1:ncol(X1)) {
            if ( ncol(X1) == 0 | ncol(X2) == 0 ) {
                message("Peak with zero length!")
                varX <- varY <- numeric(0)
            } else {
                # Compute Cov of columns: Var1 is the pooled var of all vars of columns of X
                # print(paste("+++kk = ", kk))
                # print(paste("+ ncol X1 = ", ncol(X1) ))
                var1 <- mean(diag(var(X1[, 1:kk, drop = FALSE])))
                var2 <- mean(diag(var(X2[, 1:kk, drop = FALSE])))
                X <- colMeans(X1[, 1:kk, drop = FALSE]) - colMeans(X2[, 1:kk, drop = FALSE])

                varX[kk] <- var1; varY[kk] <- var2
                X <- X/sqrt(var1/n1 + var2/n2)
                tmp <- sum(X^2 - 1)/sqrt(2 * kk)

                if ( (tmp > Tstar) & (!is.na(tmp)) ) {
                    Tstar <- tmp
                    kstar <- kk
                }
            }
        }
    }
    else {
        for (kk in 1: ncol(X1)) {
            if ( ncol(X1) == 0 | ncol(X2) == 0 ) {
                message("Peak with zero length!")
                #varX <- varY <- T_stat <- numeric(0)
                return(list(statistic = NA, p.value = NA, kstar = NA, varX = numeric(0), varY = numeric(0)))
            }
            else {
                # Compute Cov of columns: Var1 is the pooled var of all vars of columns of X
                var1 <- poolVarX[kk]
                var2 <- poolVarY[kk]
                X <- colMeans(X1[, 1:kk, drop = FALSE]) - colMeans(X2[,1:kk, drop = FALSE])

                varX[kk] <- var1; varY[kk] <- var2
                X <- X/sqrt(var1/n1 + var2/n2)
                tmp <- sum(X^2 - 1)/sqrt(2 * kk)
                if ( (tmp > Tstar) & (!is.na(tmp)) ) {
                    Tstar <- tmp
                    kstar <- kk
                }
            }
        }
    }
    T_stat <- sqrt(2 * log(log(p))) * Tstar - (2 * log(log(p)) + log(log(log(p)))/2 -
                                              log(4 * pi)/2)
    p <- 1 - exp(-exp(-T_stat))
    res <- list(statistic = T_stat, p.value = p, kstar = kstar, varX = varX, varY = varY)
    class(res) <- "htest"
    res
}


#' Compute Adaptive Neyman Test after Discrete Fourier Transform.
#'
#' @param X A matrix represents coverage for condition 1 (T x n1).
#' @param Y A matrix represents coverage for condition 1 (T x n2).
#' @param maxdim the number of the maximum frequency to be tested.
#'   maxdim = L/2 seems a good choice, but you are welcome to use maxdim = L.
#' @param poolVar
#' @param numPar
#' @param candK a vector, candidate values for the true number of components to be tested.
#' @param na.rm
#'
#' @return a list of \code{mhat} and Adaptive Neyman test \code{T_AN}.
#' @export
#'
#' @examples
aneyman <- function(X,Y, maxdim = floor(nrow(X) / 2), poolVar = FALSE, numPar = 5, candK = 1:nrow(X), na.rm = FALSE)
{
    if (na.rm) {
        X1 <- t(X)
        X2 <- t(Y)
        na1 <- apply(X1, 2, FUN = function(x) sum(is.na(x)))
        na2 <- apply(X2, 2, FUN = function(x) sum(is.na(x)))
        idxs <- which((na1 == 0) & (na2 == 0))
        X1 <- X1[, idxs]
        X2 <- X2[, idxs]
        X <- t(X1)
        Y <- t(X2)
    }

    # Avoid var = 0 when two counts at a position t are the same
    id1_same <- c()
    id2_same <- c()
    id_same <- c()
    X1 <- t(X)
    X2 <- t(Y)
    for (t in 1:ncol(X1)) {
        if (max(X1[,t]) == min(X1[,t]) ) {
            id1_same <- c(id1_same,t)
        }

        if (max(X2[,t]) == min(X2[,t]) ) {
            id2_same <- c(id2_same,t)
        }
    }
    id_same <- union(id1_same, id2_same)
    if (length(id_same) > 0 ) {
        X1 <- X1[, -id_same]
        X2 <- X2[, -id_same]
    }
    X <- t(X1)
    Y <- t(X2)

    ##
    X <- apply(X[candK, ],2,fft)                 #sum( exp(-2*pi* i*(0:(T-1))/T j)*X )
    Y <- apply(Y[candK, ],2,fft)
    # L = T, the number of sampling points
    L <- nrow(X)
    R1 <- Re(X); I1 <- Im(X)
    R2 <- Re(Y); I2 <- Im(Y)
    # store real parts of X, Y
    X <- R1; Y <- R2
    # even indices stored for real parts
    index <- seq(2,L,2)
    X[index,] <- R1[2:floor(L/2+1),]
    Y[index,] <- R2[2:floor(L/2+1),]
    # odd indices stored for imaginary parts
    index <- seq(3,L,2)
    X[index, ] <- I1[2:floor(L/2+0.5),]
    Y[index, ] <- I2[2:floor(L/2+0.5),]
    # Compute means
    R1 <- apply(X,1,mean)
    R2 <- apply(Y,1,mean)
    # Compute vars
    if (poolVar==FALSE){
        I1 <- apply(X,1,var)
        I2 <- apply(Y,1,var)
    } else {
        I1 <- poolVar(X, numPar = numPar)
        I2 <- poolVar(Y, numPar = numPar)
    }
    # standardized difference Z*(t)
    stddiff <- (R1-R2)/sqrt(I1/ncol(X)+I2/ncol(Y))
    # n = the number of components used to compute adaptive Neyman test
    n <- maxdim
    stddiff <- stddiff[1:n]
    # stddiff = T*_AN
    stddiff <- cumsum(stddiff^2-1)/sqrt(2*(1:n))
    # mhat defined in section 2.1 Adaptive Neyman Test
    mhat <- order(stddiff)[n]
    T <-  sqrt(2* log(log(n))) *stddiff - (2*log(log(n))+
                                               0.5*log(log(log(n))) - 0.5*log(4*pi) )
    T <- T[mhat]
    p <- 1 - exp(-exp(-T))
    res <- list(statistic = T, p.value = p, mhat = mhat)
    class(res) <- "htest"
    res
}

#' Eval p-values from different quantile.
#'
#' @param P a slot PvalList of class tanDb.
#' @param total Total genomic intervals.
#' @param nSamples The sample size for each conditions (n1=n2=n).
#' @param quant Quantile used for combined pvales.
#' @param BH  A logical if multiple testing correction using B-H method is applied. Defaul is FALSE.
#' @param na.rm logical; if true, any NA and NaN's are removed before computing.
#'
#' @author Duy Nguyen on February 26, 2016
#' @return Combined (Raw) p-values with the given quantile.
#' @export
#'
#' @examples
evalPvals <- function(P, total = nrow(P[['pval']]), quant = 1, nSamples, BH = FALSE, na.rm = TRUE ) {
    ## Extract original p.vals matrix p without quantile q
    p <- P[['pval']]
    p <- p[, 1:( ncol(p) - 1)]
    Between_cols <- numeric()
    Pc <- rep(NA, dim(p)[1])
    if (nSamples == 4) {
        Between_cols <- 6 * 6
        Pc <- apply(p[,1:Between_cols], 1, function(x) quantile(x, probs = quant, na.rm = na.rm))
    }
    else if (nSamples == 3) {
        Between_cols <- 3 * 3
        Pc <- apply(p[,1:Between_cols], 1, function(x) quantile(x, probs = quant, na.rm = na.rm))
    }
    else if (nSamples == 2) {
        Pc <- p
    }
    if (BH) {
        Pc <- p.adjust(Pc, method = "BH")
    }
    # fdr <- p.adjust(as.vector(p), method='BH')
    # FDR <- matrix(0,nrow=nrow(p),ncol=ncol(p)+1)
    # for (i in 1:ncol(p)){
    #     FDR[,i] <- fdr[(i-1)*nrow(p)+(1:nrow(p))]
    # }
    # if (any(is.na(FDR))){
    #     message(length(which(is.na(FDR))),' NAs found  (of ', length(FDR),')')
    #     FDR[is.na(FDR)]=min(FDR[!is.na(FDR)])
    # }
    # if (nSamples %in% c(3,4)) {
    #     FDR[, i+1] <- apply(FDR[,1:Between_cols], 1, function(x) quantile(x, probs = quant))
    # }
    # else if (nSamples == 2) {
    #     FDR[, i+1] <- FDR[, 1]
    # }
    # fdr <- matrix(NA,nrow=total, ncol=ncol(P[['FDR']]))
    # fdr[1:total,] <- FDR
    # colnames(fdr) <- colnames(P[['FDR']])
    # last_col <- dim(fdr)[2]
    # Pvals = fdr[,last_col]
    return(Pc)
}

#' Title
#'
#' @param p.list A list resulted from running \code{computePvals_batch}.
#' @param totalPeaks A total number of peaks tested in batch mode.
#' @param nSamples The sample size for each conditions (n1=n2=n).
#'
#' @return
#' @export
#'
#' @examples
create_pMat <- function(p.list, totalPeaks, nSamples) {
    pMat <- matrix()
    if (nSamples == 4) {
        pMat <- matrix(NA, nrow = totalPeaks, ncol = 6 * 6 + 3 * 2)
    }
    else if (nSamples == 3) {
        pMat <- matrix(NA, nrow = totalPeaks, ncol = 3 * 3 + 3 * 2)
    }
    else if (nSamples == 2) {
        pMat <- matrix(NA, nrow = totalPeaks, ncol = 1)
    }
    for (bin in 1:length(p.list)) {
        between <- p.list[[bin]]
        if (length(between) > 0) {
            print(paste("Combine bin :", bin, sep = " "))
            for (j in 1:length(between)) {
                pMat[between[[j]]$sites, j] <- between[[j]]$values
            }
        }
    }
    pMat
}

#' @useDynLib tan
#' @importFrom Rcpp sourceCpp
NULL
#' @import methods
NULL
#' @import foreach
NULL


