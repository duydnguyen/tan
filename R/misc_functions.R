#' Give the tag Counts and genomic positions of a given consensus interval.
#'
#' @param site Site indexed in c_ranges.
#' @param Rep Replication ChIP data which is obtained from loading RData.
#' @param c_ranges Ranges(GRanges object).
#' @param normalized A vector of normalized values for each ChIP samples.
#' @author Duy Nguyen on April 6, 2015
#' @return coord (Genomic Positions) and y.coord (normalized tag Counts)
#' @export
#'
#' @examples
getXY <- function(site, Rep, c_ranges, normalized = 1) {
    # Get boundary interval
    c_start <- c_ranges@start[site]
    c_end <- c_start + c_ranges@width[site] - 1
    x.coord = Rep[[site]][, 1]
    bd_left = which(x.coord == c_start)
    bd_right = which(x.coord == c_end)
    x.coord = Rep[[site]][ bd_left:bd_right, 1]
    y.coord = Rep[[site]][ bd_left:bd_right, 2] * normalized
    return(cbind(x.coord, y.coord))
}


#' Compute the MMD distance in MMDiff package
#'
#' @param X An 1-by-T matrix to store 2 curves
#' @param labels 1 for first curve; -1 for second curve
#' @param sigma NULL
#'
#' @return A data frame of MMD; MMD distance is stored in d$val
#' @export
#'
#' @examples
mmd <-function(X,labels,sigma=NULL) {

    poslabels <- which(labels==1);
    neglabels <- which(labels==-1);

    if (length(poslabels)==0 || length(neglabels)==0){
        stop('kmd(): each sample should have at least one observation ');
    }

    # KERNEL

    # RULE OF THUMB to determine Kernel width
    # calculates the kernel matrix and estimate sigma

    x <- as.matrix(X[poslabels,],nrow=length(poslabels));
    y <- as.matrix(X[neglabels,],nrow=length(neglabels));

    ## Kxx
    G <- x%*%t(x);
    L <- ncol(G);
    nor <- rep(G[seq(1,L^2,L+1)],L);
    Kxx <- -2*G +matrix(nor, nrow=L,byrow=TRUE)+matrix(nor, nrow=L)
    rm(ls='G', 'L', 'nor')

    ## Kyy
    G <- y%*%t(y);
    L <- ncol(G);
    nor <- rep(G[seq(1,L^2,L+1)],L);
    Kyy <- -2*G +matrix(nor, nrow=L,byrow=TRUE)+matrix(nor, nrow=L)
    rm(ls='G', 'L', 'nor')

    ## Kxy
    G <- x%*%t(y);
    Ly <- ncol(G);
    Lx <- nrow(G);
    norx <- rep(rowSums(x*x),Ly);
    nory <- rep(rowSums(y*y),Lx);
    Kxy <- t(-2*G +matrix(norx, nrow=Lx)+matrix(nory, nrow=Lx,byrow=TRUE))
    rm(ls='G', 'Lx','Ly', 'norx','nory')


    #now get the median distance
    mdist <- median(Kxy[Kxy!=0]);
    mdist <- sqrt(mdist/2);
    if (is.null(sigma)){
        sigma <- mdist;
        if (sigma ==0){
            sigma <- 1;}
    }
    #apply RBF
    Kxx <- exp(-1/2/sigma^2 * Kxx);
    Kyy <- exp(-1/2/sigma^2 * Kyy);
    Kxy <- exp(-1/2/sigma^2 * Kxy);

    rm('x' ,'y')

    # fprintf('Using RBF kernel with sigma=#1.2f\n',kernelsize);

    ############
    # calculates MMD

    m <- nrow(Kxx);
    n <- nrow(Kyy);

    N <- max(m,n);
    M <- min(m,n);


    #Kxx
    sumKxx <- sum(Kxx);
    if (m!=n){
        sumKxx_M <- sum(sum(Kxx[1:M,1:M]));
    } else {
        sumKxx_M <- sumKxx;
    }
    dgxx <- Kxx[seq(1,m^2,m+1)];

    sumKxxnd <- sumKxx - sum(dgxx); #no diags

    # R = max(dgxx); #upper bound on kernel, should be one
    # R_M = max(dgxx[1:M]);

    # h_u = colSums(Kxx[1:M,1:M]) - dgxx[1:M]; #one sided sum

    #Kyy
    sumKyy <- sum(Kyy);
    if (m!=n){
        sumKyy_M <- sum(sum(Kyy[1:M,1:M]));
    } else {
        sumKyy_M <- sumKyy;
    }

    dgyy <- Kyy[seq(1,n^2,n+1)];

    sumKyynd <- sumKyy - sum(dgyy); #no diags

    # R = max(R,max(dgyy)); #upper bound on kernel, should be one
    # R_M = max(R,max(dgyy[1:M]));

    # h_u = h_u + colSums(Kyy[1:M,1:M]) - dgyy[1:M]; #one sided sum


    #Kxy
    sumKxy <- sum(Kxy);
    if (m!=n){
        sumKxy_M <- sum(sum(Kxy[1:M,1:M]));
    } else {
        sumKxy_M <- sumKxy;
    }

    # dg = Kxy[seq(1,nrow(Kxy)*M,nrow(Kxy)+1)]; #up to M only

    # h_u = h_u -  colSums(Kxy[1:M,1:M]) - colSums(t(Kxy[1:M,1:M])) + 2*dg; #one sided sum


    #compute MMDs
    #corrolary 11 (biased)
    biased <- sqrt(sumKxx/(m*m) +  sumKyy/(n * n) - 2/m/n * sumKxy);

    #equation 14
    #only for m==n (=M) (unbiased)
    #unbiased = sum(h_u)/M/(M-1);
    #mmd=list(biased,unbiased,sigma)
    #names(mmd)=c('biased','unbiased','kernelsize')

    mmd <- list(biased,sigma,mdist)
    names(mmd) <- c('biased','kernelsize','mdist')
    return(mmd)

}

getLabel <- function(scores, cutoff, pval=FALSE) {
    length = length(scores)
    lab_pred <- rep(NA, length )
    if (pval==FALSE) {
        for (i in 1:length) {
            if (scores[i] > cutoff & scores[i]!=100) lab_pred[i] <- 1
            else lab_pred[i] <- 0
        }
    } else {
        for (i in 1:length) {
            if (scores[i] < cutoff & scores[i]!=100) lab_pred[i] <- 1
            else lab_pred[i] <- 0
        }
    }
    return(lab_pred)
}

evalPC <- function(scores, lenT = 100, FT = FALSE, coverage, pval=FALSE ) {
    labels <- coverage$labels
    index_0 <- which(labels ==0)
    index_1 <- which(labels ==1)
    idNA <- which(is.na(scores)==TRUE)
    if (length(idNA) > 0) {
        ## label NA = 100
        scores[idNA] <- 100
    }
    min <- min(scores, na.rm = TRUE)
    max <- max(scores, na.rm = TRUE)

    if (FT == TRUE) {
        T_seq <- c(seq(min, 5, by = 0.1), seq(5, max, by = 10) )
    }
    else {
        T_seq <- seq(min, max, length.out = lenT)
    }
    T_seq <- c(-Inf, T_seq, Inf)
    lenT = length(T_seq)
    Precision <- rep(0, lenT)
    Sensitivity <- rep(0, lenT)
    Specificity <- rep(0, lenT)
    for (i in 1:lenT) {
        print(paste(i,' out of', lenT))
        cutoff <- T_seq[i]
        lab_pred <- getLabel(scores, cutoff, pval = pval)
        index_1_pred <- which(lab_pred == 1)
        index_0_pred <- which(lab_pred == 0)
        ## Compute TP
        TP <- length( intersect(index_1, index_1_pred) )
        FP <- length(index_1_pred) - TP
        FN <- length(index_1) - TP
        TN <- length(index_0) - FP
        Precision[i] = TP / (TP + FP)
        Sensitivity[i] = TP / (TP + FN)
        Specificity[i] = TN / (FP + TN)
    }
    PC <- data.frame(Precision, Sensitivity, Specificity)
    colnames(PC) <- c('Precision', 'Sensitivity', 'Specificity')
    return(PC)
}

evalPC_sep <- function(scores, lenT = 100, affinity_check= TRUE, coverage, pval = FALSE) {
    labels <- coverage$labels
    index_0 <- which(labels ==0)
    idNA <- which(is.na(scores)==TRUE)
    if (length(idNA) > 0) {
        ## label NA = 100
        scores[idNA] <- 100
    }
    min <- min(scores, na.rm = TRUE)
    max <- max(scores, na.rm = TRUE)
    T_seq <- seq(min, max, length.out = lenT)
    #T_seq <- c(seq(min, 5, by = 0.1), seq(5, max, by = 10) )
    T_seq <- c(-Inf, T_seq, Inf)
    lenT = length(T_seq)
    Precision <- rep(0, lenT)
    Sensitivity <- rep(0, lenT)
    Specificity <- rep(0, lenT)
    if (affinity_check) {
        index_1 <- which(coverage$labels_new == 1)
    } else {
        index_1 <- which(coverage$labels_new == 2)
    }
    for (i in 1:lenT) {
        print(paste(i,' out of', lenT))
        cutoff <- T_seq[i]
        lab_pred <- getLabel(scores, cutoff, pval = pval)
        index_1_pred <- which(lab_pred == 1)
        index_0_pred <- which(lab_pred == 0)
        ## Compute TP
        TP <- length( intersect(index_1, index_1_pred) )
        FP <- length(index_1_pred) - TP
        FN <- length(index_1) - TP
        TN <- length(index_0) - FP
        Precision[i] = TP / (TP + FP)
        Sensitivity[i] = TP / (TP + FN)
        Specificity[i] = TN / (FP + TN)
    }
    PC <- data.frame(Precision, Sensitivity, Specificity)
    colnames(PC) <- c('Precision', 'Sensitivity', 'Specificity')
    return(PC)
}

## 02.03.2016: use this version for PR curve in folder 01.23.2016
evalPC_sep2 <- function(Pvals, lenT = 100) {
    scores <- Pvals
    idNA <- which(is.na(scores)==TRUE)
    if (length(idNA) > 0) {
        ## label NA = 100
        scores[idNA] <- 100
    }
    min <- min(scores, na.rm = TRUE)
    max <- max(scores, na.rm = TRUE)
    T_seq <- seq(min, max, length.out = lenT)
    T_seq <- c(-Inf, T_seq, Inf)
    lenT = length(T_seq)
    Precision <- rep(0, lenT)
    Sensitivity <- rep(0, lenT)
    Specificity <- rep(0, lenT)
    eFDR <- rep(0, lenT)
    for (i in 1:lenT) {
        print(paste(i,' out of', lenT))
        cutoff <- T_seq[i]
        lab_pred <- getLabel(scores, cutoff, pval = TRUE)
        index_1_pred <- which(lab_pred == 1)
        index_0_pred <- which(lab_pred == 0)
        ## Compute TP
        TP <- length( intersect(9801:9900, index_1_pred) )
        FP <- length(index_1_pred) - TP
        FN <- length(9801:9900) - TP
        TN <- length(1:9800) - FP
        Precision[i] = TP / (TP + FP)
        Sensitivity[i] = TP / (TP + FN)
        Specificity[i] = TN / (FP + TN)
        eFDR[i] = FP / (FP + TP)
    }
    PC <- data.frame(Precision, Sensitivity, Specificity, eFDR)
    colnames(PC) <- c('Precision', 'Sensitivity', 'Specificity', 'eFDR')
    return(PC)
}

### Moving Average
# x: the vector
# n: the number of samples
# centered: if FALSE, then average current sample and previous (n-1) samples
#           if TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.)
movingAverage <- function(x, n=1, centered=FALSE) {

    if (centered) {
        before <- floor  ((n-1)/2)
        after  <- ceiling((n-1)/2)
    } else {
        before <- n-1
        after  <- 0
    }

    # Track the sum and count of number of non-NA items
    s     <- rep(0, length(x))
    count <- rep(0, length(x))

    # Add the centered data
    new <- x
    # Add to count list wherever there isn't a
    count <- count + !is.na(new)
    # Now replace NA_s with 0_s and add to total
    new[is.na(new)] <- 0
    s <- s + new

    # Add the data from before
    i <- 1
    while (i <= before) {
        # This is the vector with offset values to add
        new   <- c(rep(NA, i), x[1:(length(x)-i)])

        count <- count + !is.na(new)
        new[is.na(new)] <- 0
        s <- s + new

        i <- i+1
    }

    # Add the data from after
    i <- 1
    while (i <= after) {
        # This is the vector with offset values to add
        new   <- c(x[(i+1):length(x)], rep(NA, i))

        count <- count + !is.na(new)
        new[is.na(new)] <- 0
        s <- s + new

        i <- i+1
    }

    # return sum divided by count
    s/count
}

##' <description>
##' bindDf(): return an approriate data frame to getPlot()
##' <details>
##' @title
##' @param data: data frame from plot_profiles() (package Segvis)
##' @return return an approriate data frame to getPlot()
##' @author Duy Nguyen on February 27, 2016
bindDf <- function(data) {
    df <- data.frame(matrix(ncol = 0, nrow = dim(data)[1] ))
    colx <- data[condition==1, x]
    df <- cbind(df,colx)
    for (i in 1:9) {
        coly <- data[condition==i, y]
        df <- cbind(df,coly)
    }
    names(df) <- c('Genomic_coordinates', 'Rep1.1', 'Rep1.2', 'Rep2' ,'Rep3', 'Rep4',
                   'Input_Rep1', 'Input_Rep2', 'Input_Rep3' ,'Input_Rep4')
    rownames(df) <- NULL
    return(df)
}

##' <description>
##' getPlot(): ggplot2 of samples
##' <details>
##' @title
##' @param df: a data frame from running bindDf()
##' @param condition: see Segvig's plot_profiles()
##' @param chr: chromosome number
##' @param combined: TRUE if plotting averages; FALSE if plotting separate samples
##' @return return an approriate data frame to getPlot()
##' @author Duy Nguyen on February 27, 2016
getPlot <- function(df, condition, chr, combined = FALSE, sizeChIP = 1.5) {

    if (!combined) {
        plot <-  ggplot( df, aes(Genomic_coordinates)) + geom_line(aes(y =df[,2] , colour ='Rep1.1')) +
            geom_line(aes(y =df[,3] , colour ='Rep1.2'), size = sizeChIP) +
            geom_line(aes(y =df[,4] , colour ='Rep2'), size = sizeChIP) +
            geom_line(aes(y =df[,5] , colour ='Rep3'), size = sizeChIP) +
            geom_line(aes(y =df[,6] , colour ='Rep4'), size = sizeChIP) +

            geom_line(aes(y =df[,7] , colour ='Rep1.1'), linetype="dotted", size = 1) +
            geom_line(aes(y =df[,8] , colour ='Rep2'), linetype="dotted", size = 1) +
            geom_line(aes(y =df[,9] , colour ='Rep3'), linetype="dotted", size = 1) +
            geom_line(aes(y =df[,10] , colour ='Rep4'), linetype="dotted", size = 1) +
            labs(title = paste(condition, ', ' , chr , ': ',df[,1][1], '-', df[,1][dim(df)[1]]
                               , ', width:', df[,1][dim(df)[1]] - df[,1][1] + 1  )) +
            xlab("Genomic Coordinates") + ylab("Read Counts") + theme(legend.title=element_blank(), legend.position="top")

        print(plot)
        return(plot)
    }
    else {
        df_combined <- data.frame(matrix(ncol = 3, nrow = dim(df)[1] ))
        names(df_combined) <- c("Genomic_coordinates", "mChip", "mInput")
        df_combined[, 1] <- df[,1]
        df_combined[, 2] <- apply(df[,3:6], MARGIN = 1, FUN = mean) # ChIP
        df_combined[, 3] <- apply(df[,7:10], MARGIN = 1, FUN = mean) # input
        df <- df_combined
        p1 <-  ggplot( df, aes(Genomic_coordinates)) +
            geom_line(aes(y =df[,2] , colour ='minusChIP'), size = sizeChIP) +
            geom_line(aes(y =df[,3] , colour ='minusChIP'), linetype="dotted", size = 1) +
            labs(title = paste(condition, ', ' , chr , ': ',df[,1][1], '-', df[,1][dim(df)[1]]
                               , ', width:', df[,1][dim(df)[1]] - df[,1][1] + 1  )) +
            xlab("Genomic Coordinates") + ylab("Read Counts") + theme(legend.title=element_blank(), legend.position="top")
        print(p1)
        return(p1)
    }

}
