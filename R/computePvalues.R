.computePvalues <- function(object, quant, poolQuant, movAve) {
    if (object@nSamples == 4) {
        print(paste("Computing p-values for sample size n = ", object@nSamples), sep = "")
        ## Name columns:
        # colnames(object@W1) <- c('ab vs cd', 'ac vs bd', 'ad vs bc' )
        colnames(object@W2) <- c('AB vs CD', 'AC vs BD', 'AD vs BC' )
        minLabs <- c('ab','ac','ad','bc','bd','cd')
        plusLabs <- c('AB','AC','AD','BC','BD','CD')
        plusLabs <- toupper(minLabs)
        ## Create colnames for matrix B
        labs <- ''
        for (mlab in minLabs) {
            for (plab in plusLabs) {
                labs <- c(labs, paste(mlab, 'vs', plab))
            }
        }
        labs <- labs[-c(1)]
        ### Get dictinary {a:1, b:2, c:3, d:4}
        getDict <- function(s) {
            switch(s, a = 1, b = 2, c = 3, d = 4,
                   A = 5, B = 6, C = 7, D = 8)
        }
        ### Main ### STOP
        total <-  nPeaks <- length(object@coverage)
        # number of columns between B: 6*6 columns
        Between_cols <- 36
        Within <- cbind(object@W1, object@W2)
        AvsB <- matrix(NA, nrow = total, ncol = Between_cols)
        colnames(AvsB) <- labs
        AvsB <- cbind(AvsB, Within)
        ncomps <- ncol(AvsB)
        p <- matrix(NA, nrow = nrow(AvsB), ncol = ncol(AvsB))
        colnames(p) <- colnames(AvsB)
        H0 <- as.vector(Within)
        if (ncol(object@Ns)>1){
            temp <- colnames(Within)
            sampleNames <- unlist(strsplit(temp,' vs '))
            ## For each interval, take average of all ChIP samples "ac" "bd" "ad" "bc" "AC" "BD" "AD" "BC"; 'ac' means take average between a and c for that interval and so on
            H0.Ns <- rowMeans(object@Ns[,sampleNames])
            H0.Ns <-rep(H0.Ns,ncol(Within))
        } else {
            H0.Ns <-rep(object@Ns,ncol(Within))
        }

        if (max(object@dN)<max(object@Ns)){
            object@dN <- c(object@dN, max(object@Ns))
        }

        for (i in 1:length(object@dN)){
            print(paste('+++Bin = ', i))
            if (i==1){
                if (ncol(object@Ns)==1){
                    idx <- which(object@Ns<=object@dN[[i]])
                }
                H0.idx <- which(H0.Ns<=object@dN[[i]])
            } else {
                if (ncol(object@Ns)==1){
                    idx <- which(object@Ns > object@dN[[i-1]] & object@Ns <= object@dN[[i]])
                }
                H0.idx <- which(H0.Ns > object@dN[[i-1]] & H0.Ns <= object@dN[[i]])
            }
            ##################################################################
            ###
            ### Testing between: j for column, idx for rows of matrix p
            ###
            ##################################################################
            for (j in 1:ncomps){
                print(paste('Testing', j ,colnames(p)[j]))
                # eval idx = bin (test j)
                if (ncol(object@Ns)>1){
                    comp <- colnames(AvsB)[j]
                    sampleNames <- unlist(strsplit(comp,' vs '))
                    CompNs <- rowMeans(object@Ns[,sampleNames])
                    if (i==1){
                        idx <- which(CompNs<=object@dN[[i]])
                    } else {
                        idx <- which(CompNs > object@dN[[i-1]] & CompNs<=object@dN[[i]])
                    }
                    if (length(idx)==0) next
                }
                # print(paste('+++idx = ', idx))
                ## Get label indices
                m1 <- getDict(substr(sampleNames[1], 1, 1))
                m2 <- getDict(substr(sampleNames[1], 2, 2))
                p1 <- getDict(substr(sampleNames[2], 1, 1))
                p2 <- getDict(substr(sampleNames[2], 2, 2))
                minusVar_idx <- plusVar_idx <- list()
                minGlobal <- Inf
                # compute between-TAN
                sitesUnused <- c()
                for (ii in 1:length(idx)) {
                    # print(paste("ii = ", ii))
                    site <- idx[ii]
                    # vector to stored sites in bin( test j ) whose lengths of var <= Global_lower
                    get_m1 <- object@coverage[[site]][m1,]
                    get_m2 <- object@coverage[[site]][m2,]
                    get_p1 <- object@coverage[[site]][p1,]
                    get_p2 <- object@coverage[[site]][p2,]
                    X <- rbind( get_m1, get_m2)
                    Y <- rbind( get_p1, get_p2)
                    if ( dim(X)[2] < object@s.size ) {
                        # test <-  AN.test(X, Y, na.rm=TRUE)
                        test <- tan::AN_test(X, Y, na_rm = TRUE, pool = FALSE, poolVarX = NA, poolVarY = NA)
                    }
                    else {
                        design <- object@Designs[site, ]
                        # test <- AN.test(X[, design], Y[, design], na.rm=TRUE, candK = candK)
                        test <- AN_test(X[, design], Y[, design], na_rm = TRUE, pool = FALSE, poolVarX = NA, poolVarY = NA)
                    }
                    minIndex <- min(length(test$varX), length(test$varY))
                    ## check minIndex > Global_lower (lower bound for pooled var vector of each bins)
                    if (minIndex > Global_lower) {
                        if (minGlobal > minIndex) {
                            minGlobal <- minIndex
                        }
                        minusVar_idx[[ii]] <- test$varX
                        plusVar_idx[[ii]] <- test$varY
                    } else { # case minIndex <= Global_lower
                        # print(paste(" +++ minIndex <= Global_lower: ", minIndex, sep = ""))
                        minusVar_idx[[ii]] <- c()
                        plusVar_idx[[ii]] <- c()
                        sitesUnused <- c(sitesUnused, ii)
                    }
                } # end of for (ii in 1:length(idx))
                print(paste(" +++ total sites Unused: ", length(sitesUnused), sep = ""))
                print(paste(" +++ minGlobal = : ", minGlobal, sep = ""))
                ### The following only works if 0 < minGlobal < Inf
                ##################################################################
                ### Compute pooled variances for between-TAN in bin(test j) or idx
                ##################################################################
                if (minGlobal != Inf) {
                    matVar_minus <- matVar_plus <- matrix(NA, nrow = length(idx), ncol = minGlobal)
                    # check if sites whose len of var > Global_lower: if TRUE --> use pooled, o.w. use unpool
                    for (ii in 1:length(minusVar_idx)) {
                        if (length(minusVar_idx[[ii]]) > 0 ) {
                            matVar_minus[ii,] <- minusVar_idx[[ii]][1:minGlobal]
                            matVar_plus[ii,] <- plusVar_idx[[ii]][1:minGlobal]
                        }
                    }
                    if (length(sitesUnused) > 0 ) {
                        varMinus <- apply(matVar_minus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varMinus) >= movAve ) {
                            varMinus <- tan::movingAverage(varMinus, movAve)
                        }
                        varPlus <- apply(matVar_plus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varPlus) >= movAve ) {
                            varPlus <- tan::movingAverage(varPlus, movAve)
                        }
                    } else { # case sitesUnused = c(): use unpooled
                        varMinus <- apply(matVar_minus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varMinus) >= movAve ) {
                            varMinus <- tan::movingAverage(varMinus, movAve)
                        }
                        varPlus <- apply(matVar_plus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varPlus) >= movAve ) {
                            varPlus <- tan::movingAverage(varPlus, movAve)
                        }
                    }
                } else {
                    print("minGlobal = Inf")
                }
                #########################################################
                ### Compute Adaptive Neyman tests with pooled variance
                #########################################################
                ANT <- c()
                if (minGlobal != Inf) {
                    for (ii in 1:length(idx)) {
                        site <- idx[ii]
                        get_m1 <- object@coverage[[site]][m1,]
                        get_m2 <- object@coverage[[site]][m2,]
                        get_p1 <- object@coverage[[site]][p1,]
                        get_p2 <- object@coverage[[site]][p2,]
                        X <- rbind( get_m1, get_m2)
                        Y <- rbind( get_p1, get_p2)
                        if ( dim(X)[2] < object@s.size ) {
                            #######################################################
                            ##  check if site is in sitesUnused --> use unpooled var
                            ##  if site is NOT in sitesUnused --> use pooled var
                            #######################################################
                            if (ii %in% sitesUnused ) {
                                ANT[ii] <- tan::AN_test(X, Y, na_rm=TRUE, pool= FALSE, poolVarX = NA,
                                                   poolVarY = NA)$statistic
                            } else {
                                clen <- 1:length(varMinus)
                                ANT[ii] <- tan::AN_test(X[, clen], Y[, clen], na_rm=TRUE, pool= TRUE, poolVarX = varMinus,
                                                   poolVarY = varPlus)$statistic
                            }
                        } # end of if (dim(X)[2] < s.size)
                        else {
                            if (ii %in% sitesUnused ) {
                                design <- object@Designs[site, ]
                                ANT[ii] <- tan::AN_test(X[,design], Y[,design], na_rm = TRUE, pool = FALSE, poolVarX = NA,
                                                   poolVarY = NA)$statistic
                            } else {
                                design <- object@Designs[site, ]
                                clen <- 1:length(varMinus)
                                ANT[ii] <- tan::AN_test(X[,design[clen]], Y[,design[clen]], na_rm = TRUE, pool = TRUE, poolVarX = varMinus,
                                                   poolVarY = varPlus)$statistic
                            }
                        }

                    } # end of for (ii in 1:length(idx))
                    if (length(H0.idx) > 0 ) {
                        p[idx,j] <- sapply(ANT,function(tanTest) {
                            length(which(H0[H0.idx]>=tanTest))/length(H0.idx)
                            })
                    }

                } # end if (minGlobal != Inf)
            } # end of for (j in 1:ncomps)
        } # end of for (i in 1:length(dN))
        # impute missing values
        if (any(is.na(p))){
            message(length(which(is.na(p))),' NAs found  (of ', length(p),')')
            p[is.na(p)]=min(p[!is.na(p)])
        }
        Pc <- apply(p[,1:Between_cols], 1, function(x) quantile(x, probs = quant, na.rm =TRUE))
        fdr <- p.adjust(as.vector(p),method='BH')
        FDR <- matrix(0,nrow=nrow(p),ncol=ncol(p)+1)

        for (i in 1:ncol(p)){
            FDR[,i] <- fdr[(i-1)*nrow(p)+(1:nrow(p))]
        }
        if (any(is.na(FDR))){
            message(length(which(is.na(FDR))),' NAs found  (of ', length(FDR),')')
            FDR[is.na(FDR)]=min(FDR[!is.na(FDR)])
        }
        ### CHANGE THIS TO THE RIGHT QUANTILE OF P-VALUES
        FDR[,i+1] <- apply(FDR[,1:Between_cols], 1, function(x) quantile(x, probs = quant, na.rm = TRUE))
        p <- cbind(p,Pc)
        colnames(p)[i+1] <- 'combined'
        colnames(FDR) <- colnames(p)
        P <- list(p,FDR)
        names(P) <- c('pval','FDR')
        object@PvalList <- P
    } # end of if (n=4)
    object
}

setMethod("computePvalues", signature("tanDb"), .computePvalues)
