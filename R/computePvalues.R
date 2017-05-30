evaluate_pvalues <- function(coverage, nSamples, wSites, use_cpp = TRUE,
                            s.size, Designs, sitesUnused, W1, W2,
                            Ns, dN, ignore_sitesUnused, quant,
                            poolQuant, movAve, na_impute) {
    ## Get dictionary {a:1, b:2, c:3, d:4,...}
    getDict <- function(s) {
        if (s == tolower(s)) {
            return(match(s, letters))
        }
        else {
            return(match(tolower(s), letters) + nSamples)
        }
    }

    ## create pair sample for within labels
    ## n4: "ab vs cd" "ac vs bd" "ad vs bc" "AB vs CD" "AC vs BD" "AD vs BC"
    create_labels <- function(nSamples) {
        if (nSamples > 3) {
            # create colnames for Within matrix
            matLab <- combn(letters[1:nSamples], 2)
            labs <- apply(matLab, 2, FUN = function(x) paste(x[1], x[2], sep='') )
            temp <- c()
            check_labs <- rep(FALSE, length(labs)); names(check_labs) <- labs
            for (l in labs) {
                l_exclude <- setdiff(letters[1:nSamples], c( substr(l, 1, 1), substr(l, 2, 2)))
                mat_exclude <-  apply(combn(l_exclude, 2), 2, FUN = function(x) paste(x[1], x[2], sep='') )
                for (i in 1:length(mat_exclude)) {
                    if (!check_labs[mat_exclude[i]]) {
                        temp <- c(temp, paste(l, mat_exclude[i], sep = ''))
                    }
                }
                check_labs[l] <- TRUE
            }
            temp_plus <- toupper(temp)
            temp <- unlist(lapply(temp, function(s) paste(substr(s, 1, 2), substr(s, 3, 4), sep = ' vs ' )))
            temp <- c(temp,unlist(lapply(temp_plus, function(s) paste(substr(s, 1, 2), substr(s, 3, 4), sep = ' vs ' ))))
            return(list('withinLabel' = temp, 'withinPair' = labs))
        }
        else if (nSamples == 3) {
            temp <- c("ab vs ac", "ab vs bc", "ac vs bc",
                     "AB vs AC", "AB vs BC", "AC vs BC")
            labs <- c('ab', 'ac', 'bc')
            return(list('withinLabel' = temp, 'withinPair' = labs))
        }
    }
    ## get index {1,2,...} for within's labels
    create_indexList <- function(nSamples) {
        indexList <- list()
        if (nSamples > 3) {
            withinLabel <- create_labels(nSamples)[[1]]
            numTests <- floor(length(withinLabel) / 2)
            minusLabel <- withinLabel[1:numTests]
            # sampleNames <- unlist(strsplit(minusLabel, ' vs '))
            for (i in 1:numTests) {
                pairSample <- unlist(strsplit(minusLabel[i], ' vs '))
                s <- c(substr(pairSample[1], 1, 1), substr(pairSample[1], 2, 2),
                      substr(pairSample[2], 1, 1), substr(pairSample[2], 2, 2))
                s <- sapply(s, getDict)
                names(s) <- NULL
                indexList[[i]] <- s
            }
            return(indexList)
        }
        else if (nSamples == 3) {
            indexList[[1]] <- c(1,2,1,3) # ab vs ac
            indexList[[2]] <- c(1,2,2,3) # ab vs bc
            indexList[[3]] <- c(1,3,2,3) # ac vs bc
            return(indexList)
        }
    }
    ### MAIN ###
    # n=4: "ab" "ac" "ad" "bc" "bd" "cd"
    minLabs <- create_labels(nSamples)[['withinPair']]
    plusLabs <- toupper(minLabs)
    ## Create colnames for matrix between B
    labs <- ''
    for (mlab in minLabs) {
        for (plab in plusLabs) {
            labs <- c(labs, paste(mlab, 'vs', plab))
        }
    }
    labs <- labs[-c(1)]

    total <- length(coverage)
    #n=4: number of columns between B equals 6*6 columns
    Between_cols <- length(labs)
    sitesUnused_Within <- sitesUnused
    Within <- cbind(W1, W2)
    AvsB <- matrix(NA, nrow = total, ncol = Between_cols)
    colnames(AvsB) <- labs
    AvsB <- cbind(AvsB, Within)
    ncomps <- ncol(AvsB)
    p <- matrix(NA, nrow = nrow(AvsB), ncol = ncol(AvsB))
    colnames(p) <- colnames(AvsB)
    H0 <- as.vector(Within)
    if (ncol(Ns) > 1) {
        temp <- colnames(Within)
        sampleNames <- unlist(strsplit(temp,' vs '))
        ## For each interval, take average of all ChIP samples "ac" "bd" "ad" "bc" "AC" "BD" "AD" "BC";
        # 'ac' means take average between a and c for that interval and so on
        H0.Ns <- rowMeans(Ns[, sampleNames]) # index from 1:total
        H0.Ns <- rep(H0.Ns, ncol(Within))
    } else {
        H0.Ns <- rep(Ns, ncol(Within))
    }

    if (max(dN) < max(Ns)){
        dN <- c(dN, max(Ns))
    }
    # evaluate H0.idx
    for (i in 1:length(dN)) {
        print(paste('+++Bin = ', i))
        if (i == 1) {
            if (ncol(Ns) == 1){
                idx <- which(Ns <= dN[[i]])
            }
            H0.idx <- which(H0.Ns <= dN[[i]])
        } else {
            if (ncol(Ns) == 1){
                idx <- which(Ns > dN[[i-1]] & Ns <= dN[[i]])
            }
            H0.idx <- which(H0.Ns > dN[[i-1]] & H0.Ns <= dN[[i]])
        }
        ##################################################################
        ###
        ### Testing between: j for column, idx for rows of matrix p
        ###
        ##################################################################
        for (j in 1:ncomps){
            print(paste('Testing', j, colnames(p)[j]))
            # eval idx = bin (test j)
            if (ncol(Ns) > 1){
                comp <- colnames(AvsB)[j]
                sampleNames <- unlist(strsplit(comp,' vs '))
                CompNs <- rowMeans(Ns[, sampleNames])
                if (i == 1){
                    idx <- which(CompNs <= dN[[i]])
                } else {
                    idx <- which(CompNs > dN[[i-1]] & CompNs <= dN[[i]])
                }
                if (length(idx) == 0) next
            }
            ## Get label indices
            m1 <- getDict(substr(sampleNames[1], 1, 1))
            m2 <- getDict(substr(sampleNames[1], 2, 2))
            p1 <- getDict(substr(sampleNames[2], 1, 1))
            p2 <- getDict(substr(sampleNames[2], 2, 2))
            minusVar_idx <- plusVar_idx <- list()
            minGlobal <- Inf
            ## compute between-TAN variance
            sitesUnused <- c()
            for (ii in 1:length(idx)) {
                # print(paste("ii = ", ii))
                site <- idx[ii]
                ## vector to store sites in bin( test j )
                ## whose lengths of var <= Global_lower
                get_m1 <- coverage[[site]][m1,]
                get_m2 <- coverage[[site]][m2,]
                get_p1 <- coverage[[site]][p1,]
                get_p2 <- coverage[[site]][p2,]
                X <- rbind(get_m1, get_m2)
                Y <- rbind(get_p1, get_p2)
                if ( dim(X)[2] < s.size ) {
                    if (use_cpp) {
                        test <- tan::AN_test(X, Y, na_rm = TRUE,
                                            pool = FALSE, poolVarX = NA, poolVarY = NA)
                    }
                    else {
                        test <-  tan::AN.test(X, Y, na.rm=TRUE)
                    }
                }
                else {
                    design <- Designs[site, ]
                    if (use_cpp) {
                        test <- tan::AN_test(X[, design], Y[, design], na_rm = TRUE,
                                            pool = FALSE, poolVarX = NA, poolVarY = NA)
                    }
                    else {
                        test <- tan::AN.test(X[, design], Y[, design], na.rm=TRUE)
                    }
                }
                minIndex <- min(length(test$varX), length(test$varY))
                ## check minIndex > Global_lower (lower bound for pooled var vector of each bins)
                if (minIndex > Global_lower) {
                    if (minGlobal > minIndex) {
                        minGlobal <- minIndex
                    }
                    minusVar_idx[[ii]] <- test$varX
                    plusVar_idx[[ii]] <- test$varY
                }
                ## store all unused sites (minIndex <= Global_lower):
                ## minIndex = 0 -> var empty due to flat peak,
                ## or many repeated counts, or peakLength too small.
                else {
                    # print(paste(" +++ minIndex <= Global_lower: ", minIndex, sep = ""))
                    minusVar_idx[[ii]] <- c()
                    plusVar_idx[[ii]] <- c()
                    sitesUnused <- c(sitesUnused, ii)
                }
            } # end of for (ii in 1:length(idx))
            print(paste(" +++ total sites Unused: ", length(sitesUnused), sep = ""))
            print(paste(" +++ minGlobal = : ", minGlobal, sep = ""))
            ##################################################################
            ### Compute pooled variances for between-TAN in bin(test j) or idx
            ##################################################################
            ### The following only works if 0 < minGlobal < Inf
            if (minGlobal != Inf) {
                matVar_minus <- matVar_plus <- matrix(NA, nrow = length(idx), ncol = minGlobal)
                # check if sites whose len of var > Global_lower: if TRUE --> use pooled, o.w. use unpool
                for (ii in 1:length(minusVar_idx)) {
                    if (length(minusVar_idx[[ii]]) > 0 ) {
                        matVar_minus[ii, ] <- minusVar_idx[[ii]][1:minGlobal]
                        matVar_plus[ii, ] <- plusVar_idx[[ii]][1:minGlobal]
                    }
                }
                # TODO: COMBINE these two cases length(sitesUnused) > 0 or else
                if (length(sitesUnused) > 0 ) {
                    varMinus <- apply(matVar_minus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                    if ( length(varMinus) >= movAve ) {
                        varMinus <- tan::movingAverage(varMinus, movAve)
                    }
                    varPlus <- apply(matVar_plus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                    if ( length(varPlus) >= movAve ) {
                        varPlus <- tan::movingAverage(varPlus, movAve)
                    }
                }
                # case sitesUnused = c(): use unpooled.
                else {
                    varMinus <- apply(matVar_minus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                    if ( length(varMinus) >= movAve ) {
                        varMinus <- tan::movingAverage(varMinus, movAve)
                    }
                    varPlus <- apply(matVar_plus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                    if ( length(varPlus) >= movAve ) {
                        varPlus <- tan::movingAverage(varPlus, movAve)
                    }
                }
            }
            # case minGlobal = Inf
            else {
                message("minGlobal = Inf. Please consider larger bin sizes
                   and/or filtering sites with zero variances and small peakLengths!")
            }
            # add siteUnused_between to slot:
            sitesUnused <- unique(c(sitesUnused, idx[sitesUnused]))
            #########################################################
            ### Compute Adaptive Neyman tests with pooled variance
            #########################################################
            ANT <- c()
            if (minGlobal != Inf) {
                for (ii in 1:length(idx)) {
                    site <- idx[ii]
                    get_m1 <- coverage[[site]][m1,]
                    get_m2 <- coverage[[site]][m2,]
                    get_p1 <- coverage[[site]][p1,]
                    get_p2 <- coverage[[site]][p2,]
                    X <- rbind( get_m1, get_m2)
                    Y <- rbind( get_p1, get_p2)
                    if ( dim(X)[2] < s.size ) {
                        #######################################################
                        ##  check if site is in sitesUnused --> use unpooled var
                        ##  if site is NOT in sitesUnused --> use pooled var
                        #######################################################
                        if (ii %in% sitesUnused ) {
                            if (use_cpp) {
                                if (ignore_sitesUnused == FALSE) {
                                    ANT[ii] <- tan::AN_test(X, Y, na_rm=TRUE,
                                                           pool= FALSE, poolVarX = NA, poolVarY = NA)$statistic
                                    #@: return NA when ignore_SitesUnused = TRUE ?
                                }
                                else {
                                    ANT[ii] <- NA
                                }
                            }
                            else {
                                if (ignore_sitesUnused == FALSE) {
                                    ANT[ii] <- tan::AN.test(X, Y, na.rm=TRUE)$statistic
                                }
                                else {
                                    ANT[ii] <- NA
                                }
                            }
                        }
                        else {
                            clen <- 1:length(varMinus)
                            if (use_cpp) {
                                ANT[ii] <- tan::AN_test(X[, clen], Y[, clen], na_rm=TRUE,
                                                       pool= TRUE, poolVarX = varMinus,
                                                       poolVarY = varPlus)$statistic
                            }
                            else {
                                ANT[ii] <- tan::AN.test(X[, clen], Y[, clen], na.rm=TRUE,
                                                       pool= TRUE, poolVarX = varMinus,
                                                       poolVarY = varPlus)$statistic
                            }
                        }
                    } # end of if (dim(X)[2] < s.size)
                    else {
                        if (ii %in% sitesUnused ) {
                            design <- Designs[site, ]
                            if (use_cpp) {
                                if (ignore_sitesUnused == FALSE) {
                                    ANT[ii] <- tan::AN_test(X[,design], Y[,design], na_rm = TRUE, pool = FALSE, poolVarX = NA,
                                                           poolVarY = NA)$statistic
                                }
                                else {
                                    ANT[ii] <- NA
                                }
                            }
                            else {
                                if (ignore_sitesUnused == FALSE) {
                                    ANT[ii] <- tan::AN.test(X[,design], Y[,design], na.rm = TRUE)$statistic
                                }
                                else {
                                    ANT[ii] <- NA
                                }
                            }
                        } else {
                            design <- Designs[site, ]
                            clen <- 1:length(varMinus)
                            if (use_cpp) {
                                ANT[ii] <- tan::AN_test(X[,design[clen]], Y[,design[clen]], na_rm = TRUE,
                                                       pool = TRUE, poolVarX = varMinus,
                                                       poolVarY = varPlus)$statistic
                            }
                            else {
                                ANT[ii] <- tan::AN.test(X[,design[clen]], Y[,design[clen]], na.rm = TRUE,
                                                       pool = TRUE, poolVarX = varMinus,
                                                       poolVarY = varPlus)$statistic
                            }
                        }
                    }
                } # end of for (ii in 1:length(idx))
                if (length(H0.idx) > 0 ) {
                    if (ignore_sitesUnused == FALSE) {
                        p[idx, j] <- sapply(ANT, function(tanTest) {
                            length(which(H0[H0.idx] >= tanTest )) / length(H0.idx)
                        })
                    }
                    # new: 03/24/2017
                    else {
                        # ignore unused sites from testing between
                        na_indices <- which(is.na(ANT) == TRUE)
                        # ignore unused sites from testing within
                        H0.idx_ <- setdiff(H0.idx, sitesUnused_Within)
                        H0.idx_ <- rep(H0.idx_, ncol(Within))
                        p[idx, j] <- sapply(ANT, function(tanTest) {
                            length(which(H0[H0.idx_] >= tanTest )) / length(H0.idx_)
                        })
                        p[idx[na_indices], j] <- NA
                    }
                }
            } # end if (minGlobal != Inf)
        } # end of for (j in 1:ncomps)
    } # end of for (i in 1:length(dN))
    # impute missing values
    if (any(is.na(p))) {
        message(length(which(is.na(p))),' NAs found  (of ', length(p),')')
        if (na_impute) {
            p[is.na(p)] = min(p[!is.na(p)])
        }
    }
    Pc <- apply(p[, 1:Between_cols], 1, function(x) quantile(x, probs = quant, na.rm =TRUE))
    fdr <- p.adjust(as.vector(p),method='BH')
    FDR <- matrix(0,nrow = nrow(p), ncol = ncol(p)+1)
    for (i in 1:ncol(p)){
        FDR[, i] <- fdr[ (i-1) * nrow(p) + (1:nrow(p)) ]
    }
    if (any(is.na(FDR))){
        message(length(which(is.na(FDR))),' NAs found  (of ', length(FDR),')')
        if (na_impute) {
            FDR[is.na(FDR)]=min(FDR[!is.na(FDR)])
        }
    }
    FDR[, i+1] <- apply(FDR[,1:Between_cols], 1, function(x) quantile(x, probs = quant, na.rm = TRUE))
    p <- cbind(p,Pc)
    colnames(p)[i+1] <- 'combined'
    colnames(FDR) <- colnames(p)
    P <- list(p,FDR)
    names(P) <- c('pval','FDR')
    return(list('PvalList' = P, 'sitesUnused' = sitesUnused))
}

.computePvalues <- function(object, quant, poolQuant, movAve, Global_lower, ...) {
    if (object@nSamples > 3) {
        print(paste("Computing p-values for sample size n = ", object@nSamples), sep = "")

        resultList <- evaluate_pvalues(coverage = object@coverage, nSamples = object@nSamples,
                              wSites = object@wSites, use_cpp = use_cpp, s.size = object@s.size,
                              Designs = object@Designs, sitesUnused = object@sitesUnused,
                              W1 = object@W1, W2 = object@W2, Ns = object@Ns, dN = object@dN,
                              ignore_sitesUnused = ignore_sitesUnused, quant = quant,
                              poolQuant = poolQuant, movAve = movAve, na_impute = na_impute)

        ## ## ## Name columns:
        ## ## colnames(object@W1) <- c('ab vs cd', 'ac vs bd', 'ad vs bc' )
        ## ## colnames(object@W2) <- c('AB vs CD', 'AC vs BD', 'AD vs BC' )
        ## minLabs <- c('ab','ac','ad','bc','bd','cd')
        ## plusLabs <- c('AB','AC','AD','BC','BD','CD')
        ## plusLabs <- toupper(minLabs)
        ## ## Create colnames for matrix B
        ## labs <- ''
        ## for (mlab in minLabs) {
        ##     for (plab in plusLabs) {
        ##         labs <- c(labs, paste(mlab, 'vs', plab))
        ##     }
        ## }
        ## labs <- labs[-c(1)]
        ## ### Get dictionary {a:1, b:2, c:3, d:4}
        ## getDict <- function(s) {
        ##     switch(s, a = 1, b = 2, c = 3, d = 4,
        ##            A = 5, B = 6, C = 7, D = 8)
        ## }
        ## ### Main ###
        ## total <-  nPeaks <- length(object@coverage)
        ## # number of columns between B: 6*6 columns
        ## Between_cols <- 36
        ## sitesUnused_Within <- object@sitesUnused # @
        ## Within <- cbind(object@W1, object@W2)
        ## AvsB <- matrix(NA, nrow = total, ncol = Between_cols)
        ## colnames(AvsB) <- labs
        ## AvsB <- cbind(AvsB, Within)
        ## ncomps <- ncol(AvsB)
        ## p <- matrix(NA, nrow = nrow(AvsB), ncol = ncol(AvsB))
        ## colnames(p) <- colnames(AvsB)
        ## H0 <- as.vector(Within)
        ## if (ncol(object@Ns) > 1) {
        ##     temp <- colnames(Within)
        ##     sampleNames <- unlist(strsplit(temp,' vs '))
        ##     ## For each interval, take average of all ChIP samples "ac" "bd" "ad" "bc" "AC" "BD" "AD" "BC";
        ##     # 'ac' means take average between a and c for that interval and so on
        ##     H0.Ns <- rowMeans(object@Ns[, sampleNames]) # index from 1:total
        ##     H0.Ns <- rep(H0.Ns, ncol(Within))
        ## } else {
        ##     H0.Ns <- rep(object@Ns, ncol(Within))
        ## }

        ## if (max(object@dN) < max(object@Ns)){
        ##     object@dN <- c(object@dN, max(object@Ns))
        ## }
        ## # evaluate H0.idx
        ## for (i in 1:length(object@dN)) {
        ##     print(paste('+++Bin = ', i))
        ##     if (i == 1) {
        ##         if (ncol(object@Ns) == 1){
        ##             idx <- which(object@Ns <= object@dN[[i]])
        ##         }
        ##         H0.idx <- which(H0.Ns <= object@dN[[i]])
        ##     } else {
        ##         if (ncol(object@Ns) == 1){
        ##             idx <- which(object@Ns > object@dN[[i-1]] & object@Ns <= object@dN[[i]])
        ##         }
        ##         H0.idx <- which(H0.Ns > object@dN[[i-1]] & H0.Ns <= object@dN[[i]])
        ##     }
        ##     ##################################################################
        ##     ###
        ##     ### Testing between: j for column, idx for rows of matrix p
        ##     ###
        ##     ##################################################################
        ##     for (j in 1:ncomps){
        ##         print(paste('Testing', j, colnames(p)[j]))
        ##         # eval idx = bin (test j)
        ##         if (ncol(object@Ns) > 1){
        ##             comp <- colnames(AvsB)[j]
        ##             sampleNames <- unlist(strsplit(comp,' vs '))
        ##             CompNs <- rowMeans(object@Ns[, sampleNames])
        ##             if (i == 1){
        ##                 idx <- which(CompNs <= object@dN[[i]])
        ##             } else {
        ##                 idx <- which(CompNs > object@dN[[i-1]] & CompNs<=object@dN[[i]])
        ##             }
        ##             if (length(idx) == 0) next
        ##         }
        ##         ## Get label indices
        ##         m1 <- getDict(substr(sampleNames[1], 1, 1))
        ##         m2 <- getDict(substr(sampleNames[1], 2, 2))
        ##         p1 <- getDict(substr(sampleNames[2], 1, 1))
        ##         p2 <- getDict(substr(sampleNames[2], 2, 2))
        ##         minusVar_idx <- plusVar_idx <- list()
        ##         minGlobal <- Inf
        ##         ## compute between-TAN variance
        ##         sitesUnused <- c()
        ##         for (ii in 1:length(idx)) {
        ##             # print(paste("ii = ", ii))
        ##             site <- idx[ii]
        ##             # vector to stored sites in bin( test j ) whose lengths of var <= Global_lower
        ##             get_m1 <- object@coverage[[site]][m1,]
        ##             get_m2 <- object@coverage[[site]][m2,]
        ##             get_p1 <- object@coverage[[site]][p1,]
        ##             get_p2 <- object@coverage[[site]][p2,]
        ##             X <- rbind( get_m1, get_m2)
        ##             Y <- rbind( get_p1, get_p2)
        ##             if ( dim(X)[2] < object@s.size ) {
        ##                 if (use_cpp) {
        ##                     test <- tan::AN_test(X, Y, na_rm = TRUE, pool = FALSE, poolVarX = NA, poolVarY = NA)
        ##                 }
        ##                 else {
        ##                     test <-  tan::AN.test(X, Y, na.rm=TRUE)
        ##                 }
        ##             }
        ##             else {
        ##                 design <- object@Designs[site, ]
        ##                 if (use_cpp) {
        ##                     test <- tan::AN_test(X[, design], Y[, design], na_rm = TRUE, pool = FALSE, poolVarX = NA, poolVarY = NA)
        ##                 }
        ##                 else {
        ##                     test <- tan::AN.test(X[, design], Y[, design], na.rm=TRUE)
        ##                 }
        ##             }
        ##             minIndex <- min(length(test$varX), length(test$varY))
        ##             ## check minIndex > Global_lower (lower bound for pooled var vector of each bins)
        ##             if (minIndex > Global_lower) {
        ##                 if (minGlobal > minIndex) {
        ##                     minGlobal <- minIndex
        ##                 }
        ##                 minusVar_idx[[ii]] <- test$varX
        ##                 plusVar_idx[[ii]] <- test$varY
        ##             }
        ##             # store all unused sites (minIndex <= Global_lower): minIndex = 0 -> var empty due to flat peak,
        ##             # or many repeated counts, or peakLength too small. 03/24/17
        ##             else {
        ##                 # print(paste(" +++ minIndex <= Global_lower: ", minIndex, sep = ""))
        ##                 minusVar_idx[[ii]] <- c()
        ##                 plusVar_idx[[ii]] <- c()
        ##                 sitesUnused <- c(sitesUnused, ii)
        ##             }
        ##         } # end of for (ii in 1:length(idx))
        ##         print(paste(" +++ total sites Unused: ", length(sitesUnused), sep = ""))
        ##         print(paste(" +++ minGlobal = : ", minGlobal, sep = ""))
        ##         ##################################################################
        ##         ### Compute pooled variances for between-TAN in bin(test j) or idx
        ##         ##################################################################
        ##         ### The following only works if 0 < minGlobal < Inf
        ##         if (minGlobal != Inf) {
        ##             matVar_minus <- matVar_plus <- matrix(NA, nrow = length(idx), ncol = minGlobal)
        ##             # check if sites whose len of var > Global_lower: if TRUE --> use pooled, o.w. use unpool
        ##             for (ii in 1:length(minusVar_idx)) {
        ##                 if (length(minusVar_idx[[ii]]) > 0 ) {
        ##                     matVar_minus[ii, ] <- minusVar_idx[[ii]][1:minGlobal]
        ##                     matVar_plus[ii, ] <- plusVar_idx[[ii]][1:minGlobal]
        ##                 }
        ##             }
        ##             # TODO: COMBINE these two cases length(sitesUnused) > 0 or else
        ##             if (length(sitesUnused) > 0 ) {
        ##                 varMinus <- apply(matVar_minus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
        ##                 if ( length(varMinus) >= movAve ) {
        ##                     varMinus <- tan::movingAverage(varMinus, movAve)
        ##                 }
        ##                 varPlus <- apply(matVar_plus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
        ##                 if ( length(varPlus) >= movAve ) {
        ##                     varPlus <- tan::movingAverage(varPlus, movAve)
        ##                 }
        ##             }
        ##             # case sitesUnused = c(): use unpooled.
        ##             else {
        ##                 varMinus <- apply(matVar_minus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
        ##                 if ( length(varMinus) >= movAve ) {
        ##                     varMinus <- tan::movingAverage(varMinus, movAve)
        ##                 }
        ##                 varPlus <- apply(matVar_plus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
        ##                 if ( length(varPlus) >= movAve ) {
        ##                     varPlus <- tan::movingAverage(varPlus, movAve)
        ##                 }
        ##             }
        ##         }
        ##         # case minGlobal = Inf
        ##         else {
        ##             message("minGlobal = Inf. Please consider larger bin sizes and/or filtering sites with zero variances and small peakLengths!")
        ##         }
        ##         # add siteUnused_between to slot:
        ##         object@sitesUnused <- unique(c(object@sitesUnused, idx[sitesUnused] ))
        ##         #########################################################
        ##         ### Compute Adaptive Neyman tests with pooled variance
        ##         #########################################################
        ##         ANT <- c()
        ##         if (minGlobal != Inf) {
        ##             for (ii in 1:length(idx)) {
        ##                 site <- idx[ii]
        ##                 get_m1 <- object@coverage[[site]][m1,]
        ##                 get_m2 <- object@coverage[[site]][m2,]
        ##                 get_p1 <- object@coverage[[site]][p1,]
        ##                 get_p2 <- object@coverage[[site]][p2,]
        ##                 X <- rbind( get_m1, get_m2)
        ##                 Y <- rbind( get_p1, get_p2)
        ##                 if ( dim(X)[2] < object@s.size ) {
        ##                     #######################################################
        ##                     ##  check if site is in sitesUnused --> use unpooled var
        ##                     ##  if site is NOT in sitesUnused --> use pooled var
        ##                     #######################################################
        ##                     if (ii %in% sitesUnused ) {
        ##                         if (use_cpp) {
        ##                             if (ignore_sitesUnused == FALSE) {
        ##                                 ANT[ii] <- tan::AN_test(X, Y, na_rm=TRUE, pool= FALSE, poolVarX = NA, poolVarY = NA)$statistic #@: return NA when ignore_SitesUnused = TRUE ???
        ##                             }
        ##                             else {
        ##                                 ANT[ii] <- NA
        ##                             }
        ##                         }
        ##                         else {
        ##                             if (ignore_sitesUnused == FALSE) {
        ##                                 ANT[ii] <- tan::AN.test(X, Y, na.rm=TRUE)$statistic
        ##                             }
        ##                             else {
        ##                                 ANT[ii] <- NA
        ##                             }
        ##                         }
        ##                     }
        ##                     else {
        ##                         clen <- 1:length(varMinus)
        ##                         if (use_cpp) {
        ##                             ANT[ii] <- tan::AN_test(X[, clen], Y[, clen], na_rm=TRUE, pool= TRUE, poolVarX = varMinus,
        ##                                                         poolVarY = varPlus)$statistic
        ##                         }
        ##                         else {
        ##                            ANT[ii] <- tan::AN.test(X[, clen], Y[, clen], na.rm=TRUE, pool= TRUE, poolVarX = varMinus,
        ##                                                         poolVarY = varPlus)$statistic
        ##                         }
        ##                     }
        ##                 } # end of if (dim(X)[2] < s.size)
        ##                 else {
        ##                     if (ii %in% sitesUnused ) {
        ##                         design <- object@Designs[site, ]
        ##                         if (use_cpp) {
        ##                             if (ignore_sitesUnused == FALSE) {
        ##                                 ANT[ii] <- tan::AN_test(X[,design], Y[,design], na_rm = TRUE, pool = FALSE, poolVarX = NA,
        ##                                                     poolVarY = NA)$statistic
        ##                             }
        ##                             else {
        ##                                 ANT[ii] <- NA
        ##                             }
        ##                         }
        ##                         else {
        ##                             if (ignore_sitesUnused == FALSE) {
        ##                                 ANT[ii] <- tan::AN.test(X[,design], Y[,design], na.rm = TRUE)$statistic
        ##                             }
        ##                             else {
        ##                                 ANT[ii] <- NA
        ##                             }
        ##                         }
        ##                     } else {
        ##                         design <- object@Designs[site, ]
        ##                         clen <- 1:length(varMinus)
        ##                         if (use_cpp) {
        ##                             ANT[ii] <- tan::AN_test(X[,design[clen]], Y[,design[clen]], na_rm = TRUE, pool = TRUE, poolVarX = varMinus,
        ##                                                     poolVarY = varPlus)$statistic
        ##                         }
        ##                         else {
        ##                             ANT[ii] <- tan::AN.test(X[,design[clen]], Y[,design[clen]], na.rm = TRUE, pool = TRUE, poolVarX = varMinus,
        ##                                                     poolVarY = varPlus)$statistic
        ##                         }
        ##                     }
        ##                 }
        ##             } # end of for (ii in 1:length(idx))
        ##             if (length(H0.idx) > 0 ) {
        ##                 if (ignore_sitesUnused == FALSE) {
        ##                     p[idx, j] <- sapply(ANT, function(tanTest) {
        ##                         length(which(H0[H0.idx] >= tanTest )) / length(H0.idx)
        ##                     })
        ##                 }
        ##                 # new: 03/24/2017
        ##                 else {
        ##                     # ignore unused sites from testing between
        ##                     na_indices <- which(is.na(ANT) == TRUE)
        ##                     # ignore unused sites from testing within
        ##                     H0.idx_ <- setdiff(H0.idx, sitesUnused_Within)
        ##                     H0.idx_ <- rep(H0.idx_, ncol(Within))
        ##                     p[idx, j] <- sapply(ANT, function(tanTest) {
        ##                         length(which(H0[H0.idx_] >= tanTest )) / length(H0.idx_)
        ##                     })
        ##                     p[idx[na_indices], j] <- NA
        ##                 }
        ##             }
        ##         } # end if (minGlobal != Inf)
        ##     } # end of for (j in 1:ncomps)
        ## } # end of for (i in 1:length(dN))
        ## # impute missing values
        ## if (any(is.na(p))) {
        ##     message(length(which(is.na(p))),' NAs found  (of ', length(p),')')
        ##     if (na_impute) {
        ##         p[is.na(p)] = min(p[!is.na(p)])
        ##     }
        ## }
        ## Pc <- apply(p[, 1:Between_cols], 1, function(x) quantile(x, probs = quant, na.rm =TRUE))
        ## fdr <- p.adjust(as.vector(p),method='BH')
        ## FDR <- matrix(0,nrow = nrow(p), ncol = ncol(p)+1)

        ## for (i in 1:ncol(p)){
        ##     FDR[, i] <- fdr[ (i-1) * nrow(p) + (1:nrow(p)) ]
        ## }
        ## if (any(is.na(FDR))){
        ##     message(length(which(is.na(FDR))),' NAs found  (of ', length(FDR),')')
        ##     if (na_impute) {
        ##         FDR[is.na(FDR)]=min(FDR[!is.na(FDR)])
        ##     }
        ## }
        ## FDR[, i+1] <- apply(FDR[,1:Between_cols], 1, function(x) quantile(x, probs = quant, na.rm = TRUE))
        ## p <- cbind(p,Pc)
        ## colnames(p)[i+1] <- 'combined'
        ## colnames(FDR) <- colnames(p)
        ## P <- list(p,FDR)
        ## names(P) <- c('pval','FDR')
        ## object@PvalList <- P

        object@PvalList <- resultList[['PvalList']]
        object@sitesUnused <- resultList[['sitesUnused']]

    } # end of if (n=4)
    else if (object@nSamples == 3) {
        print(paste("Computing p-values for sample size n = ", object@nSamples), sep = "")
        ## Name columns:
        colnames(object@W1) <- c('ab vs ac', 'ab vs bc', 'ac vs bc')
        colnames(object@W2) <- c('AB vs AC', 'AB vs BC', 'AC vs BC')
        minLabs <- c('ab', 'ac', 'bc')
        plusLabs <- toupper(minLabs)
        ## Create colnames for matrix B
        labs <- ''
        for (mlab in minLabs) {
            for (plab in plusLabs) {
                labs <- c(labs, paste(mlab, 'vs', plab))
            }
        }
        labs <- labs[-c(1)]
        ### Get dictinary {a:1, b:2, c:3, ...}
        getDict <- function(s) {
            switch(s, a = 1, b = 2, c = 3,
                   A = 4, B = 5, C = 6)
        }
        ### Main ###
        total <-  nPeaks <- length(object@coverage)
        # number of columns between B: 3*3 columns
        Between_cols <- 9
        sitesUnused_Within <- object@sitesUnused # @
        Within <- cbind(object@W1, object@W2)
        AvsB <- matrix(NA, nrow = total, ncol = Between_cols)
        colnames(AvsB) <- labs
        AvsB <- cbind(AvsB, Within)
        ncomps <- ncol(AvsB)
        p <- matrix(NA, nrow = nrow(AvsB), ncol = ncol(AvsB))
        colnames(p) <- colnames(AvsB)
        H0 <- as.vector(Within)
        if (ncol(object@Ns) > 1) {
            temp <- colnames(Within)
            sampleNames <- unlist(strsplit(temp,' vs '))
            ## For each interval, take average of all ChIP samples "ac" "bd" "ad" "bc" "AC" "BD" "AD" "BC";
            ## 'ac' means take average between a and c for that interval and so on
            H0.Ns <- rowMeans(object@Ns[, sampleNames])
            H0.Ns <- rep(H0.Ns, ncol(Within))
        } else {
            H0.Ns <- rep(object@Ns, ncol(Within))
        }
        if (max(object@dN) < max(object@Ns)){
            object@dN <- c(object@dN, max(object@Ns))
        }
        # evaluate H0.idx
        for (i in 1:length(object@dN)){
            print(paste('+++Bin = ', i))
            if (i==1){
                if (ncol(object@Ns) == 1){
                    idx <- which(object@Ns<=object@dN[[i]])
                }
                H0.idx <- which(H0.Ns <= object@dN[[i]])
            } else {
                if (ncol(object@Ns) == 1){
                    idx <- which(object@Ns > object@dN[[i-1]] & object@Ns<= object@dN[[i]])
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
                if (ncol(object@Ns) > 1){
                    comp <- colnames(AvsB)[j]
                    sampleNames <- unlist(strsplit(comp,' vs '))
                    CompNs <- rowMeans(object@Ns[,sampleNames])
                    if (i==1) {
                        idx <- which(CompNs <= object@dN[[i]])
                    } else {
                        idx <- which(CompNs > object@dN[[i-1]] & CompNs<=object@dN[[i]])
                    }
                    if (length(idx)==0) next
                }
                ## Get label indices
                m1 <- getDict(substr(sampleNames[1], 1, 1))
                m2 <- getDict(substr(sampleNames[1], 2, 2))
                p1 <- getDict(substr(sampleNames[2], 1, 1))
                p2 <- getDict(substr(sampleNames[2], 2, 2))
                minusVar_idx <- plusVar_idx <- list()
                minGlobal <- Inf
                # compute between-TAN variance
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
                        if (use_cpp) {
                            test <- tan::AN_test(X, Y, na_rm = TRUE, pool = FALSE, poolVarX = NA, poolVarY = NA)
                        }
                        else {
                            test <-  tan::AN.test(X, Y, na.rm=TRUE)
                        }
                    }
                    else {
                        design <- object@Designs[site, ]
                        if (use_cpp) {
                            test <- tan::AN_test(X[, design], Y[, design], na_rm = TRUE, pool = FALSE, poolVarX = NA, poolVarY = NA)
                        }
                        else {
                            test <- tan::AN.test(X[, design], Y[, design], na.rm=TRUE)
                        }
                    }
                    minIndex <- min(length(test$varX), length(test$varY))
                    ## check minIndex > Global_lower (lower bound for pooled var vector of each bins)
                    if (minIndex > Global_lower) {
                        if (minGlobal > minIndex) {
                            minGlobal <- minIndex
                        }
                        minusVar_idx[[ii]] <- test$varX
                        plusVar_idx[[ii]] <- test$varY
                    }
                    # store all unused sites (minIndex <= Global_lower): minIndex = 0 -> var empty due to flat peak,
                    # or many repeated counts, or peakLength too small. 03/24/17
                    else {
                        # print(paste(" +++ minIndex <= Global_lower: ", minIndex, sep = ""))
                        minusVar_idx[[ii]] <- c()
                        plusVar_idx[[ii]] <- c()
                        sitesUnused <- c(sitesUnused, ii)
                    }
                } # end of for (ii in 1:length(idx))
                print(paste(" +++ total sites Unused: ", length(sitesUnused), sep = ""))
                print(paste(" +++ minGlobal = : ", minGlobal, sep = ""))
                ##################################################################
                ### Compute pooled variances for between-TAN in bin(test j) or idx
                ##################################################################
                ### The following only works if 0 < minGlobal < Inf
                if (minGlobal != Inf) {
                    matVar_minus <- matVar_plus <- matrix(NA, nrow = length(idx), ncol = minGlobal)
                    # check if sites whose len of var > Global_lower: if TRUE --> use pooled, otherwise use unpooled
                    for (ii in 1:length(minusVar_idx)) {
                        if (length(minusVar_idx[[ii]]) > 0 ) {
                            matVar_minus[ii,] <- minusVar_idx[[ii]][1:minGlobal]
                            matVar_plus[ii,] <- plusVar_idx[[ii]][1:minGlobal]
                        }
                    }
                    # TODO: COMBINE these two cases length(sitesUnused) > 0 or else
                    if (length(sitesUnused) > 0 ) {
                        varMinus <- apply(matVar_minus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varMinus) >= movAve ) {
                            varMinus <- tan::movingAverage(varMinus, movAve)
                        }
                        varPlus <- apply(matVar_plus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varPlus) >= movAve ) {
                            varPlus <- tan::movingAverage(varPlus, movAve)
                        }
                    }
                    # case sitesUnused = c(): use unpooled
                    else {
                        varMinus <- apply(matVar_minus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varMinus) >= movAve ) {
                            varMinus <- tan::movingAverage(varMinus, movAve)
                        }
                        varPlus <- apply(matVar_plus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varPlus) >= movAve ) {
                            varPlus <- tan::movingAverage(varPlus, movAve)
                        }
                    }
                }
                # case minGlobal = Inf
                else {
                    message("minGlobal = Inf. Please consider larger bin sizes and/or filtering sites with zero variances and small peakLengths!")
                }
                # add siteUnused_between to slot:
                object@sitesUnused <- unique(c(object@sitesUnused, idx[sitesUnused] ))
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
                            ##  check if site is in sitesUnused --> use unpooled var or return NA
                            ##  if site is NOT in sitesUnused --> use pooled var
                            #######################################################
                            if (ii %in% sitesUnused ) {
                                if (use_cpp) {
                                    if (ignore_sitesUnused == FALSE) {
                                        ANT[ii] <- tan::AN_test(X, Y, na_rm=TRUE, pool= FALSE, poolVarX = NA,
                                                                poolVarY = NA)$statistic
                                    }
                                    else {
                                        ANT[ii] <- NA
                                    }
                                }
                                # use_cpp = FALSE
                                else {
                                    if (ignore_sitesUnused == FALSE) {
                                        ANT[ii] <- tan::AN.test(X, Y, na.rm=TRUE)$statistic
                                    }
                                    else {
                                        ANT[ii] <- NA
                                    }
                                }
                            }
                            else {
                                clen <- 1:length(varMinus)
                                if (use_cpp) {
                                    ANT[ii] <- tan::AN_test(X[, clen], Y[, clen], na_rm=TRUE, pool= TRUE, poolVarX = varMinus,
                                                            poolVarY = varPlus)$statistic
                                }
                                # use_cpp = FALSE
                                else {
                                    ANT[ii] <- tan::AN.test(X[, clen], Y[, clen], na.rm=TRUE, pool= TRUE, poolVarX = varMinus,
                                                            poolVarY = varPlus)$statistic
                                }
                            }
                        } # end of if (dim(X)[2] < s.size)
                        else {
                            if (ii %in% sitesUnused ) {
                                design <- object@Designs[site, ]
                                if (use_cpp) {
                                    if (ignore_sitesUnused == FALSE) {
                                        ANT[ii] <- tan::AN_test(X[,design], Y[,design], na_rm = TRUE, pool = FALSE, poolVarX = NA,
                                                                poolVarY = NA)$statistic
                                    }
                                    else {
                                        ANT[ii] <- NA
                                    }
                                }
                                # use_cpp = FALSE
                                else {
                                    if (ignore_sitesUnused == FALSE) {
                                        ANT[ii] <- tan::AN.test(X[,design], Y[,design], na.rm = TRUE)$statistic
                                    }
                                    else {
                                        ANT[ii] <- NA
                                    }
                                }
                            }
                            else {
                                design <- object@Designs[site, ]
                                clen <- 1:length(varMinus)
                                if (use_cpp) {
                                    ANT[ii] <- tan::AN_test(X[, design[clen]], Y[, design[clen]], na_rm = TRUE, pool = TRUE, poolVarX = varMinus,
                                                            poolVarY = varPlus)$statistic
                                }
                                # use_cpp = FALSE
                                else {
                                    ANT[ii] <- tan::AN.test(X[,design[clen]], Y[,design[clen]], na.rm = TRUE, pool = TRUE, poolVarX = varMinus,
                                                            poolVarY = varPlus)$statistic
                                }
                            }
                        }
                    } # end of for (ii in 1:length(idx))
                    if (length(H0.idx) > 0 ) {
                        if (ignore_sitesUnused == FALSE) {
                            p[idx, j] <- sapply(ANT, function(tanTest) {
                                length(which(H0[H0.idx] >= tanTest )) / length(H0.idx)
                            })
                        }
                        # new: 03/24/2017
                        else {
                            # ignore unused sites from testing between
                            na_indices <- which(is.na(ANT) == TRUE)
                            # ignore unused sites from testing within
                            H0.idx_ <- setdiff(H0.idx, sitesUnused_Within)
                            H0.idx_ <- rep(H0.idx_, ncol(Within))
                            p[idx, j] <- sapply(ANT, function(tanTest) {
                                length(which(H0[H0.idx_] >= tanTest )) / length(H0.idx_)
                            })
                            p[idx[na_indices], j] <- NA
                        }
                    }
                } # end if (minGlobal != Inf)
            } # end of for (j in 1:ncomps)
        } # end of for (i in 1:length(dN))
        # impute missing values
        if (any(is.na(p))) {
            message(length(which(is.na(p))),' NAs found  (of ', length(p),')')
            if (na_impute) {
                p[is.na(p)] = min(p[!is.na(p)])
            }
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
        FDR[,i+1] <- apply(FDR[,1:Between_cols], 1, function(x) quantile(x, probs = quant, na.rm = TRUE))
        p <- cbind(p,Pc)
        colnames(p)[i+1] <- 'combined'
        colnames(FDR) <- colnames(p)
        P <- list(p,FDR)
        names(P) <- c('pval','FDR')
        object@PvalList <- P
    } # end of if (n=3)
    else if (object@nSamples == 2) {
        print(paste("Computing p-values for sample size n = ", object@nSamples), sep = "")
        ### Get dictinary {a:1, b:2, c:3, ...}
        getDict <- function(s) {
            switch(s, a = 1, b = 2,
                   A = 3, B = 4)
        }
        ### Main ###
        total <-  nPeaks <- length(object@coverage)
        sitesUnused_Within <- object@sitesUnused # @
        ncomps <- 1 # only compute p-value for ab vs AB
        p <- matrix(0, nrow = total, ncol = 1)
        colnames(p) <- colnames(object@W)[1]
        H0 <- as.vector(object@W)
        if (ncol(object@Ns)>1){
            #temp <- colnames(W)
            sampleNames <- c('ab', 'AB')
            ## For each interval, take average of all ChIP samples "ac" "bd" "ad" "bc" "AC" "BD" "AD" "BC";
            # 'ac' means take average between a and c for that interval and so on
            H0.Ns <- rowMeans(object@Ns[ ,sampleNames])
            H0.Ns <-rep(H0.Ns, ncol(object@W)) # ncol(W) = number of total performed tests
        } else {
            H0.Ns <-rep(object@Ns,ncol(object@W))
        }
        if (max(object@dN)<max(object@Ns)){
            object@dN <- c(object@dN, max(object@Ns))
        }
        # evaluate H0.idx
        for (i in 1:length(object@dN)){
            print(paste('+++Bin = ', i))
            if (i == 1){
                if (ncol(object@Ns) == 1){
                    idx <- which(object@Ns <= object@dN[[i]])
                }
                H0.idx <- which(H0.Ns <= object@dN[[i]])
            } else {
                if (ncol(object@Ns) == 1) {
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
                if (ncol(object@Ns) > 1) {
                    comp <- "ab vs AB"
                    sampleNames <- unlist(strsplit(comp,' vs '))
                    CompNs <- rowMeans(object@Ns[, sampleNames])
                    if (i==1){
                        idx <- which(CompNs <= object@dN[[i]])
                    } else {
                        idx <- which(CompNs > object@dN[[i-1]] & CompNs <= object@dN[[i]])
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
                # compute between-TAN variance
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
                        if (use_cpp) {
                            test <- tan::AN_test(X, Y, na_rm = TRUE, pool = FALSE, poolVarX = NA, poolVarY = NA)
                        }
                        # use_cpp = FALSE
                        else {
                            test <-  tan::AN.test(X, Y, na.rm = TRUE)
                        }
                    }
                    else {
                        design <- object@Designs[site, ]
                        if (use_cpp) {
                            test <- tan::AN_test(X[, design], Y[, design], na_rm = TRUE, pool = FALSE, poolVarX = NA, poolVarY = NA)
                        }
                        # use_cpp = FALSE
                        else {
                            test <- tan::AN.test(X[, design], Y[, design], na.rm = TRUE)
                        }
                    }
                    minIndex <- min(length(test$varX), length(test$varY))
                    ## check minIndex > Global_lower (lower bound for pooled var vector of each bins)
                    if (minIndex > Global_lower) {
                        if (minGlobal > minIndex) {
                            minGlobal <- minIndex
                        }
                        minusVar_idx[[ii]] <- test$varX
                        plusVar_idx[[ii]] <- test$varY
                    }
                    # store all unused sites (minIndex <= Global_lower): minIndex = 0 -> var empty due to flat peak,
                    # or many repeated counts, or peakLength too small. 03/24/17
                    else {
                        # print(paste(" +++ minIndex <= Global_lower: ", minIndex, sep = ""))
                        minusVar_idx[[ii]] <- c()
                        plusVar_idx[[ii]] <- c()
                        sitesUnused <- c(sitesUnused, ii)
                    }
                } # end of for (ii in 1:length(idx))
                print(paste(" +++ total sites Unused: ", length(sitesUnused), sep = ""))
                print(paste(" +++ minGlobal = : ", minGlobal, sep = ""))
                ##################################################################
                ### Compute pooled variances for between-TAN in bin(test j) or idx
                ##################################################################
                ### The following only works if 0 < minGlobal < Inf
                if (minGlobal != Inf) {
                    matVar_minus <- matVar_plus <- matrix(NA, nrow = length(idx), ncol = minGlobal)
                    # check if sites whose len of var > Global_lower: if TRUE --> use pooled, o.w. use unpool
                    for (ii in 1:length(minusVar_idx)) {
                        if (length(minusVar_idx[[ii]]) > 0 ) {
                            matVar_minus[ii,] <- minusVar_idx[[ii]][1:minGlobal]
                            matVar_plus[ii,] <- plusVar_idx[[ii]][1:minGlobal]
                        }
                    }
                    # TODO: COMBINE these two cases length(sitesUnused) > 0 or else
                    if (length(sitesUnused) > 0 ) {
                        varMinus <- apply(matVar_minus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varMinus) >= movAve ) {
                            varMinus <- tan::movingAverage(varMinus, movAve)
                        }
                        varPlus <- apply(matVar_plus[-sitesUnused, ], 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varPlus) >= movAve ) {
                            varPlus <- tan::movingAverage(varPlus, movAve)
                        }
                    }
                    # case sitesUnused = c(): use unpooled
                    else {
                        varMinus <- apply(matVar_minus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varMinus) >= movAve ) {
                            varMinus <- tan::movingAverage(varMinus, movAve)
                        }
                        varPlus <- apply(matVar_plus, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(varPlus) >= movAve ) {
                            varPlus <- tan::movingAverage(varPlus, movAve)
                        }
                    }
                }
                # case minGlobal = Inf
                else {
                    message("minGlobal = Inf. Please consider larger bin sizes and/or filtering sites with zero variances and small peakLengths!")
                }
                # add siteUnused_between to slot:
                object@sitesUnused <- unique(c(object@sitesUnused, idx[sitesUnused] ))
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
                                if (use_cpp) {
                                    if (ignore_sitesUnused == FALSE) {
                                        ANT[ii] <- tan::AN_test(X, Y, na_rm=TRUE, pool= FALSE, poolVarX = NA, poolVarY = NA)$statistic
                                    }
                                    else {
                                        ANT[ii] <- NA
                                    }
                                }
                                # use_cpp = FALSE
                                else {
                                    if (ignore_sitesUnused == FALSE) {
                                        ANT[ii] <- tan::AN.test(X, Y, na.rm=TRUE)$statistic
                                    }
                                    else {
                                        ANT[ii] <- NA
                                    }
                                }
                            }
                            ####
                            else {
                                clen <- 1:length(varMinus)
                                if (use_cpp) {
                                    ANT[ii] <- tan::AN_test(X[, clen], Y[, clen], na_rm=TRUE, pool= TRUE, poolVarX = varMinus,
                                                            poolVarY = varPlus)$statistic
                                }
                                else {
                                    ANT[ii] <- tan::AN.test(X[, clen], Y[, clen], na.rm=TRUE, pool= TRUE, poolVarX = varMinus,
                                                            poolVarY = varPlus)$statistic
                                }
                            }
                        } # end of if (dim(X)[2] < s.size)
                        else {
                            if (ii %in% sitesUnused ) {
                                design <- object@Designs[site, ]
                                if (use_cpp) {
                                    if (ignore_sitesUnused == FALSE) {
                                        ANT[ii] <- tan::AN_test(X[, design], Y[, design], na_rm = TRUE, pool = FALSE, poolVarX = NA,
                                                                poolVarY = NA)$statistic
                                    }
                                    else {
                                        ANT[ii] <- NA
                                    }
                                }
                                # use_cpp = FALSE
                                else {
                                    if (ignore_sitesUnused == FALSE) {
                                        ANT[ii] <- tan::AN.test(X[, design], Y[, design], na.rm = TRUE)$statistic
                                    }
                                    else {
                                        ANT[ii] <- NA
                                    }
                                }
                            } else {
                                design <- object@Designs[site, ]
                                clen <- 1:length(varMinus)
                                if (use_cpp) {
                                    ANT[ii] <- tan::AN_test(X[,design[clen]], Y[,design[clen]], na_rm = TRUE, pool = TRUE, poolVarX = varMinus,
                                                            poolVarY = varPlus)$statistic
                                }
                                else {
                                    ANT[ii] <- tan::AN.test(X[,design[clen]], Y[,design[clen]], na.rm = TRUE, pool = TRUE, poolVarX = varMinus,
                                                            poolVarY = varPlus)$statistic
                                }
                            }
                        }
                    } # end of for (ii in 1:length(idx))
                    if (length(H0.idx) > 0 ) {
                        if (ignore_sitesUnused == FALSE) {
                            p[idx, j] <- sapply(ANT, function(tanTest) {
                                length(which(H0[H0.idx] >= tanTest )) / length(H0.idx)
                            })
                        }
                        # new: 03/24/2017
                        else {
                            # ignore unused sites from testing between
                            na_indices <- which(is.na(ANT) == TRUE)
                            # ignore unused sites from testing within
                            H0.idx_ <- setdiff(H0.idx, sitesUnused_Within)
                            # H0.idx_ <- rep(H0.idx_, ncol(Within))
                            H0.idx_ <- rep(H0.idx_, ncol(object@W))
                            p[idx, j] <- sapply(ANT, function(tanTest) {
                                length(which(H0[H0.idx_] >= tanTest )) / length(H0.idx_)
                            })
                            p[idx[na_indices], j] <- NA
                        }
                    }
                } # end if (minGlobal != Inf)
            } # end of for (j in 1:ncomps)
        } # end of for (i in 1:length(dN))
        # impute missing values
        if (any(is.na(p))) {
            message(length(which(is.na(p))),' NAs found  (of ', length(p),')')
            if (na_impute) {
                p[is.na(p)] = min(p[!is.na(p)])
            }
        }
        Pc <- p
        fdr <- p.adjust(as.vector(p),method='BH')
        FDR <- matrix(0,nrow=nrow(p),ncol=ncol(p)+1)
        for (i in 1:ncol(p)){
            FDR[,i] <- fdr[(i-1)*nrow(p)+(1:nrow(p))]
        }
        if (any(is.na(FDR))){
            message(length(which(is.na(FDR))),' NAs found  (of ', length(FDR),')')
            FDR[is.na(FDR)]=min(FDR[!is.na(FDR)])
        }
        FDR[, i+1] <- FDR[, 1]
        p <- cbind(p,Pc)
        colnames(p)[i+1] <- 'combined'
        colnames(FDR) <- colnames(p)
        P <- list(p,FDR)
        names(P) <- c('pval','FDR')
        object@PvalList <- P
    }
    object
}

setMethod("computePvalues", signature("tanDb"), .computePvalues)
