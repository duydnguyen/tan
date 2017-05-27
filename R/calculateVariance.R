evaluate_variance <- function(coverage, nSamples, wSites, lab_pool, minus_condition,
                             use_cpp = TRUE, s.size, Designs, Global_lower) {
    ## Get dictionary {a:1, b:2, c:3, d:4,...}
    getDict <- function(s) {
        # only for lower case mapping
        return(match(tolower(s), letters))
    }
    ## create pair sample for within labels
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
            return(list('withinLabel' = temp))
        }
        else if (nSamples == 3) {
            temp <- c("ab vs ac", "ab vs bc", "ac vs bc",
                     "AB vs AC", "AB vs BC", "AC vs BC")
            return(list('withinLabel' = temp))
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
    Var <- list()
    sitesUnused <- c()
    ## n = 4: ab vs cd, ac vs bd, ad vs bc, and other three for second cond
    ## n = 3: ab vs ac, ab vs bc, ac vs bc, and other three for second cond
    withinLabel <- create_labels(nSamples)[['withinLabel']]
    for (bin in 1:length(wSites)) {
        print(paste('bin = ', bin))
        if (length(wSites[[bin]]) > 0 ) {
            sites <- wSites[[bin]]
            df <- data.frame()
            varList <- list()
            minGlobal <- Inf
            count <- 1
            if (length(sites)>0) {
                for (site in sites) {
                    # print(paste(' -------------------- site = ', site))
                    testList <- list()
                    withinX <- withinY <- list()
                    numTests <-  c()
                    indexList <- list()
                    if (minus_condition == TRUE) {
                        ## X1 <- coverage[[site]][1:2,]
                        ## Y1 <- coverage[[site]][3:4,]
                        ## X2 <- coverage[[site]][c(1,3),]
                        ## Y2 <- coverage[[site]][c(2,4),]
                        ## X3 <- coverage[[site]][c(1,4),]
                        ## Y3 <- coverage[[site]][c(2,3),]

                        # new TODO
                        indexList <- create_indexList(nSamples)
                        numTests <- length(indexList)
                        for (tt in 1:numTests) {
                            ids <- indexList[[tt]]
                            withinX[[tt]] <- coverage[[site]][ids[1:2], ]
                            withinY[[tt]] <- coverage[[site]][ids[3:4],]
                        }
                    }
                    else {
                        ## X1 <- coverage[[site]][5:6,]
                        ## Y1 <- coverage[[site]][7:8,]
                        ## X2 <- coverage[[site]][c(5,7),]
                        ## Y2 <- coverage[[site]][c(6,8),]
                        ## X3 <- coverage[[site]][c(5,8),]
                        ## Y3 <- coverage[[site]][c(6,7),]

                        # new TODO
                        indexList <- create_indexList(nSamples)
                        numTests <- length(indexList)
                        for (tt in 1:numTests) {
                            ids <- indexList[[tt]] + nSamples
                            withinX[[tt]] <- coverage[[site]][ids[1:2], ]
                            withinY[[tt]] <- coverage[[site]][ids[3:4],]
                        }
                    }
                    ## if ( dim(X1)[2] < s.size ) {
                    if ( dim(withinX[[1]])[2] < s.size ) {
                        if (use_cpp) {
                            ## test1 <- tan::compute_Var(X1, Y1, na_rm = TRUE, pool = FALSE)
                            ## test2 <- tan::compute_Var(X2, Y2, na_rm = TRUE, pool = FALSE)
                            ## test3 <- tan::compute_Var(X3, Y3, na_rm = TRUE, pool = FALSE)

                            for (tt in 1:numTests) {
                                X <- withinX[[tt]]
                                Y <- withinY[[tt]]
                                testList[[tt]] <- tan::compute_Var(X, Y, na_rm = TRUE, pool = FALSE)
                            }
                        }
                        else {
                            ## test1 <-  tan::AN.test(X1, Y1, na.rm=TRUE)
                            ## test2 <-  tan::AN.test(X2, Y2, na.rm=TRUE)
                            ## test3 <-  tan::AN.test(X3, Y3, na.rm=TRUE)

                            for (tt in 1:numTests) {
                                X <- withinX[[tt]]
                                Y <- withinY[[tt]]
                                testList[[tt]] <- tan::AN.test(X, Y, na_rm = TRUE)
                            }

                        }
                    }
                    else {
                        design <- Designs[site, ]
                        if (use_cpp) {
                            ## test1 <- tan::compute_Var(X1[, design], Y1[, design], na_rm = TRUE, pool = FALSE)
                            ## test2 <- tan::compute_Var(X2[, design], Y2[, design], na_rm = TRUE, pool = FALSE)
                            ## test3 <- tan::compute_Var(X3[, design], Y3[, design], na_rm = TRUE, pool = FALSE)

                            for (tt in 1:numTests) {
                                X <- withinX[[tt]]
                                Y <- withinY[[tt]]
                                testList[[tt]] <- tan::compute_Var(X[, design], Y[, design], na_rm = TRUE, pool = FALSE)
                            }
                        }
                        else {
                            ## test1 <- tan::AN.test(X1[, design], Y1[, design], na.rm=TRUE)
                            ## test2 <- tan::AN.test(X2[, design], Y2[, design], na.rm=TRUE)
                            ## test3 <- tan::AN.test(X3[, design], Y3[, design], na.rm=TRUE)

                            for (tt in 1:numTests) {
                                X <- withinX[[tt]]
                                Y <- withinY[[tt]]
                                testList[[tt]] <- tan::AN.test(X[, design], Y[, design], na_rm = TRUE)
                            }
                        }
                    }
                    ## minIndex <- min(c(length(test1$varX), length(test2$varX),length(test3$varX),
                    ##                  length(test1$varY), length(test2$varY), length(test3$varY)))

                    lenIndices <- c()
                    for (tt in 1:numTests) {
                        test_ <- testList[[tt]]
                        lenIndices <- c(lenIndices, length(test_$varX), length(test_$varY))
                    }
                    minIndex <- min(lenIndices)

                    ## check minIndex > Global_lower (lower bound for pooled var vector of each bins)
                    if (minIndex > Global_lower) {
                        if (minGlobal > minIndex) {
                            minGlobal <- minIndex
                        }
                        ## df <- data.frame('ab' = test1$varX[1:minIndex], 'ac' = test2$varX[1:minIndex],
                        ##                 'ad' = test3$varX[1:minIndex], 'bc' = test3$varY[1:minIndex],
                        ##                 'bd' = test2$varY[1:minIndex], 'cd' = test1$varY[1:minIndex])

                        df <- data.frame(matrix(NA, nrow = minIndex, ncol = length(lab_pool)))
                        if (nSamples > 3) {
                            col_id <- 1
                            for (tt in 1:numTests) {
                                test_ <- testList[[tt]]
                                ids <- indexList[[tt]]
                                df[, col_id] <- test_$varX[1:minIndex]
                                df[, col_id + 1] <- test_$varY[1:minIndex]
                                colnames(df)[col_id:(col_id+1)] <- c(paste(letters[ids[1:2]], collapse = ""),
                                                                    paste(letters[ids[3:4]], collapse = ""))
                                col_id <- col_id + 2
                            }
                        } else if (nSamples == 3) {
                            df <- data.frame('ab' = testList[[1]]$varX[1:minIndex],
                                            'ac' = testList[[1]]$varY[1:minIndex],
                                            'bc' = testList[[2]]$varY[1:minIndex])
                        }

                        varList[[count]] <- df
                        count <- count + 1
                    }
                    # store all unused sites (minIndex < Global_lower): minIndex = 0 -> var empty due to flat peak,
                    # or many repeated counts, or peakLength too small. 03/24/17
                    else {
                        sitesUnused <- c(sitesUnused, site)
                    }
                } # end of for (site in sites)
                ## Pooling variances across sites in bin
                poolVar <- list()
                print(paste(" +++ minGlobal = ", minGlobal, sep = ""))
                ## Case: minGlobal < Inf
                if (minGlobal < Inf) {
                    matVar <- matrix(NA, nrow = length(varList), ncol = minGlobal) # @: Case minGlobal = Inf
                    for (pair in lab_pool) {
                        for (i in 1:length(varList)) {
                            matVar[i, ] <- varList[[i]][1:minGlobal, pair]
                        }
                        var <- apply(matVar, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(var) >= movAve ) {
                            var <- tan::movingAverage(var, movAve)
                        }
                        poolVar[[pair]] <- var
                    }
                    Var[[bin]] <- poolVar
                }
                ## Case: minGlobal = Inf
                else {
                    message("minGlobal = Inf: 1. Variance vector for this bin returned NA,
                               and/or 2. Sites in this bin stored in sitesUnused slot")
                    Var[[bin]] <- NA
                }
            }
        }
    }
    return(list('Var' = Var, 'sitesUnused' = sitesUnused))
}

.calculateVariance <- function(object, minus_condition, Global_lower, poolQuant, movAve, ...) {
    if (object@nSamples > 2 ) {
        if (minus_condition == TRUE) {
            print("Calculating Variance for first condition")
        } else {
            print("Calculating Variance for second condition")
        }
    }
    else if (object@nSamples == 2) {
        print("Calculating pool Variance for both conditions")
    }
    ### MAIN ###
    Var <- list()
    if (object@nSamples > 3) {
        print(paste("Calculating pooled variance for sample size n = ", object@nSamples), sep = "")
        # lab_pool <- c('ab', 'ac', 'ad', 'bc', 'bd', 'cd')
        lab_pool <- colnames(object@Ns)[1:(dim(object@Ns)[2] / 2)]
        sitesUnused <- c()
        resultList <- evaluate_variance(coverage = object@coverage, nSamples = object@nSamples,
                          wSites = object@wSites, lab_pool = lab_pool,
                          minus_condition = minus_condition, use_cpp = use_cpp,
                          s.size = object@s.size, Designs = object@Designs,
                          Global_lower = Global_lower)
        Var <- resultList[['Var']]
        sitesUnused <- resultList[['sitesUnused']]
        object@sitesUnused <- unique(c(object@sitesUnused, sitesUnused))
    } # end of if (n > 3)
    else if (object@nSamples == 3) {
        print(paste("Calculating pooled variance for sample size n = ", object@nSamples), sep = "")
        # lab_pool <- c('ab', 'ac', 'bc') # label of S_
        lab_pool <- colnames(object@Ns)[1:(dim(object@Ns)[2] / 2)]
        sitesUnused <- c()

        resultList <- evaluate_variance(coverage = object@coverage, nSamples = object@nSamples,
                          wSites = object@wSites, lab_pool = lab_pool,
                          minus_condition = minus_condition, use_cpp = use_cpp,
                          s.size = object@s.size, Designs = object@Designs,
                          Global_lower = Global_lower)
        Var <- resultList[['Var']]
        sitesUnused <- resultList[['sitesUnused']]

        ## for (bin in 1:length(object@wSites)) {
        ##     print(paste('bin = ', bin))
        ##     if (length(object@wSites[[bin]]) > 0 ) {
        ##         sites <- object@wSites[[bin]]
        ##         df <- data.frame()
        ##         varList <- list()
        ##         minGlobal <- Inf
        ##         count <- 1
        ##         if (length(sites) > 0) {
        ##             for (site in sites) {
        ##                 if (minus_condition == TRUE) {
        ##                     X1 <- object@coverage[[site]][1:2,] #ab
        ##                     Y1 <- object@coverage[[site]][c(1,3),] #ac
        ##                     X2 <- object@coverage[[site]][1:2,] #ab
        ##                     Y2 <- object@coverage[[site]][2:3,] #bc
        ##                     X3 <- object@coverage[[site]][c(1,3),] #ac
        ##                     Y3 <- object@coverage[[site]][2:3,] #bc
        ##                 }
        ##                 else {
        ##                     X1 <- object@coverage[[site]][4:5,]
        ##                     Y1 <- object@coverage[[site]][c(4,6),]
        ##                     X2 <- object@coverage[[site]][4:5,]
        ##                     Y2 <- object@coverage[[site]][5:6,]
        ##                     X3 <- object@coverage[[site]][c(4,6),]
        ##                     Y3 <- object@coverage[[site]][5:6,]

        ##                 }
        ##                 if ( dim(X1)[2] < object@s.size ) {
        ##                     if (use_cpp) {
        ##                         test1 <- tan::compute_Var(X1, Y1, na_rm = TRUE, pool = FALSE)
        ##                         test2 <- tan::compute_Var(X2, Y2, na_rm = TRUE, pool = FALSE)
        ##                         test3 <- tan::compute_Var(X3, Y3, na_rm = TRUE, pool = FALSE)
        ##                     }
        ##                     else {
        ##                         test1 <-  tan::AN.test(X1, Y1, na.rm=TRUE)
        ##                         test2 <-  tan::AN.test(X2, Y2, na.rm=TRUE)
        ##                         test3 <-  tan::AN.test(X3, Y3, na.rm=TRUE)
        ##                     }

        ##                 }
        ##                 else {
        ##                     design <- object@Designs[site, ]
        ##                     if (use_cpp) {
        ##                         test1 <- tan::compute_Var(X1[, design], Y1[, design], na_rm = TRUE, pool = FALSE)
        ##                         test2 <- tan::compute_Var(X2[, design], Y2[, design], na_rm = TRUE, pool = FALSE)
        ##                         test3 <- tan::compute_Var(X3[, design], Y3[, design], na_rm = TRUE, pool = FALSE)
        ##                     }
        ##                     else {
        ##                         test1 <- tan::AN.test(X1[, design], Y1[, design], na.rm=TRUE)
        ##                         test2 <- tan::AN.test(X2[, design], Y2[, design], na.rm=TRUE)
        ##                         test3 <- tan::AN.test(X3[, design], Y3[, design], na.rm=TRUE)
        ##                     }
        ##                 }
        ##                 minIndex <- min(c(length(test1$varX), length(test2$varX),length(test3$varX),
        ##                                   length(test1$varY), length(test2$varY), length(test3$varY)))
        ##                 ## check minIndex > Global_lower (lower bound for pooled var vector of each bins)
        ##                 if (minIndex > Global_lower) {
        ##                     if (minGlobal > minIndex) {
        ##                         minGlobal <- minIndex
        ##                     }
        ##                     df <- data.frame('ab' = test1$varX[1:minIndex], 'ac' = test2$varX[1:minIndex], 'ad' = test3$varX[1:minIndex],
        ##                                      'bc' = test3$varY[1:minIndex], 'bd' = test2$varY[1:minIndex], 'cd' = test1$varY[1:minIndex])
        ##                     varList[[count]] <- df
        ##                     count <- count + 1
        ##                 }
        ##                 # store all unused sites (minIndex < Global_lower): minIndex = 0 -> var empty due to flat peak,
        ##                 # or manyrepeated counts: 03/24/17
        ##                 else {
        ##                     sitesUnused <- c(sitesUnused, site)
        ##                 }
        ##             } # end of for (site in sites)
        ##             ## Pooling variances across sites in bin
        ##             poolVar <- list()
        ##             print(paste(" +++ minGlobal = ", minGlobal, sep = ""))
        ##             ## Case: minGlobal < Inf
        ##             if (minGlobal < Inf) {
        ##                 matVar <- matrix(NA, nrow = length(varList), ncol = minGlobal)
        ##                 for (pair in lab_pool) {
        ##                     for (i in 1:length(varList)) {
        ##                         matVar[i, ] <- varList[[i]][1:minGlobal, pair]
        ##                     }
        ##                     var <- apply(matVar, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
        ##                     if ( length(var) >= movAve ) {
        ##                         var <- tan::movingAverage(var, movAve)
        ##                     }
        ##                     poolVar[[pair]] <- var
        ##                 }
        ##                 Var[[bin]] <- poolVar
        ##             }
        ##             ## Case: minGlobal = Inf
        ##             else {
        ##                 message("minGlobal = Inf: 1. Variance vector for this bin returned NA,
        ##                             and/or 2. Sites in this bin stored in sitesUnused slot")
        ##                 Var[[bin]] <- NA
        ##             }
        ##         }
        ##     }
        ## } # end of bin
        object@sitesUnused <- unique(c(object@sitesUnused, sitesUnused))
    }
    else if (object@nSamples == 2) {
        print(paste("Calculating pooled variance for sample size n = ", object@nSamples), sep = "")
        Var <- list()
        lab_pool <- c('ab', 'aA', 'aB', 'AB', 'bB', 'Ab')
        sitesUnused <- c()
        for (bin in 1:length(object@wSites)) {
            print(paste('bin = ', bin))
            sites <- object@wSites[[bin]]
            df <- data.frame()
            varList <- list()
            minGlobal <- Inf
            count <- 1
            if (length(sites) > 0) {
                for (site in sites) {
                    # print(paste('site = ', site))
                    geta <- object@coverage[[site]][1,]
                    getb <- object@coverage[[site]][2,]
                    getA <- object@coverage[[site]][3,]
                    getB <- object@coverage[[site]][4,]
                    X1 <- rbind( geta, getb)
                    Y1 <- rbind( getA, getB)
                    X2 <- rbind( geta, getA)
                    Y2 <- rbind( getb, getB)
                    X3 <- rbind( geta, getB)
                    Y3 <- rbind( getA, getb)
                    if ( dim(X1)[2] < object@s.size ) {
                        if (use_cpp) {
                            test1 <- tan::compute_Var(X1, Y1, na_rm = TRUE, pool = FALSE)
                            test2 <- tan::compute_Var(X2, Y2, na_rm = TRUE, pool = FALSE)
                            test3 <- tan::compute_Var(X3, Y3, na_rm = TRUE, pool = FALSE)
                        }
                        else {
                            test1 <-  tan::AN.test(X1, Y1, na.rm=TRUE)
                            test2 <-  tan::AN.test(X2, Y2, na.rm=TRUE)
                            test3 <-  tan::AN.test(X3, Y3, na.rm=TRUE)
                        }
                    }
                    else {
                        design <- object@Designs[site, ]
                        if (use_cpp) {
                            test1 <- tan::compute_Var(X1[, design], Y1[, design], na_rm = TRUE, pool = FALSE)
                            test2 <- tan::compute_Var(X2[, design], Y2[, design], na_rm = TRUE, pool = FALSE)
                            test3 <- tan::compute_Var(X3[, design], Y3[, design], na_rm = TRUE, pool = FALSE)
                        }
                        else {
                            test1 <- tan::AN.test(X1[, design], Y1[, design], na.rm=TRUE)
                            test2 <- tan::AN.test(X2[, design], Y2[, design], na.rm=TRUE)
                            test3 <- tan::AN.test(X3[, design], Y3[, design], na.rm=TRUE)
                        }
                    }
                    minIndex <- min(c(length(test1$varX), length(test2$varX),length(test3$varX),
                                      length(test1$varY), length(test2$varY), length(test3$varY)))
                    ## check minIndex > Global_lower (lower bound for pooled var vector of each bins)
                    if (minIndex > Global_lower) {
                        if (minGlobal > minIndex) {
                            minGlobal <- minIndex
                        }
                        df <- data.frame('ab' = test1$varX[1:minIndex], 'aA' = test2$varX[1:minIndex], 'aB' = test3$varX[1:minIndex],
                                         'AB' = test1$varY[1:minIndex], 'bB' = test2$varY[1:minIndex], 'Ab' = test3$varY[1:minIndex])
                        varList[[count]] <- df
                        count <- count + 1 # keep track of all sites in bin i
                    } # end of if (minIndex > Global_lower)
                    # store all unused sites (minIndex < Global_lower): minIndex = 0 -> var empty due to flat peak,
                    # or manyrepeated counts: 03/24/17
                    else {
                        sitesUnused <- c(sitesUnused, site)
                    }
                } # end of for (site in sites)
                ## Pooling variances across sites in bin
                poolVar <- list()
                print(paste(" +++ minGlobal = ", minGlobal, sep = ""))
                ## Case: minGlobal < Inf
                if (minGlobal < Inf) {
                    matVar <- matrix(NA, nrow = length(varList), ncol = minGlobal)
                    for (pair in lab_pool) {
                        for (i in 1:length(varList)) {
                            matVar[i,] <- varList[[i]][1:minGlobal, pair]
                        }
                        var <- apply(matVar, 2, function(x) quantile(x, probs = poolQuant, na.rm = TRUE))
                        if ( length(var) >= movAve ) {
                            var <- tan::movingAverage(var, movAve)
                        }
                        poolVar[[pair]] <- var
                    }
                    Var[[bin]] <- poolVar
                }
                ## Case: minGlobal = Inf
                else {
                    message("minGlobal = Inf: 1. Variance vector for this bin returned NA, and 2. Sites in this bin stored in sitesUnused slot")
                    Var[[bin]] <- NA
                }
            }
        } # end of bin
        object@sitesUnused <- unique(c(object@sitesUnused, sitesUnused))
    }
    # return results
    if (object@nSamples == 2) {
        object@poolVar <- Var
    }
    else if (object@nSamples > 2 ) {
        if (minus_condition) {
            object@minusVar <- Var
        } else {
            object@plusVar <- Var
        }
    }
    object
}

setMethod("calculateVariance", signature("tanDb"), .calculateVariance)
