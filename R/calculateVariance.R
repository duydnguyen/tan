.calculateVariance <- function(object, minus_condition, Global_lower, poolQuant, movAve, ...) {
    if (object@nSamples %in% c(3,4)) {
        if (minus_condition == TRUE) {
            print("Calculating Variance for first condition")
        } else {
            print("Calculating Variance for second condition")
        }
    }
    else if (object@nSamples == 2) {
        print("Calculating pool Variance for both conditions")
    }
    Var <- list()
    ### nSamples = 4 ###
    if (object@nSamples == 4) {
        print(paste("Calculating pooled variance for sample size n = ", object@nSamples), sep = "")
        lab_pool <- c('ab', 'ac', 'ad', 'bc', 'bd', 'cd')
        sitesUnused <- c()
        for (bin in 1:length(object@wSites)) {
            print(paste('bin = ', bin))
            if (length(object@wSites[[bin]]) > 0 ) {
                sites <- object@wSites[[bin]]
                df <- data.frame()
                varList <- list()
                minGlobal <- Inf
                count <- 1
                if (length(sites)>0) {
                    for (site in sites) {
                        # print(paste(' -------------------- site = ', site))
                        if (minus_condition == TRUE) {
                            X1 <- object@coverage[[site]][1:2,]
                            Y1 <- object@coverage[[site]][3:4,]
                            X2 <- object@coverage[[site]][c(1,3),]
                            Y2 <- object@coverage[[site]][c(2,4),]
                            X3 <- object@coverage[[site]][c(1,4),]
                            Y3 <- object@coverage[[site]][c(2,3),]
                        }
                        else {
                            X1 <- object@coverage[[site]][5:6,]
                            Y1 <- object@coverage[[site]][7:8,]
                            X2 <- object@coverage[[site]][c(5,7),]
                            Y2 <- object@coverage[[site]][c(6,8),]
                            X3 <- object@coverage[[site]][c(5,8),]
                            Y3 <- object@coverage[[site]][c(6,7),]
                        }
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
                            df <- data.frame('ab' = test1$varX[1:minIndex], 'ac' = test2$varX[1:minIndex], 'ad' = test3$varX[1:minIndex],
                                             'bc' = test3$varY[1:minIndex], 'bd' = test2$varY[1:minIndex], 'cd' = test1$varY[1:minIndex])
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
                        Var[[bin]] <- NA
                    }
                }
            }
        } # end of bin
        object@sitesUnused <- unique(c(object@sitesUnused, sitesUnused))
    } # end of if (n = 4)
    else if (object@nSamples == 3) {
        print(paste("Calculating pooled variance for sample size n = ", object@nSamples), sep = "")
        lab_pool <- c('ab', 'ac', 'bc') # label of S_
        sitesUnused <- c()
        for (bin in 1:length(object@wSites)) {
            print(paste('bin = ', bin))
            if (length(object@wSites[[bin]]) > 0 ) {
                sites <- object@wSites[[bin]]
                df <- data.frame()
                varList <- list()
                minGlobal <- Inf
                count <- 1
                if (length(sites)>0) {
                    for (site in sites) {
                        # print(paste('site = ', site))
                        if (minus_condition == TRUE) {
                            X1 <- object@coverage[[site]][1:2,] #ab
                            Y1 <- object@coverage[[site]][c(1,3),] #ac
                            X2 <- object@coverage[[site]][1:2,] #ab
                            Y2 <- object@coverage[[site]][2:3,] #bc
                            X3 <- object@coverage[[site]][c(1,3),] #ac
                            Y3 <- object@coverage[[site]][2:3,] #bc
                        }
                        else {
                            X1 <- object@coverage[[site]][4:5,]
                            Y1 <- object@coverage[[site]][c(4,6),]
                            X2 <- object@coverage[[site]][4:5,]
                            Y2 <- object@coverage[[site]][5:6,]
                            X3 <- object@coverage[[site]][c(4,6),]
                            Y3 <- object@coverage[[site]][5:6,]

                        }
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
                            df <- data.frame('ab' = test1$varX[1:minIndex], 'ac' = test2$varX[1:minIndex], 'ad' = test3$varX[1:minIndex],
                                             'bc' = test3$varY[1:minIndex], 'bd' = test2$varY[1:minIndex], 'cd' = test1$varY[1:minIndex])
                            varList[[count]] <- df
                            count <- count + 1
                        }
                        # store all unused sites (minIndex < Global_lower): minIndex = 0 -> var empty due to flat peak,
                        # or manyrepeated counts: 03/24/17
                        else {
                            sitesUnused <- c(sitesUnused, site)
                        }
                    } # end of for (site in sites)
                    ## Pooling variances across sites in bin
                    poolVar <- list()
                    print(paste(" +++ minGlobal = ", minGlobal, sep = ""))
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
            }
        } # end of bin
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
        } # end of bin
        object@sitesUnused <- unique(c(object@sitesUnused, sitesUnused))
    }
    # return results
    if (object@nSamples == 2) {
        object@poolVar <- Var
    }
    else if (object@nSamples %in% c(3,4)) {
        if (minus_condition) {
            object@minusVar <- Var
        } else {
            object@plusVar <- Var
        }
    }
    object
}

setMethod("calculateVariance", signature("tanDb"), .calculateVariance)
