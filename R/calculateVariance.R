.calculateVariance <- function(object, minus_condition, Global_lower, poolQuant, movAve) {
    if (minus_condition == TRUE) {
        print("Calculating Variance for first condition")
    } else {
        print("Calculating Variance for second condition")
    }
    Var <- list()
    ### nSamples = 4 ###
    if (object@nSamples == 4) {
        print(paste("Calculating pooled variance for sample size n = ", object@nSamples), sep = "")
        lab_pool <- c('ab', 'ac', 'ad', 'bc', 'bd', 'cd')
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
                            test1 <- tan::compute_Var(X1, Y1, na_rm = TRUE, pool = FALSE)
                            test2 <- tan::compute_Var(X2, Y2, na_rm = TRUE, pool = FALSE)
                            test3 <- tan::compute_Var(X3, Y3, na_rm = TRUE, pool = FALSE)
                        }
                        else {
                            design <- object@Designs[site, ]
                            test1 <- tan::compute_Var(X1[, design], Y1[, design], na_rm = TRUE, pool = FALSE)
                            test2 <- tan::compute_Var(X2[, design], Y2[, design], na_rm = TRUE, pool = FALSE)
                            test3 <- tan::compute_Var(X3[, design], Y3[, design], na_rm = TRUE, pool = FALSE)
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

    } # end of if (n = 4)
    else if (object@nSamples == 3) {
        print(paste("Calculating pooled variance for sample size n = ", object@nSamples), sep = "")
        lab_pool <- c('ab', 'ac', 'bc') # label of S_
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
                        if (minus == TRUE) {
                            X1 <- object@coverage[[site]][1:2,] #ab
                            Y1 <- object@coverage[[site]][c(1,3),] #ac
                            X2 <- object@coverage[[site]][1:2,] #ab
                            Y2 <- object@coverage[[site]][2:3,] #bc
                            X3 <- object@coverage[[site]][c(1,3),] #ac
                            Y3 <- object@coverage[[site]][2:3,] #bc
                        }
                        else {
                            X1 <- object@coverage[[site]][5:6,]
                            Y1 <- object@coverage[[site]][c(5,7),]
                            X2 <- object@coverage[[site]][5:6,]
                            Y2 <- object@coverage[[site]][6:7,]
                            X3 <- object@coverage[[site]][c(5,7),]
                            Y3 <- object@coverage[[site]][6:7,]

                        }
                        if ( dim(X1)[2] < object@s.size ) {
                            test1 <- tan::compute_Var(X1, Y1, na_rm = TRUE, pool = FALSE)
                            test2 <- tan::compute_Var(X2, Y2, na_rm = TRUE, pool = FALSE)
                            test3 <- tan::compute_Var(X3, Y3, na_rm = TRUE, pool = FALSE)
                        }
                        else {
                            design <- object@Designs[site, ]
                            test1 <- tan::compute_Var(X1[, design], Y1[, design], na_rm = TRUE, pool = FALSE)
                            test2 <- tan::compute_Var(X2[, design], Y2[, design], na_rm = TRUE, pool = FALSE)
                            test3 <- tan::compute_Var(X3[, design], Y3[, design], na_rm = TRUE, pool = FALSE)
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
    } else if (object@nSamples == 2) {
        #TODO
    }
    # return results
    if (minus_condition) {
        object@minusVar <- Var
    } else {
        object@plusVar <- Var
    }
    object
}

setMethod("calculateVariance", signature("tanDb"), .calculateVariance)
