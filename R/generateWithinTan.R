.generateWithinTan <- function(object, minus_condition, ...) {
    if (object@nSamples %in% c(3,4)) {
        if (minus_condition == TRUE) {
            print("Generating the Within Adaptive test for first condition")
        } else {
            print("Generating the Within Adaptive test for second condition")
        }
    }
    else if (object@nSamples == 2) {
        print("Calculating the Within Adaptive tests for both conditions")
    }
    W <- matrix()
    if (object@nSamples == 4) {
        print(paste("Generating the Within Adaptive test for sample size n = ", object@nSamples), sep = "")
        total <- length(object@coverage)
        # number of columns within W1, W2: (ab vs. cd, ac vs. bd, ad vs. bc; sim. for plus)
        Within_cols <- 3
        Sites <- 1:total
        # Create a bin label for each sites from wSites
        binLabs <- c()
        for (bin in 1:length(object@wSites)) {
            sites <- object@wSites[[bin]]
            if (length(sites) > 0) {
                binLabs[sites] <- bin
            }
        }
        # Initilize W
        W <- matrix(NA, nrow = total, ncol = Within_cols)
        index <- 1 # index for W matrix
        for (site in Sites) {
            bin <- binLabs[site]
            if (minus_condition == TRUE) {
                pooled <- object@minusVar[[bin]] # could returned NA if Bin's sites are "low quality"
                colnames(W) <- c('ab vs cd', 'ac vs bd', 'ad vs bc' )
                X1 <- object@coverage[[site]][1:2,]
                Y1 <- object@coverage[[site]][3:4,]
                X2 <- object@coverage[[site]][c(1,3),]
                Y2 <- object@coverage[[site]][c(2,4),]
                X3 <- object@coverage[[site]][c(1,4),]
                Y3 <- object@coverage[[site]][c(2,3),]
            }
            else {
                pooled <- object@plusVar[[bin]]
                colnames(W) <- toupper(c('ab vs cd', 'ac vs bd', 'ad vs bc' ))
                X1 <- object@coverage[[site]][5:6,]
                Y1 <- object@coverage[[site]][7:8,]
                X2 <- object@coverage[[site]][c(5,7),]
                Y2 <- object@coverage[[site]][c(6,8),]
                X3 <- object@coverage[[site]][c(5,8),]
                Y3 <- object@coverage[[site]][c(6,7),]
            }
            if (site %% 1000 == 0) {
                print(paste('W, site: ', site))
            }
            if ( dim(X1)[2] < object@s.size ) {
                if (use_cpp) {
                    # pooled is NOT NA
                    if (is.na(pooled) == FALSE) {
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[site,1] <- tan::AN_test(X1[, clen], Y1[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                  poolVarY = poolVarY)$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[site,2] <- tan::AN_test(X2[, clen], Y2[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                  poolVarY = poolVarY)$statistic
                        poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[site,3] <- tan::AN_test(X3[, clen], Y3[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                  poolVarY = poolVarY)$statistic
                        index <- index + 1
                    }
                    # pooled is NA
                    else {
                        W[site, 1] <- W[site, 2] <- W[site, 3] <- NA
                    }
                }
                # if use_cpp = FALSE
                else {
                    if (is.na(pooled) == FALSE) {
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[site,1] <- tan::AN.test(X1[, clen], Y1[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                  poolVarY = poolVarY)$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[site,2] <- tan::AN.test(X2[, clen], Y2[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                  poolVarY = poolVarY)$statistic
                        poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[site,3] <- tan::AN.test(X3[, clen], Y3[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                  poolVarY = poolVarY)$statistic
                        index <- index + 1
                    }
                    # pooled is NA
                    else {
                        W[site, 1] <- W[site, 2] <- W[site, 3] <- NA
                    }
                }
            }
            # peakLength > s.size
            else {
                if (use_cpp) {
                    if (is.na(pooled) == FALSE) {
                        design <- object@Designs[site, ]
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN_test(X1[, design[clen]], Y1[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[site,1] <- rAN$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN_test(X2[, design[clen]], Y2[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[site,2] <- rAN$statistic
                        poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN_test(X3[, design[clen]], Y3[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[site,3] <- rAN$statistic
                        index <- index +1
                    }
                    # pooled = NA
                    else {
                        W[site, 1] <- W[site, 2] <- W[site, 3] <- NA
                    }
                }
                # use_cpp = FALSE
                else {
                    if (is.na(pooled) == FALSE) {
                        design <- object@Designs[site, ]
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN.test(X1[, design[clen]], Y1[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[site,1] <- rAN$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN.test(X2[, design[clen]], Y2[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[site,2] <- rAN$statistic
                        poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN.test(X3[, design[clen]], Y3[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[site,3] <- rAN$statistic
                        index <- index +1
                    }
                    # pooled = NA
                    else {
                        W[site, 1] <- W[site, 2] <- W[site, 3] <- NA
                    }
                }
            }
        } # end of for(Sites)
    } # end of if (n=4)
    else if (object@nSamples == 3) {
        print(paste("Generating the Within Adaptive test for sample size n = ", object@nSamples), sep = "")
        total <- length(object@coverage)
        # number of columns within W1, W2: (ab vs. ac, ab vs bc, ac vs. bc; sim. for plus)
        Within_cols <- 3
        Sites <- 1:total
        # Create a bin label for each sites from wSites
        binLabs <- c()
        for (bin in 1:length(object@wSites)) {
            sites <- object@wSites[[bin]]
            if (length(sites) > 0) {
                binLabs[sites] <- bin
            }
        }
        # Initilize W
        W <- matrix(NA, nrow = total , ncol = Within_cols)
        index <- 1 # index for W matrix
        for (site in Sites) {
            bin <- binLabs[site]
            if (minus_condition == TRUE) {
                pooled <- object@minusVar[[bin]]
                colnames(W) <- c('ab vs ac', 'ab vs bc', 'ac vs bc')
                X1 <- object@coverage[[site]][1:2,] #ab
                Y1 <- object@coverage[[site]][c(1,3),] #ac
                X2 <- object@coverage[[site]][1:2,] #ab
                Y2 <- object@coverage[[site]][2:3,] #bc
                X3 <- object@coverage[[site]][c(1,3),] #ac
                Y3 <- object@coverage[[site]][2:3,] #bc
            }
            else {
                pooled <- object@plusVar[[bin]]
                colnames(W) <- toupper(c('ab vs ac', 'ab vs bc', 'ac vs bc'))
                X1 <- object@coverage[[site]][4:5,]
                Y1 <- object@coverage[[site]][c(4,6),]
                X2 <- object@coverage[[site]][4:5,]
                Y2 <- object@coverage[[site]][5:6,]
                X3 <- object@coverage[[site]][c(4,6),]
                Y3 <- object@coverage[[site]][5:6,]
            }
            if (site %% 1000 == 0) {
                print(paste('W, site: ', site))
            }
            if ( dim(X1)[2] < object@s.size ) {
                if (use_cpp) {
                    # pooled is NOT NA
                    if (is.na(pooled) == FALSE) {
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[2]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index, 1] <- tan::AN_test(X1[, clen], Y1[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[3]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index, 2] <- tan::AN_test(X2[, clen], Y2[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[3]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index, 3] <- tan::AN_test(X3[, clen], Y3[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        index <- index + 1
                    }
                    # pooled is NA
                    else {
                        W[index, 1] <- W[index, 2] <- W[index, 3] <- NA
                        index <- index + 1
                    }
                }
                # if use_cpp = FALSE
                else {
                    if (is.na(pooled) == FALSE) {
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[2]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index,1] <- tan::AN.test(X1[, clen], Y1[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[3]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index,2] <- tan::AN.test(X2[, clen], Y2[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[3]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index,3] <- tan::AN.test(X3[, clen], Y3[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        index <- index + 1
                    }
                    # pooled is NA
                    else {
                        W[index, 1] <- W[index, 2] <- W[index, 3] <- NA
                        index <- index + 1
                    }
                }
            }
            # peakLength > s.size
            else {
                if (use_cpp) {
                    if (is.na(pooled) == FALSE) {
                        design <- object@Designs[site, ]
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[2]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN_test(X1[, design[clen]], Y1[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index,1] <- rAN$statistic
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[3]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN_test(X2[, design[clen]], Y2[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index,2] <- rAN$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[3]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN_test(X3[, design[clen]], Y3[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index,3] <- rAN$statistic
                        index <- index +1
                    }
                    # pooled = NA
                    else {
                        W[index, 1] <- W[index, 2] <- W[index, 3] <- NA
                        index <- index + 1
                    }
                }
                # use_cpp = FALSE
                else {
                    if (is.na(pooled) == FALSE) {
                        design <- object@Designs[site, ]
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[2]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN.test(X1[, design[clen]], Y1[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index,1] <- rAN$statistic
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[3]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN.test(X2[, design[clen]], Y2[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index,2] <- rAN$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[3]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN.test(X3[, design[clen]], Y3[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index,3] <- rAN$statistic
                        index <- index +1
                    }
                    # pooled = NA
                    else {
                        W[index, 1] <- W[index, 2] <- W[index, 3] <- NA
                        index <- index + 1
                    }
                }
            }
        } # end of for(Sites)
    }
    else if (object@nSamples == 2) {
        print(paste("Generating the Within Adaptive test for sample size n = ", object@nSamples), sep = "")
        total <- length(object@coverage)
        # Create a bin label for each sites from wSites
        binLabs <- c()
        for (bin in 1:length(object@wSites)) {
            sites <- object@wSites[[bin]]
            if (length(sites) > 0) {
                binLabs[sites] <- bin
            }
        }
        # Initilize W
        W <- matrix(NA, nrow = total , ncol = 3)
        index <- 1 # index for W matrix
        for (site in 1:total) {
            bin <- binLabs[site]
            pooled <- object@poolVar[[bin]]
            geta <- object@coverage[[site]][1,]
            getb <- object@coverage[[site]][2,]
            getA <- object@coverage[[site]][3,]
            getB <- object@coverage[[site]][4,]
            colnames(W) <- c('ab vs AB', 'aA vs bB', 'aB vs Ab')
            X1 <- rbind( geta, getb)
            Y1 <- rbind( getA, getB)
            X2 <- rbind( geta, getA)
            Y2 <- rbind( getb, getB)
            X3 <- rbind( geta, getB)
            Y3 <- rbind( getA, getb)
            if (site %% 1000 == 0) {
                print(paste('W, site: ', site))
            }
            if ( dim(X1)[2] < object@s.size ) {
                if (use_cpp) {
                    # pooled is NOT NA
                    if (is.na(pooled) == FALSE) {
                        # ab vs AB
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[4]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index, 1] <- tan::AN_test(X1[, clen], Y1[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        # aA vs bB
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index, 2] <- tan::AN_test(X2[, clen], Y2[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        # aB vs Ab
                        poolVarX <- pooled[[3]]; poolVarY <- pooled[[6]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index, 3] <- tan::AN_test(X3[, clen], Y3[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        index <- index + 1
                    }
                    # pooled is NA
                    else {
                        W[index, 1] <- W[index, 2] <- W[index, 3] <- NA
                        index <- index + 1
                    }
                }
                # if use_cpp = FALSE
                else {
                    if (is.na(pooled) == FALSE) {
                        # ab vs AB
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[4]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index, 1] <- tan::AN.test(X1[, clen], Y1[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        # aA vs bB
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index, 2] <- tan::AN.test(X2[, clen], Y2[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        # aB vs Ab
                        poolVarX <- pooled[[3]]; poolVarY <- pooled[[6]]
                        # clen <- 1:length(poolVarX)
                        clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                        W[index, 3] <- tan::AN.test(X3[, clen], Y3[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                        index <- index + 1
                    }
                    # pooled is NA
                    else {
                        W[index, 1] <- W[index, 2] <- W[index, 3] <- NA
                        index <- index + 1
                    }
                }
            }
            else {
                if (use_cpp) {
                    if (is.na(pooled) == FALSE) {
                        design <- object@Designs[site, ]
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[4]]
                        clen <- 1:length(poolVarX)
                        rAN <-tan::AN_test(X1[, design[clen]], Y1[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                           poolVarY=poolVarY)
                        W[index, 1] <- rAN$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN_test(X2[, design[clen]], Y2[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index, 2] <- rAN$statistic
                        poolVarX <- pooled[[3]]; poolVarY <- pooled[[6]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN_test(X3[, design[clen]], Y3[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index, 3] <- rAN$statistic
                        index <- index + 1
                    }
                    # pooled is NA
                    else {
                        W[index, 1] <- W[index, 2] <- W[index, 3] <- NA
                        index <- index + 1
                    }
                }
                else {
                    if (is.na(pooled) == FALSE) {
                        design <- object@Designs[site, ]
                        poolVarX <- pooled[[1]]; poolVarY <- pooled[[4]]
                        clen <- 1:length(poolVarX)
                        rAN <-tan::AN.test(X1[, design[clen]], Y1[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                           poolVarY=poolVarY)
                        W[index,1] <- rAN$statistic
                        poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN.test(X2[, design[clen]], Y2[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index,2] <- rAN$statistic
                        poolVarX <- pooled[[3]]; poolVarY <- pooled[[6]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN.test(X3[, design[clen]], Y3[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                                            poolVarY=poolVarY)
                        W[index,3] <- rAN$statistic
                        index <- index + 1
                    }
                    # pooled is NA
                    else {
                        W[index, 1] <- W[index, 2] <- W[index, 3] <- NA
                        index <- index + 1
                    }
                }
            }
        }
    } # end if (n=2)
    # return results
    if (object@nSamples == 2) {
        object@W <- W
    }
    else if (object@nSamples %in% c(3,4)) {
        if (minus_condition) {
            object@W1 <- W
        } else {
            object@W2 <- W
        }
    }
    object
}

setMethod("generateWithinTan", signature("tanDb"), .generateWithinTan)
