evaluate_within <- function(coverage, nSamples, wSites, minus_condition, use_cpp = TRUE,
                           s.size, Designs, minusVar, plusVar) {
    ## Get dictionary {a:1, b:2, c:3, d:4,...}
    getDict <- function(s) {
                                        # only for lower case mapping
        return(match(tolower(s), letters))
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
    total <- length(coverage)
    withinLabel <- create_labels(nSamples)[[1]]
    # Within_cols <- 3 (for nSamples = 4: "ab vs cd" "ac vs bd" "ad vs bc")
    Within_cols <- floor(length(withinLabel) / 2)
    # Initilize Within matrix W
    W <- matrix(NA, nrow = total, ncol = Within_cols)
    # index <- 1 # index for W matrix
    # Create a bin label for each sites from wSites
    binLabs <- c()
    for (bin in 1:length(wSites)) {
        sites <- wSites[[bin]]
        if (length(sites) > 0) {
            binLabs[sites] <- bin
        }
    }
    for (site in 1:total) {
        bin <- binLabs[site]

        # testList <- list()
        withinX <- withinY <- list()
        numTests <-  c()
        indexList <- list()

        if (minus_condition == TRUE) {
            pooled <- minusVar[[bin]] # could returned NA if Bin's sites are "low quality"
            # colnames(W) <- c('ab vs cd', 'ac vs bd', 'ad vs bc' )
            colnames(W) <- withinLabel[1:Within_cols]
            X1 <- coverage[[site]][1:2,]
            Y1 <- coverage[[site]][3:4,]
            X2 <- coverage[[site]][c(1,3),]
            Y2 <- coverage[[site]][c(2,4),]
            X3 <- coverage[[site]][c(1,4),]
            Y3 <- coverage[[site]][c(2,3),]

            indexList <- create_indexList(nSamples)
            numTests <- length(indexList)
            for (tt in 1:numTests) {
                ids <- indexList[[tt]]
                withinX[[tt]] <- coverage[[site]][ids[1:2], ]
                withinY[[tt]] <- coverage[[site]][ids[3:4],]
            }

        }
        else {
            pooled <- plusVar[[bin]]
            # colnames(W) <- toupper(c('ab vs cd', 'ac vs bd', 'ad vs bc' ))
            colnames(W) <- withinLabel[(Within_cols + 1):(2 * Within_cols)]
            X1 <- coverage[[site]][5:6,]
            Y1 <- coverage[[site]][7:8,]
            X2 <- coverage[[site]][c(5,7),]
            Y2 <- coverage[[site]][c(6,8),]
            X3 <- coverage[[site]][c(5,8),]
            Y3 <- coverage[[site]][c(6,7),]

            indexList <- create_indexList(nSamples)
            numTests <- length(indexList)
            for (tt in 1:numTests) {
                ids <- indexList[[tt]] + nSamples
                withinX[[tt]] <- coverage[[site]][ids[1:2], ]
                withinY[[tt]] <- coverage[[site]][ids[3:4],]
            }

        }
        if (site %% 1000 == 0) {
            print(paste('W, site: ', site))
        }
        ## if ( dim(X1)[2] < s.size ) {

        if ( dim(withinX[[1]])[2] < s.size ) {

            if (use_cpp) {
                # pooled is NOT NA
                if (is.na(pooled)[1] == FALSE) {
                    ## poolVarX <- pooled[['ab']]; poolVarY <- pooled[['cd']]
                    ## # clen <- 1:length(poolVarX)
                    ## clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                    ## W[site,1] <- tan::AN_test(X1[, clen], Y1[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                           poolVarY = poolVarY)$statistic
                    ## poolVarX <- pooled[['ac']]; poolVarY <- pooled[['bd']]
                    ## # clen <- 1:length(poolVarX)
                    ## clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                    ## W[site,2] <- tan::AN_test(X2[, clen], Y2[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                           poolVarY = poolVarY)$statistic
                    ## poolVarX <- pooled[['ad']]; poolVarY <- pooled[['bc']]
                    ## # clen <- 1:length(poolVarX)
                    ## clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                    ## W[site,3] <- tan::AN_test(X3[, clen], Y3[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                           poolVarY = poolVarY)$statistic
                    ## # index <- index + 1

                    for (tt in 1:numTests) {
                        X <- withinX[[tt]]
                        Y <- withinY[[tt]]
                        pairNames <- tolower(unlist(strsplit(colnames(W)[tt], ' vs '))) # e.g, c('ab', 'cd')
                        poolVarX <- pooled[[pairNames[1]]]; poolVarY <- pooled[[pairNames[2]]]
                        clen <- 1:min(length(poolVarX), dim(X)[2]) # modified on 04/12/16
                        W[site, tt] <- tan::AN_test(X[, clen], Y[, clen], na_rm=TRUE,
                                                   pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                    }

                }
                # pooled is NA
                else {
                    ## W[site, 1] <- W[site, 2] <- W[site, 3] <- NA

                    W[site, 1:numTests] <- rep(NA, numTests)

                }
            }
            # if use_cpp = FALSE
            else {
                if (is.na(pooled)[1] == FALSE) {
                    ## poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
                    ## # clen <- 1:length(poolVarX)
                    ## clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                    ## W[site,1] <- tan::AN.test(X1[, clen], Y1[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                           poolVarY = poolVarY)$statistic
                    ## poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                    ## # clen <- 1:length(poolVarX)
                    ## clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                    ## W[site,2] <- tan::AN.test(X2[, clen], Y2[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                           poolVarY = poolVarY)$statistic
                    ## poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
                    ## # clen <- 1:length(poolVarX)
                    ## clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
                    ## W[site,3] <- tan::AN.test(X3[, clen], Y3[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                           poolVarY = poolVarY)$statistic
                    ## # index <- index + 1

                    for (tt in 1:numTests) {
                        X <- withinX[[tt]]
                        Y <- withinY[[tt]]
                        pairNames <- tolower(unlist(strsplit(colnames(W)[tt], ' vs '))) # e.g, c('ab', 'cd')
                        poolVarX <- pooled[[pairNames[1]]]; poolVarY <- pooled[[pairNames[2]]]
                        clen <- 1:min(length(poolVarX), dim(X)[2]) # modified on 04/12/16
                        W[site, tt] <- tan::AN.test(X[, clen], Y[, clen], na.rm=TRUE,
                                                   pool=TRUE, poolVarX = poolVarX,
                                                   poolVarY = poolVarY)$statistic
                    }

                }
                # pooled is NA
                else {
                    ## W[site, 1] <- W[site, 2] <- W[site, 3] <- NA

                    W[site, 1:numTests] <- rep(NA, numTests)

                }
            }
        }
        # peakLength > s.size
        else {
            if (use_cpp) {
                if (is.na(pooled)[1] == FALSE) {
                    ## design <- Designs[site, ]
                    ## poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
                    ## clen <- 1:length(poolVarX)
                    ## rAN <- tan::AN_test(X1[, design[clen]], Y1[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                     poolVarY=poolVarY)
                    ## W[site,1] <- rAN$statistic
                    ## poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                    ## clen <- 1:length(poolVarX)
                    ## rAN <- tan::AN_test(X2[, design[clen]], Y2[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                     poolVarY=poolVarY)
                    ## W[site,2] <- rAN$statistic
                    ## poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
                    ## clen <- 1:length(poolVarX)
                    ## rAN <- tan::AN_test(X3[, design[clen]], Y3[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                     poolVarY=poolVarY)
                    ## W[site,3] <- rAN$statistic
                    ## # index <- index +1

                    for (tt in 1:numTests) {
                        design <- Designs[site, ]
                        X <- withinX[[tt]]
                        Y <- withinY[[tt]]
                        pairNames <- tolower(unlist(strsplit(colnames(W)[tt], ' vs '))) # e.g, c('ab', 'cd')
                        poolVarX <- pooled[[pairNames[1]]]; poolVarY <- pooled[[pairNames[2]]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN_test(X[, design[clen]], Y[, design[clen]], na_rm=TRUE,
                                           pool=TRUE, poolVarX = poolVarX,
                                           poolVarY = poolVarY)
                        W[site, tt] <- rAN$statistic
                    }

                }
                # pooled = NA
                else {
                    ## W[site, 1] <- W[site, 2] <- W[site, 3] <- NA

                    W[site, 1:numTests] <- rep(NA, numTests)

                }
            }
            # use_cpp = FALSE
            else {
                if (is.na(pooled)[1] == FALSE) {
                    ## design <- Designs[site, ]
                    ## poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
                    ## clen <- 1:length(poolVarX)
                    ## rAN <- tan::AN.test(X1[, design[clen]], Y1[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                     poolVarY=poolVarY)
                    ## W[site,1] <- rAN$statistic
                    ## poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
                    ## clen <- 1:length(poolVarX)
                    ## rAN <- tan::AN.test(X2[, design[clen]], Y2[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                     poolVarY=poolVarY)
                    ## W[site,2] <- rAN$statistic
                    ## poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
                    ## clen <- 1:length(poolVarX)
                    ## rAN <- tan::AN.test(X3[, design[clen]], Y3[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
                    ##                     poolVarY=poolVarY)
                    ## W[site,3] <- rAN$statistic
                    ## # index <- index +1

                    for (tt in 1:numTests) {
                        design <- Designs[site, ]
                        X <- withinX[[tt]]
                        Y <- withinY[[tt]]
                        pairNames <- tolower(unlist(strsplit(colnames(W)[tt], ' vs '))) # e.g, c('ab', 'cd')
                        poolVarX <- pooled[[pairNames[1]]]; poolVarY <- pooled[[pairNames[2]]]
                        clen <- 1:length(poolVarX)
                        rAN <- tan::AN.test(X[, design[clen]], Y[, design[clen]], na.rm=TRUE,
                                           pool=TRUE, poolVarX = poolVarX,
                                           poolVarY = poolVarY)
                        W[site, tt] <- rAN$statistic
                    }

                }
                # pooled = NA
                else {
                    ## W[site, 1] <- W[site, 2] <- W[site, 3] <- NA

                    W[site, 1:numTests] <- rep(NA, numTests)

                }
            }
        }
    } # end of for(Sites)
    return(W)
}

.generateWithinTan <- function(object, minus_condition, ...) {
    if (object@nSamples > 2) {
        if (minus_condition == TRUE) {
            print("Generating the Within Adaptive test for first condition")
        } else {
            print("Generating the Within Adaptive test for second condition")
        }
    }
    else if (object@nSamples == 2) {
        print("Calculating the Within Adaptive tests for both conditions")
    }
    ### MAIN ###
    W <- matrix()
    if (object@nSamples > 3) {
        print(paste("Generating the Within Adaptive test for sample size n = ", object@nSamples), sep = "")
        # total <- length(object@coverage)
        # withinLabel <- create_labels(object@nSamples)[[1]]
        ## number of columns within W1, W2: (ab vs. cd, ac vs. bd, ad vs. bc; sim. for plus)
        #Within_cols <- 3
        #Within_cols <- floor(length(withinLabel) / 2)
        #Sites <- 1:total

        W <-  evaluate_within(coverage = object@coverage, nSamples = object@nSamples,
                              wSites = object@wSites, minus_condition = minus_condition,
                              use_cpp = use_cpp, s.size = object@s.size,
                              Designs = object@Designs, minusVar = object@minusVar,
                              plusVar = object@plusVar)


        ## # Create a bin label for each sites from wSites
        ## binLabs <- c()
        ## for (bin in 1:length(object@wSites)) {
        ##     sites <- object@wSites[[bin]]
        ##     if (length(sites) > 0) {
        ##         binLabs[sites] <- bin
        ##     }
        ## }
        ## # Initilize W
        ## W <- matrix(NA, nrow = total, ncol = Within_cols)
        ## index <- 1 # index for W matrix
        ## for (site in Sites) {
        ##     bin <- binLabs[site]
        ##     if (minus_condition == TRUE) {
        ##         pooled <- object@minusVar[[bin]] # could returned NA if Bin's sites are "low quality"
        ##         colnames(W) <- c('ab vs cd', 'ac vs bd', 'ad vs bc' )
        ##         X1 <- object@coverage[[site]][1:2,]
        ##         Y1 <- object@coverage[[site]][3:4,]
        ##         X2 <- object@coverage[[site]][c(1,3),]
        ##         Y2 <- object@coverage[[site]][c(2,4),]
        ##         X3 <- object@coverage[[site]][c(1,4),]
        ##         Y3 <- object@coverage[[site]][c(2,3),]
        ##     }
        ##     else {
        ##         pooled <- object@plusVar[[bin]]
        ##         colnames(W) <- toupper(c('ab vs cd', 'ac vs bd', 'ad vs bc' ))
        ##         X1 <- object@coverage[[site]][5:6,]
        ##         Y1 <- object@coverage[[site]][7:8,]
        ##         X2 <- object@coverage[[site]][c(5,7),]
        ##         Y2 <- object@coverage[[site]][c(6,8),]
        ##         X3 <- object@coverage[[site]][c(5,8),]
        ##         Y3 <- object@coverage[[site]][c(6,7),]
        ##     }
        ##     if (site %% 1000 == 0) {
        ##         print(paste('W, site: ', site))
        ##     }
        ##     if ( dim(X1)[2] < object@s.size ) {
        ##         if (use_cpp) {
        ##             # pooled is NOT NA
        ##             if (is.na(pooled)[1] == FALSE) {
        ##                 poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
        ##                 # clen <- 1:length(poolVarX)
        ##                 clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
        ##                 W[site,1] <- tan::AN_test(X1[, clen], Y1[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                           poolVarY = poolVarY)$statistic
        ##                 poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
        ##                 # clen <- 1:length(poolVarX)
        ##                 clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
        ##                 W[site,2] <- tan::AN_test(X2[, clen], Y2[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                           poolVarY = poolVarY)$statistic
        ##                 poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
        ##                 # clen <- 1:length(poolVarX)
        ##                 clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
        ##                 W[site,3] <- tan::AN_test(X3[, clen], Y3[, clen], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                           poolVarY = poolVarY)$statistic
        ##                 index <- index + 1
        ##             }
        ##             # pooled is NA
        ##             else {
        ##                 W[site, 1] <- W[site, 2] <- W[site, 3] <- NA
        ##             }
        ##         }
        ##         # if use_cpp = FALSE
        ##         else {
        ##             if (is.na(pooled)[1] == FALSE) {
        ##                 poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
        ##                 # clen <- 1:length(poolVarX)
        ##                 clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
        ##                 W[site,1] <- tan::AN.test(X1[, clen], Y1[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                           poolVarY = poolVarY)$statistic
        ##                 poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
        ##                 # clen <- 1:length(poolVarX)
        ##                 clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
        ##                 W[site,2] <- tan::AN.test(X2[, clen], Y2[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                           poolVarY = poolVarY)$statistic
        ##                 poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
        ##                 # clen <- 1:length(poolVarX)
        ##                 clen <- 1:min(length(poolVarX), dim(X1)[2]) # modified on 04/12/16
        ##                 W[site,3] <- tan::AN.test(X3[, clen], Y3[, clen], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                           poolVarY = poolVarY)$statistic
        ##                 index <- index + 1
        ##             }
        ##             # pooled is NA
        ##             else {
        ##                 W[site, 1] <- W[site, 2] <- W[site, 3] <- NA
        ##             }
        ##         }
        ##     }
        ##     # peakLength > s.size
        ##     else {
        ##         if (use_cpp) {
        ##             if (is.na(pooled)[1] == FALSE) {
        ##                 design <- object@Designs[site, ]
        ##                 poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
        ##                 clen <- 1:length(poolVarX)
        ##                 rAN <- tan::AN_test(X1[, design[clen]], Y1[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                     poolVarY=poolVarY)
        ##                 W[site,1] <- rAN$statistic
        ##                 poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
        ##                 clen <- 1:length(poolVarX)
        ##                 rAN <- tan::AN_test(X2[, design[clen]], Y2[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                     poolVarY=poolVarY)
        ##                 W[site,2] <- rAN$statistic
        ##                 poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
        ##                 clen <- 1:length(poolVarX)
        ##                 rAN <- tan::AN_test(X3[, design[clen]], Y3[, design[clen]], na_rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                     poolVarY=poolVarY)
        ##                 W[site,3] <- rAN$statistic
        ##                 index <- index +1
        ##             }
        ##             # pooled = NA
        ##             else {
        ##                 W[site, 1] <- W[site, 2] <- W[site, 3] <- NA
        ##             }
        ##         }
        ##         # use_cpp = FALSE
        ##         else {
        ##             if (is.na(pooled)[1] == FALSE) {
        ##                 design <- object@Designs[site, ]
        ##                 poolVarX <- pooled[[1]]; poolVarY <- pooled[[6]]
        ##                 clen <- 1:length(poolVarX)
        ##                 rAN <- tan::AN.test(X1[, design[clen]], Y1[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                     poolVarY=poolVarY)
        ##                 W[site,1] <- rAN$statistic
        ##                 poolVarX <- pooled[[2]]; poolVarY <- pooled[[5]]
        ##                 clen <- 1:length(poolVarX)
        ##                 rAN <- tan::AN.test(X2[, design[clen]], Y2[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                     poolVarY=poolVarY)
        ##                 W[site,2] <- rAN$statistic
        ##                 poolVarX <- pooled[[3]]; poolVarY <- pooled[[4]]
        ##                 clen <- 1:length(poolVarX)
        ##                 rAN <- tan::AN.test(X3[, design[clen]], Y3[, design[clen]], na.rm=TRUE, pool=TRUE, poolVarX = poolVarX,
        ##                                     poolVarY=poolVarY)
        ##                 W[site,3] <- rAN$statistic
        ##                 index <- index +1
        ##             }
        ##             # pooled = NA
        ##             else {
        ##                 W[site, 1] <- W[site, 2] <- W[site, 3] <- NA
        ##             }
        ##         }
        ##     }
        ## } # end of for(Sites)
    } # end of if (n > 3)
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
                # colnames(W) <- toupper(c('ab vs ac', 'ab vs bc', 'ac vs bc'))
                colnames(W) <- c('AB vs AC', 'AB vs BC', 'AC vs BC')
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
                    if (is.na(pooled)[1] == FALSE) {
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
                    if (is.na(pooled)[1] == FALSE) {
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
                    if (is.na(pooled)[1] == FALSE) {
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
                    if (is.na(pooled)[1] == FALSE) {
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
                    if (is.na(pooled)[1] == FALSE) {
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
                    if (is.na(pooled)[1] == FALSE) {
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
                    if (is.na(pooled)[1] == FALSE) {
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
                    if (is.na(pooled)[1] == FALSE) {
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
    else if (object@nSamples > 2) {
        if (minus_condition) {
            object@W1 <- W
        } else {
            object@W2 <- W
        }
    }
    object
}

setMethod("generateWithinTan", signature("tanDb"), .generateWithinTan)
