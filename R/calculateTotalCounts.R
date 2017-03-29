# n = 4. Ns_cols is a function of nSamples
.calculateTotalCounts <- function(object, nSamples, bNormWidth, bSampleMean) {
    ## parallel backend
    total <- length(object@coverage)
    Ns_cols <- c()
    n_cores <- parallel::detectCores()
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    sizeGrid <- floor(total / n_cores)
    Ncount <- matrix()
    ### nSamples = 4 ###
    if (nSamples == 4) {
        print(paste("Calculating Total Counts for sample size n = ", nSamples), sep = "")
        Ns_cols <- 12
        Ncount  <- foreach::foreach(grid = c(1:n_cores), .combine = rbind, .multicombine = T) %do% {
            if (grid != n_cores) {
                start <- (grid-1)*sizeGrid + 1
                end <- grid*sizeGrid
                Sites <- start:end
            } else {
                start <- (grid-1)*sizeGrid + 1
                end <- total
                Sites <- start:end
            }
            Counts <- matrix(NA, nrow = length(Sites), ncol = nSamples * 2)
            index <- 1 # index for Counts
            for (site in Sites) {
                if (site %% 1000 == 0) {
                    print(paste(site,' out of', total))
                }
                # NA values:
                allNA <- is.na(object@coverage[[site]])
                idNA <- c()
                for (i in 1:(nSamples * 2)) {
                    if (sum(allNA[i,]) > 0) {
                        id <- which(allNA[i,]==TRUE)
                        idNA <- union(idNA, id)
                        print(paste('+++length of NAs =',length(idNA)))
                    }
                }
                # dim X <- nxT (number of reps x width of interval)
                X <- object@coverage[[site]][1:4,]
                X[, idNA] <- 0
                # dim Y <- nxT (number of reps x width of interval)
                Y <- object@coverage[[site]][5:8,]
                Y[, idNA] <- 0
                if ( dim(X)[2] < object@s.size) {
                    XY <- rbind(X,Y)
                    if (bNormWidth) {
                        Counts[index, ] <- apply(XY, 1, sum)/dim(X)[2]
                        index <- index + 1
                    } else {
                        Counts[index, ] <- apply(XY, 1, sum)
                        index <- index + 1
                    }
                }
                else {
                    # Get Latin Hypercube design sampling:
                    design <- object@Designs[site, ]
                    XY <- rbind(X,Y)
                    if (bNormWidth) {
                        Counts[index, ] <- apply(XY[, design], 1, sum)/ object@s.size
                        index <- index + 1
                    } else {
                        Counts[index, ] <-apply(XY[, design], 1, sum)
                        index <- index + 1
                    }
                }
                if (bSampleMean) {
                    Ns <- as.matrix(rowMeans(Counts),ncol=1,nrow=nrow(Counts))
                } else {
                    Ns <- Counts
                    colnames(Ns)= c('a','b','c','d','A','B','C','D')
                }
            } # end 'for' of site
            return(Ns)
        }

    } # end if (nSamples == 4)
    else if (nSamples == 3) {
        print(paste("Calculating Total Counts for sample size n = ", nSamples), sep = "")
        # number of columns for Ns =  (ab,ac,bc; similar for plus condition)
        Ns_cols <- 3*2
        Ncount  <- foreach::foreach(grid = c(1:n_cores), .combine = rbind, .multicombine = T) %do% {
            if (grid != n_cores) {
                start <- (grid-1)*sizeGrid + 1
                end <- grid*sizeGrid
                Sites <- start:end
            } else {
                start <- (grid-1)*sizeGrid + 1
                end <- total
                Sites <- start:end
            }
            Counts <- matrix(NA, nrow = length(Sites), ncol = nSamples*2)
            index <- 1 # index for Counts
            for (site in Sites) {
                if (site %% 1000 == 0) {
                    print(paste(site,' out of', total))
                }
                # NA values:
                allNA <- is.na(object@coverage[[site]])
                idNA <- c()
                for (i in 1:(nSamples * 2)) {
                    if (sum(allNA[i,]) > 0) {
                        id <- which(allNA[i,]==TRUE)
                        idNA <- union(idNA, id)
                        print(paste('+++length of NAs =',length(idNA)))
                    }
                }
                X <- object@coverage[[site]][1:3,]
                X[,idNA] <- 0
                #Y <- object@coverage[[site]][5:7, ]
                Y <- object@coverage[[site]][4:6, ]
                Y[, idNA] <- 0
                if ( dim(X)[2] < object@s.size) {
                    XY <- rbind(X,Y)
                    if (bNormWidth) {
                        Counts[index, ] <- apply(XY, 1, sum)/dim(X)[2]
                        index <- index + 1
                    } else {
                        Counts[index, ] <- apply(XY, 1, sum)
                        index <- index + 1
                    }
                }
                else {
                    # Get Latin Hypercube design sampling:
                    design <- object@Designs[site, ]
                    XY <- rbind(X,Y)
                    if (bNormWidth) {
                        Counts[index, ] <- apply(XY[, design], 1, sum)/ s.size
                        index <- index + 1
                    } else {
                        Counts[index, ] <-apply(XY[, design], 1, sum)
                        index <- index + 1
                    }
                }
                if (bSampleMean) {
                    Ns <- as.matrix(rowMeans(Counts),ncol=1,nrow=nrow(Counts))
                } else {
                    Ns <- Counts
                    colnames(Ns)= c('a','b','c','A','B','C')
                }
            } # end 'for' of site
            return(Ns)
        }
    } # end if (nSamples == 3)
    else if (nSamples == 2) {
        print(paste("Calculating Total Counts for sample size n = ", nSamples), sep = "")
        # number of columns for Ns =  (ab, AB)
        Ns_cols <- 1*2
        Ncount  <- foreach::foreach(grid = c(1:n_cores), .combine = rbind, .multicombine = T) %do% {
            if (grid != n_cores) {
                start <- (grid-1)*sizeGrid + 1
                end <- grid*sizeGrid
                Sites <- start:end
            } else {
                start <- (grid-1)*sizeGrid + 1
                end <- total
                Sites <- start:end
            }
            Counts <- matrix(NA, nrow = length(Sites), ncol = nSamples*2)
            index <- 1 # index for Counts
            for (site in Sites) {
                if (site %% 1000 == 0) {
                    print(paste(site,' out of', total))
                }
                # NA values:
                allNA <- is.na(object@coverage[[site]])
                idNA <- c()
                for (i in 1:(nSamples * 2)) {
                    if (sum(allNA[i,])>0) {
                        id <- which(allNA[i,]==TRUE)
                        idNA <- union(idNA, id)
                        print(paste('+++length of NAs =',length(idNA)))
                    }
                }
                X <- object@coverage[[site]][1:2,]
                X[,idNA] <- 0
                # dim Y <- nxT (number of reps x width of interval)
                Y <- object@coverage[[site]][3:4,]
                Y[, idNA] <- 0
                if ( dim(X)[2] < object@s.size) {
                    XY <- rbind(X,Y)
                    if (bNormWidth) {
                        Counts[index, ] <- apply(XY, 1, sum)/dim(X)[2]
                        index <- index + 1
                    } else {
                        Counts[index, ] <- apply(XY, 1, sum)
                        index <- index + 1
                    }
                }
                else {
                    # Get Latin Hypercube design sampling:
                    design <- object@Designs[site, ]
                    XY <- rbind(X,Y)
                    if (bNormWidth) {
                        Counts[index, ] <- apply(XY[, design], 1, sum)/ s.size
                        index <- index + 1
                    } else {
                        Counts[index, ] <-apply(XY[, design], 1, sum)
                        index <- index + 1
                    }
                }
                if (bSampleMean) {
                    Ns <- as.matrix(rowMeans(Counts),ncol=1,nrow=nrow(Counts))
                } else {
                    Ns <- Counts
                    colnames(Ns)= c('a','b','A','B')
                }
            } # end 'for' of site
            return(Ns)
        }
    }
    parallel::stopCluster(cl)
    # Add labels for Ns
    Ns <- matrix(NA, nrow = total, ncol = Ns_cols)
    #colnames(Ncount) <- c(letters[1:4], toupper(letters[1:4]))
    colnames(Ncount) <- c(letters[1:nSamples], toupper(letters[1:nSamples]))
    if (nSamples == 4) {
        minLabs <- c('ab','ac','ad','bc','bd','cd')
        plusLabs <- toupper(minLabs)
        colnames(Ns) <- c(minLabs, plusLabs)
        for (st in colnames(Ns)) {
            Ns[, st] <- rowMeans(Ncount[,c( substr(st,1,1), substr(st,2,2))])
        }
    } else if (nSamples == 3) {
        minLabs <- c('ab','ac','bc')
        plusLabs <- toupper(minLabs)
        colnames(Ns) <- c(minLabs, plusLabs)
        for (st in colnames(Ns)) {
            Ns[, st] <- rowMeans(Ncount[,c( substr(st,1,1), substr(st,2,2))])
        }
    } else if (nSamples == 2) {
        minLabs <- c('ab')
        plusLabs <- toupper(minLabs)
        colnames(Ns) <- c(minLabs, plusLabs)
        for (st in colnames(Ns)) {
            Ns[, st] <- rowMeans(Ncount[,c( substr(st,1,1), substr(st,2,2))])
        }
    }

    object@Ns <- Ns
    object@nSamples <- nSamples
    object
}

setMethod("calculateTotalCounts", signature("tanDb"), .calculateTotalCounts)
