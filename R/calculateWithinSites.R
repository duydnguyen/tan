.calculateWithinSites <- function(object, quantprobs) {
    ## create lables for pair samples for within, and between
    create_labels <- function(nSamples, Ns) {
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
        ## ## Create colnames for Between matrix
        ## minLabs <- colnames(Ns)[1:(dim(Ns)[2] / 2)]
        ## plusLabs <- toupper(minLabs)
        ## labs <- ''
        ## for (mlab in minLabs) {
        ##     for (plab in plusLabs) {
        ##         labs <- c(labs, paste(mlab, 'vs', plab))
        ##     }
        ## }
        ## labs <- labs[-c(1)]
        ## allCol <- c(labs, temp)
        return(list('withinLabel' = temp))
    }

    ## create sites for each quantile bins under H0
    createWithinSites <- function(Ns, dN, withinLabels) {
        wSites <- list()
        if (ncol(Ns)>1){
            sampleNames <- unlist(strsplit(withinLabels,' vs '))
            ## For each interval, take average of all ChIP samples "ac" "bd" "ad" "bc" "AC" "BD" "AD" "BC";
            ##'ac' means take average between a and c for that interval and so on
            H0.Ns <- rowMeans(Ns[, sampleNames])
            ## totColumn = ncol(Within) or length(temp)
            # totColumn <- length(temp)
            totColumn <- 1
            H0.Ns <- rep(H0.Ns, totColumn)
        } else {
            H0.Ns <- rep(Ns, totColumn)
        }
        if (max(dN) < max(Ns)){
            dN <- c(dN, max(Ns))
        }
        for (i in 1:length(dN)){
            if (i==1){
                if (ncol(Ns)==1){
                    idx <- which(Ns<=dN[[i]])
                }
                H0.idx <- which(H0.Ns<=dN[[i]])
                wSites[[i]] <- H0.idx
            } else {
                if (ncol(Ns)==1){
                    idx <- which(Ns>dN[[i-1]] & Ns<=dN[[i]])
                }
                H0.idx <- which(H0.Ns>dN[[i-1]] & H0.Ns<=dN[[i]])
                wSites[[i]] <- H0.idx
            }
        }
        return(wSites)
    }
    ### MAIN ###
    wSites <- list()
    dN <- quantile(object@Ns, probs = quantprobs, na.rm = TRUE)
    ### nSamples > 3 ###
    if (object@nSamples > 3) {
        print(paste("Calculating within sites for sample size n = ", object@nSamples), sep = "")
        allColList <- create_labels(object@nSamples, object@Ns)
        temp <- allColList[['withinLabel']]
        wSites <- createWithinSites(Ns = object@Ns, dN = dN, withinLabels = temp)
    }
    else if (object@nSamples == 3) {
        print(paste("Calculating within sites for sample size n = ", object@nSamples), sep = "")
        temp <- c('ab vs ac', 'ab vs bc', 'ac vs bc', 'AB vs AC', 'AB vs BC', 'AC vs BC')
        wSites <- createWithinSites(Ns = object@Ns, dN = dN, withinLabels = temp)
    }
    else if (object@nSamples == 2) {
        print(paste("Calculating within sites for sample size n = ", object@nSamples), sep = "")
        temp <- allCol <- c('ab vs AB')
        wSites <- createWithinSites(Ns = object@Ns, dN = dN, withinLabels = temp)
    }
    # return results
    object@wSites <- wSites
    object@dN <- dN
    object
}

setMethod("calculateWithinSites", signature("tanDb"), .calculateWithinSites)
