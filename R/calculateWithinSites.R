.calculateWithinSites <- function(object, quantprobs) {
    wSites <- list()
    dN <- quantile(object@Ns, probs = quantprobs, na.rm = TRUE)
    ### nSamples = 4 ###
    if (object@nSamples == 4) {
        print(paste("Calculating within sites for sample size n = ", object@nSamples), sep = "")
        temp <- c('ab vs cd', 'ac vs bd', 'ad vs bc','AB vs CD', 'AC vs BD', 'AD vs BC' )
        ## Create colnames for matrix B
        minLabs <- c('ab','ac','ad','bc','bd','cd')
        plusLabs <- toupper(minLabs)
        labs <- ''
        for (mlab in minLabs) {
            for (plab in plusLabs) {
                labs <- c(labs, paste(mlab, 'vs', plab))
            }
        }
        labs <- labs[-c(1)]
        allCol <- c(labs, temp)
        # ncomps = ncol(AvsB) : number of columns in matB and Within
        ncomps <- length(allCol)
        if (ncol(object@Ns)>1){
            sampleNames <- unlist(strsplit(temp,' vs '))
            ## For each interval, take average of all ChIP samples "ac" "bd" "ad" "bc" "AC" "BD" "AD" "BC";
            ##'ac' means take average between a and c for that interval and so on
            H0.Ns <- rowMeans(object@Ns[,sampleNames])
            ## totColumn = ncol(Within) or length(temp)
            # totColumn <- length(temp)
            totColumn <- 1
            H0.Ns <-rep(H0.Ns, totColumn)
        } else {
            H0.Ns <-rep(object@Ns, totColumn)
        }
        if (max(dN) < max(object@Ns)){
            dN <- c(dN, max(object@Ns))
        }

        for (i in 1:length(dN)){
            # print(paste('+++Bin = ', i))
            if (i==1){
                if (ncol(object@Ns)==1){
                    idx <- which(object@Ns<=dN[[i]])
                }
                H0.idx <- which(H0.Ns<=dN[[i]])
                # print(paste('@@@H0.idx = ', H0.idx))
                wSites[[i]] <- H0.idx
            } else {
                if (ncol(object@Ns)==1){
                    idx <- which(object@Ns>dN[[i-1]] & object@Ns<=dN[[i]])
                }
                H0.idx <- which(H0.Ns>dN[[i-1]]&H0.Ns<=dN[[i]])
                # print(paste('@@@H0.idx = ', H0.idx))
                wSites[[i]] <- H0.idx
            }
        }
        ### Create a bin label for each sites from wSites
        binLabs <- c()
        for (bin in 1: length(wSites)) {
            sites <- wSites[[bin]]
            if (length(sites) > 0) {
                binLabs[sites] <- bin
            }
        }
    } # end of n=4
    else if (object@nSamples == 3) {
        print(paste("Calculating within sites for sample size n = ", object@nSamples), sep = "")
        temp <- c('ab vs ac', 'ab vs bc', 'ac vs bc', 'AB vs AC', 'AB vs BC', 'AC vs BC')
        ## Create colnames for matrix B
        minLabs <- c('ab','ac','bc')
        plusLabs <- toupper(minLabs)
        labs <- ''
        for (mlab in minLabs) {
            for (plab in plusLabs) {
                labs <- c(labs, paste(mlab, 'vs', plab))
            }
        }
        labs <- labs[-c(1)]
        allCol <- c(labs, temp)
        # ncomps = ncol(AvsB) : number of columns in matB and Within
        ncomps <- length(allCol)

        if (ncol(object@Ns)>1){
            sampleNames <- unlist(strsplit(temp,' vs '))
            ## For each interval, take average of all ChIP samples "ac" "bd" "ad" "bc" "AC" "BD" "AD" "BC";
            ##'ac' means take average between a and c for that interval and so on
            H0.Ns <- rowMeans(object@Ns[, sampleNames])
            ## totColumn = ncol(Within) or length(temp)
            # totColumn <- length(temp)
            totColumn <- 1
            H0.Ns <-rep(H0.Ns, totColumn)
        } else {
            H0.Ns <-rep(object@Ns,totColumn)
        }
        if (max(dN)<max(object@Ns)){
            dN <- c(dN, max(object@Ns))
        }

        for (i in 1:length(dN)){
            # print(paste('+++Bin = ', i))
            if (i==1){
                if (ncol(object@Ns)==1){
                    idx <- which(object@Ns<=dN[[i]])
                }
                H0.idx <- which(H0.Ns<=dN[[i]])
                # print(paste('@@@H0.idx = ', H0.idx))
                wSites[[i]] <- H0.idx
            } else {
                if (ncol(object@Ns)==1){
                    idx <- which(object@Ns>dN[[i-1]] & object@Ns<=dN[[i]])
                }
                H0.idx <- which(H0.Ns>dN[[i-1]]&H0.Ns<=dN[[i]])
                # print(paste('@@@H0.idx = ', H0.idx))
                wSites[[i]] <- H0.idx
            }
        }
        ### Create a bin label for each sites from wSites
        binLabs <- c()
        for (bin in 1: length(wSites)) {
            sites <- wSites[[bin]]
            if (length(sites) > 0) {
                binLabs[sites] <- bin
            }
        }
    }
    else if (object@nSamples == 2) {
        print(paste("Calculating within sites for sample size n = ", object@nSamples), sep = "")
        temp <- allCol <- c('ab vs AB')
        # ncomps = ncol(AvsB) : number of columns in matB and Within
        ncomps <- length(allCol)

        if (ncol(object@Ns)>1){
            sampleNames <- unlist(strsplit(temp,' vs '))
            ## For each interval (or row), take average of all ChIP samples
            H0.Ns <- rowMeans(object@Ns[,sampleNames])
            ## totColumn = ncol(Within) or length(temp)
            # totColumn <- length(temp)
            totColumn <- 1
            H0.Ns <-rep(H0.Ns, totColumn)
        } else {
            H0.Ns <-rep(object@Ns,totColumn)
        }
        if (max(dN)<max(object@Ns)){
            dN <- c(dN, max(object@Ns))
        }

        for (i in 1:length(dN)){
            # print(paste('+++Bin = ', i))
            if (i==1){
                if (ncol(object@Ns)==1){
                    idx <- which(object@Ns<=dN[[i]])
                }
                H0.idx <- which(H0.Ns<=dN[[i]])
                # print(paste('@@@H0.idx = ', H0.idx))
                wSites[[i]] <- H0.idx
            } else {
                if (ncol(object@Ns)==1){
                    idx <- which(object@Ns>dN[[i-1]] & object@Ns<=dN[[i]])
                }
                H0.idx <- which(H0.Ns>dN[[i-1]]&H0.Ns<=dN[[i]])
                # print(paste('@@@H0.idx = ', H0.idx))
                wSites[[i]] <- H0.idx
            }
        }
        ### Create a bin label for each sites from wSites
        binLabs <- c()
        for (bin in 1: length(wSites)) {
            sites <- wSites[[bin]]
            if (length(sites) > 0) {
                binLabs[sites] <- bin
            }
        }
    }
    # return results
    object@wSites <- wSites
    object@dN <- dN
    object
}

setMethod("calculateWithinSites", signature("tanDb"), .calculateWithinSites)
