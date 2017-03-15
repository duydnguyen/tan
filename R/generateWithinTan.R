.generateWithinTan <- function(object, minus_condition) {
    if (minus_condition == TRUE) {
        print("Generating the Within Adaptive test for first condition")
    } else {
        print("Generating the Within Adaptive test for second condition")
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
                pooled <- object@minusVar[[bin]]
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
            else {
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
        } # end of for(Sites)
    } # end of if (n=4)
    # return results
    if (minus_condition) {
        object@W1 <- W
    } else {
        object@W2 <- W
    }
    object
}

setMethod("generateWithinTan", signature("tanDb"), .generateWithinTan)
