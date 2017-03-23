.createDesigns <- function(object, s.size, LHD, Uniform) {
    total <- length(object@coverage)
    Designs <- matrix(NA, nrow = total, ncol = s.size)
    for (site in 1:total) {
        if (site %%1000 == 0) {
            print(paste(site,' out of', total))
        }
        peakWidth <- length(object@coverage[[site]][1,])
        if (s.size <= peakWidth) {
            Designs[site, ] <- tan::Sampling(s.size = s.size, start = 1 , end = peakWidth, LHD = LHD, Uniform = Uniform)
        } else {
            message("s.size > width of genomic peak: No Sampling required! ")
            Designs[site, ] <- tan::Sampling(s.size = s.size, start = 1 , end = peakWidth, LHD = LHD, Uniform = Uniform)
        }

    }
    object@Designs <- Designs
    object@s.size <- s.size
    object
}

setMethod("createDesigns", signature("tanDb"), .createDesigns)

## or write function directly to setMethod()
# setMethod("createDesigns", signature("tanDb"),
#           function(object, s.size, LHD, Uniform) {
#             total <- length(object@coverage)
#             Designs <- matrix(NA, nrow = total, ncol = s.size)
#             for (site in 1:total) {
#                 if (site %%1000 == 0) {
#                     print(paste(site,' out of', total))
#                 }
#                 Designs[site, ] <- tan::Sampling(s.size = s.size, start = 1 , end = length(object@coverage[[site]][1,]), LHD = LHD, Uniform = Uniform)
#             }
#             object@Designs <- Designs
#             object
# })
