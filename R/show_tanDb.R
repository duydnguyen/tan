.show <- function(object) {
    cat("class: tanDb", "\n")
    cat(2*object@nSamples, " Samples, 2 conditions,",
        length(object@coverage), "regions", "\n")
}

setMethod("show", signature("tanDb"), .show)
