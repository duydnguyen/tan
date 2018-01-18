.show <- function(object) {
    cat("class: tanDb", "\n")
    cat(2*object@nSamples, paste(", 2 conditions, "),
        length(object@coverage), paste("regions"), "\n")
}

setMethod("show", signature("tanDb"), .show)
