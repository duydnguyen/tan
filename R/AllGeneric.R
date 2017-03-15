#' Create a matrix to store sampling locations for genomic regions.
#'
#' @param object A \code{tanDb} object
#' @param s.size The number of sampling points allowed in genomic regions.
#' @param LHD TRUE if using even spacing for sampling points; always set to TRUE.
#' @param Uniform TRUE if allowing uniform point within each grid of continuous sampling points;
#'    always set to FALSE
#'
#' @return Designs, a matrix (total regions x s.size) to store sampling locations for genomic regions
#' @export
#'
#' @examples
setGeneric("createDesigns", function(object, s.size, LHD, Uniform) {
    standardGeneric("createDesigns")
})

#' Calculate the matrix of total counts for each genomic regions.
#'
#' @param object A \code{tanDb} object
#' @param nSamples Sample sizes for each conditions (n1=n2=n)
#' @param bNormWidth Parameter for Binning TODO
#' @param bSampleMean Parameter for Binning TODO
#'
#' @return the matrix of total counts \code{Ns}
#' @export
#'
#' @examples
setGeneric("calculateTotalCounts", function(object, nSamples, bNormWidth, bSampleMean) {
    standardGeneric("calculateTotalCounts")
})

#' Get a list of within-site indices for each bin based
#'   on total counts \code{Ns} and \code{quantprobs}.
#'
#' @param object A \code{tanDb} object
#' @param quantprobs A quantile vector for binning.
#'
#' @return wSites, a list of within-site indices for each bin;
#'    and dN, a vector of count ranges for each bins.
#' @export
#'
#' @examples
setGeneric("calculateWithinSites", function(object, quantprobs) {
    standardGeneric("calculateWithinSites")
})

#' Calculuate pooled variance
#'
#' @param object A \code{tanDb} object
#' @param minus_condition TRUE if calculating the first (or minus ) condition;
#'    FALSE if calculating the second (or plus) condition.
#' @param Global_lower set lower bound for minGlobal (length of pooled var vector for each bins) TODO.
#'    Recommend using \code{Global_lower <- floor(s.size/2)}
#' @param poolQuant A pooled quantile at each genomic position.
#' @param movAve A number of points for moving average.
#'
#' @return  A list of pooled variances for the given condition
#' @export
#'
#' @examples
setGeneric("calculateVariance", function(object, minus_condition, Global_lower, poolQuant, movAve ) {
    standardGeneric("calculateVariance")
})


#' Generate the within (or null) Adaptive Neyman tests
#'
#' @param object A \code{tanDb} object
#' @param minus_condition TRUE if calculating the first (or minus ) condition;
#'    FALSE if calculating the second (or plus) condition.
#'
#' @return
#' @export
#'
#' @examples
setGeneric("generateWithinTan", function(object, minus_condition) {
    standardGeneric("generateWithinTan")
})

#' Compute p values for adaptive tests
#'
#' @param object A \code{tanDb} object
#' @param quant A quantile to obtain the combined p-values.
#' @param poolQuant A quantile to pool the variances; DEFAULT is median.
#' @param movAve A parameter to smooth the variance.
#'
#' @return
#' @export
#'
#' @examples
setGeneric("computePvalues", function(object, quant, poolQuant, movAve) {
    standardGeneric("computePvalues")
})

