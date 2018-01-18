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
#' @param use_cpp TRUE if the implementation uses cpp functions;
#' FALSE if the implementation uses R functions. R version is under development.
#'
#' @return  A list of pooled variances for the given condition.
#' @export
#'
#' @examples
setGeneric("calculateVariance", function(object, minus_condition, Global_lower, poolQuant, movAve, use_cpp = TRUE, ... ) {
    standardGeneric("calculateVariance")
})


#' Generate the within (or null) Adaptive Neyman tests
#'
#' @param object A \code{tanDb} object
#' @param minus_condition TRUE if calculating the first (or minus ) condition;
#'    FALSE if calculating the second (or plus) condition.
#' @param use_cpp TRUE if the implementation uses cpp functions;
#' FALSE if the implementation uses R functions. R version is under development.
#'
#' @return
#' @export
#'
#' @examples
setGeneric("generateWithinTan", function(object, minus_condition, use_cpp = TRUE, ...) {
    standardGeneric("generateWithinTan")
})

#' Compute p values for adaptive tests
#'
#' @param object A \code{tanDb} object
#' @param quant A quantile to obtain the combined p-values.
#' @param poolQuant A quantile to pool the variances; DEFAULT is median.
#' @param movAve A parameter to smooth the variance.
#' @param Global_lower set lower bound for minGlobal (length of pooled var vector for each bins)
#'    Recommend using \code{Global_lower <- floor(s.size/2)}
#' @param use_cpp TRUE if the implementation uses cpp functions;
#' FALSE if the implementation uses R functions. R version is under development.
#' @param ignore_sitesUnused FALSE by default if considering sitesUnused when computing p-values of testing between;
#' recommend using TRUE since sitesUnused caused pvalues = 0. #TODO: clean up
#' @param na_impute Impute missing value of \code{p-values} and \code{FDR} matrices by taking global mean for all p;TRUE by default. Recommend using FALSE
#'
#' @return
#' @export
#'
#' @examples
setGeneric("computePvalues", function(object, quant, poolQuant, movAve, Global_lower, use_cpp = TRUE, ignore_sitesUnused = FALSE, na_impute = TRUE, ...) {
    standardGeneric("computePvalues")
})

#' Compute p values for adaptive tests in batch mode
#'
#' @param object A \code{tanDb} object
#' @param quant A quantile to obtain the combined p-values.
#' @param poolQuant A quantile to pool the variances; DEFAULT is median.
#' @param movAve A parameter to smooth the variance.
#' @param Global_lower set lower bound for minGlobal (length of pooled var vector for each bins)
#'    Recommend using \code{Global_lower <- floor(s.size/2)}
#' @param use_cpp TRUE if the implementation uses cpp functions;
#' FALSE if the implementation uses R functions. R version is under development.
#' @param ignore_sitesUnused FALSE by default if considering sitesUnused when computing p-values of testing between;
#' recommend using TRUE since sitesUnused caused pvalues = 0. #TODO: clean up
#' @param na_impute Impute missing value of \code{p-values} and \code{FDR} matrices by taking global mean for all p;TRUE by default. Recommend using FALSE
#' @param bins A vector stores bin indices to be run in batch mode.
#' @param create_pMat A logical value. TRUE to generate p-value matrix \code{pMat} by combining p.list in \code{binsCompleted}.
#'    Default is FALSE.
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
setGeneric("computePvalues_batch", function(object, quant, poolQuant, movAve, Global_lower, use_cpp = TRUE, ignore_sitesUnused = FALSE, na_impute = TRUE,
                                            bins, create_pMat = FALSE,...) {
    standardGeneric("computePvalues_batch")
})

##'  show \code{tanDb} object in interactive mode
##'
##' content for details
##' @title
##' @param object A \code{tanDb} object
##' @return output an overview of \code{tanDb} object
##' @export
##' @author Duy Nguyen
setGeneric("show", function(object) {
    standardGeneric("show")
})
