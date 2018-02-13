#' A class to store ChIP data and all related objects for differential analysis.
#'
#' @slot coverage A list stores all coverage data. The length of this list
#'   is equal to the number of pre-defined genomic regions to be tested.
#' @slot Designs A matrix stores the sampling locations for genomic regions.
#' @slot Ns A matrix stores the total counts for each genomic regions and
#'   across library
#' @slot s.size The number of sampling points allowed in genomic regions.
#' @slot wSites A list of within-site indices for each bin based
#'   on total counts \code{Ns} and \code{quantprobs}.
#' @slot quantprobs  A quantile vector for binning.
#' @slot nSamples The sample size for each conditions (n1=n2=n)
#' @slot dN a vector of count ranges for each bins based on \code{quantprobs}.
#' @slot minusVar A list of pooled variances for the first condition. Pooling method for
#'    variance is disscussed in our paper.
#' @slot plusVar A list of pooled variances for the second condition. Pooling method for
#'    variance is disscussed in our paper.
#' @slot poolVar A list of pooled variances for case sanple size \code{n = 2}.
#' @slot W1 A matrix stores within Adaptive Neyman tests for first condition.
#' @slot W2 A matrix stores within Adaptive Neyman tests for second condition.
#' @slot W A matrix stores Adaptive Neyman tests under the null hypothesis. This is only for case sample size \code{n=2}.
#' @slot PvalList A list of results (\code{('pval','FDR')})for P-values.
#' @slot pMat A matrix with number of rows equal to number of peaks \code{x}. \code{pMat} stores all the p-values
#'    from running computePvalues in batch mode (useful when testing very large number of peaks).
#' @slot p.list A list to store results in data frame format when running batch mode.
#' @slot binsCompleted A vector stores all bins that are succesfullly run in batch mode.
#' @slot sitesUnused A vector of indices for unused sites due to bad quality peaks during testing within and between TAN.
#'
#' @return A prediction of differential regions
#' @export
#'
#' @examples
setClass("tanDb",
         representation = representation(coverage = "list", Designs = "matrix", Ns = "matrix",
                                         s.size = "numeric", wSites = "list", quantprobs = "numeric",
                                         dN = "numeric", nSamples = "numeric",
                                         minusVar = "list", plusVar = "list", poolVar = "list",
                                         W1 = "matrix", W2 = "matrix", W = "matrix",
                                         PvalList = "list", pMat = "matrix", p.list = "list", binsCompleted = "numeric", sitesUnused = "numeric"),
         prototype = prototype(coverage = list(), Designs = matrix(), pMat = matrix())
         )

