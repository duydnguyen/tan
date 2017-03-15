#' Perform Latin Hypercube sampling or uniform sampling on a genomic interval.
#'
#' @param s.size size of sampling points.
#' @param start start of the genomic interval.
#' @param end end of the genomic interval.
#' @param LHD True if performing LH sampling; False if performing uniform sampling.
#' @param Uniform for case LHD = TRUE, if Uniform = TRUE, we
#'   sample points uniformly for each bin; FALSE if sampling points are end points.
#'
#' @return An one dimenstional array storing the sampling genomic positions
#' @export
#' @author Duy Nguyen on April 6, 2015; modified on February 24, 2016
#' @examples
Sampling <- function(s.size, start = 1, end =1000, LHD = TRUE, Uniform = TRUE) {
    design <- matrix(NA, nrow = 1, ncol = s.size)
    if (LHD) {
        par <- round(seq(start, end, length.out=s.size+1))

        if (Uniform == TRUE) {
            for (i in 1:( length(par) - 1) ) {
                design[i] <- round(runif(n = 1, min = par[i], max = par[i+1]), digits = 0)
            }
        } else {
            design <- par[1:s.size]
        }

    }
    else design <- sort(sample.int(n = end, size = s.size, replace = FALSE))
    return(design)
}
