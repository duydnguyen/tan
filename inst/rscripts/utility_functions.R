### the following functions used to create coverage from Segvis blocks
## .filter_sb <- function (object, cond)
## {
##     if (any(is.na(cond))) {
##         warning("There are regions impossible to evaluate")
##         cond[is.na(cond)] = FALSE
##     }
##     V1 <- NULL
##     coverage <- copy(cover_table(object))
##     lengths <- coverage[, length(coord), by = list(chr, match)][,
##                                                                 (V1)]
##     extended_cond <- unlist(mapply(function(x, l) rep(x, l),
##                                    cond, lengths, SIMPLIFY = FALSE))
##     out_regions <- regions(object)[cond]
##     rm(cond)
##     coverage[, `:=`(cond, extended_cond)]
##     coverage <- coverage[cond == TRUE]
##     coverage[, `:=`(cond, NULL)]
##     out <- new("segvis_block", name = name(object), regions = out_regions,
##                bandwidth = bandwidth(object), normConst = normConst(object),
##                cover_table = coverage, .isScaled = object@.isScaled)
##     return(out)
## }

## create_plot_data <- function (counts, name, coord)
## {
##     dt <- data.table(x = coord, y = counts, condition = name)
##     return(dt)
## }

## .subset_logical <- function (object, condition_call)
## {
##     cond <- as.logical(eval(condition_call, as(regions(object),
##                                                "data.frame"), parent.frame()))
##     return(cond)
## }

create_profile <- function (object, condition, coord = NULL, mc, FUN, ...)
{
    .filter_sb <- function (object, cond)
    {
        if (any(is.na(cond))) {
            warning("There are regions impossible to evaluate")
            cond[is.na(cond)] = FALSE
        }
        V1 <- NULL
        coverage <- copy(cover_table(object))
        lengths <- coverage[, length(coord), by = list(chr, match)][,
            (V1)]
        extended_cond <- unlist(mapply(function(x, l) rep(x, l),
                                      cond, lengths, SIMPLIFY = FALSE))
        out_regions <- regions(object)[cond]
        rm(cond)
        coverage[, `:=`(cond, extended_cond)]
        coverage <- coverage[cond == TRUE]
        coverage[, `:=`(cond, NULL)]
        out <- new("segvis_block", name = name(object), regions = out_regions,
                  bandwidth = bandwidth(object), normConst = normConst(object),
                  cover_table = coverage, .isScaled = object@.isScaled)
        return(out)
    }
    create_plot_data <- function (counts, name, coord)
    {
        dt <- data.table(x = coord, y = counts, condition = name)
        return(dt)
    }
    .subset_logical <- function (object, condition_call)
    {
        cond <- as.logical(eval(condition_call, as(regions(object),
                                                  "data.frame"), parent.frame()))
        return(cond)
    }
    ### MAIN ###
    if (is.null(names(object))) {
        nms <- 1:length(object)
    }
    else {
        nms <- names(object)
    }
    if (!missing(condition)) {
        conds <- mclapply(object, .subset_logical, substitute(condition),
                          mc.cores = mc)
    }
    else {
        nregions <- lapply(object, function(x) length(regions(x)))
        conds <- lapply(nregions, function(x) rep(TRUE, x))
    }
    subsets <- mcmapply(.filter_sb, object, conds, SIMPLIFY = FALSE,
                        mc.cores = mc)
    widths <- mclapply(subsets, function(x) width(regions(x)),
                       mc.cores = mc)
    if (length(unique(unlist(widths))) > 1) {
        stop("The supplied regions doesn't have the same width")
    }
    else plot_width <- unique(unlist(widths))
    if (is.null(coord)) {
        coord <- 1:plot_width
    }
    profiles <- mclapply(subsets, summarize, FUN, ..., mc.cores = mc)
    plot_data <- mcmapply(create_plot_data, profiles, nms, MoreArgs = list(coord),
                          SIMPLIFY = FALSE, mc.silent = TRUE, mc.cores = mc, mc.preschedule = TRUE)
    # plot_data <- do.call(rbind, plot_data)
    return(plot_data)
}

######
## toCoverage <- function(p_minus, p_plus) {
##     df <- data.frame(matrix(ncol = 0, nrow = dim(p_minus[[1]])[1] ))
##     for (i in 1:3) {
##         df <- cbind(df, p_minus[[i]][, y])
##     }
##     for (i in 1:3) {
##         df <- cbind(df, p_plus[[i]][, y])
##     }
##     names(df) <- c('mRep1.1', 'mRep2' ,'mRep3',
##                    'pRep1.1', 'pRep2' ,'pRep3')
##     rownames(df) <- NULL
##     return(df)
## }


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title get read coverage
##' @param plot_profile : a list returned by running \code{plot_profiles2}
##' @return a data frame contains read coverage
##' @author Duy Nguyen
get_coverage <- function(plot_profile) {
    dfCoverage <- data.frame(matrix(ncol = 0, nrow = dim(plot_profile[[1]])[1] ))
    nSample <- length(plot_profile)
    for (i in 1:nSample) {
        dfCoverage <- cbind(dfCoverage, plot_profile[[i]][, y])
    }
    dfCoverage <- t(dfCoverage)
    rownames(dfCoverage) <- NULL
    return(dfCoverage)
}


#' Extract reads from a Segvis object
#'
#' @param object A Segvis_block_list object.
#' @param condition
#' @param coord a vector storing the genomic coordinates to be extracted.
#' @param mc number of cores for parallel pipeline.
#' @param FUN
#' @param ...
#'
#' @return a data frame contains read coverage.
#' @export
#'
#' @examples
get_coverage2 <- function (object, condition, coord = NULL, mc, FUN, ...)
{
    .filter_sb <- function (object, cond)
    {
        if (any(is.na(cond))) {
            warning("There are regions impossible to evaluate")
            cond[is.na(cond)] = FALSE
        }
        V1 <- NULL
        coverage <- copy(cover_table(object))
        lengths <- coverage[, length(coord), by = list(chr, match)][,
            (V1)]
        extended_cond <- unlist(mapply(function(x, l) rep(x, l),
                                      cond, lengths, SIMPLIFY = FALSE))
        out_regions <- regions(object)[cond]
        rm(cond)
        coverage[, `:=`(cond, extended_cond)]
        coverage <- coverage[cond == TRUE]
        coverage[, `:=`(cond, NULL)]
        out <- new("segvis_block", name = name(object), regions = out_regions,
                  bandwidth = bandwidth(object), normConst = normConst(object),
                  cover_table = coverage, .isScaled = object@.isScaled)
        return(out)
    }
    create_plot_data <- function (counts, name, coord)
    {
        dt <- data.table(x = coord, y = counts, condition = name)
        return(dt)
    }
    .subset_logical <- function (object, condition_call)
    {
        cond <- as.logical(eval(condition_call, as(regions(object),
                                                  "data.frame"), parent.frame()))
        return(cond)
    }
    ### MAIN ###
    if (is.null(names(object))) {
        nms <- 1:length(object)
    }
    else {
        nms <- names(object)
    }
    if (!missing(condition)) {
        conds <- mclapply(object, .subset_logical, substitute(condition),
                          mc.cores = mc)
    }
    else {
        nregions <- lapply(object, function(x) length(regions(x)))
        conds <- lapply(nregions, function(x) rep(TRUE, x))
    }
    subsets <- mcmapply(.filter_sb, object, conds, SIMPLIFY = FALSE,
                        mc.cores = mc)
    widths <- mclapply(subsets, function(x) width(regions(x)),
                       mc.cores = mc)
    if (length(unique(unlist(widths))) > 1) {
        stop("The supplied regions doesn't have the same width")
    }
    else plot_width <- unique(unlist(widths))
    if (is.null(coord)) {
        coord <- 1:plot_width
    }
    profiles <- mclapply(subsets, summarize, FUN, ..., mc.cores = mc)
    plot_data <- mcmapply(create_plot_data, profiles, nms, MoreArgs = list(coord),
                          SIMPLIFY = FALSE, mc.silent = TRUE, mc.cores = mc, mc.preschedule = TRUE)
    # plot_data <- do.call(rbind, plot_data)

    dfCoverage <- data.frame(matrix(ncol = 0, nrow = dim(plot_data[[1]])[1] ))
    nSample <- length(plot_data)
    for (i in 1:nSample) {
        dfCoverage <- cbind(dfCoverage, plot_data[[i]][, y])
    }
    dfCoverage <- t(dfCoverage)
    rownames(dfCoverage) <- NULL
    return(dfCoverage)
}







##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title coverage plot for GenePeak
##' @param coverage : a coverage resulting from Segvis block (see extractCoverage.R)
##' @param geneNames : name of (target) gene
##' @param size
##' @param showLegend : show legend or not
##' @return ggplot2 plot
##' @author Duy Nguyen
plotCoverage <- function(coverage, geneNames, size = 0.5, showLegend = TRUE) {
    geneChIP <- t(coverage)
    minus <- apply(geneChIP[, 1:3], 1, mean)
    plus <- apply(geneChIP[, 4:6], 1, mean)
    x_axis <- rep(seq(1,length(minus)), 2)
    samples <- c(rep('Condition 1', length(minus)), rep('Condition 2', length(minus)))
    df <- data.frame('RPM' = c(minus,plus), 'samples' = samples, 'x' = x_axis)
    #p <- qplot(x = x, y = RPM, color = samples, data = df, geom = "line", main = geneNames, xlab = 'Genomic Region (5 -> 3)') + theme(legend.position = "top")
    p <- ggplot(df, aes(x, RPM, color = samples)) + geom_line(size = size) +
        ylab("Mean normalized counts") + xlab("Genomic region (5 -> 3)") +
        labs(title= geneNames) + theme_bw() +
        scale_colour_manual(values = c('#619CFF','#CC6666'))
    if (showLegend) {
        return(p + theme(legend.position="top", legend.direction="horizontal", legend.title = element_blank(),legend.key.size = unit(0.65, "cm")))
        print(p)
    } else {
        return(p + theme(legend.position="none", legend.direction="horizontal",
                         legend.title = element_blank(),legend.key.size = unit(0.65, "cm"),
                         axis.text.x = element_blank(),
                         axis.ticks.x = element_blank()))
    }

}
