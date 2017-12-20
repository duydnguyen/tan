## ----style, echo = FALSE, results = 'asis'---------------------------------

    library(BiocStyle)
    markdown(css.files = c('custom.css'))


## ----extraload, include = FALSE, echo = FALSE, eval = TRUE-----------------
    library(tan)
    library(knitr)
    library(ggplot2)
    opts_chunk$set(fig.align = "center")


## ----load,include=TRUE,echo=TRUE,eval=FALSE--------------------------------
#      library(tan)
#      library(tanExample)
#  

## ----loadSegvis_show, include = TRUE, echo = TRUE, eval = FALSE------------
#      library(GenomicRanges)
#      library(GenomicAlignments)
#      library(data.table)
#  

## ----loadSegvis_noshow, include = FALSE, echo = FALSE, eval = TRUE---------
    library(GenomicRanges)
    library(GenomicAlignments)
    library(data.table)
    

## ----extractBam, include=TRUE, echo=TRUE, eval=TRUE------------------------
files = list.files(system.file("extdata/bzw2",
        package = "tanExample"),full.names = TRUE)
basename(files[c(1,3,6,9,12)])
bam_files <- files[c(3,6, 9, 12)]
bed_files <- files[1]
# checking if there is an index 
file.exists(gsub(".bam$", ".bam.bai", files[c(1,3,6,9,12)]))

## ----params, include=TRUE, echo=TRUE, eval=TRUE----------------------------
binsize <- 150
sm <- 1 # smooth parameter
mc_cores <- 3
bed_content <- read.table(file = files[1], stringsAsFactors = FALSE)
gr <- GRanges(seqnames = bed_content[, 1],
              ranges = IRanges(bed_content[,2], bed_content[,3]), 
              strand = "*")
chromosomes <- c("chr7")
gr

## ----bamCoverage, include=TRUE, echo=TRUE, eval=TRUE-----------------------
coverage <- tan::bamCoverage(bam_files = bam_files, bed_files = bed_files,
                         mc_cores = mc_cores, sm = sm, binsize = binsize)

class(coverage)
length(coverage)

## ----bamCoverage_plot, include=TRUE, echo=TRUE, eval=TRUE------------------

tan::plotCoverage(coverage[[1]], title = toString(gr[1]))

## ----GR_info, include = TRUE, echo = TRUE, eval = TRUE---------------------
files = list.files(system.file("extdata",
        package = "tanExample"),full.names = TRUE)
basename(files[2:3])
load(files[3])
gr_sitesSelect

## ----tanDb_construct, include = TRUE, echo = TRUE, eval = FALSE------------
#  load(files[2])
#  tanDb <- new("tanDb", coverage = coverage)

## ----design, include = TRUE, echo = TRUE, eval = FALSE---------------------
#  tanDb <- createDesigns(tanDb, s.size = 500, LHD = TRUE, Uniform = FALSE )

## ----totalCounts, include = FALSE, echo = FALSE, eval = TRUE---------------
load(files[4])

## ----totalCounts2, include = TRUE, echo = TRUE, eval = FALSE---------------
#  tanDb <- calculateTotalCounts(tanDb, nSamples = 3, bNormWidth = FALSE,
#                                bSampleMean = FALSE)
#  head(tanDb@Ns)

## ----totalCounts3, include = TRUE, echo = FALSE, eval = TRUE---------------
print(head(tanDb@Ns))

## ----calculateWithinSites, include = TRUE, echo = TRUE, eval = FALSE-------
#  # quantile vector for binning
#  quantprobs <- seq(0, 1, 0.05)
#  tanDb <- calculateWithinSites(tanDb, quantprobs = quantprobs)

## ----calculateVariance, include = TRUE, echo = TRUE, eval = FALSE----------
#  Global_lower <- 100
#  ## pooled quantile at each genomic position:
#  poolQuant <- 0.5
#  # number of points for moving average
#  movAve <- 20
#  tanDb <- calculateVariance(tanDb, minus_condition = TRUE,
#                             Global_lower = Global_lower,
#                             poolQuant = poolQuant, movAve = movAve )
#  tanDb <- calculateVariance(tanDb, minus_condition = FALSE,
#                             Global_lower = Global_lower,
#                             poolQuant = poolQuant, movAve = movAve )

## ----generateWithinTan, include = TRUE, echo = TRUE, eval = FALSE----------
#  tanDb <- generateWithinTan(tanDb, minus_condition = TRUE)
#  tanDb <- generateWithinTan(tanDb, minus_condition = FALSE)

## ----computePvalues, include = TRUE, echo = TRUE, eval = FALSE-------------
#  tanDb <- computePvalues(tanDb, quant = 1, poolQuant = poolQuant,
#                          movAve = movAve, Global_lower = Global_lower,
#                          ignore_sitesUnused = TRUE,
#                          na_impute = FALSE)

