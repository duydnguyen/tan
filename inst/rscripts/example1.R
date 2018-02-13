rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))

# devtools::install_github("duydnguyen/tan-coverage", auth_token = "a37d127a798d298abb2b3a1ea306e3ab08babc0b")

library(Segvis)
library(GenomicRanges)
library(ggplot2)
library(data.table)

# parameters
filedir <- "../extdata/H3K27me3"
gnames <- "BZW2"
mc_cores <- 4
# prefix <- paste('../DeriveData/', gnames, sep ='')

# Build a GRanges from bed file
bed_content <- read.table(file = file.path(filedir, "BZW2.bed"), stringsAsFactors = FALSE)
gr <- GRanges(seqnames = bed_content[, 1],
              ranges = IRanges(bed_content[,2], bed_content[,3]), strand = "*")

chromosomes <- c("chr1","chr2", "chr3","chr4","chr5","chr6","chr7","chr8",
                 "chr9","chr10","chr11","chr12","chr13","chr14", "chr15","chr16",
                 "chr17","chr18","chr19","chr20","chr21","chr22","chrX")

### REP1 ###
### Minus, ChIP
# Create a segvis object
rep1.1_minus_ChIP <- buildSegvis(name = "rep1.1_minus_ChIP",
                       file = file.path(filedir, "BZW2_rep1_minus_sorted.bam"),
                       maxBandwidth = 101, fragLen = 300, isPET = FALSE,
                       chr = chromosomes)
regions(rep1.1_minus_ChIP) <- gr
rep1.1_minus_ChIP
# Create segvis_block object
rep1.1_minus_ChIP <- loadReads(rep1.1_minus_ChIP, mc = mc_cores)
rep1.1_minus_ChIP <- matchReads(rep1.1_minus_ChIP, mc = mc_cores)
rep1.1_minus_ChIP <- getCoverage(rep1.1_minus_ChIP, mc = mc_cores)
rep1.1_minus_ChIP_block <- Segvis_block(rep1.1_minus_ChIP, bw = 1, mc = mc_cores)
normConst(rep1.1_minus_ChIP_block) <- 87626218
rep1.1_minus_ChIP_block <- normalize(rep1.1_minus_ChIP_block, base = 10^8)
### Plus, ChIP
# Create a segvis object
rep1.1_plus_ChIP <- buildSegvis(name = "rep1.1_plus_ChIP",
                     file = file.path(filedir, "BZW2_rep1_plus_sorted.bam"),
                           maxBandwidth = 101, fragLen = 300, isPET = FALSE,
                           chr = chromosomes)
regions(rep1.1_plus_ChIP) <- gr
rep1.1_plus_ChIP
# Create segvis_block object
rep1.1_plus_ChIP <- loadReads(rep1.1_plus_ChIP, mc = mc_cores)
rep1.1_plus_ChIP <- matchReads(rep1.1_plus_ChIP, mc = mc_cores)
rep1.1_plus_ChIP <- getCoverage(rep1.1_plus_ChIP, mc = mc_cores)
rep1.1_plus_ChIP_block <- Segvis_block(rep1.1_plus_ChIP, bw = 1, mc = mc_cores)
normConst(rep1.1_plus_ChIP_block) <- 90287636
rep1.1_plus_ChIP_block <- normalize(rep1.1_plus_ChIP_block, base = 10^8)

### REP2 ###
### Minus, ChIP
# Create a segvis object
rep2_minus_ChIP <- buildSegvis(name = "rep2_minus_ChIP",
                          file = file.path(filedir, "BZW2_rep2_minus_sorted.bam"),
                          maxBandwidth = 101, fragLen = 200, isPET = FALSE,
                          chr = chromosomes)
regions(rep2_minus_ChIP) <- gr
rep2_minus_ChIP
# Create segvis_block object
rep2_minus_ChIP <- loadReads(rep2_minus_ChIP, mc = mc_cores)
rep2_minus_ChIP <- matchReads(rep2_minus_ChIP, mc = mc_cores)
rep2_minus_ChIP <- getCoverage(rep2_minus_ChIP, mc = mc_cores)
rep2_minus_ChIP_block <- Segvis_block(rep2_minus_ChIP, bw = 1, mc = mc_cores)
normConst(rep2_minus_ChIP_block) <- 63142317
rep2_minus_ChIP_block <- normalize(rep2_minus_ChIP_block, base = 10^8)
### Plus, ChIP
# Create a segvis object
rep2_plus_ChIP <- buildSegvis(name = "rep2_plus_ChIP",
                         file = file.path(filedir, "BZW2_rep2_plus_sorted.bam"), maxBandwidth = 101, fragLen = 200, isPET = FALSE, chr = chromosomes)
regions(rep2_plus_ChIP) <- gr
rep2_plus_ChIP
# Create segvis_block object
rep2_plus_ChIP <- loadReads(rep2_plus_ChIP, mc = mc_cores)
rep2_plus_ChIP <- matchReads(rep2_plus_ChIP, mc = mc_cores)
rep2_plus_ChIP <- getCoverage(rep2_plus_ChIP, mc = mc_cores)
rep2_plus_ChIP_block <- Segvis_block(rep2_plus_ChIP, bw = 1, mc = mc_cores)
normConst(rep2_plus_ChIP_block) <- 83789159
rep2_plus_ChIP_block <- normalize(rep2_plus_ChIP_block, base = 10^8)


### REP3 ###
### Minus, ChIP
# Create a segvis object
rep3_minus_ChIP <- buildSegvis(name = "rep3_minus_ChIP",
                          file = file.path(filedir, "BZW2_rep3_minus_sorted.bam"),
                          maxBandwidth = 101, fragLen = 200, isPET = FALSE,
                          chr = chromosomes)
regions(rep3_minus_ChIP) <- gr
rep3_minus_ChIP
# Create segvis_block object
rep3_minus_ChIP <- loadReads(rep3_minus_ChIP, mc = mc_cores)
rep3_minus_ChIP <- matchReads(rep3_minus_ChIP, mc = mc_cores)
rep3_minus_ChIP <- getCoverage(rep3_minus_ChIP, mc = mc_cores)
rep3_minus_ChIP_block <- Segvis_block(rep3_minus_ChIP, bw = 1, mc = mc_cores)
normConst(rep3_minus_ChIP_block) <- 58779937
rep3_minus_ChIP_block <- normalize(rep3_minus_ChIP_block, base = 10^8)
### Plus, ChIP
# Create a segvis object
rep3_plus_ChIP <- buildSegvis(name = "rep3_plus_ChIP",
                         file = file.path(filedir, "BZW2_rep3_plus_sorted.bam"),
                         maxBandwidth = 101, fragLen = 200, isPET = FALSE,
                         chr = chromosomes)
regions(rep3_plus_ChIP) <- gr
rep3_plus_ChIP
# Create segvis_block object
rep3_plus_ChIP <- loadReads(rep3_plus_ChIP, mc = mc_cores)
rep3_plus_ChIP <- matchReads(rep3_plus_ChIP, mc = mc_cores)
rep3_plus_ChIP <- getCoverage(rep3_plus_ChIP, mc = mc_cores)
rep3_plus_ChIP_block <- Segvis_block(rep3_plus_ChIP, bw = 1, mc = mc_cores)
normConst(rep3_plus_ChIP_block) <- 56388356
rep3_plus_ChIP_block <- normalize(rep3_plus_ChIP_block, base = 10^8)

##########################
#### Extract Coverage ####
##########################
#### minus Block: ChIP, Input
block_minus_list <- Segvis_block_list(rep1.1_minus_ChIP_block, rep2_minus_ChIP_block, rep3_minus_ChIP_block)
#### plus Block: ChIP, Input
block_plus_list <- Segvis_block_list(rep1.1_plus_ChIP_block,rep2_plus_ChIP_block, rep3_plus_ChIP_block)

source("utility_functions.R")

## coverage <- list()
## for (site in 1:length(gr)) {
##     print(paste("site = ", site))
##     rstart <- start(gr)[site]
##     rend <- end(gr)[site]
##     chr <- as.character(seqnames(gr)[site])
##     iden <- function(x)x
##     p_minus <- plot_profiles2(block_minus_list, condition = seqnames == chr & start == rstart,
##                              coord = rstart:rend, FUN = iden, mc = mc_cores)
##     p_plus <- plot_profiles2(block_plus_list, condition = seqnames == chr & start == rstart,
##                             coord = rstart:rend, FUN = iden, mc= mc_cores)
##     df <- toCoverage(p_minus, p_plus) # fix toCoverage for general sample sizes
##     coverage[[site]] <- t(df)
##     ## df <- rbind(get_coverage(p_minus), get_coverage(p_plus))
##     ## coverage[[site]] <- df
## }
## old <- coverage

coverage <- list()
for (site in 1:length(gr)) {
    print(paste("site = ", site))
    rstart <- start(gr)[site]
    rend <- end(gr)[site]
    chr <- as.character(seqnames(gr)[site])
    iden <- function(x)x
    p_minus <- create_profile(block_minus_list, condition = seqnames == chr & start == rstart,
                             coord = rstart:rend, FUN = iden, mc = mc_cores)
    p_plus <- create_profile(block_plus_list, condition = seqnames == chr & start == rstart,
                            coord = rstart:rend, FUN = iden, mc= mc_cores)
    df <- rbind(get_coverage(p_minus), get_coverage(p_plus))
    coverage[[site]] <- df
}



## # plot
sampleSize <- 100
sizePlot <- 0.5

coverage <- coverage[[1]]

design <- seq(1, dim(coverage)[2], length = sampleSize)

p <- plotCoverage(coverage[, design], geneNames = gnames, size = sizePlot, showLegend = TRUE)
p
