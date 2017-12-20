rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))

# devtools::install_github("duydnguyen/tan-coverage", auth_token = "a37d127a798d298abb2b3a1ea306e3ab08babc0b")

library(Segvis)
library(GenomicRanges)
library(ggplot2)
library(data.table)
library(tan)

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

block_list <- Segvis_block_list(rep1.1_minus_ChIP_block)
rstart <- start(gr)[1]
rend <- end(gr)[1]
chr <- as.character(seqnames(gr)[1])
iden <- function(x) x

### Visualize
p1 <- plot_profiles(block_list,
                   condition = seqnames == chr & start == rstart,
                   coord = rstart:rend,FUN = iden,mc=mc_cores)
p1

### Extract read coverage
source("utility_functions.R")

p <- create_profile(block_list, condition = seqnames == chr & start == rstart,
                         coord = rstart:rend, FUN = iden, mc = mc_cores)

coverage <- get_coverage(p)

coverage2 <- get_coverage2(block_list, condition = seqnames == chr & start == rstart,
                         coord = rstart:rend, FUN = iden, mc = mc_cores)

coverage3 <- tan::generate_coverage(block_list, condition = seqnames == chr & start == rstart,
                         coord = rstart:rend, FUN = iden, mc = mc_cores)

identical(coverage, coverage2)
