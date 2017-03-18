context("calculateTotalCounts")

test_that(" down-sampling data with n = 4", {
    ### parameters ###
    s.size <- 500; LHD <- TRUE; Uniform <- FALSE
    bNormWidth <- FALSE; bSampleMean <- FALSE
    # quantile vector for binning
    quantprobs <- seq(0, 1, 0.05) # this is the default binning for down sampling
    ## set lower bound for minGlobal (length of pooled var vector for each bins)
    Global_lower <- floor(s.size/2) # defaul value: Global_lower <- 0
    ## pooled quantile at each genomic position: (for calculateVariance())
    poolQuant <- 0.5
    # number of points for moving average
    movAve <- 20
    ## pooled quantile at each genomic position:
    poolQuant <- 0.5
    ### coverage_sample
    load("/p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/coverage_downsample.RData")
    fooDb <- new("tanDb", coverage = coverage)
    fooDb <- createDesigns(fooDb, s.size = s.size, LHD = LHD, Uniform = Uniform )
    fooDb <- calculateTotalCounts(fooDb, nSamples = 4, bNormWidth = bNormWidth, bSampleMean = bSampleMean)
    # test:
    load("p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/DiffAnalysis_n4/Ns.RData")
    expect_equal(fooDb@Ns, Ns)
})