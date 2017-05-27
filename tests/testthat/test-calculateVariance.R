context("calculate Variance")

# test_that(" down-sampling data with n = 4", {
#     ### parameters ###
#     s.size <- 500; LHD <- TRUE; Uniform <- FALSE
#     bNormWidth <- FALSE; bSampleMean <- FALSE
#     # quantile vector for binning
#     quantprobs <- seq(0, 1, 0.05) # this is the default binning for down sampling
#     ## set lower bound for minGlobal (length of pooled var vector for each bins)
#     Global_lower <- floor(s.size/2) # defaul value: Global_lower <- 0
#     ## pooled quantile at each genomic position: (for calculateVariance())
#     poolQuant <- 0.5
#     # number of points for moving average
#     movAve <- 20
#     ## pooled quantile at each genomic position:
#     poolQuant <- 0.5
#     ### coverage_sample
#     load("/p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/coverage_downsample.RData")
#     fooDb <- new("tanDb", coverage = coverage)
#     fooDb <- createDesigns(fooDb, s.size = s.size, LHD = LHD, Uniform = Uniform )
#     fooDb <- calculateTotalCounts(fooDb, nSamples = 4, bNormWidth = bNormWidth, bSampleMean = bSampleMean)
#     fooDb <- calculateWithinSites(fooDb, quantprobs = quantprobs)
#     fooDb <- calculateVariance(fooDb, minus_condition = TRUE, Global_lower = Global_lower, poolQuant = poolQuant, movAve = movAve )
#     fooDb <- calculateVariance(fooDb, minus_condition = FALSE, Global_lower = Global_lower, poolQuant = poolQuant, movAve = movAve )
#     # test:
#     load("/p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/DiffAnalysis_n4/minusVar.RData")
#     load("/p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/DiffAnalysis_n4/plusVar.RData")
#     expect_equal(fooDb@minusVar, minusVar)
#     expect_equal(fooDb@plusVar, plusVar)
# })
#
# test_that(" down-sampling data with n = 3", {
#     ### parameters ###
#     s.size <- 500; LHD <- TRUE; Uniform <- FALSE
#     bNormWidth <- FALSE; bSampleMean <- FALSE
#     # quantile vector for binning
#     quantprobs <- seq(0, 1, 0.05) # this is the default binning for down sampling
#     ## set lower bound for minGlobal (length of pooled var vector for each bins)
#     Global_lower <- floor(s.size/2) # defaul value: Global_lower <- 0
#     ## pooled quantile at each genomic position: (for calculateVariance())
#     poolQuant <- 0.5
#     # number of points for moving average
#     movAve <- 20
#     ## pooled quantile at each genomic position:
#     poolQuant <- 0.5
#     ### coverage_sample
#     load("/p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/coverage_downsample.RData")
#     # get coverage for samples 1,2,3
#     coverage123 <- list()
#     for (site in 1: length(coverage)) {
#         coverage123[[site]] <- coverage[[site]][c(1:3, 5:7),]
#     }
#     fooDb <- new("tanDb", coverage = coverage123)
#     fooDb <- createDesigns(fooDb, s.size = s.size, LHD = LHD, Uniform = Uniform )
#     fooDb <- calculateTotalCounts(fooDb, nSamples = 3, bNormWidth = bNormWidth, bSampleMean = bSampleMean)
#     fooDb <- calculateWithinSites(fooDb, quantprobs = quantprobs)
#     fooDb <- calculateVariance(fooDb, minus_condition = TRUE, Global_lower = Global_lower, poolQuant = poolQuant, movAve = movAve )
#     fooDb <- calculateVariance(fooDb, minus_condition = FALSE, Global_lower = Global_lower, poolQuant = poolQuant, movAve = movAve )
#     # test:
#     load("/p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/DiffAnalysis_n3/Samples_123/minusVar.RData")
#     load("/p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/DiffAnalysis_n3/Samples_123/plusVar.RData")
#     expect_equal(fooDb@minusVar, minusVar)
#     expect_equal(fooDb@plusVar, plusVar)
# })
#
# test_that(" down-sampling data with n = 2", {
#     ### parameters ###
#     s.size <- 500; LHD <- TRUE; Uniform <- FALSE
#     bNormWidth <- FALSE; bSampleMean <- FALSE
#     # quantile vector for binning
#     quantprobs <- seq(0, 1, 0.05) # this is the default binning for down sampling
#     ## set lower bound for minGlobal (length of pooled var vector for each bins)
#     Global_lower <- floor(s.size/2) # defaul value: Global_lower <- 0
#     ## pooled quantile at each genomic position: (for calculateVariance())
#     poolQuant <- 0.5
#     # number of points for moving average
#     movAve <- 20
#     ## pooled quantile at each genomic position:
#     poolQuant <- 0.5
#     ### coverage_sample
#     load("/p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/coverage_downsample.RData")
#     # get coverage for samples 1,2
#     coverage12 <- list()
#     for (site in 1: length(coverage)) {
#         coverage12[[site]] <- coverage[[site]][c(1:2,5:6),]
#     }
#     fooDb <- new("tanDb", coverage = coverage12)
#     print("create Designs")
#     fooDb <- createDesigns(fooDb, s.size = s.size, LHD = LHD, Uniform = Uniform )
#     fooDb <- calculateTotalCounts(fooDb, nSamples = 2, bNormWidth = bNormWidth, bSampleMean = bSampleMean)
#     fooDb <- calculateWithinSites(fooDb, quantprobs = quantprobs)
#     fooDb <- calculateVariance(fooDb, minus_condition = NA, Global_lower = Global_lower, poolQuant = poolQuant, movAve = movAve )
#     # test:
#     load("/p/keles/DBChIP/volume4/Summaries/2016_06_18/Generated/Simulation1/DiffAnalysis_n2/Samples_12/poolVar.RData")
#     expect_equal(fooDb@poolVar, poolVar)
# })

# test_that(" Simulation with n = 4", {
#     ### parameters ###
#     s.size <- 500; LHD <- TRUE; Uniform <- FALSE
#     bNormWidth <- FALSE; bSampleMean <- FALSE
#     # quantile vector for binning
#     quantprobs <- seq(0, 1, 0.05) # this is the default binning for down sampling
#     ## set lower bound for minGlobal (length of pooled var vector for each bins)
#     Global_lower <- 0 # defaul value: Global_lower <- 0
#     ## pooled quantile at each genomic position: (for calculateVariance())
#     poolQuant <- 0.5
#     # number of points for moving average
#     movAve <- 20
#     ## pooled quantile at each genomic position:
#     poolQuant <- 0.5
#     ### coverage from Simulation
#     path <- "/p/keles/DBChIP/volume2/Summaries/2016_01_23/Generated/EmpPvalPool/"
#     load(file = "/u/d/n/dnguyen/Projects/Histone_Modifications/2015_10_27/Generated/CoverageList.RData")
#     sim <- 1
#     coverage <- coverageList[[sim]]
#     coverage[["labels"]] <- coverage[["labels_new"]] <- NULL
#     fooDb <- new("tanDb", coverage = coverage)
#     print("create Designs")
#     fooDb <- createDesigns(fooDb, s.size = s.size, LHD = LHD, Uniform = Uniform )
#     fooDb <- calculateTotalCounts(fooDb, nSamples = 4, bNormWidth = bNormWidth, bSampleMean = bSampleMean)
#     fooDb <- calculateWithinSites(fooDb, quantprobs = quantprobs)
#     fooDb <- calculateVariance(fooDb, minus_condition = TRUE, Global_lower = Global_lower, poolQuant = poolQuant, movAve = movAve )
#     fooDb <- calculateVariance(fooDb, minus_condition = FALSE, Global_lower = Global_lower, poolQuant = poolQuant, movAve = movAve )
#     # test:
#     load(file = paste(path, 'minusVar.RData', sep =''))
#     load(file = paste(path, 'plusVar.RData', sep =''))
#     expect_equal(fooDb@minusVar, minusVar)
#     expect_equal(fooDb@plusVar, plusVar)
# })

test_that(" Simulation with n = 4", {
    source("../../../../Simulation/R/params_DA.R")
    ### coverage from Simulation
    path <- "../../../../Simulation/DeriveData/"
    load(file = paste(path, "coverageList.RData", sep = ""))
    sim <- 1
    coverage <- coverageList[[sim]]
    coverage[["labels"]] <- coverage[["labels_new"]] <- NULL
    tanDb_test <- new("tanDb", coverage = coverage)
    print("create Designs")
    tanDb_test <- createDesigns(tanDb_test, s.size = s.size, LHD = LHD, Uniform = Uniform )
    tanDb_test <- calculateTotalCounts(tanDb_test, nSamples = 4, bNormWidth = bNormWidth, bSampleMean = bSampleMean)

    tanDb_test <- calculateWithinSites(tanDb_test, quantprobs = quantprobs)

    tanDb_test <- calculateVariance(tanDb_test, minus_condition = TRUE,
                                   Global_lower = Global_lower, poolQuant = poolQuant,
                                   movAve = movAve )

    tanDb_test <- calculateVariance(tanDb_test, minus_condition = FALSE,
                                   Global_lower = Global_lower, poolQuant = poolQuant,
                                   movAve = movAve )

    # test:
    load(file = paste(path, 'tanDb_n4.RData', sep =''))
    expect_equal(tanDb_test@minusVar, tanDb@minusVar)
    expect_equal(tanDb_test@plusVar, tanDb@plusVar)
})

test_that(" Simulation with n = 3", {
    source("../../../../Simulation/R/params_DA.R")
    ### coverage from Simulation
    path <- "../../../../Simulation/DeriveData/"
    load(file = paste(path, "coverageList.RData", sep = ""))
    sim <- 1
    coverage <- coverageList[[sim]]
    coverage[["labels"]] <- coverage[["labels_new"]] <- NULL
    tanDb_test <- new("tanDb", coverage = coverage)
    print("create Designs")
    tanDb_test <- createDesigns(tanDb_test, s.size = s.size, LHD = LHD, Uniform = Uniform )
    tanDb_test <- calculateTotalCounts(tanDb_test, nSamples = 3, bNormWidth = bNormWidth, bSampleMean = bSampleMean)

    tanDb_test <- calculateWithinSites(tanDb_test, quantprobs = quantprobs)

    tanDb_test <- calculateVariance(tanDb_test, minus_condition = TRUE,
                                   Global_lower = Global_lower, poolQuant = poolQuant,
                                   movAve = movAve )

    ## tanDb_test <- calculateVariance(tanDb_test, minus_condition = FALSE,
    ##                                Global_lower = Global_lower, poolQuant = poolQuant,
    ##                                movAve = movAve )

    # test:
    load(file = paste(path, 'tanDb_n3.RData', sep =''))
    # expect_equal(tanDb_test@minusVar, tanDb@minusVar)
    expect_equal(tanDb_test@minusVar[[2]]$ac, tanDb@minusVar[[2]]$ac)

    ## expect_equal(tanDb_test@plusVar, tanDb@plusVar)

})

## test_that(" Simulation with n = 2", {
##     source("../../../../Simulation/R/params_DA.R")
##     ### coverage from Simulation
##     path <- "../../../../Simulation/DeriveData/"
##     load(file = paste(path, "coverageList.RData", sep = ""))
##     sim <- 1
##     coverage <- coverageList[[sim]]
##     coverage[["labels"]] <- coverage[["labels_new"]] <- NULL
##     tanDb_test <- new("tanDb", coverage = coverage)
##     print("create Designs")
##     tanDb_test <- createDesigns(tanDb_test, s.size = s.size, LHD = LHD, Uniform = Uniform )
##     tanDb_test <- calculateTotalCounts(tanDb_test, nSamples = 2, bNormWidth = bNormWidth, bSampleMean = bSampleMean)

##     tanDb_test <- calculateWithinSites(tanDb_test, quantprobs = quantprobs)

##     tanDb_test <- calculateVariance(tanDb_test, minus_condition = TRUE,
##                                    Global_lower = Global_lower, poolQuant = poolQuant,
##                                    movAve = movAve )

##     tanDb_test <- calculateVariance(tanDb_test, minus_condition = FALSE,
##                                    Global_lower = Global_lower, poolQuant = poolQuant,
##                                    movAve = movAve )

##     # test:
##     load(file = paste(path, 'tanDb_n2.RData', sep =''))
##     expect_equal(tanDb_test@minusVar, tanDb@minusVar)
##     expect_equal(tanDb_test@plusVar, tanDb@plusVar)

## })
