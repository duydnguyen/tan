context("calculateWithinSites")

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

    # test:
    load(file = paste(path, 'tanDb_n4.RData', sep =''))
    expect_equal(tanDb_test@wSites, tanDb@wSites)
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

    # test:
    load(file = paste(path, 'tanDb_n3.RData', sep =''))
    expect_equal(tanDb_test@wSites, tanDb@wSites)
})

test_that(" Simulation with n = 2", {
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
    tanDb_test <- calculateTotalCounts(tanDb_test, nSamples = 2, bNormWidth = bNormWidth, bSampleMean = bSampleMean)

    tanDb_test <- calculateWithinSites(tanDb_test, quantprobs = quantprobs)

    # test:
    load(file = paste(path, 'tanDb_n2.RData', sep =''))
    expect_equal(tanDb_test@wSites, tanDb@wSites)
})
