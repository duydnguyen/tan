context("create Designs")

#createDesigns
test_that(" down-sampling data with n = 4", {
    s.size <- 500; LHD <- TRUE; Uniform <- FALSE
    ### coverage_sample
    load("~/Google Drive/DiffAna_sample/Generated/coverage_downsample.RData")
    load("~/Google Drive/DiffAna_sample/Generated/Designs.RData")
    fooDb <- new("tanDb", coverage = coverage)
    fooDb <- createDesigns(fooDb, s.size = s.size, LHD = LHD, Uniform = Uniform )
    print("Testing createDesigns(): ")
    expect_equal(fooDb@Designs, Designs)
})
