library(testthat)
library(CytoPipeline)

data(OMIP021Samples)
runDumpDistributions <- FALSE # activate for full mounty testing (but slow)
test_check("CytoPipeline")
