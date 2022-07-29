test_that("singletsGate works", {
  
  
  mySingletsGate <- singletsGate(OMIP021Samples[[1]])
  
  selectedSinglets <- flowCore::filter(OMIP021Samples[[1]], 
                                       mySingletsGate)
  ff_l <- OMIP021Samples[[1]][selectedSinglets@subSet, ]
  ff_l <- appendCellID(ff_l, which(selectedSinglets@subSet))
  
  linRange <- c(0,250000)
  p <- ggplotFilterEvents(
    ffPre = OMIP021Samples[[1]],
    ffPost = ff_l,
    seed = 1,
    xChannel = "FSC-A", xLinearRange = linRange,
    yChannel = "FSC-H", yLinearRange = linRange)
  
  vdiffr::expect_doppelganger("singletsGate default channels", fig = p)
  
  mySingletsGate <- singletsGate(OMIP021Samples[[1]], nmad = 3)
  
  selectedSinglets <- flowCore::filter(OMIP021Samples[[1]], 
                                       mySingletsGate)
  ff_l <- OMIP021Samples[[1]][selectedSinglets@subSet, ]
  ff_l <- appendCellID(ff_l, which(selectedSinglets@subSet))
  
  p <- ggplotFilterEvents(
    ffPre = OMIP021Samples[[1]],
    ffPost = ff_l,
    seed = 1,
    xChannel = "FSC-A", xLinearRange = linRange,
    yChannel = "FSC-H", yLinearRange = linRange)
  
  vdiffr::expect_doppelganger("singletsGate default channels with fixed nmad",
                              fig = p)
  
  mySingletsGate <- singletsGate(OMIP021Samples[[1]],
                                 channel1 = "SSC-A",
                                 channel2 = "SSC-H")
  
  selectedSinglets <- flowCore::filter(OMIP021Samples[[1]], 
                                       mySingletsGate)
  ff_l <- OMIP021Samples[[1]][selectedSinglets@subSet, ]
  ff_l <- appendCellID(ff_l, which(selectedSinglets@subSet))
  
  p <- ggplotFilterEvents(
    ffPre = OMIP021Samples[[1]],
    ffPost = ff_l,
    seed = 1,
    xChannel = "SSC-A", xLinearRange = linRange,
    yChannel = "SSC-H", yLinearRange = linRange)
  
  vdiffr::expect_doppelganger("singletsGate selected channels", fig = p)
  
  # test application of two gates one after the other
  singletsGate1 <- singletsGate(OMIP021Samples[[1]], nmad = 3)
  singletsGate2 <- singletsGate(OMIP021Samples[[1]],
                                channel1 = "SSC-A",
                                channel2 = "SSC-H",
                                filterId = "Singlets2")
  
  singletCombinedGate <- singletsGate1 & singletsGate2
  
  selectedSinglets <- flowCore::filter(OMIP021Samples[[1]], 
                                       singletCombinedGate)
  
  ff_l <- OMIP021Samples[[1]][selectedSinglets@subSet, ]
  ff_l <- appendCellID(ff_l, which(selectedSinglets@subSet))
  
  p1 <- ggplotFilterEvents(
    ffPre = OMIP021Samples[[1]],
    ffPost = ff_l,
    seed = 1,
    xChannel = "FSC-A", xLinearRange = linRange,
    yChannel = "FSC-H", yLinearRange = linRange)
  
  p2 <- ggplotFilterEvents(
    ffPre = OMIP021Samples[[1]],
    ffPost = ff_l,
    seed = 1,
    xChannel = "SSC-A", xLinearRange = linRange,
    yChannel = "SSC-H", yLinearRange = linRange)
  
  
  vdiffr::expect_doppelganger("singletsGates one after the other - fig1", 
                              fig = p1)
  vdiffr::expect_doppelganger("singletsGates one after the other - fig2", 
                              fig = p2)
  
})
