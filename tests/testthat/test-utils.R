# CytoPipeline - Copyright (C) <2022>
# <Université catholique de Louvain (UCLouvain), Belgique>
#
#   Description and complete License: see LICENSE file.
#
# This program (CytoPipeline) is free software:
#   you can redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (<http://www.gnu.org/licenses/>).


test_that("areSignalCols works", {
    ff <- OMIP021Samples[[1]]

    expectedRes <- c(rep(TRUE, times = 20), FALSE, FALSE)

    res <- areSignalCols(ff)
    expect_equal(res, expectedRes, ignore_attr = TRUE)
})

test_that("areFluoCols works", {
    ff <- OMIP021Samples[[1]]

    expectedRes <- c(
        rep(FALSE, times = 4),
        rep(TRUE, times = 16),
        FALSE, FALSE
    )

    res <- areFluoCols(ff)
    expect_equal(res, expectedRes, ignore_attr = TRUE)
})



test_that("subsample works", {
    ff <- OMIP021Samples[[1]]

    nSamples <- 50
    ffSub <- subsample(ff, nSamples)
    expect_equal(flowCore::nrow(ffSub), nSamples)

    # subsample with more samples than original nrow (5000)
    ffSub <- subsample(ff, 10000)

    nOriginal <- flowCore::nrow(ff)
    expect_equal(flowCore::nrow(ffSub), nOriginal)
})


test_that(".addCompensation2FluoChannelNames works", {
    ff <- OMIP021Samples[[1]]

    ff2 <- .addCompensation2FluoChannelNames(ff)

    expect_equal(
        flowCore::colnames(ff2),
        c(
            "FSC-A", "FSC-H", "SSC-A", "SSC-H",
            "Comp-450/50Violet-A", "Comp-525/50Violet-A", "Comp-540/30Violet-A",
            "Comp-585/15Violet-A", "Comp-610/20Violet-A", "Comp-670/30Violet-A",
            "Comp-670/14Red-A", "Comp-730//45Red-A", "Comp-780/60Red-A",
            "Comp-530/30Blue-A", "Comp-710/50Blue-A", "Comp-582/15Yellow-A",
            "Comp-610/20Yellow-A", "Comp-670/30Yellow-A", "Comp-710/50Yellow-A",
            "Comp-780/60Yellow-A", "Time", "Original_ID"
        )
    )

    expect_error(.addCompensation2FluoChannelNames(OMIP021Samples),
        regexp = "type not recognized"
    )
})

test_that("runCompensation works", {
    ff <- OMIP021Samples[[1]]

    compMatrix <- flowCore::spillover(ff)$SPILL
    ff1 <- flowCore::compensate(ff, spillover = compMatrix)
    ff2 <- runCompensation(ff,
        spillover = compMatrix,
        updateChannelNames = FALSE
    )

    expect_equal(ff2, ff1)

    ff2 <- runCompensation(ff, spillover = compMatrix)
    # the following avoids comparing name attributes
    expect_true(all(flowCore::exprs(ff2) == flowCore::exprs(ff1)))
    expect_equal(
        flowCore::colnames(ff2),
        c(
            "FSC-A", "FSC-H", "SSC-A", "SSC-H",
            "Comp-450/50Violet-A", "Comp-525/50Violet-A", "Comp-540/30Violet-A",
            "Comp-585/15Violet-A", "Comp-610/20Violet-A", "Comp-670/30Violet-A",
            "Comp-670/14Red-A", "Comp-730//45Red-A", "Comp-780/60Red-A",
            "Comp-530/30Blue-A", "Comp-710/50Blue-A", "Comp-582/15Yellow-A",
            "Comp-610/20Yellow-A", "Comp-670/30Yellow-A", "Comp-710/50Yellow-A",
            "Comp-780/60Yellow-A", "Time", "Original_ID"
        )
    )

    fs <- runCompensation(OMIP021Samples,
        spillover = compMatrix
    )
    ff3 <- flowCore::compensate(OMIP021Samples[[1]],
        spillover = compMatrix
    )
    # the following avoids comparing name attributes
    expect_true(all(flowCore::exprs(fs[[1]]) == flowCore::exprs(ff3)))
    expect_equal(
        flowCore::colnames(fs),
        c(
            "FSC-A", "FSC-H", "SSC-A", "SSC-H",
            "Comp-450/50Violet-A", "Comp-525/50Violet-A", "Comp-540/30Violet-A",
            "Comp-585/15Violet-A", "Comp-610/20Violet-A", "Comp-670/30Violet-A",
            "Comp-670/14Red-A", "Comp-730//45Red-A", "Comp-780/60Red-A",
            "Comp-530/30Blue-A", "Comp-710/50Blue-A", "Comp-582/15Yellow-A",
            "Comp-610/20Yellow-A", "Comp-670/30Yellow-A", "Comp-710/50Yellow-A",
            "Comp-780/60Yellow-A", "Time", "Original_ID"
        )
    )
})

test_that("aggregateAndSample works", {
    nCells <- 100
    agg <- aggregateAndSample(
        fs = OMIP021Samples,
        nTotalEvents = nCells,
        seed = 1
    )

    expect_equal(nrow(flowCore::exprs(agg)), nCells)

    ind1 <- which(flowCore::exprs(agg)[, "File"] == 1)
    expect_equal(nrow(flowCore::exprs(agg)[ind1, ]), nCells / 2)
    ind2 <- which(flowCore::exprs(agg)[, "File"] == 2)
    expect_equal(nrow(flowCore::exprs(agg)[ind2, ]), nCells / 2)
})

test_that("aggregateAndSample still works with flow frame", {
    nCells <- 100
    agg <- aggregateAndSample(fs = OMIP021Samples[[1]],
                              nTotalEvents = nCells,
                              seed = 1)
    expect_equal(nrow(flowCore::exprs(agg)), nCells)
})

test_that("getTransfoParams works", {
    # hybrid transformation list :
    # - two channels are logicle-ly transformed with automatic param estimates
    # - one channel has explicit logicle transfo with default parameters
    # - one channel has linear transformation
    # - other channels have no transformation
    translist <- flowCore::estimateLogicle(
        OMIP021Samples[[1]],
        c("450/50Violet-A", "525/50Violet-A")
    )
    translist <- c(
        translist,
        flowCore::transformList(
            "FSC-A",
            flowCore::linearTransform(
                a = 0.1,
                b = 0
            )
        ),
        flowCore::transformList(
            "540/30Violet-A",
            flowCore::logicleTransform()
        )
    )

    ret <- getTransfoParams(translist, channel = "SSC-A")
    expect_equal(is.null(ret), TRUE)

    ret <- getTransfoParams(translist, channel = "FSC-A")
    expect_equal(ret$type, "linear")
    expect_equal(ret$paramsList$a, 0.1)
    expect_equal(ret$paramsList$b, 0.)

    ret <- getTransfoParams(translist, channel = "525/50Violet-A")
    expect_equal(ret$type, "logicle")
    expect_equal(ret$paramsList$a, 0.)
    myW <- 0.26281575
    expect_equal(ret$paramsList$w, myW)
    expect_equal(ret$paramsList$m, 4.5)
    expect_equal(ret$paramsList$t, 262143)

    ret <- getTransfoParams(translist, channel = "540/30Violet-A")
    expect_equal(ret$type, "logicle")
    expect_equal(ret$paramsList$a, 0.)
    expect_equal(ret$paramsList$w, 0.5)
    expect_equal(ret$paramsList$m, 4.5)
    expect_equal(ret$paramsList$t, 262144)

    tf <- flowCore::linearTransform(a = 1.1, b = 0.2)
    otherTransList <- flowCore::transformList(
        from = "FSC-A",
        tfun = tf
    )
    ret <- getTransfoParams(otherTransList,
        channel = "FSC-A"
    )
    expect_equal(ret$type, "linear")
    expect_equal(ret$paramsList$a, 1.1)
    expect_equal(ret$paramsList$b, 0.2)
})

test_that("computeScatterChannelsLinearScale works", {
    ff <- OMIP021Samples[[1]]
    refMarker <- "APCCy7 - CD4"
    refChannel <- "780/60Red-A"

    targetFSCA <- list()
    targetFSCA$type <- "linear"
    targetFSCA$paramsList <- list()
    targetFSCA$paramsList$a <- 1.3504601e-05
    targetFSCA$paramsList$b <- 0.41438898

    targetSSCA <- list()
    targetSSCA$type <- "linear"
    targetSSCA$paramsList <- list()
    targetSSCA$paramsList$a <- 1.07290554e-05
    targetSSCA$paramsList$b <- 0.38032844

    transList <- flowCore::estimateLogicle(ff,
        channels = refChannel
    )

    # base case, reference channel exists, transList pre-filled with logicle
    # transfo for ref channel
    retTransList <-
        computeScatterChannelsLinearScale(ff,
            transList = transList,
            referenceChannel = refMarker,
            silent = TRUE
        )

    retFSCA <- getTransfoParams(retTransList, channel = "FSC-A")
    expect_equal(retFSCA, targetFSCA)

    retFSCH <- getTransfoParams(retTransList, channel = "FSC-H")
    expect_equal(retFSCH, targetFSCA)

    retSSCA <- getTransfoParams(retTransList, channel = "SSC-A")
    expect_equal(retSSCA, targetSSCA)

    retSSCH <- getTransfoParams(retTransList, channel = "SSC-H")
    expect_equal(retSSCH, targetSSCA)

    # test with a reference channel that does not exist
    expect_error(computeScatterChannelsLinearScale(
        ff,
        transList = transList,
        referenceChannel = "Yipee",
        silent = TRUE
    ), regexp = "can't find")

    # test with a reference channel that is not a fluo channel
    expect_error(computeScatterChannelsLinearScale(
        ff,
        transList = transList,
        referenceChannel = "SSC-A",
        silent = TRUE
    ), regexp = "should be a fluorochrome channel")

    # test with a NULL transList

    targetFSCA$paramsList$a <- 0.047610702
    targetFSCA$paramsList$b <- -491.62068
    targetSSCA$paramsList$a <- 0.037825469
    targetSSCA$paramsList$b <- -611.70173
    retTransList <-
        computeScatterChannelsLinearScale(ff,
            transList = NULL,
            referenceChannel = refMarker,
            silent = TRUE
        )

    retFSCA <- getTransfoParams(retTransList, channel = "FSC-A")
    expect_equal(retFSCA, targetFSCA)

    retFSCH <- getTransfoParams(retTransList, channel = "FSC-H")
    expect_equal(retFSCH, targetFSCA)

    retSSCA <- getTransfoParams(retTransList, channel = "SSC-A")
    expect_equal(retSSCA, targetSSCA)

    retSSCH <- getTransfoParams(retTransList, channel = "SSC-H")
    expect_equal(retSSCH, targetSSCA)

    # test with a transList that does not contain the reference channel,
    # but well an already existing SSC-A transformation
    linTrans <- flowCore::linearTransform()
    stupidTransList <- flowCore::transformList(
        from = "SSC-A",
        tfun = linTrans
    )
    retTransList <-
        computeScatterChannelsLinearScale(ff,
            transList = stupidTransList,
            referenceChannel = refMarker,
            silent = TRUE
        )


    retFSCA <- getTransfoParams(retTransList, channel = "FSC-A")
    expect_equal(retFSCA, targetFSCA)

    retFSCH <- getTransfoParams(retTransList, channel = "FSC-H")
    expect_equal(retFSCH, targetFSCA)

    retSSCA <- getTransfoParams(retTransList, channel = "SSC-A")
    expect_equal(retSSCA, targetSSCA)

    retSSCH <- getTransfoParams(retTransList, channel = "SSC-H")
    expect_equal(retSSCH, targetSSCA)
})

test_that("findTimeChannel works", {
    # with flow set
    ret <- findTimeChannel(OMIP021Samples)
    expect_equal(ret, "Time")

    # with flow frame
    ret2 <- findTimeChannel(OMIP021Samples[[1]])
    expect_equal(ret2, "Time")

    # test exclude channels parameter
    ret3 <- findTimeChannel(OMIP021Samples[[1]],
        excludeChannels = "Time"
    )
    expect_null(ret3)
})

test_that("getChannelNamesFromMarkers works", {
    # with existing markers
    ret <- getChannelNamesFromMarkers(
        OMIP021Samples[[1]],
        c(
            "FSC-A",
            "L/D Aqua - Viability",
            "FITC - gdTCR",
            "PECy5 - CD28"
        )
    )

    expected <- c("FSC-A", "525/50Violet-A", "530/30Blue-A", "670/30Yellow-A")
    expect_equal(ret, expected)

    # with boolean vector
    indices <- c(1, 6, 14, 18)
    boolInput <- rep(FALSE, 21)
    boolInput[indices] <- TRUE
    ret <- getChannelNamesFromMarkers(
        OMIP021Samples[[1]],
        boolInput
    )

    expect_equal(ret, expected)


    # with indices vector
    ret <- getChannelNamesFromMarkers(
        OMIP021Samples[[1]],
        indices
    )
    expect_equal(ret, expected)

    # missing channel
    missingInput <- c("CD4")

    expect_error(getChannelNamesFromMarkers(
        OMIP021Samples[[1]],
        missingInput
    ),
    regexp = "could not be found"
    )
})

test_that("updateMarkerName works", {
    retFF <- updateMarkerName(OMIP021Samples[[1]],
                            channel = "FSC-A",
                            newMarkerName = "Fwd Scatter-A")
    
    checkMkName <- getChannelNamesFromMarkers(retFF, markers = "Fwd Scatter-A")
    expect_equal(checkMkName, "FSC-A")
    
    expect_error(updateMarkerName(OMIP021Samples[[1]], channel = "FFF-A", 
                                  newMarkerName = "whatever"), 
                 regexp = "channel not found")
})

test_that("removeChannels works", {
    retFF <- removeChannels(OMIP021Samples[[1]],
                            channels = c("FSC-A", "SSC-A"))
    nRemainingChannels <- length(flowCore::colnames(retFF))
    expect_equal(nRemainingChannels, 20)
    
    expect_warning(removeChannels(OMIP021Samples[[1]],
                                  channels = c("FSC-A", "SSC-S")),
                   regexp = "not found in flowFrame")
    
    
})