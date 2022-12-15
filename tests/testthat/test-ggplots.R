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



# obtain OMIP021UTSamples, light-weight version used specifically for these 
# unit tests
path <- system.file("scripts",
                    package = "CytoPipeline"
)

source(paste0(path,"/MakeOMIP021UTSamples.R"))


test_that("ggplotFlowRate works", {
    p <- ggplotFlowRate(OMIP021UTSamples[[1]])
    vdiffr::expect_doppelganger("ggplotFlowRate single", fig = p)

    p <- ggplotFlowRate(OMIP021UTSamples)
    vdiffr::expect_doppelganger("ggplotFlowRate double", fig = p)

    p <- ggplotFlowRate(OMIP021UTSamples[[1]], title = "Test Flow Rate plot")
    vdiffr::expect_doppelganger("ggplotFlowRate single with title", fig = p)

    p <- ggplotFlowRate(OMIP021UTSamples[[1]], timeUnit = 50)
    vdiffr::expect_doppelganger(
        "ggplotFlowRate single with time unit",
        fig = p
    )
})



test_that("ggplotEvents with 1D works", {
    expect_error(ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "450/50Violet-A",
        xScale = "bi-exponential"
    ),
    regexp = "should be one of"
    )

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "FSC-A",
        xScale = "linear"
    )
    vdiffr::expect_doppelganger("ggplotEvents 1D linear single", fig = p)

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "FSC-A",
        xScale = "linear",
        xLinearRange = c(0, 250000)
    )
    vdiffr::expect_doppelganger(
        "ggplotEvents 1D linear explicit range - single",
        fig = p
    )

    p <- ggplotEvents(OMIP021UTSamples, xChannel = "FSC-A", xScale = "linear")
    vdiffr::expect_doppelganger("ggplotEvents 1D linear double", fig = p)

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "450/50Violet-A",
        xScale = "logicle"
    )
    vdiffr::expect_doppelganger("ggplotEvents 1D logicle single", fig = p)

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "450/50Violet-A",
        xScale = "logicle", xLogicleParams = list(
            a = 1,
            w = 2,
            m = 7,
            t = 270000
        )
    )

    vdiffr::expect_doppelganger(
        "ggplotEvents 1D logicle with explicit params - single",
        fig = p
    )

    p <- ggplotEvents(OMIP021UTSamples[[2]],
        xChannel = "450/50Violet-A",
        xScale = "logicle", nDisplayCells = 500, seed = 1
    )

    vdiffr::expect_doppelganger("ggplotEvents 1D sub-sampling", fig = p)

    p <- ggplotEvents(OMIP021UTSamples[[2]],
        xChannel = "450/50Violet-A",
        xScale = "logicle", alpha = 0.5, fill = "red"
    )

    vdiffr::expect_doppelganger("ggplotEvents 1D fill and color", fig = p)

    compensationMatrix <- flowCore::spillover(OMIP021UTSamples[[1]])$SPILL

    ffC <- runCompensation(OMIP021UTSamples[[1]],
        spillover = compensationMatrix,
        updateChannelNames = FALSE
    )

    transList <- flowCore::estimateLogicle(
        ffC,
        colnames(compensationMatrix)
    )

    transList <-
        c(
            transList,
            flowCore::transformList(
                "FSC-A",
                flowCore::linearTransform(a = 0.00001)
            )
        )

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "450/50Violet-A",
        xScale = "linear", transList = transList,
        runTransforms = FALSE
    )
    vdiffr::expect_doppelganger(
        "ggplotEvents 1D transformList logicle not run",
        fig = p
    )

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "450/50Violet-A",
        xScale = "linear", transList = transList,
        runTransforms = TRUE
    )
    vdiffr::expect_doppelganger("ggplotEvents 1D transformList logicle run",
        fig = p
    )

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "FSC-A",
        xScale = "logicle", transList = transList,
        runTransforms = FALSE
    )
    vdiffr::expect_doppelganger(
        "ggplotEvents 1D transformList linear not run",
        fig = p
    )

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "FSC-A",
        xScale = "logicle", transList = transList,
        runTransforms = TRUE
    )
    vdiffr::expect_doppelganger("ggplotEvents 1D transformList linear run",
        fig = p
    )
})

test_that("ggplotEvents with 2D works", {
    p <- ggplotEvents(OMIP021UTSamples,
        xChannel = "FSC-A",
        xScale = "linear",
        yChannel = "SSC-A",
        yScale = "linear"
    )
    vdiffr::expect_doppelganger("ggplotEvents 2D linear double", fig = p)

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "FSC-A",
        xScale = "linear",
        yChannel = "610/20Violet-A",
        yScale = "logicle"
    )
    vdiffr::expect_doppelganger("ggplotEvents 2D x linear y logicle", fig = p)

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "450/50Violet-A",
        xScale = "logicle",
        yChannel = "SSC-A",
        yScale = "linear"
    )
    vdiffr::expect_doppelganger("ggplotEvents 2D x logicle y linear", fig = p)

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "TETaGC",
        xScale = "logicle",
        yChannel = "CD27",
        yScale = "logicle"
    )
    vdiffr::expect_doppelganger(
        "ggplotEvents 2D x logicle y logicle by markers",
        fig = p
    )

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "TETaGC",
        xScale = "logicle",
        yChannel = "CD27",
        yScale = "logicle",
        bins = 128
    )
    vdiffr::expect_doppelganger(
        "ggplotEvents 2D x logicle y logicle bins",
        fig = p
    )


    compensationMatrix <- flowCore::spillover(OMIP021UTSamples[[1]])$SPILL

    ffC <- runCompensation(OMIP021UTSamples[[1]],
        spillover = compensationMatrix,
        updateChannelNames = FALSE
    )

    transList <- flowCore::estimateLogicle(
        ffC,
        colnames(compensationMatrix)
    )

    transList <-
        c(
            transList,
            flowCore::transformList(
                "FSC-A",
                flowCore::linearTransform(a = 0.00001)
            )
        )

    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "TETaGC",
        xScale = "logicle",
        yChannel = "CD27",
        yScale = "logicle",
        transList = transList,
        runTransforms = FALSE
    )
    vdiffr::expect_doppelganger(
        "ggplotEvents 2D x logicle y logicle transList not run",
        fig = p
    )
    p <- ggplotEvents(OMIP021UTSamples[[1]],
        xChannel = "TETaGC",
        xScale = "logicle",
        yChannel = "CD27",
        yScale = "logicle",
        transList = transList,
        runTransforms = TRUE
    )
    vdiffr::expect_doppelganger(
        "ggplotEvents 2D x logicle y logicle transList run",
        fig = p
    )
})


test_that("ggplotFilterEvents works", {
    ffPre <- OMIP021UTSamples[[1]]

    LDMarker <- "L/D Aqua - Viability"

    LDChannel <- getChannelNamesFromMarkers(ffPre, markers = LDMarker)
    liveGateMatrix <- matrix(
        data = c(
            50000, 50000, 100000, 200000, 200000,
            100, 1000, 2000, 2000, 1
        ),
        ncol = 2,
        dimnames = list(
            c(),
            c("FSC-A", LDChannel)
        )
    )

    liveGate <- flowCore::polygonGate(
        filterId = "Live",
        .gate = liveGateMatrix
    )

    selectedLive <- flowCore::filter(ffPre, liveGate)
    ffL <- flowCore::Subset(ffPre, selectedLive)
    #ffL <- .appendCellID(ffL, which(selectedLive@subSet))

    p <- ggplotFilterEvents(
        ffPre = ffPre,
        ffPost = ffL,
        nDisplayCells = 1000,
        xChannel = "FSC-A", xScale = "linear",
        yChannel = LDMarker, yScale = "logicle"
    ) +
        ggtitle("Live gate filter")
    vdiffr::expect_doppelganger(
        "ggplotFilterEvents 2D - all points",
        fig = p
    )

    p <- ggplotFilterEvents(
        ffPre = ffPre,
        ffPost = ffL,
        nDisplayCells = 1000,
        xChannel = "FSC-A", xScale = "linear",
        yChannel = LDMarker, yScale = "logicle",
        interactive = TRUE
    ) +
        ggtitle("Live gate filter")
    vdiffr::expect_doppelganger(
        "ggplotFilterEvents 2D - all points - interactive",
        fig = p
    )

    p <- ggplotFilterEvents(
        ffPre = ffPre,
        ffPost = ffL,
        nDisplayCells = 500,
        seed = 1,
        xChannel = "FSC-A", xScale = "linear",
        yChannel = LDMarker, yScale = "logicle"
    ) +
        ggtitle("Live gate filter")
    vdiffr::expect_doppelganger(
        "ggplotFilterEvents 2D - 5000 points",
        fig = p
    )

    p <- ggplotFilterEvents(
        ffPre = ffPre,
        ffPost = ffL,
        nDisplayCells = Inf,
        xChannel = "FSC-A", xScale = "linear",
        yChannel = LDMarker, yScale = "logicle"
    ) +
        ggtitle("Live gate filter")
    vdiffr::expect_doppelganger(
        "ggplotFilterEvents 2D - Infinite nb of points",
        fig = p
    )

    p <- ggplotFilterEvents(
        size = 0.1,
        ffPre = ffPre,
        ffPost = ffL,
        seed = 1,
        xChannel = "FSC-A", xScale = "linear",
        yChannel = LDMarker, yScale = "logicle"
    ) +
        ggtitle("Live gate filter")
    vdiffr::expect_doppelganger(
        "ggplotFilterEvents 2D - small size",
        fig = p
    )
    
    # calculate transformation list for the next steps
    compensationMatrix <- flowCore::spillover(OMIP021UTSamples[[1]])$SPILL
    
    ffC <- runCompensation(OMIP021UTSamples[[1]],
                           spillover = compensationMatrix,
                           updateChannelNames = FALSE
    )
    
    transList <- flowCore::estimateLogicle(
        ffC,
        colnames(compensationMatrix)
    )
    
    transList <-
        c(
            transList,
            flowCore::transformList(
                "FSC-A",
                flowCore::linearTransform(a = 0.00001)
            )
        )
    
    p <- ggplotFilterEvents(
        ffPre = ffPre,
        ffPost = ffL,
        seed = 1,
        xChannel = "FSC-A", 
        yChannel = LDMarker, 
        transList = transList,
        runTransforms = FALSE
    ) +
        ggtitle("Live gate filter")
    
    vdiffr::expect_doppelganger(
        "ggplotFilterEvents 2D - transList - not run",
        fig = p
    )
    
    p <- ggplotFilterEvents(
        ffPre = ffPre,
        ffPost = ffL,
        seed = 1,
        xChannel = "FSC-A", 
        yChannel = LDMarker, 
        transList = transList,
        runTransforms = TRUE
    ) +
        ggtitle("Live gate filter")
    
    vdiffr::expect_doppelganger(
        "ggplotFilterEvents 2D - transList - run",
        fig = p
    )
})
