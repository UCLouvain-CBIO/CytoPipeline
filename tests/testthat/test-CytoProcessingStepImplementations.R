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


test_that("estimateScaleTransforms work", {
    compMatrix <- flowCore::spillover(OMIP021UTSamples[[1]])$SPILL
    ff_c <- runCompensation(OMIP021UTSamples[[1]], spillover = compMatrix)

    transList <-
        suppressMessages(estimateScaleTransforms(
            ff = ff_c,
            fluoMethod = "estimateLogicle",
            scatterMethod = "linear",
            scatterRefMarker = "BV785 - CD3"
        ))

    transListPath <- paste0(system.file("extdata", 
                                        package = "CytoPipeline"),
                            "/OMIP021_TransList.rds")
    refTransList <- readRDS(transListPath)

    # saveRDS(transList, transListPath)

    refFF <- flowCore::transform(ff_c, refTransList)
    thisFF <- flowCore::transform(ff_c, transList)

    expect_equal(
        flowCore::exprs(thisFF),
        flowCore::exprs(refFF)
    )
})




test_that("readSampleFiles works", {
    rawDataDir <-
        paste0(system.file("extdata", package = "CytoPipeline"), "/")
    sampleFiles <-
        paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))

    truncateMaxRange <- FALSE
    minLimit <- NULL

    fs_raw <-
        flowCore::read.flowSet(sampleFiles,
            truncate_max_range = truncateMaxRange,
            min.limit = minLimit
        )
    fs_raw <- flowCore::fsApply(fs_raw, FUN = .appendCellID)

    res <- readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = "all",
        truncate_max_range = truncateMaxRange,
        min.limit = minLimit
    )
    expect_equal(
        flowCore::exprs(res[[1]]),
        flowCore::exprs(fs_raw[[1]])
    )
    expect_equal(
        flowCore::exprs(res[[2]]),
        flowCore::exprs(fs_raw[[2]])
    )

    res2 <- readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = 2,
        truncate_max_range = truncateMaxRange,
        min.limit = minLimit
    )

    expect_equal(
        flowCore::exprs(res2),
        flowCore::exprs(fs_raw[[2]])
    )
})

test_that("removeMarginsPeacoQC works", {
    fs_raw <- OMIP021UTSamples

    ff_m <-
        suppressWarnings(removeMarginsPeacoQC(x = fs_raw[[1]]))

    ref_ff_m <- readRDS(test_path("fixtures", "ff_m.rds"))

    # saveRDS(ff_m, test_path("fixtures", "ff_m.rds"))

    expect_equal(
        flowCore::exprs(ff_m),
        flowCore::exprs(ref_ff_m)
    )
})


test_that("compensateFromMatrix works", {
    ref_ff_m <- readRDS(test_path("fixtures", "ff_m.rds"))

    ff_c <-
        compensateFromMatrix(ref_ff_m,
            matrixSource = "fcs"
        )

    ref_ff_c <- readRDS(test_path("fixtures", "ff_c.rds"))

    # saveRDS(ff_c, test_path("fixtures", "ff_c.rds"))

    expect_equal(
        flowCore::exprs(ff_c),
        flowCore::exprs(ref_ff_c)
    )
})

test_that("removeDoubletsFlowStats works", {
    ref_ff_c <- readRDS(test_path("fixtures", "ff_c.rds"))


    ff_s <-
        removeDoubletsFlowStats(ref_ff_c,
            areaChannels = c("FSC-A", "SSC-A"),
            heightChannels = c("FSC-H", "SSC-H"),
            widerGate = TRUE
        )

    ref_ff_s <- readRDS(test_path("fixtures", "ff_s.rds"))

    # saveRDS(ff_s, test_path("fixtures", "ff_s.rds"))

    expect_equal(
        flowCore::exprs(ff_s),
        flowCore::exprs(ref_ff_s)
    )
})


test_that("removeDoubletsPeacoQC works", {
    ref_ff_c <- readRDS(test_path("fixtures", "ff_c.rds"))

    ff_s2 <-
        removeDoubletsPeacoQC(ref_ff_c,
            areaChannels = c("FSC-A", "SSC-A"),
            heightChannels = c("FSC-H", "SSC-H"),
            nmads = c(3, 5)
        )

    ref_ff_s2 <- readRDS(test_path("fixtures", "ff_s2.rds"))

    # saveRDS(ff_s2, test_path("fixtures", "ff_s2.rds"))

    expect_equal(
        flowCore::exprs(ff_s2),
        flowCore::exprs(ref_ff_s2)
    )
})

test_that("removeDoubletsCytoPipeline works", {
    ref_ff_c <- readRDS(test_path("fixtures", "ff_c.rds"))

    ff_s3 <-
        removeDoubletsCytoPipeline(ref_ff_c,
            areaChannels = c("FSC-A", "SSC-A"),
            heightChannels = c("FSC-H", "SSC-H"),
            nmads = c(3, 5)
        )


    ref_ff_s3 <- readRDS(test_path("fixtures", "ff_s3.rds"))

    # saveRDS(ff_s3, test_path("fixtures", "ff_s3.rds"))

    expect_equal(
        flowCore::exprs(ff_s3),
        flowCore::exprs(ref_ff_s3)
    )
})

test_that("removeDebrisFlowClustTmix works", {
    ref_ff_s <- readRDS(test_path("fixtures", "ff_s.rds"))

    ff_cells <-
        removeDebrisFlowClustTmix(ref_ff_s,
            FSCChannel = "FSC-A",
            SSCChannel = "SSC-A",
            nClust = 3,
            level = 0.97,
            B = 100
        )

    ref_ff_cells <- readRDS(test_path("fixtures", "ff_cells.rds"))

    # saveRDS(ff_cells, test_path("fixtures", "ff_cells.rds"))

    expect_equal(
        flowCore::exprs(ff_cells),
        flowCore::exprs(ref_ff_cells)
    )
})

test_that("removeDebrisManualGate works", {
    ref_ff_s <- readRDS(test_path("fixtures", "ff_s.rds"))
    
    remDebrisGateData <- c(73615, 110174, 213000, 201000, 126000,
                           47679, 260500, 260500, 113000, 35000)
    
    ff_cells <-
        removeDebrisManualGate(ref_ff_s,
                               FSCChannel = "FSC-A",
                               SSCChannel = "SSC-A",
                               gateData = remDebrisGateData
        )
    
    ref_ff_cells <- readRDS(test_path("fixtures", "ff_cells_manual_gate.rds"))
    
    #saveRDS(ff_cells, test_path("fixtures", "ff_cells_manual_gate.rds"))
    
    expect_equal(
        flowCore::exprs(ff_cells),
        flowCore::exprs(ref_ff_cells)
    )
})

test_that("removeDeadCellsGateTail works", {
    ref_ff_cells <- readRDS(test_path("fixtures", "ff_cells.rds"))

    transListPath <- paste0(system.file("extdata", 
                                        package = "CytoPipeline"),
                            "/OMIP021_TransList.rds")
    refTransList <- readRDS(transListPath)

    ff_lcells <-
        removeDeadCellsGateTail(ref_ff_cells,
            preTransform = TRUE,
            transList = refTransList,
            LDMarker = "L/D Aqua - Viability",
            num_peaks = 2,
            ref_peak = 2,
            strict = FALSE,
            positive = FALSE
        )

    ref_ff_lcells <- readRDS(test_path("fixtures", "ff_lcells.rds"))

    # saveRDS(ff_lcells, test_path("fixtures", "ff_lcells.rds"))

    expect_equal(
        flowCore::exprs(ff_lcells),
        flowCore::exprs(ref_ff_lcells)
    )
})

test_that("removeDeadCellsManualGate works", {
    ref_ff_cells <- readRDS(test_path("fixtures", "ff_cells.rds"))
    
    remDeadCellsGateData <- c(0, 0, 250000, 250000,
                              0, 650, 650, 0)
    
    ff_lcells <-
        removeDeadCellsManualGate(ref_ff_cells,
                                  FSCChannel = "FSC-A",
                                  LDMarker = "L/D Aqua - Viability",
                                  gateData = remDeadCellsGateData)
    
    ref_ff_lcells <- readRDS(test_path("fixtures", "ff_lcells_manual_gate.rds"))
    
    # saveRDS(ff_lcells, test_path("fixtures", "ff_lcells_manual_gate.rds"))
    
    expect_equal(
        flowCore::exprs(ff_lcells),
        flowCore::exprs(ref_ff_lcells)
    )
})

test_that("qualityControlPeacoQC", {
    ref_ff_lcells <- readRDS(test_path("fixtures", "ff_lcells.rds"))

    transListPath <- paste0(system.file("extdata", 
                                        package = "CytoPipeline"),
                            "/OMIP021_TransList.rds")
    refTransList <- readRDS(transListPath)

    expect_error(suppressWarnings(qualityControlPeacoQC(
            ref_ff_lcells,
            preTransform = TRUE,
            transList = refTransList,
            min_cells = 30, # decreased from the 150 default
            max_bins = 500, # default
            MAD = 6, # default
            IT_limit = 0.55, # default
            force_IT = 150, # default
            peak_removal = (1 / 3), # default
            min_nr_bins_peakdetection = 10 # default
    )), regexp = "must be strictly positive and finite")


    # ref_ff_qualityControl <-
    #     readRDS(test_path("fixtures", "ff_QC_PeacoQC.rds"))
    # 
    # # saveRDS(ff_QualityControl, test_path("fixtures", "ff_QC_PeacoQC.rds"))
    # 
    # expect_equal(
    #     flowCore::exprs(ff_QualityControl),
    #     flowCore::exprs(ref_ff_qualityControl)
    # )
})




test_that("qualityControlFlowAI works", {
    fs_raw <- OMIP021UTSamples

    ff_QualityControl <- suppressWarnings(
        qualityControlFlowAI(fs_raw[[1]],
            remove_from = "all", # all default
            second_fractionFR = 0.1,
            deviationFR = "MAD",
            alphaFR = 0.01,
            decompFR = TRUE,
            outlier_binsFS = FALSE,
            pen_valueFS = 500,
            max_cptFS = 3,
            sideFM = "both",
            neg_valuesFM = 1
        )
    )

    ref_ff_qualityControl_flowAI <-
        readRDS(test_path("fixtures", "ff_QC_flowAI.rds"))

    # saveRDS(ff_QualityControl, test_path("fixtures", "ff_QC_flowAI.rds"))

    expect_equal(
        flowCore::exprs(ff_QualityControl),
        flowCore::exprs(ref_ff_qualityControl_flowAI)
    )
})

test_that("qualityControlFlowCut works", {
    fs_raw <- OMIP021UTSamples

    ff_QualityControl <- suppressMessages(
        qualityControlFlowCut(fs_raw[[1]],
            MaxContin = 0.1,
            MeanOfMeans = 0.13,
            MaxOfMeans = 0.15,
            MaxValleyHgt = 0.1,
            MaxPercCut = 0.3,
            LowDensityRemoval = 0.1,
            RemoveMultiSD = 7,
            AlwaysClean = FALSE,
            IgnoreMonotonic = FALSE,
            MonotonicFix = NULL,
            Measures = c(1:8)
        )
    )


    ref_ff_qualityControl_flowCut <-
        readRDS(test_path("fixtures", "ff_QC_flowCut.rds"))

    # saveRDS(ff_QualityControl, test_path("fixtures", "ff_QC_flowCut.rds"))

    expect_equal(
        flowCore::exprs(ff_QualityControl),
        flowCore::exprs(ref_ff_qualityControl_flowCut)
    )
})

test_that("qualityControlFlowClean works", {
    fs_raw <- OMIP021UTSamples

    ff_QualityControl <- suppressWarnings(
        qualityControlFlowClean(fs_raw[[1]],
            binSize = 0.01, # default
            nCellCutoff = 500, # default
            cutoff = "median", # default
            fcMax = 1.3, # default
            nstable = 5
        )
    )

    ref_ff_qualityControl_flowClean <-
        readRDS(test_path("fixtures", "ff_QC_flowClean.rds"))

    # saveRDS(ff_QualityControl, test_path("fixtures", "ff_QC_flowClean.rds"))

    expect_equal(
        flowCore::exprs(ff_QualityControl),
        flowCore::exprs(ref_ff_qualityControl_flowClean)
    )
})
