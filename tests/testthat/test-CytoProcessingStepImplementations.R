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

source(file.path(path,"MakeOMIP021UTSamples.R"))


test_that("estimateScaleTransforms work", {
    compMatrix <- flowCore::spillover(OMIP021UTSamples[[1]])$SPILL
    ff_c <- runCompensation(OMIP021UTSamples[[1]], spillover = compMatrix)

    transList <-
        suppressMessages(estimateScaleTransforms(
            ff = ff_c,
            fluoMethod = "estimateLogicle",
            scatterMethod = "linearQuantile",
            scatterRefMarker = "BV785 - CD3"
        ))

    transListPath <- file.path(system.file("extdata", 
                                           package = "CytoPipeline"),
                               "OMIP021_TransList.rds")
    refTransList <- readRDS(transListPath)

    # saveRDS(transList, transListPath)

    refFF <- flowCore::transform(ff_c, refTransList)
    thisFF <- flowCore::transform(ff_c, transList)

    expect_equal(
        flowCore::exprs(thisFF),
        flowCore::exprs(refFF)
    )
    
    transList2 <- estimateScaleTransforms(
        ff = ff_c,
        fluoMethod = "estimateLogicle",
        scatterMethod = "linearQuantile",
        scatterRefMarker = "BV785 - CD3",
        specificScatterChannels = c("SSC-A", "SSC-H")
    )
})

# test_that("selectRandomSamples works", {
#     rawDataDir <-
#         system.file("extdata", package = "CytoPipeline")
#     sampleFiles <-
#         file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
#     
#     seed <- 2
#     nSamples <- 1
#     
#     newSampleFiles <- selectRandomSamples(sampleFiles, 
#                                           nSamples = nSamples,
#                                           seed = seed)
#     expected <- sampleFiles[1]
#     expect_equal(newSampleFiles,
#                  expected)
#     
#     newSampleFiles <- selectRandomSamples(sampleFiles, 
#                                           nSamples = 3,
#                                           seed = seed)
#     
#     expect_equal(newSampleFiles,
#                  sampleFiles)
#     
#     expect_error(selectRandomSamples(sampleFiles, 
#                                      nSamples = 0,
#                                      seed = seed),
#                  regexp = "should be a numeric >= 1")
#                                           
#         
# })

test_that("readSampleFiles with no sample files gives an error", {
    rawDataDir <-
        system.file("extdata", package = "CytoPipeline")
    sampleFiles <-
        file.path(rawDataDir, list.files(rawDataDir, pattern = "NoGood"))
    
    expect_error(
        readSampleFiles(
            sampleFiles = sampleFiles,
            whichSamples = "all",
            truncate_max_range = truncateMaxRange,
            min.limit = minLimit
        ), regexp = "no sample files to read")
})



test_that("readSampleFiles works", {
    rawDataDir <-
        system.file("extdata", package = "CytoPipeline")
    sampleFiles <-
        file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))

    truncateMaxRange <- FALSE
    minLimit <- NULL

    fs_raw <-
        flowCore::read.flowSet(sampleFiles,
            truncate_max_range = truncateMaxRange,
            min.limit = minLimit
        )
    fs_raw <- flowCore::fsApply(fs_raw, FUN = appendCellID)

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

test_that("readSampleFiles with random samples works", {
    rawDataDir <-
        system.file("extdata", package = "CytoPipeline")
    sampleFiles <-
        file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
    
    truncateMaxRange <- FALSE
    minLimit <- NULL
    
    fs_raw <-
        flowCore::read.flowSet(sampleFiles,
                               truncate_max_range = truncateMaxRange,
                               min.limit = minLimit
        )
    fs_raw <- flowCore::fsApply(fs_raw, FUN = appendCellID)
    
    seed <- 2
    nSamples <- 1
    
    res <- readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = "random",
        nSamples = nSamples,
        seed = seed,
        truncate_max_range = truncateMaxRange,
        min.limit = minLimit
    )
    
    expect_equal(
        flowCore::exprs(res),
        flowCore::exprs(fs_raw[[1]])
    )
    
    res <- readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = "random",
        nSamples = 3,
        seed = seed,
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
    
 
    expect_error(readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = "random",
        nSamples = 0,
        seed = seed,
        truncate_max_range = truncateMaxRange,
        min.limit = minLimit),
        regexp = "should be a numeric >= 1")
})

test_that("readSampleFiles with post-processing works", {
    rawDataDir <-
        system.file("extdata", package = "CytoPipeline")
    sampleFiles <-
        file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
    
    channelMarkerFile <-
        system.file("extdata/ChannelMarkerUsed.csv", package = "CytoPipeline")
    
    truncateMaxRange <- FALSE
    minLimit <- NULL
    
    fs_raw <-
        flowCore::read.flowSet(sampleFiles,
                               truncate_max_range = truncateMaxRange,
                               min.limit = minLimit
        )
    fs_raw <- flowCore::fsApply(fs_raw, FUN = appendCellID)
    
    res <- readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = "all",
        truncate_max_range = truncateMaxRange,
        min.limit = minLimit,
        channelMarkerFile = channelMarkerFile)
    
    expectedMarkerNames <- 
        flowCore::pData(flowCore::parameters(fs_raw[[1]]))$desc
    expectedMarkerNames[6] <- "Viability"
    expectedMarkerNames <- expectedMarkerNames[-5]
    
    expect_equal(
        flowCore::pData(flowCore::parameters(res[[1]]))$desc,
        expectedMarkerNames)
    expect_equal(
        flowCore::pData(flowCore::parameters(res[[2]]))$desc,
        expectedMarkerNames)
    
    res2 <- readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = 2,
        truncate_max_range = truncateMaxRange,
        min.limit = minLimit,
        channelMarkerFile = channelMarkerFile)
    
    expect_equal(
        flowCore::pData(flowCore::parameters(res2))$desc,
        expectedMarkerNames
    )
})

test_that("readSampleFiles with phenoData works", {
    rawDataDir <-
        system.file("extdata", package = "CytoPipeline")
    sampleFiles <-
        file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
    
    truncateMaxRange <- FALSE
    minLimit <- NULL
    
    phenoData <- data.frame(Donor = c(1, 2),
                            Group = c("G1", "G2"),
                            row.names = basename(sampleFiles))
    
    res <- readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = "all",
        truncate_max_range = truncateMaxRange,
        min.limit = minLimit,
        pData = phenoData)
    
    retPData <- flowCore::pData(res)
    expectedPData <- phenoData
    expectedPData$name <- c("Donor1.fcs", "Donor2.fcs")
    
    expect_equal(retPData, expectedPData)
    
    expect_equal(
        flowCore::keyword(res[[1]], "CytoPipeline_Donor")[["CytoPipeline_Donor"]],
        1)
    
    expect_equal(
        flowCore::keyword(res[[1]], "CytoPipeline_Group")[["CytoPipeline_Group"]],
        "G1")
    
    expect_equal(
        flowCore::keyword(res[[2]], "CytoPipeline_Donor")[["CytoPipeline_Donor"]],
        2)
    
    expect_equal(
        flowCore::keyword(res[[2]], "CytoPipeline_Group")[["CytoPipeline_Group"]],
        "G2")
    
    res2 <- readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = 2,
        truncate_max_range = truncateMaxRange,
        min.limit = minLimit,
        pData = phenoData)
    
    expect_equal(
        flowCore::keyword(res2, "CytoPipeline_Donor")[["CytoPipeline_Donor"]],
        2)
    
    expect_equal(
        flowCore::keyword(res2, "CytoPipeline_Group")[["CytoPipeline_Group"]],
        "G2")
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
    
    channelSpecifications <-
        list( "450/50Violet-A" = c(0, 262143),
              "PECy5 - CD28" = c(0, 20000))
    
    ff_m2 <-
        suppressWarnings(removeMarginsPeacoQC(x = fs_raw[[1]],
                         channelSpecifications = channelSpecifications))

    ref_ff_m2 <- readRDS(test_path("fixtures", "ff_m2.rds"))

    # saveRDS(ff_m2, test_path("fixtures", "ff_m2.rds"))

    expect_equal(
        flowCore::exprs(ff_m2),
        flowCore::exprs(ref_ff_m2)
    )

    channelSpecifications <-
        list( "AllFluoChannels" = c(0, 40000))

    ff_m3 <-
        suppressWarnings(removeMarginsPeacoQC(x = fs_raw[[1]],
                         channelSpecifications = channelSpecifications))

    ref_ff_m3 <- readRDS(test_path("fixtures", "ff_m3.rds"))

    # saveRDS(ff_m3, test_path("fixtures", "ff_m3.rds"))

    expect_equal(
        flowCore::exprs(ff_m3),
        flowCore::exprs(ref_ff_m3)
    )

    channelSpecifications <-
        list( "AllFluoChannels" = c(0, 40000),
              "450/50Violet-A" = c(0, 262143),
              "PECy5 - CD28" = c(0, 20000))

    ff_m4 <-
        suppressWarnings(removeMarginsPeacoQC(x = fs_raw[[1]],
                         channelSpecifications = channelSpecifications))

    ref_ff_m4 <- readRDS(test_path("fixtures", "ff_m4.rds"))

    # saveRDS(ff_m4, test_path("fixtures", "ff_m4.rds"))

    expect_equal(
        flowCore::exprs(ff_m4),
        flowCore::exprs(ref_ff_m4)
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

test_that("compensateFromMatrix works with pData", {
    rawDataDir <-
        system.file("extdata", package = "CytoPipeline")
    sampleFiles <-
        file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
    
    truncateMaxRange <- FALSE
    minLimit <- NULL
    
    phenoData <- data.frame(Donor = c(1, 2),
                            Group = c("G1", "G2"),
                            row.names = basename(sampleFiles))
    
    fs <- readSampleFiles(
        sampleFiles = sampleFiles,
        whichSamples = "all",
        truncate_max_range = truncateMaxRange,
        min.limit = minLimit,
        pData = phenoData)
    
    pDataVar <- "Group"
    pDataPathMapping <-
        list(G1 = "compMatForG1.csv",
             G2 = "compMatForG2.csv")
    
    expect_error(
        compensateFromMatrix(fs[[1]],
                             matrixSource = "pData",
                             pDataVar = pDataVar,
                             pDataPathMapping = pDataPathMapping),
        regexp = "compMatForG1.csv")
    
    expect_error(
        compensateFromMatrix(fs[[2]],
                             matrixSource = "pData",
                             pDataVar = pDataVar,
                             pDataPathMapping = pDataPathMapping),
        regexp = "compMatForG2.csv")
    
    expect_error(
        compensateFromMatrix(fs,
                             matrixSource = "pData",
                             pDataVar = pDataVar,
                             pDataPathMapping = pDataPathMapping),
        regexp = "compMatForG1.csv")
    
    pDataPathMapping <- 
        list(W1 = "compMatForG1.csv",
             G2 = "compMatForG2.csv")
    
    expect_error(
        compensateFromMatrix(fs[[1]],
                             matrixSource = "pData",
                             pDataVar = pDataVar,
                             pDataPathMapping = pDataPathMapping),
        regexp = "No mapping found for variable")
    
    pDataVar <- "NonExistentVar"
    expect_error(
        compensateFromMatrix(fs[[1]],
                             matrixSource = "pData",
                             pDataVar = pDataVar,
                             pDataPathMapping = pDataPathMapping),
        regexp = "not found in flowFrame")
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


test_that("removeDebrisManualGate works", {
    ref_ff_s <- readRDS(test_path("fixtures", "ff_s3.rds"))
    
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



test_that("removeDeadCellsManualGate works", {
    ref_ff_cells <- readRDS(test_path("fixtures", "ff_cells_manual_gate.rds"))
    
    remDeadCellsGateData <- c(0, 0, 250000, 250000,
                              0, 650, 650, 0)
    
    ff_lcells <-
        removeDeadCellsManualGate(ref_ff_cells,
                                  FSCChannel = "FSC-A",
                                  LDMarker = "L/D Aqua - Viability",
                                  gateData = remDeadCellsGateData)
    
    ref_ff_lcells <- readRDS(test_path("fixtures", "ff_lcells_manual_gate.rds"))
    
    #saveRDS(ff_lcells, test_path("fixtures", "ff_lcells_manual_gate.rds"))
    
    expect_equal(
        flowCore::exprs(ff_lcells),
        flowCore::exprs(ref_ff_lcells)
    )
    
    # same with channel name instead of marker for L/D
    ff_lcells2 <-
        removeDeadCellsManualGate(ref_ff_cells,
                                  FSCChannel = "FSC-A",
                                  LDMarker = "Comp-525/50Violet-A",
                                  gateData = remDeadCellsGateData)
    
    expect_equal(
        flowCore::exprs(ff_lcells2),
        flowCore::exprs(ref_ff_lcells)
    )
    
})

test_that("qualityControlPeacoQC", {
    ref_ff_lcells <- readRDS(test_path("fixtures", "ff_lcells_manual_gate.rds"))

    transListPath <- file.path(system.file("extdata", 
                                           package = "CytoPipeline"),
                               "OMIP021_TransList.rds")
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


test_that("readRDSObject works", {
    expect_error(readRDSObject("dummyPath.rds"),
                 regexp = "file dummyPath.rds does not exist")
    transListPath <- file.path(system.file("extdata", 
                                           package = "CytoPipeline"),
                               "OMIP021_TransList.rds")  
                                         
    obj <- readRDSObject(transListPath)
    
    refTransList <- readRDS(transListPath)
    
    expect_identical(obj, refTransList)
})

test_that("applyScaleTransform works", {
    transListPath <- file.path(system.file("extdata", 
                                           package = "CytoPipeline"),
                               "OMIP021_TransList.rds") 
                                        
    ff_c <- readRDS(test_path("fixtures", "ff_c.rds"))
    transList <- readRDS(transListPath)
    
    refFF <- flowCore::transform(ff_c, transList)
    thisFF <- applyScaleTransforms(ff_c, transList)
    
    expect_equal(
        flowCore::exprs(thisFF),
        flowCore::exprs(refFF)
    )
})








