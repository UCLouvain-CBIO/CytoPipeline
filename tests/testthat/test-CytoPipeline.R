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


outputDir <- withr::local_tempdir()
if (!interactive()) pdf(NULL)

test_that("CytoPipeline default creation raises no error", {
    expect_error(pipL0 <- CytoPipeline(), NA)
})

test_that("Cytopipeline add/remove/clean processing step works", {
    rawDataDir <- paste0(system.file("extdata", package = "CytoPipeline"), "/")
    experimentName <- "OMIP021_PeacoQC"
    sampleFiles <- paste0(rawDataDir, list.files(rawDataDir,
        pattern = "sample_"
    ))
    transListPath <- paste0(system.file("extdata", 
                                        package = "CytoPipeline"),
                            "/OMIP021_TransList.rds")

    # main parameters : sample files and experiment name
    pipelineParams <- list()
    pipelineParams$experimentName <- experimentName
    
    pipL <- CytoPipeline(pipelineParams)
    expect_error(show(pipL), NA)
    
    pipelineParams$sampleFiles <- sampleFiles

    pipL <- CytoPipeline(pipelineParams)
    expect_error(show(pipL), NA)

    pipL <- addProcessingStep(pipL,
        whichQueue = "scale transform",
        CytoProcessingStep(
            name = "scale_transform_read",
            FUN = "readRDS",
            ARGS = list(file = transListPath)
        )
    )

    expect_equal(getNbProcessingSteps(pipL, "scale transform"), 1)

    pipL <- addProcessingStep(pipL,
        whichQueue = "scale transform",
        CytoProcessingStep(
            name = "scale_transform_sum",
            FUN = "sum",
            ARGS = list()
        )
    )

    expect_equal(getNbProcessingSteps(pipL, "scale transform"), 2)

    pipL <- removeProcessingStep(pipL,
        whichQueue = "scale transform",
        index = 2
    )
    expect_equal(getNbProcessingSteps(pipL, "scale transform"), 1)
    pS <- getProcessingStep(pipL, whichQueue = "scale transform", index = 1)
    expect_equal(pS@FUN, "readRDS")

    pipL <- addProcessingStep(pipL,
        whichQueue = "pre-processing",
        CytoProcessingStep(
            name = "pre-processing_sum",
            FUN = "sum",
            ARGS = list()
        )
    )
    expect_equal(getNbProcessingSteps(pipL, "scale transform"), 1)
    expect_equal(getNbProcessingSteps(pipL, "pre-processing"), 1)

    expect_error(pipL <- addProcessingStep(pipL,
        whichQueue = "pre-processing",
        CytoProcessingStep(
            name = "pre-processing_sum",
            FUN = "mean",
            ARGS = list()
        )
    ), regexp = "There already exist a step")

    pipL <- cleanProcessingSteps(pipL)
    expect_equal(getNbProcessingSteps(pipL, "scale transform"), 0)
    expect_equal(getNbProcessingSteps(pipL, "pre-processing"), 0)
    
    newExp <- "newExperiment"
    experimentName(pipL) <- newExp
    expect_equal(experimentName(pipL), newExp)
})

test_that("CytoPipeline with reading scale transfo only raises no error", {
    expect_error(
        {
            rawDataDir <-
                paste0(system.file("extdata", package = "CytoPipeline"), "/")
            experimentName <- "OMIP021_PeacoQC"
            sampleFiles <- paste0(rawDataDir, list.files(rawDataDir,
                pattern = "sample_"
            ))
            transListPath <- paste0(system.file("extdata", 
                                                package = "CytoPipeline"),
                                    "/OMIP021_TransList.rds")

            # main parameters : sample files and output files
            pipelineParams <- list()
            pipelineParams$experimentName <- experimentName
            pipelineParams$sampleFiles <- sampleFiles

            pipL <- CytoPipeline(pipelineParams)

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "scale transform",
                    CytoProcessingStep(
                        name = "scale_transform_read",
                        FUN = "readRDS",
                        ARGS = list(file = transListPath)
                    )
                )

            suppressWarnings(execute(pipL,
                rmCache = TRUE,
                path = outputDir
            ))
        },
        NA
    )
})


test_that("CytoPipeline with complex flows raises no error", {
    expect_error(
        {
            rawDataDir <-
                paste0(system.file("extdata", package = "CytoPipeline"), "/")
            experimentName <- "OMIP021_PeacoQC"
            sampleFiles <- paste0(rawDataDir, list.files(rawDataDir,
                pattern = "sample_"
            ))

            # main parameters : sample files and output files
            pipelineParams <- list()
            pipelineParams$experimentName <- experimentName
            pipelineParams$sampleFiles <- sampleFiles
            pipL <- CytoPipeline(pipelineParams)

            ### SCALE TRANSFORMATION STEPS ###

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "scale transform",
                    CytoProcessingStep(
                        name = "flowframe_read",
                        FUN = "readSampleFiles",
                        ARGS = list(
                            whichSamples = "all",
                            truncate_max_range = FALSE,
                            min.limit = NULL
                        )
                    )
                )

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "scale transform",
                    CytoProcessingStep(
                        name = "remove_margins",
                        FUN = "removeMarginsPeacoQC",
                        ARGS = list()
                    )
                )

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "scale transform",
                    CytoProcessingStep(
                        name = "compensate",
                        FUN = "compensateFromMatrix",
                        ARGS = list(matrixSource = "fcs")
                    )
                )

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "scale transform",
                    CytoProcessingStep(
                        name = "flowframe_aggregate",
                        FUN = "aggregateAndSample",
                        ARGS = list(
                            nTotalEvents = 10000,
                            seed = 0
                        )
                    )
                )

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "scale transform",
                    CytoProcessingStep(
                        name = "scale_transform_estimate",
                        FUN = "estimateScaleTransforms",
                        ARGS = list(
                            fluoMethod = "estimateLogicle",
                            scatterMethod = "linear",
                            scatterRefMarker = "BV785 - CD3"
                        )
                    )
                )

            ### PRE-PROCESSING STEPS ###

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "pre-processing",
                    CytoProcessingStep(
                        name = "flowframe_read",
                        FUN = "readSampleFiles",
                        ARGS = list(
                            truncate_max_range = FALSE,
                            min.limit = NULL
                        )
                    )
                )


            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "pre-processing",
                    CytoProcessingStep(
                        name = "remove_margins",
                        FUN = "removeMarginsPeacoQC",
                        ARGS = list()
                    )
                )

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "pre-processing",
                    CytoProcessingStep(
                        name = "compensate",
                        FUN = "compensateFromMatrix",
                        ARGS = list(matrixSource = "fcs")
                    )
                )

            pipL <-
                addProcessingStep(
                    pipL,
                    whichQueue = "pre-processing",
                    CytoProcessingStep(
                        name = "remove_doublets",
                        FUN = "removeDoubletsFlowStats",
                        ARGS = list(
                            areaChannels = c("FSC-A", "SSC-A"),
                            heightChannels = c("FSC-H", "SSC-H"),
                            wider_gate = TRUE
                        )
                    )
                )

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "pre-processing",
                    CytoProcessingStep(
                        name = "remove_debris",
                        FUN = "removeDebrisFlowClustTmix",
                        ARGS = list(
                            FSCChannel = c("FSC-A"),
                            SSCChannel = c("SSC-A"),
                            nClust = 3,
                            level = 0.97,
                            B = 100,
                            verbose = TRUE
                        )
                    )
                )

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "pre-processing",
                    CytoProcessingStep(
                        name = "remove_dead_cells",
                        FUN = "removeDeadCellsGateTail",
                        ARGS = list(
                            FSCChannel = c("FSC-A"),
                            LDMarker = "L/D Aqua - Viability",
                            num_peaks = 2,
                            ref_peak = 2,
                            strict = FALSE,
                            positive = FALSE
                        )
                    )
                )

            pipL <-
                addProcessingStep(
                    pipL,
                    whichQueue = "pre-processing",
                    CytoProcessingStep(
                        name = "perform_QC",
                        FUN = "qualityControlPeacoQC",
                        ARGS = list(
                            preTransform = TRUE,
                            min_cells = 150, # default
                            max_bins = 500, # default
                            step = 500, # default,
                            MAD = 6, # default
                            IT_limit = 0.55, # default
                            force_IT = 150, # default
                            peak_removal = 0.3333, # default
                            min_nr_bins_peakdetection = 10 # default
                        )
                    )
                )

            pipL <-
                addProcessingStep(pipL,
                    whichQueue = "pre-processing",
                    CytoProcessingStep(
                        name = "transform",
                        FUN = "applyScaleTransforms",
                        ARGS = list()
                    )
                )

            suppressWarnings(execute(pipL,
                rmCache = TRUE,
                path = outputDir
            ))
        },
        NA
    )
})

test_that("CytoPipeline with json input raises no error", {
    expect_error(
        {
            jsonDir <- system.file("extdata", package = "CytoPipeline")
            jsonPath <- paste0(jsonDir, "/pipelineParams.json")

            pipL2 <- CytoPipeline(jsonPath)
            withr::with_dir(new = jsonDir, {
                suppressWarnings(execute(pipL2, 
                                         rmCache = FALSE, 
                                         path = outputDir))})
        },
        NA
    )
})

test_that("CytoPipeline export as list works", {
    jsonDir <- system.file("extdata", package = "CytoPipeline")
    jsonPath <- paste0(jsonDir, "/pipelineParams.json")
    
    pipL1 <- CytoPipeline(jsonPath)
    pipList <- as.list(pipL1)
    
    pipL2 <- CytoPipeline(pipList)
    expect_identical(pipL1, pipL2)
    
})

test_that("CytoPipeline rebuilt from cache raises no error", {
    expect_error(
        {
            experimentName <- "OMIP021_PeacoQC"
            pipL3 <- buildCytoPipelineFromCache(
                experimentName = experimentName,
                path = outputDir
            )
            suppressWarnings(execute(pipL3,
                rmCache = FALSE,
                path = outputDir
            ))
        },
        NA
    )
})


test_that("CytoPipeline not in cache with warning", {
    expect_warning(
        pipL4 <- buildCytoPipelineFromCache(
            experimentName = "non_existent",
            path = outputDir
        ),
        regexp = "no cache directory found"
    )
})


test_that("Check consistency with cache works", {
    rawDataDir <- paste0(system.file("extdata", package = "CytoPipeline"), "/")
    sampleFiles <- paste0(rawDataDir, list.files(rawDataDir,
        pattern = "sample_"
    ))

    pipL5 <- CytoPipeline(experimentName = "DummyExperiment")

    sampleFiles(pipL5) <- sampleFiles
    deleteCytoPipelineCache(pipL5, path = outputDir)

    pipL5 <- addProcessingStep(
        pipL5,
        "scale transform",
        CytoProcessingStep(
            name = "flowframe_read",
            FUN = "readSampleFiles",
            ARGS = list(
                whichSamples = "all",
                truncate_max_range = FALSE,
                min.limit = NULL
            )
        )
    )
    res <- checkCytoPipelineConsistencyWithCache(pipL5, path = outputDir)

    expect_error(suppressWarnings(execute(pipL5,
        rmCache = TRUE,
        path = outputDir
    )), NA)

    res <- checkCytoPipelineConsistencyWithCache(pipL5,
        path = outputDir
    )
    expect_equal(res$isConsistent, TRUE)
    expect_equal(unname(res$scaleTransformStepStatus[1]), "run")

    pipL5 <-
        addProcessingStep(pipL5,
            whichQueue = "pre-processing",
            CytoProcessingStep(
                name = "flowframe_read",
                FUN = "readSampleFiles",
                ARGS = list(
                    truncate_max_range = FALSE,
                    min.limit = NULL
                )
            )
        )
    res <- checkCytoPipelineConsistencyWithCache(pipL5,
        path = outputDir
    )
    expect_equal(res$isConsistent, TRUE)
    expect_equal(unname(res$scaleTransformStepStatus[1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[1, 1]), "not_run")
    expect_equal(unname(res$preProcessingStepStatus[1, 2]), "not_run")
    expect_equal(res$scaleTransformStepOutputObjNames, c("flowframe_read_obj"))
    expect_equal(res$scaleTransformStepOutputClasses, c("flowSet"))
    expect_equal(res$preProcessingStepOutputObjNames, c("unknown"))
    expect_equal(res$preProcessingStepOutputClasses, c("unknown"))

    expect_error(execute(pipL5,
        rmCache = FALSE,
        path = outputDir
    ), NA)
    res <- checkCytoPipelineConsistencyWithCache(pipL5,
        path = outputDir
    )
    expect_equal(res$isConsistent, TRUE)
    expect_equal(unname(res$scaleTransformStepStatus[1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[1, 1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[1, 2]), "run")
    expect_equal(res$scaleTransformStepOutputObjNames, c("flowframe_read_obj"))
    expect_equal(res$scaleTransformStepOutputClasses, c("flowSet"))
    expect_equal(res$preProcessingStepOutputObjNames, c("flowframe_read_obj"))
    expect_equal(res$preProcessingStepOutputClasses, c("flowFrame"))

    pipL5 <-
        addProcessingStep(pipL5,
            whichQueue = "pre-processing",
            CytoProcessingStep(
                name = "remove_margins",
                FUN = "removeMarginsPeacoQC",
                ARGS = list()
            )
        )

    res <- checkCytoPipelineConsistencyWithCache(pipL5,
        path = outputDir
    )
    expect_equal(res$isConsistent, TRUE)
    expect_equal(unname(res$scaleTransformStepStatus[1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[1, 1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[1, 2]), "run")
    expect_equal(unname(res$preProcessingStepStatus[2, 1]), "not_run")
    expect_equal(unname(res$preProcessingStepStatus[2, 2]), "not_run")
    expect_equal(res$scaleTransformStepOutputObjNames, c("flowframe_read_obj"))
    expect_equal(res$scaleTransformStepOutputClasses, c("flowSet"))
    expect_equal(
        res$preProcessingStepOutputObjNames,
        c("flowframe_read_obj", "unknown")
    )
    expect_equal(
        res$preProcessingStepOutputClasses,
        c("flowFrame", "unknown")
    )


    expect_error(suppressWarnings(execute(pipL5,
        rmCache = FALSE,
        path = outputDir
    )), NA)
    res <- checkCytoPipelineConsistencyWithCache(pipL5,
        path = outputDir
    )
    expect_equal(res$isConsistent, TRUE)
    expect_equal(unname(res$scaleTransformStepStatus[1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[1, 1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[1, 2]), "run")
    expect_equal(unname(res$preProcessingStepStatus[2, 1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[2, 2]), "run")
    expect_equal(res$scaleTransformStepOutputObjNames, c("flowframe_read_obj"))
    expect_equal(res$scaleTransformStepOutputClasses, c("flowSet"))
    expect_equal(
        res$preProcessingStepOutputObjNames,
        c("flowframe_read_obj", "remove_margins_obj")
    )
    expect_equal(
        res$preProcessingStepOutputClasses,
        c("flowFrame", "flowFrame")
    )

    pipL5_bad <- pipL5

    pipL5_bad@flowFramesPreProcessingQueue[[2]]@name <- "aaaaaa"
    res <- checkCytoPipelineConsistencyWithCache(pipL5_bad,
        path = outputDir
    )
    expect_equal(res$isConsistent, FALSE)
    expect_equal(
        res$inconsistencyMsg,
        paste0(
            "inconsistent pre-processing step #2 for sample file ",
            "sample_Donor1.fcs (different in cache)"
        )
    )
    expect_equal(unname(res$scaleTransformStepStatus[1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[1, 1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[2, 1]), "inconsistent")

    expect_error(suppressWarnings(execute(pipL5_bad,
        rmCache = FALSE,
        path = outputDir
    )),
    regexp = "inconsistent pre-processing step"
    )
    
    pipL5_bad2 <- pipL5 
    pipL5_bad2@flowFramesPreProcessingQueue[[1]]@ARGS$truncate_max_range <- 
        TRUE
    res <- checkCytoPipelineConsistencyWithCache(pipL5_bad2,
                                                 path = outputDir
    )
    expect_equal(res$isConsistent, FALSE)
    expect_equal(
        res$inconsistencyMsg,
        paste0(
            "inconsistent pre-processing step #1 for sample file ",
            "sample_Donor1.fcs (different in cache)"
        )
    )
    expect_equal(unname(res$scaleTransformStepStatus[1]), "run")
    expect_equal(unname(res$preProcessingStepStatus[1, 1]), "inconsistent")
    expect_equal(unname(res$preProcessingStepStatus[2, 1]), "not_run")
    
    
    pipL5 <- removeProcessingStep(pipL5,
        whichQueue = "pre-processing",
        index = 2
    )
    
    res <- checkCytoPipelineConsistencyWithCache(pipL5,
                                                 path = outputDir
    )

    expect_equal(res$isConsistent, FALSE)
    expect_equal(
        res$inconsistencyMsg,
        paste0(
            "more pre-processing steps in cache than in CytoPipeline object"
        )
    )
})


test_that("plotCytoPipelineProcessingQueue works", {
    pipL6 <- CytoPipeline(experimentName = "DummyExperiment")

    deleteCytoPipelineCache(pipL6, path = outputDir)

    rawDataDir <- paste0(system.file("extdata", package = "CytoPipeline"), "/")
    sampleFiles <- paste0(rawDataDir, list.files(rawDataDir,
        pattern = "sample_"
    ))

    # put only second sample file for the time being
    sampleFiles(pipL6) <- sampleFiles[2]
    pipL6 <- addProcessingStep(
        pipL6,
        "scale transform",
        CytoProcessingStep(
            name = "flowframe_read",
            FUN = "readSampleFiles",
            ARGS = list(
                whichSamples = "all",
                truncate_max_range = FALSE,
                min.limit = NULL
            )
        )
    )

    expect_error(plotCytoPipelineProcessingQueue(
        pipL6,
        whichQueue = "scale transform",
        path = outputDir
    ), NA)

    expect_error(suppressWarnings(execute(pipL6,
        rmCache = TRUE,
        path = outputDir
    )), NA)

    expect_error(plotCytoPipelineProcessingQueue(
        pipL6,
        whichQueue = "scale transform",
        path = outputDir
    ), NA)
    expect_message(plotCytoPipelineProcessingQueue(
        pipL6,
        whichQueue = "pre-processing",
        path = outputDir
    ),
    regexp = "no sample file passed"
    )

    pipL6 <-
        addProcessingStep(pipL6,
            whichQueue = "pre-processing",
            CytoProcessingStep(
                name = "flowframe_read",
                FUN = "readSampleFiles",
                ARGS = list(
                    truncate_max_range = FALSE,
                    min.limit = NULL
                )
            )
        )

    expect_error(plotCytoPipelineProcessingQueue(
        pipL6,
        whichQueue = "scale transform",
        path = outputDir
    ), NA)
    expect_message(plotCytoPipelineProcessingQueue(
        pipL6,
        whichQueue = "pre-processing",
        path = outputDir
    ),
    regexp = "no sample file passed"
    )

    execute(pipL6, rmCache = FALSE, path = outputDir)

    plotCytoPipelineProcessingQueue(pipL6,
        whichQueue = "scale transform",
        path = outputDir
    )
    expect_error(plotCytoPipelineProcessingQueue(pipL6,
        sampleFile = 1,
        whichQueue = "pre-processing",
        path = outputDir
    ), NA)
    expect_error(plotCytoPipelineProcessingQueue(pipL6,
        sampleFile = 2,
        whichQueue = "pre-processing",
        path = outputDir
    ),
    regexp = "out of bounds"
    )


    pipL6 <-
        addProcessingStep(pipL6,
            whichQueue = "pre-processing",
            CytoProcessingStep(
                name = "remove_margins",
                FUN = "removeMarginsPeacoQC",
                ARGS = list()
            )
        )
    plotCytoPipelineProcessingQueue(pipL6,
        whichQueue = "scale transform",
        path = outputDir
    )
    plotCytoPipelineProcessingQueue(pipL6,
        sampleFile = 1,
        whichQueue = "pre-processing",
        path = outputDir
    )

    suppressWarnings(execute(pipL6,
        rmCache = FALSE,
        path = outputDir
    ))

    plotCytoPipelineProcessingQueue(pipL6,
        whichQueue = "scale transform",
        path = outputDir
    )
    plotCytoPipelineProcessingQueue(pipL6,
        sampleFile = 1,
        whichQueue = "pre-processing",
        path = outputDir
    )

    # add first sample file to see the impact
    sampleFiles(pipL6) <- sampleFiles

    plotCytoPipelineProcessingQueue(pipL6,
        whichQueue = "scale transform",
        path = outputDir
    )
    # following should show yellow boxes
    plotCytoPipelineProcessingQueue(pipL6,
        whichQueue = "pre-processing",
        path = outputDir
    )
    # following should show bow in green
    plotCytoPipelineProcessingQueue(pipL6,
        sampleFile = 2,
        whichQueue = "pre-processing",
        path = outputDir
    )

    suppressWarnings(execute(pipL6,
        rmCache = FALSE,
        path = outputDir
    ))

    plotCytoPipelineProcessingQueue(pipL6,
        whichQueue = "scale transform",
        path = outputDir
    )
    # following should now show the green box
    plotCytoPipelineProcessingQueue(pipL6,
        whichQueue = "pre-processing",
        path = outputDir
    )
    # following as well
    plotCytoPipelineProcessingQueue(pipL6,
        sampleFile = sampleFiles[2],
        whichQueue = "pre-processing",
        path = outputDir
    )


    pipL6@flowFramesPreProcessingQueue[[2]]@name <- "aaaaaa"

    expect_warning(plotCytoPipelineProcessingQueue(
        pipL6,
        whichQueue = "scale transform",
        path = outputDir
    ),
    regexp = "CytoPipeline object not consistent with cache"
    )

    expect_warning(plotCytoPipelineProcessingQueue(
        pipL6,
        sampleFile = 1,
        whichQueue = "pre-processing",
        path = outputDir
    ),
    regexp = "CytoPipeline object not consistent with cache"
    )


    expect_error(execute(pipL6,
        rmCache = FALSE,
        path = outputDir
    ), regexp = "inconsistent")
})

test_that("getCytoPipelineObject works", {
    expect_error(
        {
            experimentName <- "OMIP021_PeacoQC"

            rawDataDir <-
                paste0(system.file("extdata", package = "CytoPipeline"), "/")
            sampleFiles <- paste0(rawDataDir, list.files(rawDataDir,
                pattern = "sample_"
            ))

            pipL7 <- buildCytoPipelineFromCache(
                experimentName = experimentName,
                path = outputDir
            )

            plotCytoPipelineProcessingQueue(pipL7,
                sampleFile = 1,
                whichQueue = "pre-processing",
                path = outputDir
            )

            getCytoPipelineObjectInfos(pipL7,
                whichQueue = "pre-processing",
                sampleFile = sampleFiles[1],
                path = outputDir
            )
            getCytoPipelineObjectInfos(pipL7,
                whichQueue = "scale transform",
                sampleFile = sampleFiles[1],
                path = outputDir
            )

            ffFrom <- getCytoPipelineFlowFrame(pipL7,
                whichQueue = "pre-processing",
                sampleFile = sampleFiles[1],
                objectName = "compensate_obj",
                path = outputDir
            )

            ffTo <- getCytoPipelineFlowFrame(pipL7,
                whichQueue = "pre-processing",
                sampleFile = sampleFiles[1],
                objectName = "remove_doublets_obj",
                path = outputDir
            )

            ggplotFilterEvents(ffFrom, ffTo,
                xChannel = "FSC-A", yChannel = "FSC-H"
            )

            plotCytoPipelineProcessingQueue(pipL7,
                whichQueue = "scale transform",
                path = outputDir
            )
        },
        NA
    )


    expect_error(
        getCytoPipelineScaleTransform(
            pipL7,
            whichQueue = "scale transform",
            objectName = "flowframe_aggregate_obj",
            path = outputDir
        ),
        regexp = "does not appear to be a transformList"
    )


    expect_error(
        transList <-
            getCytoPipelineScaleTransform(
                pipL7,
                whichQueue = "scale transform",
                objectName = "scale_transform_estimate_obj",
                path = outputDir
            ),
        NA
    )

    expect_error(
        ffAgg <-
            getCytoPipelineFlowFrame(pipL7,
                whichQueue = "scale transform",
                objectName = "flowframe_aggregate_obj",
                path = outputDir
            ),
        NA
    )
})

