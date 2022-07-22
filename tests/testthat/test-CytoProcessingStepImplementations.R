# # flags controlling execution flow of this file
# runPipelineStepImplementations <- FALSE
#
# # main parameters : sample files and output files
#
# rawDataDir <- paste0(system.file("extdata", package = "CytoPipeline"), "/")
# rdsDir <- paste0(system.file(package = "CytoPipeline", "/tests/testthat/rds"), "/")
# sampleFiles <- paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#
# # reference compensated fcs file
# comp_matrix <- flowCore::spillover(OMIP021Samples[[1]])$SPILL
# ff_c <- FF_compensate(OMIP021Samples[[1]], spillover = comp_matrix)
#
# # reference scale transformation list
# refTransListPath <- paste0(rdsDir, "OMIP021_TransList.rds")
# refTransList <- readRDS(file = refTransListPath)
#
# # parameters for compensation
#
# #compensationMethod <- "import" # read matrix from an external file
# compensationMethod <- "fcs" #uses matrix stored in respective fcs files
# compensationParams <- list()
#
# # parameters for margin removal
#
# marginMethod <- "PeacoQC"
# marginParams <- list()
#
# transformMethod <- "compute"
# transformParams <- list()
# transformParams$aggNCells <- 100000
# transformParams$aggSeed <- 1
# transformParams$fcsReadParams <- list()
# transformParams$fcsReadParams$truncate_max_range <- FALSE
# transformParams$fcsReadParams$min.limit <- NULL
# transformParams$fluoMethod <- "estimateLogicle"
# transformParams$scatterMethod <- "linear"
# referenceMarker <- "BV785 - CD3" # Scatter values will be scaled to have the same range
# transformParams$scatterRefMarker <- referenceMarker
# transformParams$manualUpdate <- FALSE
# transformParams$saveTransfo <- FALSE
# #transformParams$rdsOutputDir <- rdsDir
#
# # parameters for import of scale transformations
# importTransformParams <- list()
# importTransformParams$manualUpdate <- FALSE
# importTransformParams$saveTransfo <- FALSE
# importTransformParams$rdsInputPath <- paste0(rdsDir, "OMIP021_TransList.rds")
#
#
# if (runPipelineStepImplementations)
#   test_that("pipeline transformation computation step works", {
#
#     transList <- suppressWarnings(
#       computeChannelsTransformations(files = sampleFiles,
#                                      transformMethod = "compute",
#                                      transformParams = transformParams,
#                                      marginMethod = marginMethod,
#                                      marginParams = marginParams,
#                                      compensationMethod = compensationMethod,
#                                      compensationParams = compensationParams
#       ))
#
#     #saveRDS(transList, refTransListPath)
#
#     refFF <- flowCore::transform(ff_c, refTransList)
#     thisFF <- flowCore::transform(ff_c, transList)
#
#     expect_equal(flowCore::exprs(thisFF),
#                  flowCore::exprs(refFF))
#
#
#   })
#
# if (runPipelineStepImplementations) test_that("pipeline transformation import step works", {
#   transList <-
#     computeChannelsTransformations(files = NULL,
#                                    transformMethod = "import",
#                                    transformParams = importTransformParams,
#                                    marginMethod = NULL,
#                                    marginParams = NULL,
#                                    compensationMethod = NULL,
#                                    compensationParams = NULL
#     )
#
#   refFF <- flowCore::transform(ff_c, refTransList)
#   thisFF <- flowCore::transform(ff_c, transList)
#
#   expect_equal(thisFF, refFF)
#
# })
#
# # parameters for fcs file reading
# fcsReadParams <- list()
# fcsReadParams$truncate_max_range <- FALSE
# fcsReadParams$min.limit <- NULL
#
# # parameters for doublets elimination
# singletsGatingMethod <- "flowStats::singletGate"
# singletsGatingParams <- list()
# singletsGatingParams$area <- c("FSC-A", "SSC-A")
# singletsGatingParams$height <- c("FSC-H", "SSC-H")
# singletsGatingParams$widerGate = TRUE
#
# # parameters for debris elimination
# cellsGatingMethod <- "flowClust::tmixFilter" #automatic
# cellsGatingParams <- list()
# cellsGatingParams$FSC <- "FSC-A"
# cellsGatingParams$SSC <- "SSC-A"
# cellsGatingParams$K <- 3
# cellsGatingParams$level <- 0.97
# cellsGatingParams$B <- 100
#
# # parameters for dead cells elimination
# liveCellsGatingMethod <- "openCyto::gate_tail"
# liveCellsGatingPreTransform <- TRUE
# liveCellsGatingParams <- list()
# liveCellsGatingParams$FSC <- "FSC-A"
# liveCellsGatingParams$LDMarker <- "L/D Aqua - Viability"
# liveCellsGatingParams$num_peaks <- 2
# liveCellsGatingParams$ref_peak <- 2
# liveCellsGatingParams$strict <- FALSE
# liveCellsGatingParams$positive <- FALSE
#
# # parameters for quality control
# qualityControlMethod <- "PeacoQC"
# qualityControlPreTransform <- TRUE
# qualityControlParams <- list()
# qualityControlParams$min_cells <- 150 #default
# qualityControlParams$max_bins <- 500 #default
# qualityControlParams$step <- 500 #default
# qualityControlParams$MAD <- 6 #default
# qualityControlParams$IT_limit <- 0.55 #default
# qualityControlParams$force_IT <- 150 #default
# qualityControlParams$peak_removal <- (1/3) #default
# qualityControlParams$min_nr_bins_peakdetection <- 10 #default
#
# fs_raw <-
#   flowCore::read.flowSet(sampleFiles,
#                          truncate_max_range = fcsReadParams$truncate_max_range,
#                          min.limit = fcsReadParams$min.limit)
#
# ref_ff_m <- readRDS(paste0(rdsDir, "ff_m.rds"))
#
# if (runPipelineStepImplementations) test_that("Margin Removal step works", {
#
#   ff_m <-
#     suppressWarnings(runMarginsRemoval(ff = fs_raw[[1]],
#                                        marginMethod = marginMethod,
#                                        marginParams = marginParams))
#
#   #saveRDS(ff_m, paste0(rdsDir, "ff_m.rds"))
#   expect_equal(ff_m, ref_ff_m)
#
# })
#
# ref_ff_c <- readRDS(paste0(rdsDir, "ff_c.rds"))
# if (runPipelineStepImplementations) test_that("Compensation step works", {
#   ff_c <-
#     runCompensation(ref_ff_m,
#                     compensationMethod = compensationMethod,
#                     compensationParams = compensationParams)
#   #saveRDS(ff_c, paste0(rdsDir, "ff_c.rds"))
#   expect_equal(ff_c, ref_ff_c)
#
# })
#
# ref_ff_s <- readRDS(paste0(rdsDir, "ff_s.rds"))
# if (runPipelineStepImplementations) test_that("Doublets removal step works", {
#   ff_s <-
#     runDoubletsRemoval(ref_ff_c,
#                        singletsGatingMethod = singletsGatingMethod,
#                        singletsGatingParams = singletsGatingParams)
#   #saveRDS(ff_s, paste0(rdsDir, "ff_s.rds"))
#   expect_equal(ff_s, ref_ff_s)
# })
#
# ref_ff_cells <- readRDS(paste0(rdsDir, "ff_cells.rds"))
# if (runPipelineStepImplementations) test_that("Debris removal step works", {
#   ff_cells <-
#     runDebrisRemoval(ref_ff_s,
#                      cellsGatingMethod = cellsGatingMethod,
#                      cellsGatingParams = cellsGatingParams)
#   #saveRDS(ff_cells, paste0(rdsDir, "ff_cells.rds"))
#   expect_equal(ff_cells, ref_ff_cells)
# })
#
# ref_ff_lcells <- readRDS(paste0(rdsDir, "ff_lcells.rds"))
# if (runPipelineStepImplementations) test_that("Dead cells removal step works", {
#   ff_lcells <-
#     runDeadCellsRemoval(ref_ff_cells,
#                         liveCellsGatingMethod = liveCellsGatingMethod,
#                         liveCellsGatingParams = liveCellsGatingParams)
#   #saveRDS(ff_lcells, paste0(rdsDir, "ff_lcells.rds"))
#   expect_equal(ff_lcells, ref_ff_lcells)
# })
#
# ref_ff_qualityControl <- readRDS(paste0(rdsDir, "ff_QC_PeacoQC.rds"))
# if (runPipelineStepImplementations) test_that("Quality control step works with PeacoQC", {
#   suppressWarnings(ff_QualityControl <-
#                      runQualityControl(ref_ff_lcells,
#                                        qualityControlMethod = qualityControlMethod,
#                                        qualityControlPreTransform = qualityControlPreTransform,
#                                        transList = refTransList,
#                                        qualityControlParams = qualityControlParams))
#   #saveRDS(ff_QualityControl, paste0(rdsDir, "ff_QC_PeacoQC.rds"))
#   expect_equal(ff_QualityControl, ref_ff_qualityControl)
# })
#
# qualityControlParamsFlowAI <-
#   list(remove_from = "all",
#        second_fractionFR = 0.1,
#        deviationFR = "MAD",
#        alphaFR = 0.01,
#        decompFR = TRUE,
#        outlier_binsFS = FALSE,
#        pen_valueFS = 500,
#        max_cptFS = 3,
#        sideFM = "both",
#        neg_valuesFM = 1)
#
# ref_ff_qualityControl_flowAI <- readRDS(paste0(rdsDir, "ff_QC_flowAI.rds"))
#
# if (runPipelineStepImplementations) test_that("Quality control step works with flowAI", {
#   ff_QualityControl <- suppressWarnings(
#     runQualityControl(fs_raw[[1]],
#                       qualityControlMethod = "flowAI",
#                       qualityControlPreTransform = FALSE,
#                       qualityControlParams = qualityControlParamsFlowAI))
#
#   #saveRDS(ff_QualityControl, paste0(rdsDir, "ff_QC_flowAI.rds"))
#   expect_equal(flowCore::exprs(ff_QualityControl),
#                flowCore::exprs(ref_ff_qualityControl_flowAI))
# })
#
# qualityControlParamsFlowClean <-
#   list(binSize = 0.01, # default
#        nCellCutoff = 500, # default
#        cutoff = "median", # default
#        fcMax = 1.3, # default
#        nstable = 5)
#
# ref_ff_qualityControl_flowClean <- readRDS(paste0(rdsDir, "ff_QC_flowClean.rds"))
#
# if (runPipelineStepImplementations) test_that("Quality control step works with flowClean", {
#   ff_QualityControl <- #suppressWarnings(
#     runQualityControl(fs_raw[[1]],
#                       qualityControlMethod = "flowClean",
#                       qualityControlPreTransform = FALSE,
#                       qualityControlParams = qualityControlParamsFlowClean)#)
#
#   #saveRDS(ff_QualityControl, paste0(rdsDir, "ff_QC_flowClean.rds"))
#   expect_equal(flowCore::exprs(ff_QualityControl),
#                flowCore::exprs(ref_ff_qualityControl_flowClean))
# })
#
# ref_ff_qualityControl_flowCut <- readRDS(paste0(rdsDir, "ff_QC_flowCut.rds"))
#
# qualityControlParamsFlowCut <-
#   list(MaxContin = 0.1,
#        MeanOfMeans = 0.13,
#        MaxOfMeans = 0.15,
#        MaxValleyHgt = 0.1,
#        MaxPercCut = 0.3,
#        LowDensityRemoval = 0.1,
#        RemoveMultiSD = 7,
#        AlwaysClean = FALSE,
#        IgnoreMonotonic = FALSE,
#        MonotonicFix = NULL,
#        Measures = c(1:8))
#
# if (runPipelineStepImplementations) test_that("Quality control step works with flowCut", {
#   ff_QualityControl <- suppressWarnings(
#     runQualityControl(fs_raw[[1]],
#                       qualityControlMethod = "flowCut",
#                       qualityControlPreTransform = FALSE,
#                       qualityControlParams = qualityControlParamsFlowCut))
#
#   #saveRDS(ff_QualityControl, paste0(rdsDir, "ff_QC_flowCut.rds"))
#   expect_equal(flowCore::exprs(ff_QualityControl),
#                flowCore::exprs(ref_ff_qualityControl_flowCut))
# })
#
# ref_ff_t <- readRDS(paste0(rdsDir, "ff_t.rds"))
# if (runPipelineStepImplementations) test_that("Scale transformation step works", {
#   ff_t <- flowCore::transform(ref_ff_qualityControl,
#                               translist = refTransList)
#   saveRDS(ff_t, paste0(rdsDir, "ff_t.rds"))
#   expect_equal(ff_t, ref_ff_t)
# })
