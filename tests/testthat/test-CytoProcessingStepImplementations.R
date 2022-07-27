# CytoPipeline - Copyright (C) <2022> <Université catholique de Louvain (UCLouvain), Belgique>
#   
#   Description and complete License: see LICENSE file.
# 
# This program (CytoPipeline) is free software: 
#   you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (<http://www.gnu.org/licenses/>).

# main parameters : sample files and output files

rawDataDir <- paste0(system.file("extdata", package = "CytoPipeline"), "/")
rdsDir <- paste0(system.file("extdata", package = "CytoPipeline"), "/")
sampleFiles <- paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))

# reference compensated fcs file
comp_matrix <- flowCore::spillover(OMIP021Samples[[1]])$SPILL
ff_c <- runCompensation(OMIP021Samples[[1]], spillover = comp_matrix)

# reference scale transformation list
refTransListPath <- paste0(rdsDir, "OMIP021_TransList.rds")
refTransList <- readRDS(file = refTransListPath)

test_that("estimateScaleTransforms work", {
  
  transList <-
    suppressMessages(estimateScaleTransforms(ff = ff_c,
                            fluoMethod = "estimateLogicle",
                            scatterMethod = "linear",
                            scatterRefMarker = "BV785 - CD3"))

  #saveRDS(transList, refTransListPath)

  refFF <- flowCore::transform(ff_c, refTransList)
  thisFF <- flowCore::transform(ff_c, transList)

  expect_equal(flowCore::exprs(thisFF),
               flowCore::exprs(refFF))

})

truncateMaxRange <- FALSE
minLimit <- NULL

fs_raw <-
  flowCore::read.flowSet(sampleFiles,
                         truncate_max_range = truncateMaxRange,
                         min.limit = minLimit)
fs_raw <- flowCore::fsApply(fs_raw, FUN = appendCellID)

test_that("readSampleFiles works", {
  res <- readSampleFiles(sampleFiles = sampleFiles,
                         whichSamples = "all",
                         truncate_max_range = truncateMaxRange,
                         min.limit = minLimit)
  expect_equal(flowCore::exprs(res[[1]]), 
               flowCore::exprs(fs_raw[[1]]))
  expect_equal(flowCore::exprs(res[[2]]), 
               flowCore::exprs(fs_raw[[2]]))
  
  res2 <- readSampleFiles(sampleFiles = sampleFiles,
                          whichSamples = 2,
                          truncate_max_range = truncateMaxRange,
                          min.limit = minLimit)
                                                    
  expect_equal(flowCore::exprs(res2), 
               flowCore::exprs(fs_raw[[2]]))
})

ref_ff_m <- readRDS(paste0(rdsDir, "ff_m.rds"))

test_that("removeMarginsPeacoQC works", {
  
  ff_m <-
    suppressWarnings(removeMarginsPeacoQC(x = fs_raw[[1]]))
  
  #saveRDS(ff_m, paste0(rdsDir, "ff_m.rds"))
  expect_equal(flowCore::exprs(ff_m), 
               flowCore::exprs(ref_ff_m))
})

ref_ff_c <- readRDS(paste0(rdsDir, "ff_c.rds"))

test_that("compensateFromMatrix works", {
  ff_c <-
    compensateFromMatrix(ref_ff_m,
                         matrixSource = "fcs")
                    
  #saveRDS(ff_c, paste0(rdsDir, "ff_c.rds"))
  expect_equal(flowCore::exprs(ff_c), 
               flowCore::exprs(ref_ff_c))
  
})

ref_ff_s <- readRDS(paste0(rdsDir, "ff_s.rds"))
test_that("removeDoubletsFlowStats works", {
  ff_s <-
    removeDoubletsFlowStats(ref_ff_c,
                            areaChannels = c("FSC-A", "SSC-A"),
                            heightChannels = c("FSC-H", "SSC-H"),
                            widerGate = TRUE)
                       
  #saveRDS(ff_s, paste0(rdsDir, "ff_s.rds"))
  expect_equal(flowCore::exprs(ff_s), 
               flowCore::exprs(ref_ff_s))
})

ref_ff_s2 <- readRDS(paste0(rdsDir, "ff_s2.rds"))
test_that("removeDoubletsPeacoQC works", {
  ff_s2 <-
    removeDoubletsFlowStats(ref_ff_c,
                            areaChannels = c("FSC-A", "SSC-A"),
                            heightChannels = c("FSC-H", "SSC-H"),
                            nmaps = c(3,5))
  
  #saveRDS(ff_s2, paste0(rdsDir, "ff_s2.rds"))
  expect_equal(flowCore::exprs(ff_s2), 
               flowCore::exprs(ref_ff_s2))
})

ref_ff_cells <- readRDS(paste0(rdsDir, "ff_cells.rds"))

test_that("removeDebrisFlowClustTmix works", {
  ff_cells <-
    removeDebrisFlowClustTmix(ref_ff_s,
                              FSCChannel = "FSC-A",
                              SSCChannel = "SSC-A",
                              nClust = 3,
                              level = 0.97,
                              B = 100)
                     
  #saveRDS(ff_cells, paste0(rdsDir, "ff_cells.rds"))
  expect_equal(flowCore::exprs(ff_cells), 
               flowCore::exprs(ref_ff_cells))
})

ref_ff_lcells <- readRDS(paste0(rdsDir, "ff_lcells.rds"))

test_that("removeDeadCellsGateTail works", {
  ff_lcells <-
    removeDeadCellsGateTail(ref_ff_cells,
                            preTransform = TRUE,
                            transList = refTransList,
                            LDMarker = "L/D Aqua - Viability",
                            num_peaks = 2,
                            ref_peak = 2,
                            strict = FALSE,
                            positive = FALSE)
                        
  #saveRDS(ff_lcells, paste0(rdsDir, "ff_lcells.rds"))
  expect_equal(flowCore::exprs(ff_lcells), 
               flowCore::exprs(ref_ff_lcells))
})


ref_ff_qualityControl <- readRDS(paste0(rdsDir, "ff_QC_PeacoQC.rds"))
test_that("qualityControlPeacoQC", {
  suppressWarnings(ff_QualityControl <-
                     qualityControlPeacoQC(ref_ff_lcells,
                                           preTransform = TRUE,
                                           transList = refTransList,
                                           min_cells = 150, #default
                                           max_bins = 500, #default
                                           MAD = 6, #default
                                           IT_limit = 0.55, #default
                                           force_IT = 150, #default
                                           peak_removal = (1/3), #default
                                           min_nr_bins_peakdetection = 10 #default
                      ))
                                       
  #saveRDS(ff_QualityControl, paste0(rdsDir, "ff_QC_PeacoQC.rds"))
  expect_equal(flowCore::exprs(ff_QualityControl), 
               flowCore::exprs(ref_ff_qualityControl))
})



ref_ff_qualityControl_flowAI <- readRDS(paste0(rdsDir, "ff_QC_flowAI.rds"))

test_that("qualityControlFlowAI works", {
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
                         neg_valuesFM = 1))

  #saveRDS(ff_QualityControl, paste0(rdsDir, "ff_QC_flowAI.rds"))
  expect_equal(flowCore::exprs(ff_QualityControl),
               flowCore::exprs(ref_ff_qualityControl_flowAI))
})

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
# test_that("Scale transformation step works", {
#   ff_t <- flowCore::transform(ref_ff_qualityControl,
#                               translist = refTransList)
#   saveRDS(ff_t, paste0(rdsDir, "ff_t.rds"))
#   expect_equal(ff_t, ref_ff_t)
# })
