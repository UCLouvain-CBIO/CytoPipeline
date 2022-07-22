##' @export
readSampleFiles <- function(sampleFiles,
                            whichSamples = "all", ...){

  if (whichSamples == "all") {
    # do nothing : sampleFiles should contain all the input sample files already
  } else if (is.integer(whichSamples)) {
    sampleFiles <- sampleFiles[whichSamples]
  } else stop("'whichSamples' should be either 'all', or a vector of indexes")

  if (length(sampleFiles) == 1) {
    res <- flowCore::read.FCS(sampleFiles, ...)
    #Add a column with Cell ID
    res <- FF_AppendCellID(res, 1:(flowCore::nrow(res)))
  } else {
    res  <- flowCore::read.flowSet(sampleFiles, ...)
    #Add a column with Cell ID
    res <- fsApply(x = res,
                   FUN = function(ff){
                     FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))
                   })
  }
}

##' @importFrom flowCore exprs
##'
##' @export
removeMarginsPeacoQC <- function(x, ...) {

  myFunc <- function(ff){
    channel4Margins <-
      flowCore::colnames(ff)[FF_areSignalCols(ff)]
    ffOut <- PeacoQC::RemoveMargins(ff, channels = channel4Margins)
    return(ffOut)
  }
  if (inherits(x, "flowFrame")) {
    return(myFunc(x))
  } else if (inherits(x, "flowSet")) {
    fsOut <- flowCore::fsApply(x, FUN = myFunc, simplify = TRUE)
    return(fsOut)
  } else stop("x should be a flowCore::flowFrame or a flowCore::flowSet")
}

##' @export
removeMarginsFlowAI <- function(x, ...) {
  if(!inherits(x, "flowFrame") && !inherits(x, "flowSet")) {
    stop("x should be a flowCore::flowFrame or a flowCore::flowSet")
  }
  channel2Exclude <-
    flowCore::colnames(x)[!FF_areSignalCols(x)]


  xOut <-
    flowAI::flow_auto_qc(x,
                         remove_from = "FM",
                         output = 1,
                         ChExcludeFM = channel2Exclude,
                         html_report = FALSE,
                         mini_report = FALSE,
                         fcs_QC = FALSE,
                         fcs_highQ = FALSE,
                         fcs_lowQ = FALSE,
                         folder_results = FALSE,
                         ...)
  return(xOut)

}

### FUNCTIONS for pre-processing / compensation ###

##' @export
getAcquiredCompensationMatrix <- function(ff){
  res <- flowCore::spillover(ff)
  if (!is.null(res$SPILL)) {
    compensationMatrix <- res$SPILL
  } else if (!is.null(res$spillover)) {
    compensationMatrix <- res$spillover
  } else if (!is.null(res$`$SPILLOVER`)) {
    compensationMatrix <- res$`$SPILLOVER`
  } else {
    fileId <- flowCore::identifier(ff)
    msg <- ("Issue retrieving compensation matrix for file ")
    msg <- paste0(msg, fileId, " : slot is NULL!")
    stop(msg)
  }
  return(compensationMatrix)
}

##' @export
compensateFromMatrix <- function(x,
                                 matrixSource = c("fcs", "import"),
                                 matrixPath = NULL,
                                 ...
){
  myFunc <- function(ff, matrixSource = c("fcs", "import"),
                   matrixPath = NULL) {
    matrixSource = match.arg(matrixSource)
    if (matrixSource == "fcs") {
      # obtains compensation matrix
      compensationMatrix <-
        getAcquiredCompensationMatrix(ff)
    } else {
      # import matrix from file (path)
      if(is.null(matrixPath)){
        stop("No path specified for compensation matrix!")
      }
      if(!file.exists(matrixPath)) stop("Compensation matrix file not found!")
      compensationMatrix <- read.csv(matrixPath,
                                     check.names = FALSE,
                                     row.names = 1)
    }

    ffOut <- FF_compensate(ff, compensationMatrix,
                                     updateChannelNames = TRUE)
    return(ffOut)
  }

  if (inherits(x, "flowFrame")) {
    return(myFunc(x, matrixSource = matrixSource, matrixPath = matrixPath))
  } else if (inherits(x, "flowSet")) {
    fsOut <- flowCore::fsApply(x, FUN = myFunc, simplify = TRUE)
    return(fsOut)
  } else stop("x should be a flowCore::flowFrame or a flowCore::flowSet")


}

##' @export
compensateFlowStats <- function(ff, ...) {
  stop("FlowStats compensation method not yet implemented")
}

##' @export
compensateAutoSpill <- function(ff, ...) {
  stop("Autospill compensation method not yet implemented")
}

### FUNCTIONS FOR DOUBLETS REMOVAL ###

##' @export
removeDoubletsPeacoQC <- function(ff,
                                  areaChannels,
                                  heightChannels,
                                  nmads,
                                  verbose = TRUE,
                                  ...) {

  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  # validate common scatter channel parameters
  nScatterFilters <- length(areaChannels)
  if(nScatterFilters < 1 || nScatterFilters > 2) {
    stop("nb of scatter channels for doublets removal should be either 1 or 2!")
  }
  if (length(heightChannels) != nScatterFilters) {
    stop("inconsistency between length of area and height channel vectors!")
  }

  if (length(nmads) != nScatterFilters) {
    stop("inconsistency between length of area channel and nMAD vectors!")
  }
  for (i in 1:nScatterFilters){
    ff <-
      PeacoQC::RemoveDoublets(ff,
                              channel1 = areaChannels[i],
                              channel2 = heightChannels[i],
                              nmad = nmads[i],
                              verbose = verbose)
  }

  return(ff)
}

##' @export
removeDoubletsCytoPipeline <- function(ff,
                                   areaChannels,
                                   heightChannels,
                                   nmads,
                                   verbose = TRUE,
                                   ...) {
  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  # validate common scatter channel parameters
  nScatterFilters <- length(areaChannels)
  if(nScatterFilters < 1 || nScatterFilters > 2) {
    stop("nb of scatter channels for doublets removal should be either 1 or 2!")
  }
  if (length(heightChannels) != nScatterFilters) {
    stop("inconsistency between length of area and height channel vectors!")
  }

  if (length(nmads) != nScatterFilters) {
    stop("inconsistency between length of area channel and nMAD vectors!")
  }
  for (i in 1:nScatterFilters){

    currentSingletGate <-
      singletsGate(ff,
                   filterId = paste0("Singlets_",
                                     areaChannels[i]),
                   channel1 = areaChannels[i],
                   channel2 = heightChannels[i],
                   nmad = nmads[i],
                   verbose = verbose)
    
    
    if (i==1) {
      singletGateCombined <- currentSingletGate
    } else {
      singletGateCombined <- singletGateCombined & currentSingletGate
    }
  }

  fltSinglet <- flowCore::filter(ff, singletGateCombined)

  ff <- ff[fltSinglet@subSet, ]

  return(ff)
}

##' @export
##'
removeDoubletsFlowStats <- function(ff,
                                    areaChannels,
                                    heightChannels,
                                    widerGate = FALSE,
                                    verbose = TRUE,
                                    ...) {
  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  # validate common scatter channel parameters
  nScatterFilters <- length(areaChannels)
  if(nScatterFilters < 1 || nScatterFilters > 2) {
    stop("nb of scatter channels for doublets removal should be either 1 or 2!")
  }
  if (length(heightChannels) != nScatterFilters) {
    stop("inconsistency between length of area and height channel vectors!")
  }

  for (i in 1:nScatterFilters){
    if (verbose) {
      currentSingletGate <-
        flowStats::singletGate(ff,
                               filterId = paste0("Singlets_",
                                                 areaChannels[i]),
                               area = areaChannels[i],
                               height = heightChannels[i],
                               wider_gate = widerGate)
    } else {
      currentSingletGate <- suppressMessages(
        flowStats::singletGate(ff,
                               filterId = paste0("Singlets_",
                                                 areaChannels[i]),
                               area = areaChannels[i],
                               height = heightChannels[i],
                               wider_gate = widerGate)
      )
    }

    if (i==1) {
      singletGateCombined <- currentSingletGate
    } else {
      singletGateCombined <- singletGateCombined & currentSingletGate
    }
  }

  fltSinglet <- flowCore::filter(ff, singletGateCombined)

  ff <- ff[fltSinglet@subSet, ]

  return(ff)
}

### FUNCTIONS FOR DEBRIS REMOVAL ###

##' @export
removeDebrisManual <- function(ff,
                               FSCChannel,
                               SSCChannel,
                               gateData,
                               verbose = TRUE,
                               ...) {
  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  cellsGateMatrix <- matrix(data = gateData, ncol = 2,
                            dimnames = list(c(), c(FSCChannel, SSCChannel)))

  if (verbose) {
    cellsGate <- flowCore::polygonGate(filterId = "Cells",
                                       .gate = cellsGateMatrix)
    selectedCells <- suppressMessages(
      flowCore::filter(ff, cellsGate))
  } else {
    cellsGate <- suppressMessages(
      flowCore::polygonGate(filterId = "Cells",
                            .gate = cellsGateMatrix)
    )
    selectedCells <- suppressMessages(
      flowCore::filter(ff, cellsGate))
  }


  ff <- ff[selectedCells@subSet, ]
}


##' @export
removeDebrisFlowClustTmix <- function(ff,
                                      FSCChannel,
                                      SSCChannel,
                                      nClust,
                                      verbose = TRUE,
                                      ...) {

  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  # handle ellipsis arguments, as 'tmixFilter' does not accept unknown args
  passedEllipsisArgs <- list(...)
  newEllipsisArgs <- list()

  argNames <-
    c("expName", "K", "B", "tol", "nu", "lambda", "nu.est", "trans",
      "min.count", "max.count", "min", "max", "level", "u.cutoff", "z.cutoff",
      "randomStart", "B.init", "tol.init", "seed", "criterion")
  for (argN in argNames){
    if (!is.null(passedEllipsisArgs[[argN]])) {
      newEllipsisArgs[[argN]] <- passedEllipsisArgs[[argN]]
    }
  }


  if (verbose) {
    cellsFilter <-
      do.call(flowClust::tmixFilter,
              args = c(list(filterId = "tmixFilter",
                            parameters = c(FSCChannel, SSCChannel),
                            K = nClust),
                       newEllipsisArgs))

  } else {
    cellsFilter <- suppressMessages(
      do.call(flowClust::tmixFilter,
              args = c(list(filterId = "tmixFilter",
                            parameters = c(FSCChannel, SSCChannel),
                            K = nClust),
                       newEllipsisArgs))
    )
  }
  if (verbose) {
    resCellsFilter <- flowCore::filter(ff, cellsFilter)
  } else {
    resCellsFilter <- suppressMessages(flowCore::filter(ff, cellsFilter))
  }


  FSCMedians <- sapply(X = 1:nClust,
                       FUN = function(x, ff, flt){
                         resCellsFltr <- flt[[x]]
                         median(flowCore::exprs(ff)[resCellsFltr@subSet, FSCChannel],
                                na.rm = TRUE)

                       },
                       ff = ff, flt = resCellsFilter)

  debrisIndex <- which.min(FSCMedians)
  keptClustersIndexes <- setdiff(1:nClust,debrisIndex)
  tokeepFilter <- resCellsFilter[[keptClustersIndexes[1]]]
  if(nClust > 2){
    for(i in keptClustersIndexes[-1]){
      tokeepFilter <- tokeepFilter | resCellsFilter[[i]]
    }
  }
  selectedCells <- flowCore::filter(ff, tokeepFilter)
  ff <- ff[selectedCells@subSet, ]

  return(ff)

}

### FUNCTIONS FOR DEAD CELLS REMOVAL ###

##' @export
removeDeadCellsManual <- function(ff,
                                  preTransform = FALSE,
                                  transList = NULL,
                                  FSCChannel,
                                  LDMarker,
                                  gateData,
                                  verbose = TRUE,
                                  ...) {
  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  if (preTransform){
    if (is.null(transList)) {
      stop("tranformation list needs to be provided if preTransform = TRUE!")
    }
    ffIn <- flowCore::transform(ff, transList)
  } else {
    ffIn <- ff
  }

  LDChannel <- FlowSOM::GetChannels(ffIn,
                                    markers = LDMarker,
                                    exact = TRUE)

  liveGateMatrix <- matrix(data = gateData, ncol = 2,
                           dimnames = list(c(), c(FSCChannel,
                                                  LDChannel)))

  if (verbose) {
    liveGate <- flowCore::polygonGate(filterId = "Live_Cells",
                                      .gate = liveGateMatrix)
  } else {
    liveGate <- suppressMessages(
      flowCore::polygonGate(filterId = "Live_Cells",
                                      .gate = liveGateMatrix))
  }

  if (verbose) {
    selectedLive <- flowCore::filter(ffIn, liveGate)
  } else {
    selectedLive <- suppressMessages(flowCore::filter(ffIn, liveGate))
  }


  ff <- ff[selectedLive@subSet, ] # note we take ff and not ffIn (no transfo)
}

##' @export
removeDeadCellsGateTail <- function(ff,
                                    preTransform = FALSE,
                                    transList = NULL,
                                    LDMarker,
                                    verbose = TRUE,
                                    ...){

  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  # handle ellipsis arguments, as 'openCyto::gate_tail' does not accept unknown
  #args
  passedEllipsisArgs <- list(...)
  newEllipsisArgs <- list()

  argNames <-
    c("num_peaks", "ref_peak", "strict", "tol", "side", "min", "max", "bias",
      "positive", "deriv", "bandwidth", "adjust", "num_points", "range.x",
      "binned", "se", "w")
  for (argN in argNames){
    if (!is.null(passedEllipsisArgs[[argN]])) {
      newEllipsisArgs[[argN]] <- passedEllipsisArgs[[argN]]
    }
  }

  if (preTransform){
    if (is.null(transList)) {
      stop("tranformation list needs to be provided if preTransform = TRUE!")
    }
    ffIn <- flowCore::transform(ff, transList)
  } else {
    ffIn <- ff
  }

  LDChannel <- FlowSOM::GetChannels(ffIn,
                                    markers = LDMarker,
                                    exact = TRUE)

  if (verbose) {
    liveGate <-
      do.call(openCyto::gate_tail,
              args = c(list(ffIn, channel = LDChannel, filterId = "Live_Cells"),
                       newEllipsisArgs))
  } else {
    liveGate <- suppressMessage(
      do.call(openCyto::gate_tail,
              args = c(list(ffIn, channel = LDChannel, filterId = "Live_Cells"),
                       newEllipsisArgs))
    )
  }
  if (verbose) {
    selectedLive <- flowCore::filter(ffIn, liveGate)
  } else {
    selectedLive <- suppressMessage(flowCore::filter(ffIn, liveGate))
  }

  ff <- ff[selectedLive@subSet, ] # note we take ff and not ffIn (no transfo)
  return(ff)
}

### FUNCTIONS for Quality Control ###

##' @export
qualityControlFlowAI <- function(ff,
                                 preTransform = FALSE,
                                 transList = NULL,
                                 outputDiagnostic = FALSE,
                                 outputDir = NULL,
                                 verbose = TRUE,
                                 ...) {
  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  if (preTransform) {
    if (is.null(transList)) {
      stop("tranformation list needs to be provided if preTransform = TRUE!")
    }
    ffIn <- flowCore::transform(ff, transList)
  } else {
    ffIn <- ff
  }

  # # check if time channel is present, and monotonically increasing
  # timeCh <- FF_findTimeChannel(ffIn)
  #
  # if (is.null(timeCh)) {
  #   stop(paste0("Time channel not found in flow frame ",
  #               flowCore::identifier(ffIn)))
  # }

  channel2Exclude <-
    flowCore::colnames(ffIn)[!FF_areSignalCols(ffIn)]
  message("Applying flowAI method...")
  if (outputDiagnostic) {
    html_report = "_QC"
    mini_report = "QCmini"
    if (!is.null(outputDir))
      folder_results = outputDir
    else folder_results = "resultsQC"
  } else {
    html_report = FALSE
    mini_report = FALSE
    folder_results = FALSE
  }

  if (verbose) {
    badEventIDs <-
      flowAI::flow_auto_qc(fcsfiles = ffIn,
                           output = 3,
                           timeCh = NULL,
                           timestep = NULL,
                           ChExcludeFS = channel2Exclude,
                           ChExcludeFM = channel2Exclude,
                           html_report = html_report,
                           mini_report = mini_report,
                           fcs_QC = FALSE,
                           fcs_highQ = FALSE,
                           fcs_lowQ = FALSE,
                           folder_results = folder_results,
                           ...)
  } else {
    badEventIDs <-
      flowAI::flow_auto_qc(fcsfiles = ffIn,
                           output = 3,
                           timeCh = NULL,
                           timestep = NULL,
                           ChExcludeFS = channel2Exclude,
                           ChExcludeFM = channel2Exclude,
                           html_report = html_report,
                           mini_report = mini_report,
                           fcs_QC = FALSE,
                           fcs_highQ = FALSE,
                           fcs_lowQ = FALSE,
                           folder_results = folder_results,
                           ...)
  }
  goodEvents <- !(1:(flowCore::nrow(ffIn)) %in% badEventIDs)
  ff <- ff[goodEvents,] # note we take ff and not ffIn (no transfo)

  return(ff)
}

##' @export
qualityControlPeacoQC <- function(ff,
                                  preTransform = FALSE,
                                  transList = NULL,
                                  outputDiagnostic = FALSE,
                                  outputDir = NULL,
                                  verbose = TRUE,
                                  ...) {
  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  if (preTransform) {
    if (is.null(transList)) {
      stop("tranformation list needs to be provided if preTransform = TRUE!")
    }
    ffIn <- flowCore::transform(ff, transList)
  } else {
    ffIn <- ff
  }

  # # check if time channel is present, and monotonically increasing
  # timeCh <- FF_findTimeChannel(ffIn)
  #
  # if (is.null(timeCh)) {
  #   stop(paste0("Time channel not found in flow frame ",
  #               flowCore::identifier(ffIn)))
  # }

  # qualityControl with PeacoQC
  message("Applying PeacoQC method...")
  channel4QualityControl <-
    flowCore::colnames(ffIn)[FF_areSignalCols(ffIn)]
  if (outputDiagnostic) {
    plot = TRUE
    report = TRUE
    if (!is.null(outputDir))
      output_directory = outputDir
    else output_directory = "."
  } else {
    plot = FALSE
    report = FALSE
    output_directory = NULL #not used
  }
  if (verbose) {
    res <- PeacoQC::PeacoQC(ff = ffIn,
                            channels = channel4QualityControl,
                            report = report,
                            plot = plot,
                            save_fcs = FALSE,
                            output_directory = output_directory,
                            ...)
  } else {
    res <- suppressMessages(
      PeacoQC::PeacoQC(ff = ffIn,
                       channels = channel4QualityControl,
                       report = report,
                       plot = plot,
                       save_fcs = FALSE,
                       output_directory = output_directory,
                       ...)
    )
  }

  ff <- ff[res$GoodCells,] # note we take ff and not ffIn (no transfo)
  return(ff)
}

##' @export
qualityControlFlowCut <- function(ff,
                                  preTransform = FALSE,
                                  transList = NULL,
                                  outputDiagnostic = FALSE,
                                  outputDir = NULL,
                                  verbose = TRUE,
                                  ...) {
  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  if (preTransform) {
    if (is.null(transList)) {
      stop("tranformation list needs to be provided if preTransform = TRUE!")
    }
    ffIn <- flowCore::transform(ff, transList)
  } else {
    ffIn <- ff
  }



  # # check if time channel is present, and monotonically increasing
  # timeCh <- FF_findTimeChannel(ffIn)
  #
  # if (is.null(timeCh)) {
  #   stop(paste0("Time channel not found in flow frame ",
  #               flowCore::identifier(ffIn)))
  # }

  # qualityControl with flowCut
  message("Applying flowCut method...")
  channelsIndices <- which(FF_areSignalCols(ffIn))
  if (outputDiagnostic) {
    Plot = "All"
    if (!is.null(outputDir))
      Directory = outputDir
    else filePrefixWithDir = "resultsQC"
  } else {
    Directory = NULL #not used
    Plot = "None"
  }

  res <-
    flowCut::flowCut(f = ffIn,
                     Channels = channelsIndices,
                     Directory = Directory,
                     FileID = NULL,
                     Plot = Plot,
                     Verbose = verbose,
                     ...)
  #browser()
  badEventIDs <- res$ind

  goodEvents <- !(1:(flowCore::nrow(ffIn)) %in% badEventIDs)
  ff <- ff[goodEvents,] # note we take ff and not ffIn (no transfo)

  return(ff)
}

##' @export
qualityControlFlowClean <- function(ff,
                                    preTransform = FALSE,
                                    transList = NULL,
                                    outputDiagnostic = FALSE,
                                    outputDir = NULL,
                                    verbose = TRUE,
                                    ...) {
  # if not present already, add a column with Cell ID
  ff <- FF_AppendCellID(ff, 1:(flowCore::nrow(ff)))

  if (preTransform) {
    if (is.null(transList)) {
      stop("tranformation list needs to be provided if preTransform = TRUE!")
    }
    ffIn <- flowCore::transform(ff, transList)
  } else {
    ffIn <- ff
  }



  # # check if time channel is present, and monotonically increasing
  # timeCh <- FF_findTimeChannel(ffIn)
  #
  # if (is.null(timeCh)) {
  #   stop(paste0("Time channel not found in flow frame ",
  #               flowCore::identifier(ffIn)))
  # }

  # qualityControl with flowClean
  message("Applying flowClean method...")
  vectMarkers <- which(FF_areSignalCols(ffIn))

  if (outputDiagnostic) {
    diagnostic = TRUE
    if (!is.null(outputDir))
      filePrefixWithDir = outputDir
    else filePrefixWithDir = "resultsQC"

    # add original fcs file name in prefix, as flowClean is designed to work
    # for one fcs at the time
    filename <- basename(flowCore::keyword(ffIn, "FILENAME")$FILENAME)
    # removing extension
    filename <- sub("([^.]+)\\.[[:alnum:]]+$", "\\1", filename)
    filePrefixWithDir <- paste0(filePrefixWithDir, filename)
  } else {
    filePrefixWithDir = NULL #not used
    diagnostic = FALSE
  }

  goodVsBadVector <-
    flowClean::clean(fF = ffIn,
                     vectMarkers = vectMarkers,
                     filePrefixWithDir = filePrefixWithDir,
                     ext = ".fcs", # not used
                     diagnostic = diagnostic,
                     announce = verbose,
                     returnVector = TRUE,
                     ...)

  areGoodEvents <- goodVsBadVector < 10000
  ff <- ffIn[areGoodEvents, ]

  return(ff)
}

##' @export
applyScaleTransforms <- function(ff, transList, ...){
  ff <- flowCore::transform(ff, transList)
  return(ff)
}

##' @export
estimateScaleTransforms <- function(ff, fluoMethod = c("estimateLogicle", "none"),
                               scatterMethod = c("linear", "none"),
                               scatterRefMarker = NULL) {
  fluoMethod <- match.arg(fluoMethod)
  scatterMethod <- match.arg(scatterMethod)

  if (fluoMethod == "estimateLogicle") {
    message("estimating logicle transformations for fluorochrome channels...")
    fluoCols <- flowCore::colnames(ff)[FF_areFluoCols(ff)]
    transList <- flowCore::estimateLogicle(ff, fluoCols)
  } # else do nothing

  if (scatterMethod == "linear") {
    if (is.null(scatterRefMarker))
      stop("linear scatter method requires a scatterRefMarker")
    msg <-
      "Estimating linear transformation for scatter channels : "
    msg <- paste0(msg,
                  "reference marker = ",
                  scatterRefMarker,
                  "..."
    )
    message(msg)
    transList <-
      FF_computeScatterChannelsLinearScale(ff,
                                           trans_list = transList,
                                           reference_channel = scatterRefMarker,
                                           silent = FALSE)
  } # else do nothing

  return(transList)
}



