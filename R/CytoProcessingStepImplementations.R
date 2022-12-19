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

#' @title estimates scale tranformations
#' @description this function estimates the scale transformations to be applied
#' on a flowFrame to obtain 'good behaving' distributions, i.e. the best
#' possible separation between + population and - population.
#' It distinguishes between scatter channels, where either linear, or no
#' transform is applied, and fluo channels, where either logicle transform
#' - using flowCore::estimateLogicle - is estimated, or no transform   
#' is applied.   
#' 
#' The idea of linear transform of scatter channels is as follows: a reference
#' channel (not a scatter one) is selected and a linear transform (Y = AX + B)
#' is applied to all scatter channel, as to align their 5 and 95 percentiles to
#' those of the reference channel
#' For the estimateLogicle function, see flowCore documentation.

#' @param ff a flowCore::flowFrame
#' @param fluoMethod method to be applied to all fluo channels
#' @param scatterMethod method to be applied to all scatter channels
#' @param scatterRefMarker the reference channel that is used to align the
#' @param specificScatterChannels vector of scatter channels for which we 
#' still want to apply the fluo method (and not the scatter Method)
#'
#' @return a flowCore::flowFrame with removed low quality events from the input
#' @export
#'
#' @examples
#' 
#' compMatrix <- flowCore::spillover(OMIP021Samples[[1]])$SPILL
#' ff_c <- runCompensation(OMIP021Samples[[1]], spillover = compMatrix)
#' 
#' transList <- 
#'     estimateScaleTransforms(        
#'         ff = ff_c,
#'         fluoMethod = "estimateLogicle",
#'         scatterMethod = "linear",
#'         scatterRefMarker = "BV785 - CD3")
#' 
estimateScaleTransforms <- function(ff,
                                    fluoMethod = c("estimateLogicle",
                                                   "none"),
                                    scatterMethod = c("linearQuantile", 
                                                      "none"),
                                    scatterRefMarker = NULL,
                                    specificScatterChannels = NULL){
    fluoMethod <- match.arg(fluoMethod)
    scatterMethod <- match.arg(scatterMethod)
    
    if (fluoMethod == "estimateLogicle") {
        message(
            "estimating logicle transformations ",
            "for fluorochrome channels..."
        )
        fluoCols <- flowCore::colnames(ff)[areFluoCols(ff)]
        transList <- flowCore::estimateLogicle(ff, fluoCols)
    } # else do nothing
    
    if (scatterMethod == "linearQuantile") {
        if (is.null(scatterRefMarker)) {
            stop("linear scatter method requires a scatterRefMarker")
        }
        
        message(
            "Estimating linear transformation for scatter channels : ",
            "reference marker = ",
            scatterRefMarker,
            "..."
        )
        transList <-
            computeScatterChannelsLinearScale(
                ff,
                transList = transList,
                referenceChannel = scatterRefMarker,
                silent = FALSE
            )
    } # else do nothing
    
    # handle specific cases of scatter channels that still need fluo method
    if (! is.null(specificScatterChannels)) {
        effectiveScatterChannels <- NULL
        for (ch in specificScatterChannels) {
            scatterChannels <-
                flowCore::colnames(ff)[!areFluoCols(ff) & areSignalCols(ff)]
            if (!(ch %in% scatterChannels)) {
                message("Specific channel [", ch, "] is not a scatter channel",
                        " => no correction of scale transformation done")
            } else {
                transList@transforms[[ch]] <- NULL
                effectiveScatterChannels <- c(effectiveScatterChannels, ch)
            } 
        }
        
        transList <- c(transList, 
                       flowCore::estimateLogicle(ff, 
                                                 effectiveScatterChannels))
    }
    
    return(transList)
}

#' @title select a nb of sample files randomly
#' @description Wrapper around sample(sampleFiles, nSamples).
#' For the needs of making it a CytoProcessingStep:
#' - takes `sampleFiles` as an explicit first parameter
#' - can manage a seed
#' - allow passing ... as a parameter (not used)
#' @param sampleFiles a vector of character path to sample files
#' @param nSamples number of samples to randomly select. If `nSamples` is
#' higher than nb of available samples, the output will be all samples
#' @param seed an optional seed parameters (provided to ease reproducibility).
#' @param ... additional parameters passed (not used here).
#'
#' @return a subset of `sampleFiles`, randomly selected
#' @export
#' 
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' nRandomSamples <- 1
#' selectSampleFiles <- selectRandomSamples(sampleFiles, 
#'                                          nSamples = nRandomSamples,
#'                                          seed = 1)
#'                                          
selectRandomSamples <- function(sampleFiles,
                                nSamples,
                                seed = NULL,
                                ...) {
    
    if (!is.numeric(nSamples) || nSamples < 1) {
        stop("[nSamples] should be a numeric >= 1")
    }
    
    nAvailableSamples <- length(sampleFiles)
    
    if (nSamples > nAvailableSamples) nSamples <- nAvailableSamples
    
    if (!is.null(seed)) {
        # set the seed locally in the execution environment,
        # restore it afterward
        withr::local_seed(seed)
    }
    
    outputSamples <- sample(sampleFiles, nSamples)
    
    outputSamples
}


#' @title Read fcs sample files
#' @description Wrapper around flowCore::read.fcs() or flowCore::read.flowSet().
#' Also adds a "Cell_ID" additional column, used in flowFrames comparison
#' @param sampleFiles a vector of character path to sample files
#' @param whichSamples either 'all' if all sample files need to be read, or
#' a vector of indexes pointing to the sampleFiles vector
#' @param channelMarkerFile an optional path to a csv file which provides the 
#' mapping between channels and markers. If provided, this csv file should 
#' contain a `Channel` column, and a `Marker` column. Optionally a 'Used' 
#' column can be provided as well (TRUE/FALSE). Channels for which the 'Used' 
#' column is set to FALSE will not be incorporated in the created flowFrame. 
#' @param ... additional parameters passed to flowCore file reading functions.
#'
#' @return either a flowCore::flowSet or a flowCore::flowFrame if
#' length(sampleFiles) == 1
#' @export
#' 
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' res <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' #res
#' 
#' # create a flowCore::flowFrame with one single sample
#' res2 <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = 2,
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' #res2
#' 
readSampleFiles <- function(sampleFiles,
                            whichSamples = "all", 
                            channelMarkerFile = NULL,
                            ...) {
    if (whichSamples == "all") {
        # do nothing : sampleFiles should contain all the input sample files
        # already
    } else if (is.numeric(whichSamples)) {
        sampleFiles <- sampleFiles[whichSamples]
    } else {
        stop("'whichSamples' should be either 'all', or a vector of indexes")
    }

    if (length(sampleFiles) == 1) {
        res <- flowCore::read.FCS(sampleFiles, ...)
        # Add a column with Cell ID
        res <- .appendCellID(res)
    } else {
        res <- flowCore::read.flowSet(sampleFiles, ...)
        # Add a column with Cell ID
        res <- flowCore::fsApply(
            x = res,
            FUN = function(ff) {
                .appendCellID(ff)
            }
        )
    } 
    
    # do we need to do any post processing to the files ?
    # => remove channels or update marker names ? 
    
    if (!is.null(channelMarkerFile)) {
        if (!file.exists(channelMarkerFile)) {
            stop("channelMarkerFile [", channelMarkerFile, "] not found!")
        }
        channelMarkerMapping <- utils::read.csv(channelMarkerFile)
        
        if (!("Channel" %in% flowCore::colnames(channelMarkerMapping))) {
            stop("channelMarkersMapping should contain [Channel] column!")
        }
        if (!("Marker" %in% flowCore::colnames(channelMarkerMapping))) {
            stop("channelMarkersMapping should contain [Marker] column!")
        }
        if ("Used" %in% colnames(channelMarkerMapping)) {
            channelMarkerMapping$Used <- as.logical(channelMarkerMapping$Used)
            # if (!is.logical(channelMarkerMapping$Used)) {
            #     stop("channelMarkerMapping [Used] column ",
            #          "should be of logical type!")
            # }
            notUsedChannels <- 
                channelMarkerMapping[!channelMarkerMapping$Used, "Channel"]
            
        } 
        
        postProcessing <- function(ff, channelMarkerMapping, notUsedChannels){
            ff <- removeChannels(ff, notUsedChannels)
            
            channels <- 
                channelMarkerMapping[channelMarkerMapping$Used, "Channel"]
            markers <- 
                channelMarkerMapping[channelMarkerMapping$Used, "Marker"]
                
            for (i in seq_along(channels)) {
                if (markers[i] != "") {
                    ff <- updateMarkerName(ff, channel = channels[i], 
                                           newMarkerName = markers[i])
                }
            }
            return(ff)
        }
        
        if (length(sampleFiles) == 1) {
            res <- postProcessing(res, channelMarkerMapping, notUsedChannels)
        } else {
            res <- flowCore::fsApply(
                x = res,
                FUN = function(ff) {
                    postProcessing(ff, channelMarkerMapping, notUsedChannels)
                }
            )
        } 
    }
    
    return (res)
}


# .updateChannelsAndMarkers <- function(ff, channelMarkersMapping) {
#     if (!inherits(ff, "flowFrame")) {
#         stop("ff type not recognized, should be a flowFrame!")
#     }
#     
#     if (!is.data.frame(channelMarkersMapping)) {
#         stop("channelMarkersMapping should be a dataframe!")
#     }
#     
#     if (!("Channel" %%in%% colnames(channelMarkersMapping))) {
#         stop("channelMarkersMapping should contain [Channel] column!")
#     }
#     if (!("Marker" %%in%% colnames(channelMarkersMapping))) {
#         stop("channelMarkersMapping should contain [Channel] column!")
#     }
#     if ("Used" %%in%% colnames(channelMarkersMapping)) {
#         if (!is.logical(channelMarkersMapping$Used)) {
#             stop("channelMarkersMapping [Used] column ",
#                  "should be of logical type!")
#         }
#     } else {
#         # if no 'Used' column, assume all channels are used by default
#         channelMarkersMapping$Used <- rep(TRUE, nrow(channelMarkersMapping))
#     }
#     
#     rownames(channelMarkersMapping) <- channelMarkersMapping[,"Channel"]
#     
#     param <- flowCore::parameters(ff)
#     paramData <- flowCore::pData(param)
#     keptChannels <- flowCore::colnames(ff)
#     
#     for (ch in rownames(channelMarkersMapping)) {
#         
#         used <- channelMarkersMapping[ch, "Used"]
#         
#         if (used == "") used <- FALSE
#         channelLine <- which(paramData[, "name"] == ch)
#         if (length(channelLine) == 0)
#         {
#             stop("channel ", ch, ", mentioned in channel marker mapping file, ",
#                  "not found in flowFrame!")
#         } else {
#             channelLine <- channelLine[1]
#         }
#         
#         if (!used) {
#             keptChannels <- keptChannels[-which(keptChannels == ch)]
#         } else {
#             marker <- channelMarkersMapping[ch, "Marker"]
#             if (marker == "") marker <- NA
#             
#             paramData[channelLine, "desc"] <- marker
#         }
#     }
#     
#     flowCore::pData(param) <- paramData
#     flowCore::parameters(ff) <- param
#     
#     ff <- ff[, keptChannels]
#     
#     return(ff)
# } 

#' @title remove margin events using PeacoQC
#' @description Wrapper around PeacoQC::RemoveMargins().
#' Also pre-selects the channels to be handled (=> all signal channels)
#' If input is a flowSet, it applies removeMargins() to each flowFrame of the
#' flowSet.
#' @param x a flowCore::flowSet or a flowCore::flowFrame
#' @param ... additional parameters passed to PeacoQC::RemoveMargins()
#'
#' @return either a flowCore::flowSet or a flowCore::flowFrame depending on
#' the input.
#' 
#' @export
#'
#' @examples
#' 
#' rawDataDir <- 
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <- 
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' fsRaw <- readSampleFiles(sampleFiles, 
#'                          truncate_max_range = truncateMaxRange,
#'                          min.limit = minLimit)
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#' ggplotFilterEvents(ffPre = fsRaw[[2]],
#'                    ffPost = ff_m,
#'                    xChannel = "FSC-A",
#'                    yChannel = "SSC-A")
removeMarginsPeacoQC <- function(x, ...) {
    myFunc <- function(ff) {
        message("Removing margins from file : ", flowCore::identifier(ff))
        channel4Margins <-
            flowCore::colnames(ff)[areSignalCols(ff)]
        ffOut <- PeacoQC::RemoveMargins(ff, channels = channel4Margins)
        return(ffOut)
    }
    if (inherits(x, "flowFrame")) {
        return(myFunc(x))
    } else if (inherits(x, "flowSet")) {
        fsOut <- flowCore::fsApply(x, FUN = myFunc, simplify = TRUE)
        return(fsOut)
    } else {
        stop("x should be a flowCore::flowFrame or a flowCore::flowSet")
    }
}

### FUNCTIONS for pre-processing / compensation ###

#' @title extract compensation matrix from a flowCore::flowFrame
#' @description helper function retrieving the compensation matrix stored
#' in fcs file (if any). It scans the following keywords: $SPILL, $spillover
#' and $SPILLOVER
#' @param ff a flowCore::flowFrame
#'
#' @return the found compensation matrix
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' compensationMatrix <- getAcquiredCompensationMatrix(fsRaw[[2]])
getAcquiredCompensationMatrix <- function(ff) {
    
    stopifnot (inherits(ff, "flowFrame"))
    
    res <- flowCore::spillover(ff)
    if (!is.null(res$SPILL)) {
        compensationMatrix <- res$SPILL
    } else if (!is.null(res$spillover)) {
        compensationMatrix <- res$spillover
    } else if (!is.null(res$`$SPILLOVER`)) {
        compensationMatrix <- res$`$SPILLOVER`
    } else {
        fileId <- flowCore::identifier(ff)
        stop(
            "Issue retrieving compensation matrix for file ",
            fileId, " : slot is NULL!"
        )
    }
    return(compensationMatrix)
}

#' @title compensation of fcs file(s) from matrix
#' @description executes the classical compensation function on a flowSet or
#' flowFrame, given a compensation matrix. The matrix can be either retrieved
#' in the fcs files themselves or provided as a csv file.
#' @param x a flowCore::flowSet or a flowCore::flowFrame
#' @param matrixSource if "fcs", the compensation matrix will be fetched from
#' the fcs files (different compensation matrices can then be applied by fcs
#' file)
#' if "import", uses matrixPath to read the matrix (should be a csv file)
#' @param matrixPath if matrixSource == "import", will be used as the input csv
#' file path
#' @param updateChannelNames if TRUE, updates the fluo channel names by
#' prefixing them with "comp-"
#' @param ... additional arguments (not used)
#' @return the compensated flowSet or flowFrame
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
compensateFromMatrix <- function(x,
                                 matrixSource = c("fcs", "import"),
                                 matrixPath = NULL,
                                 updateChannelNames = TRUE,
                                 ...) {
    myFunc <- function(ff, matrixSource = c("fcs", "import"),
                       matrixPath = NULL) {
        message("Compensating file : ", flowCore::identifier(ff))
        matrixSource <- match.arg(matrixSource)
        if (matrixSource == "fcs") {
            # obtains compensation matrix
            compensationMatrix <-
                getAcquiredCompensationMatrix(ff)
        } else {
            # import matrix from file (path)
            if (is.null(matrixPath)) {
                stop("No path specified for compensation matrix!")
            }
            if (!file.exists(matrixPath)) {
                stop("Compensation matrix file not found!")
            }
            compensationMatrix <-
                utils::read.csv(matrixPath,
                    check.names = FALSE,
                    row.names = 1
                )
        }

        ffOut <- runCompensation(ff,
            compensationMatrix,
            updateChannelNames = updateChannelNames
        )

        return(ffOut)
    }

    if (inherits(x, "flowFrame")) {
        return(myFunc(x, matrixSource = matrixSource, matrixPath = matrixPath))
    } else if (inherits(x, "flowSet")) {
        fsOut <- flowCore::fsApply(x, FUN = myFunc, simplify = TRUE)
        return(fsOut)
    } else {
        stop("x should be a flowCore::flowFrame or a flowCore::flowSet")
    }
}


### FUNCTIONS FOR DOUBLETS REMOVAL ###

#' @title remove doublets from a flowFrame, using PeacoQC
#' @description wrapper around PeacoQC::RemoveDoublets().
#' Can apply the PeacoQC function subsequently on several channel pairs,
#' e.g. (FSC-A, FSC-H) and (SSC-A, SSC-H)
#' @param ff a flowCore::flowFrame
#' @param areaChannels a character vector containing the name of the 'area type'
#' channels one wants to use
#' @param heightChannels a character vector containing the name of the
#' 'height type' channels one wants to use
#' @param nmads a numeric vector with the bandwidth above the ratio allowed, per
#' channels pair (cells are kept if the ratio between -A channel\[i\] and
#' -H channel\[i\] is smaller than the median ratio + nmad\[i\] times the median
#' absolute deviation of the ratios). Default is 4, for all channel pairs.
#' @param verbose If set to TRUE, the median ratio and width will be printed.
#' @param ... additional parameters passed to PeacoQC::RemoveDoublets()
#'
#' @return a flowCore::flowFrame with removed doublets events from the input
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#'
#' ff_s <-
#'     removeDoubletsPeacoQC(ff_c,
#'                           areaChannels = c("FSC-A", "SSC-A"),
#'                           heightChannels = c("FSC-H", "SSC-H"),
#'                           nmads = c(3, 5))                            
removeDoubletsPeacoQC <- function(ff,
                                  areaChannels,
                                  heightChannels,
                                  nmads = rep(4, length(areaChannels)),
                                  verbose = TRUE,
                                  ...) {

    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    # validate common scatter channel parameters
    nScatterFilters <- length(areaChannels)
    if (nScatterFilters < 1 || nScatterFilters > 2) {
        stop(
            "nb of scatter channels for doublets removal ",
            "should be either 1 or 2!"
        )
    }
    if (length(heightChannels) != nScatterFilters) {
        stop(
            "inconsistency between length of area and ",
            "height channel vectors!"
        )
    }

    if (length(nmads) != nScatterFilters) {
        stop("inconsistency between length of area channel and nMAD vectors!")
    }
    for (i in seq_len(nScatterFilters)) {
        ff <-
            PeacoQC::RemoveDoublets(ff,
                channel1 = areaChannels[i],
                channel2 = heightChannels[i],
                nmad = nmads[i],
                verbose = verbose
            )
    }

    return(ff)
}



#' @title remove doublets from a flowFrame, using flowStats
#' @description Wrapper around flowStats::singletGate().
#' Can apply the flowStats function subsequently on several channel pairs,
#' e.g. (FSC-A, FSC-H) and (SSC-A, SSC-H)
#' @param ff a flowCore::flowFrame
#' @param areaChannels a character vector containing the name of the 'area type'
#' channels one wants to use
#' @param heightChannels a character vector containing the name of the
#' 'height type' channels one wants to use
#' @param widerGate a boolean as wider_gate parameter to
#' flowStats::singletGate()
#' @param ... additional parameters passed to flowStats::singletGate()
#'
#' @return a flowCore::flowFrame with removed doublets events from the input
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#' 
#' ff_s <-
#'     removeDoubletsFlowStats(ff_c,
#'                             areaChannels = c("FSC-A", "SSC-A"),
#'                             heightChannels = c("FSC-H", "SSC-H"),
#'                             widerGate = TRUE)
#'                             
removeDoubletsFlowStats <- function(ff,
                                    areaChannels,
                                    heightChannels,
                                    widerGate = FALSE,
                                    ...) {
    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    # validate common scatter channel parameters
    nScatterFilters <- length(areaChannels)
    if (nScatterFilters < 1 || nScatterFilters > 2) {
        stop(
            "nb of scatter channels for doublets removal ",
            "should be either 1 or 2!"
        )
    }
    if (length(heightChannels) != nScatterFilters) {
        stop(
            "inconsistency between length of area ",
            "and height channel vectors!"
        )
    }

    for (i in seq_len(nScatterFilters)) {
        currentSingletGate <-
            flowStats::singletGate(ff,
                filterId = paste0(
                    "Singlets_",
                    areaChannels[i]
                ),
                area = areaChannels[i],
                height = heightChannels[i],
                wider_gate = widerGate
            )

        if (i == 1) {
            singletGateCombined <- currentSingletGate
        } else {
            singletGateCombined <- singletGateCombined & currentSingletGate
        }
    }

    fltSinglet <- flowCore::filter(ff, singletGateCombined)

    ff <- flowCore::Subset(ff, fltSinglet)
    #ff <- ff[fltSinglet@subSet, ]

    return(ff)
}

#' @title remove doublets from a flowFrame, using CytoPipeline custom algorithm
#' @description Wrapper around CytoPipeline::singletGate().
#' Can apply the flowStats function subsequently on several channel pairs,
#' e.g. (FSC-A, FSC-H) and (SSC-A, SSC-H)
#' @param ff a flowCore::flowFrame
#' @param areaChannels a character vector containing the name of the 'area type'
#' channels one wants to use
#' @param heightChannels a character vector containing the name of the
#' 'height type' channels one wants to use
#' @param nmads a numeric vector with the bandwidth above the ratio allowed, per
#' channels pair (cells are kept if the ratio between -A channel\[i\] and
#' -H channel\[i\] is smaller than the median ratio + nmad\[i\] times the median
#' absolute deviation of the ratios). Default is 4, for all channel pairs.
#' @param ... additional parameters passed to CytoPipeline::singletGate()
#'
#' @return a flowCore::flowFrame with removed doublets events from the input
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#' 
#' ff_s <-
#'     removeDoubletsCytoPipeline(ff_c,
#'                                areaChannels = c("FSC-A", "SSC-A"),
#'                                heightChannels = c("FSC-H", "SSC-H"),
#'                                nmads = c(3, 5))
#'                             
removeDoubletsCytoPipeline <- function(ff,
                                       areaChannels,
                                       heightChannels,
                                       nmads,
                                       ...) {
    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    # validate common scatter channel parameters
    nScatterFilters <- length(areaChannels)
    if (nScatterFilters < 1 || nScatterFilters > 2) {
        stop(
            "nb of scatter channels for doublets removal ",
            "should be either 1 or 2!"
        )
    }
    if (length(heightChannels) != nScatterFilters) {
        stop(
            "inconsistency between length of area ",
            "and height channel vectors!"
        )
    }

    if (length(nmads) != nScatterFilters) {
        stop("inconsistency between length of area channel and nMAD vectors!")
    }
    for (i in seq_len(nScatterFilters)) {
        currentSingletGate <-
            singletsGate(ff,
                filterId = paste0(
                    "Singlets_",
                    areaChannels[i]
                ),
                channel1 = areaChannels[i],
                channel2 = heightChannels[i],
                nmad = nmads[i]
            )


        if (i == 1) {
            singletGateCombined <- currentSingletGate
        } else {
            singletGateCombined <- singletGateCombined & currentSingletGate
        }
    }

    fltSinglet <- flowCore::filter(ff, singletGateCombined)

    ff <- flowCore::Subset(ff, fltSinglet)
    #ff <- ff[fltSinglet@subSet, ]

    return(ff)
}


### FUNCTIONS FOR DEBRIS REMOVAL ###

#' @title remove debris from a flowFrame using manual gating
#' @description remove debris from a flowFrame, using manual gating in the
#' FSC-A, SSC-A 2D representation. The function internally uses
#' flowCore::polygonGate()
#' @param ff a flowCore::flowFrame
#' @param FSCChannel a character containing the exact name of the forward
#' scatter channel
#' @param SSCChannel a character containing the exact name of the side scatter
#' channel
#' @param gateData a numerical vector containing the polygon gate coordinates
#' first the `FSCChannel` channel coordinates   
#' of each points of the polygon gate,   
#' then the `SSCChannel` channel coordinates of each points.
#' @param ... additional parameters passed to flowCore::polygonGate()
#'
#' @return a flowCore::flowFrame with removed debris events from the input
#' @export
#' 
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#' 
#' 
#' remDebrisGateData <- c(73615, 110174, 213000, 201000, 126000,
#'                        47679, 260500, 260500, 113000, 35000)
#' 
#' ff_cells <-
#'     removeDebrisManualGate(ff_c,
#'                            FSCChannel = "FSC-A",
#'                            SSCChannel = "SSC-A",
#'                            gateData = remDebrisGateData)
#' 
#'
removeDebrisManualGate <- function(ff,
                                   FSCChannel,
                                   SSCChannel,
                                   gateData,
                                   ...) {
    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    cellsGateMatrix <- matrix(
        data = gateData, ncol = 2,
        dimnames = list(c(), c(FSCChannel, SSCChannel))
    )

    cellsGate <- flowCore::polygonGate(
        filterId = "Cells",
        .gate = cellsGateMatrix,
        ...
    )
    selectedCells <- flowCore::filter(ff, cellsGate)

    ff <- flowCore::Subset(ff, selectedCells)
    #ff <- ff[selectedCells@subSet, ]
}


#' @title remove debris from a flowFrame, using flowClust
#' @description this function removes debris from a flowFrame,
#' using clustering capabilities of flowClust::tmixFilter(). The idea is to
#' pre-select a number of clusters to be found in the (FSC,SSC) 2D view, and
#' eliminate the cluster that is the closest to the origin.

#' @param ff a flowCore::flowFrame
#' @param FSCChannel the name of the FSC channel
#' @param SSCChannel the name of the SSC channel
#' @param nClust number of clusters to identify
#' @param ... additional parameters passed to flowClust::tmixFilter()
#'
#' @return a flowCore::flowFrame with removed debris events from the input
#' @export
#' 
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#' 
#' 
#' ff_cells <-
#'     removeDebrisFlowClustTmix(ff_c,
#'                               FSCChannel = "FSC-A",
#'                               SSCChannel = "SSC-A",
#'                               nClust = 3,
#'                               level = 0.97,
#'                               B = 100)
#' 
removeDebrisFlowClustTmix <- function(ff,
                                      FSCChannel,
                                      SSCChannel,
                                      nClust,
                                      ...) {

    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    # handle ellipsis arguments, as 'tmixFilter' does not accept unknown args
    passedEllipsisArgs <- list(...)
    newEllipsisArgs <- list()

    argNames <-
        c(
            "expName", "K", "B", "tol", "nu", "lambda", "nu.est", "trans",
            "min.count", "max.count", "min", "max", "level", "u.cutoff",
            "z.cutoff", "randomStart", "B.init", "tol.init", "seed", "criterion"
        )
    for (argN in argNames) {
        if (!is.null(passedEllipsisArgs[[argN]])) {
            newEllipsisArgs[[argN]] <- passedEllipsisArgs[[argN]]
        }
    }

    cellsFilter <-
        do.call(flowClust::tmixFilter,
            args = c(
                list(
                    filterId = "tmixFilter",
                    parameters = c(FSCChannel, SSCChannel),
                    K = nClust
                ),
                newEllipsisArgs
            )
        )


    resCellsFilter <- flowCore::filter(ff, cellsFilter)

    FSCMedians <- vapply(
        X = seq_len(nClust),
        FUN.VALUE = double(1),
        FUN = function(x, ff, flt) {
            resCellsFltr <- flt[[x]]
            # stats::median(flowCore::exprs(ff)[
            #     resCellsFltr@subSet, FSCChannel
            # ],
            # na.rm = TRUE
            # )
            stats::median(
                flowCore::exprs(
                    flowCore::Subset(ff, resCellsFltr))[, FSCChannel],
                na.rm = TRUE)
        },
        ff = ff, flt = resCellsFilter
    )

    debrisIndex <- which.min(FSCMedians)
    keptClustersIndexes <- setdiff(seq_len(nClust), debrisIndex)
    tokeepFilter <- resCellsFilter[[keptClustersIndexes[1]]]
    if (nClust > 2) {
        for (i in keptClustersIndexes[-1]) {
            tokeepFilter <- tokeepFilter | resCellsFilter[[i]]
        }
    }
    selectedCells <- flowCore::filter(ff, tokeepFilter)
    ff <- flowCore::Subset(ff, selectedCells)
    #ff <- ff[selectedCells@subSet, ]

    return(ff)
}

### FUNCTIONS FOR DEAD CELLS REMOVAL ###

#' @title remove dead cells from a flowFrame using manual gating
#' @description remove dead cells from a flowFrame, using manual gating in the
#' FSC-A, '(a)Live/Dead' 2D representation. The function uses
#' flowCore::polygonGate()
#' @param ff a flowCore::flowFrame
#' @param preTransform boolean, if TRUE: the transList list of scale transforms
#' will be applied first on the LD channel.
#' @param transList applied in conjunction with preTransform == TRUE
#' @param FSCChannel a character containing the exact name of the forward
#' scatter channel
#' @param LDMarker a character containing the exact name of the marker
#' corresponding to (a)Live/Dead channel, or the Live/Dead channel name itself
#' @param gateData a numerical vector containing the polygon gate coordinates
#' first the `FSCChannel` channel coordinates    
#' of each points of the polygon gate,   
#' then the LD channel coordinates of each points   
#' (prior to scale transform)
#' @param ... additional parameters passed to flowCore::polygonGate()
#'
#' @return a flowCore::flowFrame with removed dead cells from the input
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")    
#'                          
#' remDeadCellsGateData <- c(0, 0, 250000, 250000,
#'                           0, 650, 650, 0)  
#'
#' ff_lcells <-
#'     removeDeadCellsManualGate(ff_c,
#'                               FSCChannel = "FSC-A",
#'                               LDMarker = "L/D Aqua - Viability",
#'                               gateData = remDeadCellsGateData)
#'    
removeDeadCellsManualGate <- function(ff,
                                      preTransform = FALSE,
                                      transList = NULL,
                                      FSCChannel,
                                      LDMarker,
                                      gateData,
                                      ...) {
    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    if (preTransform) {
        if (is.null(transList)) {
            stop(
                "tranformation list needs to be provided ",
                "if preTransform = TRUE!"
            )
        }
        ffIn <- flowCore::transform(ff, transList)
    } else {
        ffIn <- ff
    }

    
    if (LDMarker %in% flowCore::colnames(ff)) {
        LDChannel <- LDMarker
    } else {
        LDChannel <- getChannelNamesFromMarkers(ffIn, markers = LDMarker)
    }

    liveGateMatrix <- matrix(
        data = gateData, ncol = 2,
        dimnames = list(c(), c(
            FSCChannel,
            LDChannel
        ))
    )

    liveGate <- flowCore::polygonGate(
        filterId = "Live_Cells",
        .gate = liveGateMatrix
    )



    selectedLive <- flowCore::filter(ffIn, liveGate)

    # note we take ff and not ffIn (no transfo)
    ff <- flowCore::Subset(ff, selectedLive)
    #ff <- ff[selectedLive@subSet, ] 
}

# @title remove dead cells from a flowFrame
# @description this function removes dead cells from a flowFrame, using a
# specific '(a)live/dead' channel, and the openCyto::gate_tail() gating
# function (see doc of the openCyto package)

# @param ff a flowCore::flowFrame
# @param preTransform if TRUE, apply the transList scale transform prior to
# running the gating algorithm
# @param transList applied in conjunction with preTransform == TRUE
# @param LDMarker a character containing the exact name of the marker
# corresponding to Live/Dead channel, or the Live/Dead channel name itself
# @param ... additional parameters passed to openCyto::gate_tail()
#
# @return a flowCore::flowFrame with removed dead cells from the input
# @export
# 
# @importFrom flowCore exprs
#
# @examples
#
# rawDataDir <-
#     paste0(system.file("extdata", package = "CytoPipeline"), "/")
# sampleFiles <-
#     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
# 
# truncateMaxRange <- FALSE
# minLimit <- NULL
# 
# # create flowCore::flowSet with all samples of a dataset
# fsRaw <- readSampleFiles(
#     sampleFiles = sampleFiles,
#     whichSamples = "all",
#     truncate_max_range = truncateMaxRange,
#     min.limit = minLimit)
# 
# suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#     
# ff_c <-
#     compensateFromMatrix(ff_m,
#                          matrixSource = "fcs")        
#
# transList <- 
#     estimateScaleTransforms(        
#         ff = ff_c,
#         fluoMethod = "estimateLogicle",
#         scatterMethod = "linear",
#         scatterRefMarker = "BV785 - CD3")
# 
# ff_lcells <-
#     removeDeadCellsGateTail(ff_c,
#                             preTransform = TRUE,
#                             transList = transList,
#                             LDMarker = "L/D Aqua - Viability",
#                             num_peaks = 2,
#                             ref_peak = 2,
#                             strict = FALSE,
#                             positive = FALSE)
#                             
# removeDeadCellsGateTail <- function(ff,
#                                     preTransform = FALSE,
#                                     transList = NULL,
#                                     LDMarker,
#                                     ...) {
#     
#     # handle ellipsis arguments, as 'openCyto::gate_tail'
#     # does not accept unknown args
#     passedEllipsisArgs <- list(...)
#     
#     newArgs <- list()
#     
#     argNames <-
#         c(
#             "num_peaks", "ref_peak", "strict", "tol", "side", "min", "max",
#             "bias", "positive", "deriv", "bandwidth", "adjust", "num_points",
#             "range.x", "binned", "se", "w"
#         )
#     for (argN in argNames) {
#         if (!is.null(passedEllipsisArgs[[argN]])) {
#             newArgs[[argN]] <- passedEllipsisArgs[[argN]]
#         }
#     }
#     
#     # will be needed by openCyto::gate_tail to support 
#     # BiocParallel::SnowParams
#     #requireNamespace("flowCore")
#     #require(flowCore, include.only = c("exprs", "rectangleGate"))
#     
#     message(paste0("Removing Dead Cells Gate Tail events from file : ", 
#                    flowCore::identifier(ff)))
#     
#     # if not present already, add a column with Cell ID
#     ff <- .appendCellID(ff)
#     
#     if (preTransform) {
#         if (is.null(transList)) {
#             stop(
#                 "tranformation list needs to be provided ",
#                 "if preTransform = TRUE!"
#             )
#         }
#         ffIn <- flowCore::transform(ff, transList)
#     } else {
#         ffIn <- ff
#     }
#     
#     if (LDMarker %in% flowCore::colnames(ff)) {
#         LDChannel <- LDMarker
#     } else {
#         LDChannel <- getChannelNamesFromMarkers(ffIn, markers = LDMarker)
#     }
#     
#     
#     
#     liveGate <-
#         do.call(openCyto::gate_tail,
#                 args = c(
#                     list(ffIn,
#                          channel = LDChannel,
#                          filterId = "Live_Cells"
#                     ),
#                     newArgs
#                 )
#         )
#     
#     selectedLive <- flowCore::filter(ffIn, liveGate)
#     
#     # note we take ff and not ffIn (no transfo)
#     ff <- flowCore::Subset(ff, selectedLive)
#     
#     return(ff)
#         
#     
#     
#     
#     
# 
# }

#' @title remove dead cells from a flowFrame
#' @description this function removes dead cells from a flowFrame, using a
#' specific '(a)live/dead' channel, and the flowDensity::deGate() gating
#' function (see doc of the flowDensity package)
#'
#' @param ff a flowCore::flowFrame
#' @param preTransform if TRUE, apply the transList scale transform prior to
#' running the gating algorithm
#' @param transList applied in conjunction with preTransform == TRUE
#' @param LDMarker a character containing the exact name of the marker
#' corresponding to Live/Dead channel, or the Live/Dead channel name itself
#' @param keepPositivePop logical flag stating whether we want to keep, 
#' after gating, the population that is positive for `LDMarker` (TRUE) 
#' or negative for `LDmarker` (FALSE)  
#' @param ... additional parameters passed to flowDensity::deGate()
#'
#' @return a flowCore::flowFrame with removed dead cells from the input
#' @export
#' 
#' @importFrom flowCore exprs
#'
#'@examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#'
#' transList <- 
#'     estimateScaleTransforms(        
#'         ff = ff_c,
#'         fluoMethod = "estimateLogicle",
#'         scatterMethod = "linear",
#'         scatterRefMarker = "BV785 - CD3")
#' 
#' ff_lcells <-
#'     removeDeadCellsDeGate(ff_c,
#'                           preTransform = TRUE,
#'                           transList = transList,
#'                           LDMarker = "L/D Aqua - Viability",
#'                           keepPositivePop = FALSE)
#'                             
removeDeadCellsDeGate <- function(ff,
                                  preTransform = FALSE,
                                  transList = NULL,
                                  LDMarker, 
                                  keepPositivePop = FALSE,
                                  ...) {

    
    message("Removing Dead Cells events (using flowDensity::deGate())", 
            " from file : ", flowCore::identifier(ff))

    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    if (preTransform) {
        if (is.null(transList)) {
            stop(
                "tranformation list needs to be provided ",
                "if preTransform = TRUE!"
            )
        }
        ffIn <- flowCore::transform(ff, transList)
    } else {
        ffIn <- ff
    }

    if (LDMarker %in% flowCore::colnames(ff)) {
        LDChannel <- LDMarker
    } else {
        LDChannel <- getChannelNamesFromMarkers(ffIn, markers = LDMarker)
    }
    
    threshold <- 
        flowDensity::deGate(ffIn,
                            channel = LDChannel,
                            ...)
    
    if (keepPositivePop) {
        selectedLive <- flowCore::exprs(ffIn)[, LDChannel] >= threshold    
    } else {
        selectedLive <- flowCore::exprs(ffIn)[, LDChannel] <= threshold
    }



    # note we take ff and not ffIn (no transfo)
    ff <- ff[selectedLive, ]

    return(ff)

}

### FUNCTIONS for Quality Control ###

#' @title perform QC with flowAI
#' @description this function is a wrapper around flowAI::flow_auto_qc()
#' function.
#' It also pre-selects the channels to be handled (=> all signal channels)
#' @param ff a flowCore::flowFrame
#' @param preTransform if TRUE, apply the transList scale transform prior to
#' running the gating algorithm
#' @param transList applied in conjunction with preTransform
#' @param outputDiagnostic if TRUE, stores diagnostic files generated by
#' flowAI in outputDir directory
#' @param outputDir used in conjunction with outputDiagnostic
#' @param ... additional parameters passed to flowAI::flow_auto_qc()
#'
#' @return a flowCore::flowFrame with removed low quality events from the input
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_QualityControl <- 
#'     qualityControlFlowAI(fsRaw[[2]],
#'                          remove_from = "all", # all default
#'                          second_fractionFR = 0.1,
#'                          deviationFR = "MAD",
#'                          alphaFR = 0.01,
#'                          decompFR = TRUE,
#'                          outlier_binsFS = FALSE,
#'                          pen_valueFS = 500,
#'                          max_cptFS = 3,
#'                          sideFM = "both",
#'                          neg_valuesFM = 1))
#' 
qualityControlFlowAI <- function(ff,
                                 preTransform = FALSE,
                                 transList = NULL,
                                 outputDiagnostic = FALSE,
                                 outputDir = NULL,
                                 ...) {
    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    if (preTransform) {
        if (is.null(transList)) {
            stop(
                "tranformation list needs to be provided ",
                "if preTransform = TRUE!"
            )
        }
        ffIn <- flowCore::transform(ff, transList)
    } else {
        ffIn <- ff
    }

    channel2Exclude <-
        flowCore::colnames(ffIn)[!areSignalCols(ffIn)]
    message("Applying flowAI method...")
    if (outputDiagnostic) {
        html_report <- "_QC"
        mini_report <- "QCmini"
        if (!is.null(outputDir)) {
            folder_results <- outputDir
        } else {
            folder_results <- "resultsQC"
        }
    } else {
        html_report <- FALSE
        mini_report <- FALSE
        folder_results <- FALSE
    }

    badEventIDs <-
        flowAI::flow_auto_qc(
            fcsfiles = ffIn,
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
            ...
        )

    goodEvents <- !(seq_len(flowCore::nrow(ffIn)) %in% badEventIDs)
    ff <- ff[goodEvents, ] # note we take ff and not ffIn (no transfo)

    return(ff)
}

#' @title perform QC with PeacoQC
#' @description this function is a wrapper around PeacoQC::PeacoQC()
#' function.
#' It also pre-selects the channels to be handled (=> all signal channels)
#' @param ff a flowCore::flowFrame
#' @param preTransform if TRUE, apply the transList scale transform prior to
#' running the gating algorithm
#' @param transList applied in conjunction with preTransform
#' @param outputDiagnostic if TRUE, stores diagnostic files generated by
#' PeacoQC in outputDir directory
#' @param outputDir used in conjunction with outputDiagnostic
#' @param ... additional parameters passed to PeacoQC::PeacoQC()
#'
#' @return a flowCore::flowFrame with removed low quality events from the input
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#' 
#' suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#'     
#' ff_c <-
#'     compensateFromMatrix(ff_m,
#'                          matrixSource = "fcs")        
#'
#' transList <- 
#'     estimateScaleTransforms(        
#'         ff = ff_c,
#'         fluoMethod = "estimateLogicle",
#'         scatterMethod = "linear",
#'         scatterRefMarker = "BV785 - CD3")
#'
#'
#' ff_QualityControl <- suppressWarnings(
#'     qualityControlPeacoQC(
#'         ff_c,
#'         preTransform = TRUE,
#'         transList = transList,
#'         min_cells = 150,
#'         max_bins = 500,
#'         MAD = 6,
#'         IT_limit = 0.55,
#'         force_IT = 150, 
#'         peak_removal = (1/3),
#'         min_nr_bins_peakdetection = 10))
#'         
qualityControlPeacoQC <- function(ff,
                                  preTransform = FALSE,
                                  transList = NULL,
                                  outputDiagnostic = FALSE,
                                  outputDir = NULL,
                                  ...) {
    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    if (preTransform) {
        if (is.null(transList)) {
            stop(
                "tranformation list needs to be provided ",
                "if preTransform = TRUE!"
            )
        }
        ffIn <- flowCore::transform(ff, transList)
    } else {
        ffIn <- ff
    }

    # qualityControl with PeacoQC
    message("Applying PeacoQC method...")
    channel4QualityControl <-
        flowCore::colnames(ffIn)[areSignalCols(ffIn)]
    if (outputDiagnostic) {
        plot <- TRUE
        report <- TRUE
        if (!is.null(outputDir)) {
            output_directory <- outputDir
        } else {
            output_directory <- "."
        }
    } else {
        plot <- FALSE
        report <- FALSE
        output_directory <- NULL # not used
    }

    res <- PeacoQC::PeacoQC(
        ff = ffIn,
        channels = channel4QualityControl,
        report = report,
        plot = plot,
        save_fcs = FALSE,
        output_directory = output_directory,
        ...
    )

    ff <- ff[res$GoodCells, ] # note we take ff and not ffIn (no transfo)
    return(ff)
}

#' @title perform QC with flowCut
#' @description this function is a wrapper around flowCut::flowCut()
#' function.
#' It also pre-selects the channels to be handled (=> all signal channels)
#' @param ff a flowCore::flowFrame
#' @param preTransform if TRUE, apply the transList scale transform prior to
#' running the gating algorithm
#' @param transList applied in conjunction with preTransform
#' @param outputDiagnostic if TRUE, stores diagnostic files generated by
#' flowCut in outputDir directory
#' @param outputDir used in conjunction with outputDiagnostic
#' @param verbose if TRUE messages comments on the QC process
#' @param ... additional parameters passed to flowCut::flowCut()
#'
#' @return a flowCore::flowFrame with removed low quality events from the input
#' @export
#'
#'
#' @examples
#'
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' fsRaw <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#'
#' ff_QualityControl <-
#'     qualityControlFlowCut(
#'         fsRaw[[2]],
#'         MaxContin = 0.1,
#'         MeanOfMeans = 0.13,
#'         MaxOfMeans = 0.15,
#'         MaxValleyHgt = 0.1,
#'         MaxPercCut = 0.3,
#'         LowDensityRemoval = 0.1,
#'         RemoveMultiSD = 7,
#'         AlwaysClean = FALSE,
#'         IgnoreMonotonic = FALSE,
#'         MonotonicFix = NULL,
#'         Measures = c(1:8))
#'       
qualityControlFlowCut <- function(ff,
                                  preTransform = FALSE,
                                  transList = NULL,
                                  outputDiagnostic = FALSE,
                                  outputDir = NULL,
                                  verbose = TRUE,
                                  ...) {
    # if not present already, add a column with Cell ID
    ff <- .appendCellID(ff)

    if (preTransform) {
        if (is.null(transList)) {
            stop(
                "tranformation list needs to be provided ",
                "if preTransform = TRUE!"
            )
        }
        ffIn <- flowCore::transform(ff, transList)
    } else {
        ffIn <- ff
    }

    # qualityControl with flowCut
    message("Applying flowCut method...")
    channelsIndices <- which(CytoPipeline::areSignalCols(ffIn))
    if (outputDiagnostic) {
        Plot <- "All"
        if (!is.null(outputDir)) {
            Directory <- outputDir
        } else {
            filePrefixWithDir <- "resultsQC"
        }
    } else {
        Directory <- NULL # not used
        Plot <- "None"
    }

    res <-
        flowCut::flowCut(
            f = ffIn,
            Channels = channelsIndices,
            Directory = Directory,
            FileID = NULL,
            Plot = Plot,
            Verbose = verbose,
            ...
        )
    # browser()
    badEventIDs <- res$ind

    goodEvents <- !(seq_len(flowCore::nrow(ffIn)) %in% badEventIDs)
    ff <- ff[goodEvents, ] # note we take ff and not ffIn (no transfo)

    return(ff)
}

# @title perform QC with flowClean
# @description this function is a wrapper around flowClean::clean()
#function.
#It also pre-selects the channels to be handled (=> all signal channels)
# @param ff a flowCore::flowFrame
# @param preTransform if TRUE, apply the transList scale transform prior to
# running the gating algorithm
# @param transList applied in conjunction with preTransform
# @param outputDiagnostic if TRUE, stores diagnostic files generated by
# flowClean in outputDir directory
# @param outputDir used in conjunction with outputDiagnostic
# @param verbose if TRUE messages comments on the QC process
# @param ... additional parameters passed to flowClean::clean()
# 
# @return a flowCore::flowFrame with removed low quality events from the input
# @export
# 
# @examples
# 
# rawDataDir <-
#     paste0(system.file("extdata", package = "CytoPipeline"), "/")
# sampleFiles <-
#     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
# 
# truncateMaxRange <- FALSE
# minLimit <- NULL
# 
# # create flowCore::flowSet with all samples of a dataset
# fsRaw <- readSampleFiles(
#     sampleFiles = sampleFiles,
#     whichSamples = "all",
#     truncate_max_range = truncateMaxRange,
#     min.limit = minLimit)
# 
# ff_QualityControl <- suppressWarnings(
#     qualityControlFlowClean(fsRaw[[2]],
#                             binSize = 0.01, # default
#                             nCellCutoff = 500, # default
#                             cutoff = "median", # default
#                             fcMax = 1.3, # default
#                             nstable = 5))
# 
# qualityControlFlowClean <- function(ff,
#                                     preTransform = FALSE,
#                                     transList = NULL,
#                                     outputDiagnostic = FALSE,
#                                     outputDir = NULL,
#                                     verbose = TRUE,
#                                     ...) {
#     #browser()
#     # if not present already, add a column with Cell ID
#     ff <- .appendCellID(ff)
# 
#     if (preTransform) {
#         if (is.null(transList)) {
#             stop(
#                 "tranformation list needs to be provided ",
#                 "if preTransform = TRUE!"
#             )
#         }
#         ffIn <- flowCore::transform(ff, transList)
#     } else {
#         ffIn <- ff
#     }
# 
#     # qualityControl with flowClean
#     message("Applying flowClean method...")
#     vectMarkers <- which(CytoPipeline::areSignalCols(ffIn))
# 
#     if (outputDiagnostic) {
#         diagnostic <- TRUE
#         if (!is.null(outputDir)) {
#             filePrefixWithDir <- outputDir
#         } else {
#             filePrefixWithDir <- "resultsQC"
#         }
# 
#         # add original fcs file name in prefix,
#         # as flowClean is designed to work for one fcs at the time
#         filename <- basename(flowCore::keyword(ffIn, "FILENAME")$FILENAME)
#         # removing extension
#         filename <- sub("([^.]+)\\.[[:alnum:]]+$", "\\1", filename)
#         filePrefixWithDir <- paste0(filePrefixWithDir, filename)
#     } else {
#         filePrefixWithDir <- NULL # not used
#         diagnostic <- FALSE
#     }
# 
#     # note we call it using returnVector = FALSE and extract 
#     # the goodVsBadVector later to work around a flowClean bug (to be corrected)
#     tempDF <-
#         flowClean::clean(
#             fF = ffIn,
#             vectMarkers = vectMarkers,
#             filePrefixWithDir = filePrefixWithDir,
#             ext = ".fcs", # not used
#             diagnostic = diagnostic,
#             announce = verbose,
#             returnVector = FALSE,
#             ...
#         )
# 
#     areGoodEvents <- tempDF[, "GoodVsBad"] < 10000
#     ff <- ffIn[areGoodEvents, ]
# 
#     return(ff)
# }

##' @title read RDS object 
##' @description wrapper around readRDS, which discards any additional 
##' parameters passed in (...)
##' @param RDSFile a RDS file containing a R object
##' object
##' @param ... other arguments (not used)
##' @return the read R object
##' @export
##' 
##' @examples
##' 
##' transListPath <- paste0(system.file("extdata", 
##'                                     package = "CytoPipeline"),
##'                         "/OMIP021_TransList.rds") 
##' 
##' transList <- readRDSObject(transListPath)
##' 
##' ff_c <- compensateFromMatrix(OMIP021Samples[[1]],
##'                              matrixSource = "fcs")  
##' 
##' ff_t <- applyScaleTransforms(ff_c, transList = transList)
##' 
readRDSObject <- function(RDSFile, ...) {
    if (!file.exists(RDSFile)) {
        stop("file ", RDSFile, " does not exist")
    }
    obj <- readRDS(file = RDSFile)
}

##' @title apply scale transforms
##' @description wrapper around flowCore::transform() that discards any 
##' additional parameter passed in (...)
##' @param ff a flowCore::flowFrame
##' @param transList a flowCore::transformList
##' @param ... other arguments (not used)
##' @return the transformed flowFrame
##' @export
##' 
##' @examples
##' 
##' transListPath <- paste0(system.file("extdata", 
##'                                     package = "CytoPipeline"),
##'                         "/OMIP021_TransList.rds") 
##' 
##' transList <- readRDSObject(transListPath)
##' 
##' ff_c <- compensateFromMatrix(OMIP021Samples[[1]],
##'                              matrixSource = "fcs")  
##' 
##' ff_t <- applyScaleTransforms(ff_c, transList = transList)
##' 
applyScaleTransforms <- function(ff, transList, ...) {
    ff <- flowCore::transform(ff, transList)
    return(ff)
}

#' @title write flowFrame to disk
#' @description wrapper around flowCore::write.FCS() or utils::write.csv
#' that discards any additional parameter passed in (...)
#' @param ff a flowCore::flowFrame
#' @param dir an existing directory to store the flowFrame, 
#' @param useFCSIdentifier if TRUE filename used will be based on flowFrame
#' identifier (using flowCOre::identifier(ff))
#' @param prefix file name prefix
#' @param suffix file name suffix
#' @param format either fcs or csv
#' @param ... other arguments (not used)
#' @return nothing
#' @export
#' 
#' @examples
#' 
#' transListPath <- paste0(system.file("extdata", 
#'                                     package = "CytoPipeline"),
#'                         "/OMIP021_TransList.rds") 
#' 
#' transList <- readRDSObject(transListPath)
#' 
#' ff_c <- compensateFromMatrix(OMIP021Samples[[1]],
#'                              matrixSource = "fcs")  
#' 
#' ff_t <- applyScaleTransforms(ff_c, transList = transList)
#' 
#' 
#' rawDataDir <-
#'     paste0(system.file("extdata", package = "CytoPipeline"), "/")
#' sampleFiles <-
#'     paste0(rawDataDir, list.files(rawDataDir, pattern = "sample_"))
#' 
#' truncateMaxRange <- FALSE
#' minLimit <- NULL
#' 
#' # create flowCore::flowSet with all samples of a dataset
#' res <- readSampleFiles(
#'     sampleFiles = sampleFiles,
#'     whichSamples = "all",
#'     truncate_max_range = truncateMaxRange,
#'     min.limit = minLimit)
#'     
#' ff_c <- compensateFromMatrix(res[[2]], matrixSource = "fcs") 
#' outputDir <- withr::local_tempdir()
#' writeFlowFrame(ff_c, 
#'                dir = outputDir,
#'                suffix = "_fcs_export",
#'                format = "csv")
#' 
writeFlowFrame <- function(ff, dir, 
                           useFCSIdentifier = TRUE,
                           prefix = "", suffix ="",
                           format = c("fcs", "csv"), ...) {
    format <- match.arg(format)
    if (!useFCSIdentifier && prefix == "" && suffix == "") {
        stop ("No file name provided! If useFCSIdentifier == FALSE, ",
              "then either prefix or suffix should be provided.")
    }
    if(!dir.exists(dir)) {
        stop ("Provided directory does not exist!")
    }
    
    fileName <- paste0(dir, "/", prefix)
    if (useFCSIdentifier) {
        fileName <- paste0(fileName, flowCore::identifier(ff))
    }
    fileName <- paste0(fileName, suffix, ".", format)
    
    if (format == "fcs") {
        flowCore::write.FCS(ff, filename = fileName)
    } else {
        utils::write.csv(flowCore::exprs(ff), file = fileName, row.names = FALSE)
    }
}

