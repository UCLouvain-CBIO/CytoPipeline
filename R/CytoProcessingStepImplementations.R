# CytoPipeline - Copyright (C) <2022-2024>
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
#' @param verbose if TRUE, send messages to the user at each step
#'
#' @return a flowCore::flowFrame with removed low quality events from the input
#' @export
#'
#' @examples
#' 
#' data(OMIP021Samples)
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
                                    scatterMethod = c("none",
                                                      "linearQuantile"),
                                    scatterRefMarker = NULL,
                                    specificScatterChannels = NULL,
                                    verbose = FALSE){
    fluoMethod <- match.arg(fluoMethod)
    scatterMethod <- match.arg(scatterMethod)
    
    if (fluoMethod == "estimateLogicle") {
        if (verbose) {
            message(
                "estimating logicle transformations ",
                "for fluorochrome channels..."
            ) 
        }
        fluoCols <- flowCore::colnames(ff)[areFluoCols(ff)]
        transList <- flowCore::estimateLogicle(ff, fluoCols)
    } # else do nothing
    
    if (scatterMethod == "linearQuantile") {
        if (is.null(scatterRefMarker)) {
            stop("linear scatter method requires a scatterRefMarker")
        }
        
        if (verbose) {
            message(
                "Estimating linear transformation for scatter channels : ",
                "reference marker = ",
                scatterRefMarker,
                "..."
            )
        }
        
        transList <-
            computeScatterChannelsLinearScale(
                ff,
                transList = transList,
                referenceChannel = scatterRefMarker,
                silent = !verbose
            )
    } # else do nothing
    
    # handle specific cases of scatter channels that still need fluo method
    if (!is.null(specificScatterChannels) && fluoMethod == "estimateLogicle") {
        effectiveScatterChannels <- NULL
        for (ch in specificScatterChannels) {
            scatterChannels <-
                flowCore::colnames(ff)[!areFluoCols(ff) & areSignalCols(ff)]
            if (!(ch %in% scatterChannels)) {
                if (verbose) {
                    message("Specific channel [", ch, "] is not a scatter channel",
                            " => no correction of scale transformation done")
                }
            } else {
                transList@transforms[[ch]] <- NULL
                effectiveScatterChannels <- c(effectiveScatterChannels, ch)
                if (verbose) {
                    message("Specific scatter channel found: [", ch, "] ",
                            "=> correcting scale transformation to logicle...")
                }
            } 
        }
        
        transList <- c(transList, 
                       flowCore::estimateLogicle(ff, 
                                                 effectiveScatterChannels))
    }
    
    return(transList)
}


#' @title Read fcs sample files
#' @description Wrapper around flowCore::read.fcs() or flowCore::read.flowSet().
#' Also adds a "Cell_ID" additional column, used in flowFrames comparison
#' @param sampleFiles a vector of character path to sample files
#' @param whichSamples one of:
#' - 'all' if all sample files need to be read
#' - 'random' if some samples need to be chosen randomly 
#' (in that case, using `nSamples` and `seed`)
#' - a vector of indexes pointing to the sampleFiles vector
#' @param nSamples number of samples to randomly select 
#' (if `whichSamples == "random"`). 
#' If `nSamples` is higher than nb of available samples, 
#' the output will be all samples
#' @param seed an optional seed parameters (provided to ease reproducibility).
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
#'     system.file("extdata", package = "CytoPipeline")
#' sampleFiles <-
#'     file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
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
                            nSamples = NULL,
                            seed = NULL,
                            channelMarkerFile = NULL,
                            ...) {
    
    if (whichSamples == "all") {
        # do nothing : sampleFiles should contain all the input sample files
        # already
    } else if (whichSamples == "random") {
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
        
        whichSamples <- sample(seq_along(sampleFiles), nSamples)
        sampleFiles <- sampleFiles[whichSamples]
    } else if (is.numeric(whichSamples)) {
        sampleFiles <- sampleFiles[whichSamples]
    } else {
        stop("'whichSamples' should be either 'all', or a vector of indexes")
    }
    

    if (length(sampleFiles) == 0) {
        stop("no sample files to read")
    } else if (length(sampleFiles) == 1) {
        res <- flowCore::read.FCS(sampleFiles, ...)
        
        # Add a column with Cell ID
        res <- appendCellID(res)
    } else {
        res <- flowCore::read.flowSet(sampleFiles,
                                      ...)
        
        # Add a column with Cell ID
        res <- flowCore::fsApply(
            x = res,
            FUN = function(ff) {
                appendCellID(ff)
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
        
        message("COL NAMES MARKER MAPPING: ", 
                paste(colnames(channelMarkerMapping), collapse = ", "))
        
        if (!("Channel" %in% colnames(channelMarkerMapping))) {
            stop("channelMarkersMapping should contain [Channel] column!")
        }
        if (!("Marker" %in% colnames(channelMarkerMapping))) {
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



#' @title remove margin events using PeacoQC
#' @description Wrapper around PeacoQC::RemoveMargins().
#' Also pre-selects the channels to be handled (=> all signal channels)
#' If input is a flowSet, it applies removeMargins() to each flowFrame of the
#' flowSet.
#' @param x a flowCore::flowSet or a flowCore::flowFrame
#' @param channelSpecifications A list of lists with parameter specifications 
#' for certain channels. This parameter should only be used if the values in 
#' the internal parameters description is too strict or wrong for a number or 
#' all channels. This should be one list per channel with first a minRange 
#' and then a maxRange value. This list should have the channel name found back 
#' in colnames(flowCore::exprs(ff)), or the corresponding marker name (found in
#' flowCore::pData(flowCore::description(ff)) ) . 
#' If a channel is not listed in this parameter, its default internal values 
#' will be used. The default of this parameter is NULL.
#' If the name of one list is set to `AllFluoChannels`, then the  
#' `minRange` and `maxRange` specified there will be taken as default   
#' for all fluorescent channels (not scatter)
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
#'     system.file("extdata", package = "CytoPipeline")
#' sampleFiles <- 
#'     file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
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
removeMarginsPeacoQC <- function(x, channelSpecifications = NULL, ...) {
    PQCChannelSpecs <- channelSpecifications
    
    if (!is.null(channelSpecifications)) {
        if (!methods::is(channelSpecifications, "list")) 
            stop("channelSpecifications should be a list of lists.")
        if (!all(lengths(channelSpecifications) == 2)) 
            stop("channel_specifications should be a list of lists. \n",
                 "Every list should have the channel name and should contain\n",
                 "a minRange and maxRange value.")
    }
    
    
    myFunc <- function(ff, channelSpecifications) {
        message("Removing margins from file : ", getFCSFileName(ff))
        channel4Margins <-
            flowCore::colnames(ff)[areSignalCols(ff)]
        
        markers4Margins <- 
            flowCore::pData(flowCore::parameters(ff))$desc[areSignalCols(ff)]
        
        PQCChannelSpecs <- channelSpecifications
        
        
        if (!is.null(PQCChannelSpecs)) {
            
            #store default fluo parameters (if any)
            defaultFluoChannelList <- PQCChannelSpecs[["AllFluoChannels"]]
            if (!is.null(defaultFluoChannelList)) {
                PQCChannelSpecs[["AllFluoChannels"]] <- NULL
            }
            
            newNames <- names(PQCChannelSpecs)
            for (l in seq_along(newNames)) {
                chName <- newNames[l]
                if (!(chName %in% flowCore::colnames(ff))) {
                    # try as marker
                    whichMarker <- which(markers4Margins == chName)
                    if (length(whichMarker) == 0)
                        stop("channelSpecifications names: ",
                             "could not find [", chName, "], neither as ",
                             "channel, nor as marker")
                    else {
                        newNames[l] <- channel4Margins[whichMarker[1]]
                    }
                }
            }
            names(PQCChannelSpecs) <- newNames
            
            # apply default fluo parameters, if any
            if (!is.null(defaultFluoChannelList)) {
                for (ch in flowCore::colnames(ff)
                     [CytoPipeline::areFluoCols(ff)]) {
                    if (!(ch %in% newNames)) {
                        PQCChannelSpecs[[ch]] <- defaultFluoChannelList
                    }
                }
            }
        }
        
        ffOut <- PeacoQC::RemoveMargins(ff, channels = channel4Margins,
                                        channel_specifications = 
                                            PQCChannelSpecs)
        return(ffOut)
    }
    if (inherits(x, "flowFrame")) {
        return(myFunc(x, channelSpecifications = channelSpecifications))
    } else if (inherits(x, "flowSet")) {
        fsOut <- flowCore::fsApply(x, FUN = myFunc, simplify = TRUE,
                                   channelSpecifications = 
                                       channelSpecifications)
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
#'     system.file("extdata", package = "CytoPipeline")
#' sampleFiles <-
#'     file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
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
        fileId <- getFCSFileName(ff)
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
#' @param x a `flowCore::flowFrame` or `flowCore::flowSet`
#' @param matrixSource if "fcs", the compensation matrix will be fetched from
#' the fcs files (different compensation matrices can then be applied by fcs
#' file)
#' if "import", uses `matrixPath` to read the matrix (should be a csv file)
#' @param matrixPath if matrixSource == "import", will be used as the input csv
#' file path
#' @param updateChannelNames if TRUE, updates the fluo channel names by
#' prefixing them with "comp-"
#' @param verbose if TRUE, displays information messages
#' @param ... additional arguments (not used)
#' @return the compensated flowSet or flowFrame
#' @export
#'
#' @examples
#'
#' rawDataDir <-
#'     system.file("extdata", package = "CytoPipeline")
#' sampleFiles <-
#'     file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
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
                                 verbose = FALSE,
                                 ...) {

    #browser()
    matrixSource <- match.arg(matrixSource)
    
    if (matrixSource == "import"){
        if (is.null(matrixPath)) {
            stop("matrixPath can't be NULL if matrixSource == 'import'")
        }
    }
    
    importCompensationMatrix <- function(ff, matrixPath) {
        # import matrix from file (path)
        if (is.null(matrixPath)) {
            stop("No path specified for compensation matrix!")
        }
        if (!file.exists(matrixPath)) {
            stop("Compensation matrix file [", matrixPath,
                 "] not found!")
        }
        compensationMatrix <-
            utils::read.csv(matrixPath,
                            check.names = FALSE,
                            row.names = 1
            )
        compensationMatrix <- 
            .updateCompMatrixLabels(compensationMatrix, ff)
        
        return(compensationMatrix)
    }
    
    compensateOneFFWithSource <- function(
        ff, matrixSource, matrixPath, verbose){
        
        if (verbose) {
            message("Compensating file : ", 
                    getFCSFileName(ff),
                    "; matrixSource = ",
                    matrixSource,
                    "; matrixPath = ",
                    matrixPath)
        }
        
        #browser()
        if (matrixSource == "fcs") {
            # obtains compensation matrix
            compensationMatrix <-
                getAcquiredCompensationMatrix(ff)
        } else {
            # find correct matrix path
            compensationMatrix <- importCompensationMatrix(ff, matrixPath)
        }
        
        ffOut <- runCompensation(ff,
                                 compensationMatrix,
                                 updateChannelNames = updateChannelNames
        )
        
        return(ffOut)
    }
            
    if (inherits(x, "flowFrame")) {
        res <- compensateOneFFWithSource(
            ff = x, 
            matrixSource = matrixSource, 
            matrixPath = matrixPath,
            verbose = verbose)
    } else if (inherits(x, "flowSet")) {
        if (matrixSource == "fcs" || 
            matrixSource == "import" && length(matrixPath) <= 1) {
            res <- flowCore::fsApply(x,
                                     FUN = compensateOneFFWithSource,
                                     simplify = TRUE,
                                     matrixSource = matrixSource,
                                     matrixPath = matrixPath,
                                     verbose = verbose)
        } else {
            # matrix path is different from flowFrame to flowFrame
            # => need to use mapply() instead of fsApply()
            res <- structure(
                mapply(x, 
                       matrixPath,
                       FUN = function(ff, matrixPath, matrixSource, verbose) {
                           ff <- compensateOneFFWithSource(
                               ff = ff,
                               matrixSource = matrixSource,
                               matrixPath = matrixPath,
                               verbose = verbose
                           )
                           ff
                       },
                       MoreArgs = list(matrixSource = matrixSource,
                                       verbose = verbose)), 
                names = flowCore::sampleNames(x))
            res <- methods::as(res, "flowSet")
            flowCore::phenoData(res) <- 
                flowCore::phenoData(x)[flowCore::sampleNames(x), , drop = FALSE]
        }
        
    } else {
        stop("x should be a flowCore::flowFrame or a flowCore::flowSet")
    }
    
    res

}



### FUNCTIONS FOR DOUBLETS REMOVAL ###



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
#'     system.file("extdata", package = "CytoPipeline")
#' sampleFiles <-
#'     file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
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
    ff <- appendCellID(ff)

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
#'     system.file("extdata", package = "CytoPipeline")
#' sampleFiles <-
#'     file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
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
    ff <- appendCellID(ff)

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
#'     system.file("extdata", package = "CytoPipeline")
#' sampleFiles <-
#'     file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
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
    ff <- appendCellID(ff)

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
#'     system.file("extdata", package = "CytoPipeline")
#' sampleFiles <-
#'     file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
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
    ff <- appendCellID(ff)

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

    goodEvents <- !(seq_len(flowCore::nrow(ffIn)) %in% badEventIDs[[1]])
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
#'     system.file("extdata", package = "CytoPipeline")
#' sampleFiles <-
#'     file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))
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
    ff <- appendCellID(ff)

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
##' data(OMIP021Samples)
##' 
##' transListPath <- file.path(system.file("extdata", 
##'                                         package = "CytoPipeline"),
##'                            "OMIP021_TransList.rds") 
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
##' Additionally, some checks regarding channels correspondance is done:  
##' if `transList` contains transformations for channels that are not present 
##' in `x`, then these transformations are first removed.
##' @param x a flowCore::flowSet or a flowCore::flowFrame
##' @param transList a flowCore::transformList
##' @param verbose if TRUE, send a message per flowFrame transformed
##' @param ... other arguments (not used)
##' @return the transformed flowFrame
##' @export
##' 
##' @examples
##' 
##' data(OMIP021Samples)
##' 
##' transListPath <- file.path(system.file("extdata", 
##'                                        package = "CytoPipeline"),
##'                            "OMIP021_TransList.rds") 
##' 
##' transList <- readRDSObject(transListPath)
##' 
##' ff_c <- compensateFromMatrix(OMIP021Samples[[1]],
##'                              matrixSource = "fcs")  
##' 
##' ff_t <- applyScaleTransforms(ff_c, transList = transList)
##' 
applyScaleTransforms <- function(x, transList, verbose = FALSE, ...) {
    
    if (!inherits(transList, "transformList")) {
        stop("transList should inherit from transformList class")
    }
    
    if (!inherits(x, "flowFrame") && !inherits(x, "flowSet")) {
        stop("x type not recognized, should be a flowFrame or a flowSet")
    }
    
    # check on channels/markers present in transList
    signalChannels <- flowCore::colnames(x)[areSignalCols(x)]
    transfoChannels <- names(transList@transforms)
    commonChannels <- intersect(signalChannels, transfoChannels)
    
    transList@transforms <- transList@transforms[commonChannels]
    
    if (inherits(x, "flowFrame")) {
        if (verbose) {
            message("scale transforming ff: ", 
                    flowCore::identifier(x))
        }
        ret <- flowCore::transform(x, transList)
    } else if (inherits(x, "flowSet")) {
        ret <- flowCore::fsApply(
            x,
            FUN = function(ff){
                if (verbose) {
                    message("scale transforming ff: ", 
                            flowCore::identifier(ff))
                }
                flowCore::transform(ff, transList)
            })
    } else {
        stop("x should be a flowCore::flowFrame or a flowCore::flowSet")
    }
    
    
    
    return(ret)
}



