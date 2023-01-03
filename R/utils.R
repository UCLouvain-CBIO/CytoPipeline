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


#' @title find flow frame columns that represent true signal
#' @description : find flow frame columns that represent true signal
#' @param ff a flowCore::flowFrame
#' @param toRemovePatterns a vector of string patterns that are to be considered
#' as non signal
#'
#' @return a vector of booleans of which the dimension is equal to the number of
#' columns in ff
#' @export
#' 
#' @examples
#' areSignalCols(OMIP021Samples[[1]])
#'
areSignalCols <- function(ff,
                          toRemovePatterns = c(
                              "Time", "Original_ID",
                              "File", "SampleID"
                          )) {
    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame")
    }


    retCols <- vapply(flowCore::colnames(ff),
        FUN.VALUE = logical(1),
        FUN = function(ch, toRemovePatterns) {
            res <- TRUE
            for (pat in toRemovePatterns) {
                res <- res & !grepl(pat, ch, ignore.case = TRUE)
            }
            res
        },
        toRemovePatterns = toRemovePatterns
    )
    return(retCols)
}
#' @title find flow frame columns that represent fluorochrome channel
#' @description : find flow frame columns that represent fluorochrome channel
#' @param ff a flowCore::flowFrame
#' @param toRemovePatterns a vector of string patterns that are to be considered
#' as non fluorochrome
#'
#' @return a vector of booleans of which the dimension is equal to the number of
#' columns in ff
#' @export
#' 
#' @examples
#' areFluoCols(OMIP021Samples[[1]])
#'
areFluoCols <- function(ff,
                        toRemovePatterns = c(
                            "FSC", "SSC",
                            "Time", "Original_ID",
                            "File", "SampleID"
                        )) {
    return(areSignalCols(ff,
        toRemovePatterns = toRemovePatterns
    ))
}



#' @title sub-sampling of a flowFrame
#' @description : sub-samples a flowFrame
#' with the specified number of samples, without replacement.
#' adds also a column 'Original_ID' if not already present in flowFrame.
#' @param ff a flowCore::flowFrame
#' @param nSamples number of samples to be obtained using sub-sampling
#' @param seed  can be set for reproducibility of event sub-sampling
#'
#' @return new flowCore::flowFrame with the obtained subset of samples
#' @export
#' 
#' @examples
#' # take first sample of dataset, subsample 100 events and create new flowFrame
#' ff <- subsample(OMIP021Samples[[1]], nSamples = 100)
#' 
#'
subsample <- function(ff, nSamples, seed = NULL) {
    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame")
    }

    eventCounts <- length(flowCore::exprs(ff)[, 1])
    nSamples <- min(eventCounts, nSamples)

    if (!is.null(seed)) {
        withr::with_seed(
            seed,
            keep <- sample(seq_len(eventCounts),
                size = nSamples,
                replace = FALSE
            )
        )
    } else {
        keep <- sample(seq_len(eventCounts),
            size = nSamples,
            replace = FALSE
        )
    }

    # add Original_ID as a new column if necessary
    ff <- appendCellID(ff)

    ff[keep, ]
}

#' @title compensate with additional options
#' @description : this is a simple wrapper around the flowCore::compensate()
#' utility, allowing to trigger an update of the fluo channel names
#' with a prefix 'comp-' (as in FlowJo)
#'
#' @param obj a flowCore::flowFrame or flowCore::flowSet
#' @param spillover compensation object or spillover matrix or a list of
#' compensation objects
#' @param updateChannelNames if TRUE, add a 'comp-' prefix to all fluorochrome
#' channels (hence does not impact the columns related to FSC, SSC, or other
#' specific keyword like TIME, Original_ID, File,...)
#' Default TRUE.
#'
#' @return a new object with compensated data, and possibly updated column names
#' @export
#' 
#' @examples 
#' ff <- OMIP021Samples[[1]]
#' compMatrix <- flowCore::spillover(ff)$SPILL
#' ff <- runCompensation(ff, 
#'                       spillover = compMatrix, 
#'                       updateChannelNames = TRUE)
#'
runCompensation <- function(obj, spillover, updateChannelNames = TRUE) {
    isFlowSet <- FALSE
    if (inherits(obj, "flowSet")) {
        isFlowSet <- TRUE
    } else if (!inherits(obj, "flowFrame")) {
        stop("obj type not recognized, should be a flowFrame or flowSet")
    }

    res <- flowCore::compensate(obj, spillover)
    if (updateChannelNames) {
        if (isFlowSet) {
            res <-
                flowCore::fsApply(res, FUN = .addCompensation2FluoChannelNames)
        } else {
            res <- .addCompensation2FluoChannelNames(res)
        }
    }
    return(res)
}

#' @title Aggregate and sample multiple flow frames of a flow set together
#' @description Aggregate multiple flow frames in order to analyze them
#' simultaneously.
#' A new FF, which contains about cTotal cells, with ceiling(cTotal/nFiles)
#' cells from each file.
#' Two new columns are added:
#' a column indicating the original file by index, and
#' a noisy version of this, for better plotting opportunities,
#' This function is based on PeacoQC::AggregateFlowframes() where file names
#' inputs have been replaced by a flowSet input.
#'
#' @param fs a flowCore::flowset
#' @param nTotalEvents Total number of cells to select from the input flow
#' frames
#' @param seed seed to be set before sampling for reproducibility.
#' Default NULL does not set any seed.
#' @param channels Channels/markers to keep in the aggregate.
#' Default NULL takes all channels of the first file.
#' @param writeOutput Whether to write the resulting flowframe to a file.
#' Default FALSE
#' @param outputFile Full path to output file. Default "aggregate.fcs"
#' @param keepOrder If TRUE, the random subsample will be ordered in the same
#' way as they were originally ordered in the file. Default = FALSE.
#'
#' @return returns a new flowCore::flowFrame
#' @export
#' 
#' @examples 
#' nCells <- 1000
#' agg <- aggregateAndSample(
#'     fs = OMIP021Samples,
#'     nTotalEvents = nCells)
aggregateAndSample <- function(fs,
                               nTotalEvents,
                               seed = NULL,
                               channels = NULL,
                               writeOutput = FALSE,
                               outputFile = "aggregate.fcs",
                               keepOrder = FALSE) {
    # browser()
    if (inherits(fs, "flowFrame")) {
        fs <- flowCore::flowSet(fs)
    }
    if (!inherits(fs, "flowSet")) {
        stop("fs object type not recognized, should be flowCore::flowSet,",
             " or flowCore::flowFrame")
    }


    if (!is.null(seed)) {
        # set the seed locally in the execution environment,
        # restore it afterward
        withr::local_seed(seed)
    }

    nFrames <- length(fs)
    cFrame <- ceiling(nTotalEvents / nFrames)
    flowFrame <- NULL
    diffNumberChannels <- FALSE
    diffMarkers <- FALSE
    for (i in seq_len(nFrames)) {
        current_ff <- fs[[i]]
        ids <- sample(seq_len(nrow(current_ff)), min(nrow(current_ff), cFrame))
        if (keepOrder) {
            ids <- sort(ids)
        }
        new_col_names <- c("File", "File_scattered", "Original_ID")
        prev_agg <- length(grep("File[0-9]*$", colnames(current_ff)))
        if (prev_agg > 0) {
            new_col_names[c(1, 2)] <- paste0(new_col_names[c(1, 2)], prev_agg +
                1)
        }
        prev_ids <- length(grep(
            "Original_ID[0-9]*$",
            flowCore::colnames(current_ff)
        ))
        if (prev_ids > 0) {
            new_col_names[3] <- paste0(new_col_names[3], prev_ids + 1)
        }
        file_ids <- rep(i, min(nrow(current_ff), cFrame))
        m <- cbind(file_ids, file_ids + stats::rnorm(
            length(file_ids),
            0, 0.1
        ), ids)
        colnames(m) <- new_col_names
        current_ff <- flowCore::fr_append_cols(current_ff[ids, ], m)
        if (is.null(flowFrame)) {
            if (is.null(channels)) {
                channels <- flowCore::colnames(current_ff)
                flowFrame <- current_ff
            } else {
                channels <- getChannelNamesFromMarkers(current_ff, channels)
                flowFrame <- current_ff[, c(channels, colnames(m)),
                    drop = FALSE
                ]
            }
            flowCore::keyword(flowFrame)[["$FIL"]] <- basename(outputFile)
            flowCore::keyword(flowFrame)[["FILENAME"]] <- basename(outputFile)
            flowCore::identifier(flowFrame) <- basename(outputFile)
        } else {
            cols_f <- flowCore::colnames(current_ff)
            cols_flowFrame <- flowCore::colnames(flowFrame)
            commonCols <- intersect(cols_f, cols_flowFrame)
            if (length(commonCols) == 0) {
                stop("No common channels between flow frames")
            }
            if (!diffNumberChannels && length(cols_flowFrame) !=
                length(commonCols)) {
                diffNumberChannels <- TRUE
            }
            if (!diffMarkers &&
                any(!flowCore::markernames(current_ff)[commonCols] %in%
                    flowCore::markernames(flowFrame)[commonCols])) {
                diffMarkers <- TRUE
            }
            flowCore::exprs(flowFrame) <-
                rbind(
                    flowCore::exprs(flowFrame)[, commonCols, drop = FALSE],
                    flowCore::exprs(current_ff)[, commonCols, drop = FALSE]
                )
        }
    }

    if (diffNumberChannels) {
        warning(
            "Flow frames do not contain the same number of channels/markers"
        )
    }
    if (diffMarkers) {
        warning("Flow frames do not contain the same markers")
    }
    if (writeOutput) {
        flowCore::write.FCS(flowFrame, filename = outputFile)
    }
    return(flowFrame)
}

#' @title get tranformation parameters for a specific channel
#' @description investigates a flowCore::tranformList object to get the type
#' and parameters of the transformation applying to a specific channel
#' @param transList a flowCore::transformList
#' @param channel channel name
#'
#' @return If the transformation exists for the specified channel, and is either
#' recognized as a logicle transfo or a linear transfo, a list with two slots:
#' - $type a character containing the transfo type ('logicle' or 'linear')
#' - $params_list a list of named numeric, according to transfo type
#'
#' Otherwise, NULL is returned.
#' @export
#' 
#' @examples

#' # set-up a hybrid transformation list :
#' # - two channels are logicle-ly transformed with automatic param estimates
#' # - one channel has explicit logicle transfo with default parameters
#' # - one channel has linear transformation
#' # - other channels have no transformation
#' translist <- flowCore::estimateLogicle(
#'     OMIP021Samples[[1]],
#'     c("450/50Violet-A", "525/50Violet-A")
#' )
#' translist <- c(
#'     translist,
#'     flowCore::transformList(
#'         "FSC-A",
#'         flowCore::linearTransform(
#'             a = 0.1,
#'             b = 0
#'        )
#'     ),
#'     flowCore::transformList(
#'         "540/30Violet-A",
#'         flowCore::logicleTransform()
#'     )
#' )
#' 
#' ret1 <- getTransfoParams(translist, channel = "FSC-A")
#' ret1$type # "linear"
#' ret1$paramsList # a = 0.1, b = 0.
#' 
#' ret2 <- getTransfoParams(translist, channel = "525/50Violet-A")
#' ret2$type # "logicle"
#' ret2$paramsList # a = 0., w = 0.2834, m = 4.5, t = 262143
#' 
#' ret3 <- getTransfoParams(translist, channel = "540/30Violet-A")
#' ret3$type # "logicle
#' ret3$paramsList # a = 0., w = 0.5, m = 4.5, t = 262144
getTransfoParams <- function(transList,
                             channel) {
    if (!inherits(transList, "transformList")) {
        stop("transList parameter should be a flowCore::transformList!")
    }
    # browser()
    transMap <- transList@transforms[[channel]]
    if (is.null(transMap)) {
        return(NULL)
    } else {
        if (methods::.hasSlot(transMap, "f")) {
            tf <- methods::new("transform", .Data = transMap@f)
        } else if (methods::.hasSlot(transMap, ".Data")) {
            tf <- methods::new("transform", .Data = transMap@.Data)
        } else {
            stop(
                "transfo on channel does not have 'f' or '.Data' slot ",
                "=> not handled"
            )
        }

        sm <- flowCore::summary(tf)
        ret <- list()
        if (!is.null(sm$k) && methods::is(sm$k, "transform")) {
            ret$type <- "logicle"
            ret$paramsList <- list(
                a = unname(sm$a),
                w = unname(sm$w),
                m = unname(sm$m),
                t = unname(sm$t)
            )
        } else if (!is.null(sm$t) && methods::is(sm$t, "transform")) {
            ret$type <- "linear"
            ret$paramsList <- list(
                a = unname(sm$a),
                b = unname(sm$b)
            )
        } else {
            stop(
                "transformation type not recognized, ",
                "currently this function only ",
                "works with linear or logicle transforms"
            )
        }
        return(ret)
    }
}


#' @title compute linear transformation of scatter channels found in ff, based
#' on 5% and 95% of referenceChannel, set as target. If there is a
#' transformation defined in transList for referenceChannel, it is applied
#' first, before computing quantiles.
#' Then the computed linear transformations (or each scatter channel) are added
#' into the transfo_list. -A channels are computed, and same linear
#' transformation is then applied to corresponding -W and -H channels
#' (if they exist in ff).
#' @description based on a referenceChannel
#' @param ff a flowCore::flowFrame
#' @param transList an initial flowCore::transformList
#' @param referenceChannel the reference channel to take target quantile values
#' from. Can be defined as marker or channel name.
#' @param silent if FALSE, will output some information on the computed
#' linear transformations
#'
#' @return the transList with added linear scale transformations
#' @export
#'
#' @examples
#' 
#' ff <- OMIP021Samples[[1]]
#' refMarker <- "APCCy7 - CD4"
#' refChannel <- "780/60Red-A"

#' transList <- flowCore::estimateLogicle(ff,
#'                                        channels = refChannel)
#' retTransList <-
#'     computeScatterChannelsLinearScale(ff,
#'                                       transList = transList,
#'                                       referenceChannel = refMarker,
#'                                       silent = TRUE
#'     )
computeScatterChannelsLinearScale <- function(ff,
                                              transList = NULL,
                                              referenceChannel,
                                              silent = TRUE) {
    # check inputs

    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame")
    }
    fluoChannels <- flowCore::colnames(ff)[areFluoCols(ff)]
    scatterChannels <-
        flowCore::colnames(ff)[!areFluoCols(ff) & areSignalCols(ff)]
    if (length(scatterChannels) == 0) {
        warning("no scatter channel to scale")
        return(transList)
    }

    referenceChannel <- flowCore::getChannelMarker(ff, referenceChannel)$name

    if (!(referenceChannel %in% fluoChannels)) {
        stop("referenceChannel should be a fluorochrome channel of ff")
    }

    if (is.null(transList)) {
        if (!silent) {
            message(
                "NULL transList found...\n",
                "Continued with no transfo applied on reference channel"
            )
        }
        ff_t <- ff[, referenceChannel]
    } else {
        if (!inherits(transList, "transformList")) {
            stop("transList parameter should be a flowCore::transformList!")
        }
        if (is.null(transList@transforms[[referenceChannel]])) {
            if (!silent) {
                message(
                    "No transformation found for referenceChannel ",
                    "in transList\n",
                    "Continued with no transfo applied on reference ",
                    "channel"
                )
            }
            ff_t <- ff[, referenceChannel]
        } else {
            transfoList <-
                flowCore::transformList(
                    from = referenceChannel,
                    tfun = transList@transforms[[referenceChannel]]@f
                )
            ff_t <- flowCore::transform(ff[, referenceChannel], transfoList)
        }
    }


    q5Goal <- stats::quantile(flowCore::exprs(ff_t)[, referenceChannel], 0.05)
    q95Goal <- stats::quantile(flowCore::exprs(ff_t)[, referenceChannel], 0.95)

    # adapt scatter channels to have the same percentiles
    # -A channels are modified independently,
    #  and the SAME transfo are applied to corresponding -H and -W channels
    # (if they exist)
    foundAreaScatter <- FALSE
    if ("FSC-A" %in% scatterChannels) {
        ch <- "FSC-A"
        q5FSCA <- stats::quantile(flowCore::exprs(ff)[, ch], 0.05)
        q95FSCA <- stats::quantile(flowCore::exprs(ff)[, ch], 0.95)
        FSCAa <- (q95Goal - q5Goal) / (q95FSCA - q5FSCA)
        FSCAb <- q5Goal - q5FSCA * (q95Goal - q5Goal) / (q95FSCA - q5FSCA)
        if (!silent) {
            message(
                "applying specific linear transformation ",
                "for FSC-A channel..."
            )
            message(
                "initial quantiles : q5 = ", round(q5FSCA, 4), " ; q95 = ",
                round(q95FSCA, 4)
            )
            message(
                "target quantiles : q5 = ", round(q5Goal, 4), " ; q95 = ",
                round(q95Goal, 4)
            )
            message(
                "a = ",
                formatC(FSCAa, format = "e", digits = 2),
                " ; b = ",
                formatC(FSCAb, format = "e", digits = 2)
            )
        }
        tf <- flowCore::linearTransform(a = FSCAa, b = FSCAb)
        if (is.null(transList)) {
            transList <-
                flowCore::transformList(
                    from = ch,
                    tfun = tf
                )
        } else {
            transList@transforms[[ch]] <- NULL
            transList <- c(transList, flowCore::transformList(ch, tf))
            # transList@transforms[[ch]] <- tf
        }

        foundAreaScatter <- TRUE

        if ("FSC-W" %in% scatterChannels) {
            ch <- "FSC-W"
            if (!silent) {
                message(
                    "applying FSC-A linear transformation ",
                    "for FSC-W channel..."
                )
            }
            transList@transforms[[ch]] <- NULL
            transList <- c(transList, flowCore::transformList(ch, tf))
            # transList@transforms[[ch]] <- tf
        }
        if ("FSC-H" %in% scatterChannels) {
            ch <- "FSC-H"
            if (!silent) {
                message(
                    "applying FSC-A linear transformation ",
                    "for FSC-H channel..."
                )
            }
            transList@transforms[[ch]] <- NULL
            transList <- c(transList, flowCore::transformList(ch, tf))
            # transList@transforms[[ch]] <- tf
        }
    }
    if ("SSC-A" %in% scatterChannels) {
        ch <- "SSC-A"
        q5SSCA <- stats::quantile(flowCore::exprs(ff)[, ch], 0.05)
        q95SSCA <- stats::quantile(flowCore::exprs(ff)[, ch], 0.95)
        SSCAa <- (q95Goal - q5Goal) / (q95SSCA - q5SSCA)
        SSCAb <- q5Goal - q5SSCA * (q95Goal - q5Goal) / (q95SSCA - q5SSCA)
        if (!silent) {
            message(
                "applying specific linear transformation ",
                "for SSC-A channel..."
            )
            message(
                "initial quantiles : q5 = ", round(q5SSCA, 4), " ; q95 = ",
                round(q95SSCA, 4)
            )
            message(
                "target quantiles : q5 = ", round(q5Goal, 4), " ; q95 = ",
                round(q95Goal, 4)
            )
            message(
                "a = ",
                formatC(SSCAa, format = "e", digits = 2),
                " ; b = ",
                formatC(SSCAb, format = "e", digits = 2)
            )
        }
        tf <- flowCore::linearTransform(a = SSCAa, b = SSCAb)
        if (is.null(transList)) {
            transList <-
                flowCore::transformList(
                    from = ch,
                    tfun = tf
                )
        } else {
            transList@transforms[[ch]] <- NULL
            transList <- c(transList, flowCore::transformList(ch, tf))
            # transList@transforms[[ch]] <- tf
        }

        foundAreaScatter <- TRUE

        if ("SSC-W" %in% scatterChannels) {
            ch <- "SSC-W"
            if (!silent) {
                message(
                    "applying SSC-A linear transformation ",
                    "for SSC-W channel..."
                )
            }
            transList@transforms[[ch]] <- NULL
            transList <- c(transList, flowCore::transformList(ch, tf))
            # transList@transforms[[ch]] <- tf
        }
        if ("SSC-H" %in% scatterChannels) {
            ch <- "SSC-H"
            if (!silent) {
                message(
                    "applying SSC-A linear transformation ",
                    "for SSC-H channel..."
                )
            }
            transList@transforms[[ch]] <- NULL
            transList <- c(transList, flowCore::transformList(ch, tf))
            # transList@transforms[[ch]] <- tf
        }
    }

    if (!foundAreaScatter) {
        warning(
            "did not find any -A scatters in channels => ",
            "did not transform anything\n",
            "The following scatter channels were found: ",
            scatterChannels
        )
    }

    return(transList)
}

#' @title find time channel in flowSet/flowFrame
#' @description tries to find a channel in a flowSet/flowFrame that could
#' be the time channel. First tries to identify a channel name containing the
#' 'time' string, then tries to identify a single monotonically increasing
#' channel.
#'
#' @param obj a flowCore::flowFrame or flowCore::flowSet
#' @param excludeChannels vector of column names to exclude in the search
#'
#' @return a character, name of the found channel that should be representing
#' time. If not found, returns NULL.
#' @export
#' @examples 
#' ret <- findTimeChannel(OMIP021Samples[[1]])
#' ret # "Time"
#'
findTimeChannel <- function(obj, excludeChannels = c()) {
    isFlowSet <- FALSE
    if (inherits(obj, "flowSet")) {
        isFlowSet <- TRUE
    } else if (inherits(obj, "flowFrame")) {
    } else {
        stop("obj type not recognized, should be a flowFrame or flowSet")
    }

    includedChannels <-
        flowCore::colnames(obj)[!(flowCore::colnames(obj) %in% excludeChannels)]

    time <- grep("^Time$", includedChannels,
        value = TRUE,
        ignore.case = TRUE
    )[1]
    if (is.na(time)) {
        if (isFlowSet) {
            xx <- flowCore::exprs(obj[[1]])[, includedChannels]
        } else if (methods::is(obj, "flowFrame")) {
            xx <- flowCore::exprs(obj)[, includedChannels]
        }
        cont <- apply(xx, 2, function(y) {
            all(sign(diff(y)) >=
                0)
        })
        cont <- apply(xx, 2, function(y) all(y == cummax(y)))
        time <- names(which(cont))
    }
    if (!length(time) || length(time) > 1) {
        time <- NULL
    }
    return(time)
}

#' @title get channel names from markers
#' @description finds name of channels corresponding to user provided markers
#'
#' @param ff a flowCore::flowFrame
#' @param markers a vector of markers, either provided as :
#' - an array of booleans (referring to flowFrame columns)
#' - an array of integers (indices in flowFrame columns)
#' - an array of characters (exact markers or channel patterns)
#'
#' @return a character vector, containing the names of the corresponding
#' channels
#' @export
#'
#' @examples 
#' # with existing markers
#' ret <- getChannelNamesFromMarkers(
#'     OMIP021Samples[[1]],
#'     c(
#'         "FSC-A",
#'         "L/D Aqua - Viability",
#'         "FITC - gdTCR",
#'         "PECy5 - CD28"
#'     ))
#'     
#' ret # c("FSC-A", "525/50Violet-A", "530/30Blue-A", "670/30Yellow-A")
#' 
#' # with boolean vector
#' indices <- c(1, 6, 14, 18)
#' boolInput <- rep(FALSE, 21)
#' boolInput[indices] <- TRUE
#' ret2 <- getChannelNamesFromMarkers(
#'     OMIP021Samples[[1]],
#'     boolInput)
#'     
#' ret2 # c("FSC-A", "525/50Violet-A", "530/30Blue-A", "670/30Yellow-A")
#' 
#' # with indices vector
#' ret3 <- getChannelNamesFromMarkers(
#'     OMIP021Samples[[1]],
#'     indices
#' )
#' ret3 # c("FSC-A", "525/50Violet-A", "530/30Blue-A", "670/30Yellow-A")
#' 
#' 
getChannelNamesFromMarkers <- function(ff, markers) {
    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame")
    }

    frameChannels <- unname(flowCore::parameters(ff)@data[["name"]])
    frameMarkers <- unname(flowCore::parameters(ff)@data[["desc"]])

    if (is.logical(markers)) {
        markers <- which(markers)
    }
    channelNames <- c()
    for (marker in markers) {
        if (is.numeric(marker)) {
            iChannel <- marker
        } else {
            marker <- paste0("^\\Q", marker, "\\E$")
            iChannel <- grep(marker, frameMarkers)
        }
        if (length(iChannel) != 0) {
            for (i in iChannel) {
                channel <- frameChannels[iChannel]
                names(channel) <- frameMarkers[iChannel]
                channelNames <- c(channelNames, channel)
            }
        } else {
            iChannel <- grep(marker, frameChannels)
            if (length(iChannel) != 0) {
                channel <- frameChannels[iChannel]
                names(channel) <- channel
                channelNames <- c(channelNames, channel)
            } else {
                stop("Marker", marker, "could not be found")
            }
        }
    }
    return(unname(channelNames))
}

#' @title update marker name of a given flowFrame channel
#' @description : in a flowCore::flowFrame, update the marker name (stored in
#' 'desc' of parameters data) of a given channel.
#' @param ff a flowCore::flowFrame
#' @param channel the channel for which to update the marker name
#' @param newMarkerName the new marker name to be given to the selected channel
#'
#' @return a new flowCore::flowFrame with the updated marker name
#' @export
#' @examples
#' retFF <- updateMarkerName(OMIP021Samples[[1]],
#'                           channel = "FSC-A",
#'                           newMarkerName = "Fwd Scatter-A")
#'
updateMarkerName <- function(ff, channel, newMarkerName) {
    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame!")
    }
    
    channelIndex <- which(flowCore::colnames(ff) == channel)
    if (length(channelIndex) == 0) {
        stop("channel not found in flowFrame!")
    } else {
        channelIndex <- channelIndex[1]
    }
    
    param <- flowCore::parameters(ff)
    paramData <- flowCore::pData(param)
    paramData[channelIndex, "desc"] <- newMarkerName
    flowCore::pData(param) <- paramData
    flowCore::parameters(ff) <- param
    return(ff)
}

#' @title remove channels from a flowFrame
#' @description : in a flowCore::flowFrame, remove the channels of the given
#' names.
#' @param ff a flowCore::flowFrame
#' @param channels the channel names to be removed
#'
#' @return a new flowCore::flowFrame with the removed channels
#' @export
#' @examples
#' retFF <- removeChannels(OMIP021Samples[[1]],
#'                         channel = "FSC-A")
#'
removeChannels <- function(ff, channels) {
    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame!")
    }    
    keptCols <- rep(TRUE, length(flowCore::colnames(ff)))
    for (ch in channels) {
        channelIndex <- which(flowCore::colnames(ff) == ch)
        if (length(channelIndex) == 0) {
            warning("channel ", ch, " not found in flowFrame => ignoring...")
        } else {
            channelIndex <- channelIndex[1]
        }
        keptCols[channelIndex] <- FALSE
    }
    ff <- ff[, keptCols]
    return(ff)
}

#' @title append 'Original_ID' column to a flowframe
#' @description : on a flowCore::flowFrame, append a 'Original_ID' column.
#' This column can be used in plots comparing the events pre and post gating.
#' If the 'Original_ID' column already exists, the function does nothing
#' @param ff a flowCore::flowFrame
#' @param eventIDs an integer vector containing the values to be added in
#' as Original ID's
#' 
#' @return new flowCore::flowFrame containing the added 'Original_ID' column
#' @export
appendCellID <- function(ff, eventIDs = seq_len(flowCore::nrow(ff))) {
    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame")
    }
    if (!("Original_ID" %in% colnames(flowCore::exprs(ff)))) {
        matrixCellIds <- matrix(
            data = eventIDs, ncol = 1,
            dimnames = list(c(), list("Original_ID"))
        )
        ff <- flowCore::fr_append_cols(ff, matrixCellIds)
    }
    return(ff)
}

# #' @title add compensation to fluo column names
# #' @description : technical utility to add "Comp-" prefix to all
# #' column names for fluorochrom channels
# #' @param ff a flowCore::flowFrame
# #' @return new flowCore::flowFrame with the new column names
.addCompensation2FluoChannelNames <- function(ff) {
    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame")
    }
    areFluoCols <- areFluoCols(ff)
    newColNames <- flowCore::colnames(ff)
    newColNames <-
        mapply(
            FUN = function(theName, need2Do) {
                newName <- theName
                if (need2Do) newName <- paste0("Comp-", theName)
                newName
            },
            newColNames, areFluoCols
        )
    flowCore::colnames(ff) <- newColNames
    return(ff)
}

#' @title get fcs file name 
#' @description get basename of $FILENAME keyword if exists
#' @param ff a flowCore::flowFrame
#' @return the basename of $FILENAME keyword
#' @export
#' @examples 
#' fName <- getFCSFileName(OMIP021Samples[[1]])
getFCSFileName <- function(ff) {
    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame")
    }
    fName <- basename(flowCore::keyword(ff, keyword = "FILENAME")$FILENAME)
    if (is.null(fName)) stop("No FILENAME keyword for flowFrame")
    return(fName)
}
