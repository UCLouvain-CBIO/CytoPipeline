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

##' @include utils.R
NULL


#' @title plot flow rate as a function of time, using ggplot2
#'
#' @param obj       a flowCore::flowFrame or flowCore::flowSet
#' @param title     a title for the graph
#' @param timeUnit  which time interval is used to calculate "instant" flow rate
#'                  (default = 100 ms)
#'
#' @return a ggplot graph
#' @import ggplot2
#' @export
#'
#' @examples 
#' 
#' data(OMIP021Samples)
#' 
#' # single flowFrame plot
#' ggplotFlowRate(OMIP021Samples[[1]])
#' 
#' # two flowFrames plot 
#' ggplotFlowRate(OMIP021Samples)
#' 
#' # single plot with title
#' ggplotFlowRate(OMIP021Samples[[1]], title = "Test Flow Rate plot")
#' 
#' # explicit time unit
#' ggplotFlowRate(OMIP021Samples[[1]], timeUnit = 50)
#' 
ggplotFlowRate <- function(obj, title = "Flow Rate", timeUnit = 100) {
    isFlowSet <- FALSE
    if (inherits(obj, "flowSet")) {
        isFlowSet <- TRUE
    } else if (inherits(obj, "flowFrame")) {
    } else {
        stop("obj type not recognized, should be a flowFrame or flowSet")
    }

    flowRateDataFrame <- function(ff, timeUnit) {
        Y <- flowCore::exprs(ff)
        timeChName <- findTimeChannel(ff)
        if (is.null(timeChName)) {
            stop(
                "No Time columns in flow frame ",
                "=> impossible to calculate flow rate"
            )
        }

        h <- graphics::hist(
            flowCore::exprs(ff)[, timeChName],
            breaks = seq(min(flowCore::exprs(ff)[, timeChName]),
                max(flowCore::exprs(ff)[, timeChName]) +
                    timeUnit,
                by = timeUnit
            ),
            plot = FALSE
        )

        df <- data.frame(
            time = h$mids,
            nbEvents = h$counts * 1000 / timeUnit,
            name = flowCore::identifier(ff)
        )

        df
    }


    if (isFlowSet) {
        df <- flowCore::fsApply(obj,
            FUN = flowRateDataFrame,
            timeUnit = timeUnit
        )
        df <- do.call(rbind.data.frame, df)
    } else {
        df <- flowRateDataFrame(obj, timeUnit = timeUnit)
    }

    # following 2 statements just to allow R cmd CHECK w/o note
    time <- NULL
    nbEvents <- NULL
    pTime <- ggplot(df) +
        theme_gray() +
        geom_point(aes(x = time, y = nbEvents)) +
        geom_line(aes(x = time, y = nbEvents)) +
        ggtitle(title) +
        facet_wrap(~name) +
        xlab("Time") +
        ylab("Nr of events per second")
    return(pTime)
}

# internal function to select the appropriate logicle scale to use from ggCyto
.my_scale_logicle <- function(axis = c("x", "y"),
                              usedScale = c("logicle", "flowjo_biexp"),
                              logicleParams,
                              axisRange,
                              interactive = FALSE) {

    axis <- match.arg(axis)
    usedScale <- match.arg(usedScale)
    
    if (usedScale == "logicle") {
        theArgs <- logicleParams
    } else {
        theArgs <- list(maxValue = logicleParams$t,
                        widthBasis = -10^logicleParams$w,
                        pos = logicleParams$m - logicleParams$w,
                        neg = logicleParams$a)
    }
    theArgs$limits <- axisRange
    if (interactive) {
        theArgs$label <- scales::scientific_format()
    }
    
    if (axis == "x") {
        if (usedScale == "logicle") {
            do.call(scale_x_logicle,
                    args = theArgs)
        } else {
            do.call(scale_x_flowjo_biexp,
                    args = theArgs)
        }
    } else {
        if (usedScale == "logicle") {
            do.call(scale_y_logicle,
                    args = theArgs)
        } else {
            do.call(scale_y_flowjo_biexp,
                    args = theArgs)
        } 
    }
}

#' @title plot events in 1D or 2D, using ggplot2
#' @description  plot events of specific channels of either :
#' flowCore::flowFrame, or flowCore::flowSet
#' in 2D or 1D, mimicking FlowJo type of graph. \cr
#' if 1D : geom_density will be used \cr
#' if 2D : geom_hex will be used \cr
#'
#' @param obj a flowCore::flowFrame or flowCore::flowSet
#' @param xChannel channel (name or index) or marker name to be displayed
#' on x axis
#' @param yChannel channel (name or index) or marker name to be displayed
#' on y axis
#' @param nDisplayCells maximum number of events that will be plotted. If
#' the number of events exceed this number, a sub-sampling will be performed
#' @param seed seed used for sub-sampling (if any)
#' @param bins used in geom_hex
#' @param fill used in geom_density
#' @param alpha used in geom_density
#' @param xScale scale to be used for the x axis 
#' (note "linear" corresponds to no transformation)
#' @param yScale scale to be used for the y axis
#' (note "linear" corresponds to no transformation)
#' @param xLogicleParams if (xScale == "logicle"), the parameters of the logicle
#' transformation to be used, as a list(w = ..., m = ..., a = ..., t = ...).
#' If NULL, these parameters will be estimated by flowCore::estimateLogicle()
#' @param yLogicleParams if (yScale == "logicle"), the parameters of the logicle
#' transformation to be used, as a list(w = ..., m = ..., a = ..., t = ...).
#' If NULL, these parameters will be estimated by flowCore::estimateLogicle()
#' @param xLinearRange if (xScale == "linear"), the x axis range to be used
#' @param yLinearRange if (yScale == "linear"), the y axis range to be used
#' @param transList optional list of scale transformations to be applied to each
#' channel. If it is non null, 'x/yScale', 'x/yLogicleParams' and
#' 'x/yLinear_range' will be discarded.
#' @param runTransforms (TRUE/FALSE) Will the application of non linear scale
#' result in data being effectively transformed ?
#' - If TRUE, than the data will undergo transformations prior to
#' visualization.
#' - If FALSE, the axis will be scaled but the data themselves will not be
#' transformed.
#' @return a list of ggplot objects
#' @import ggplot2
#' @importFrom ggcyto scale_x_logicle
#' @importFrom ggcyto scale_y_logicle
#' @importFrom ggcyto scale_x_flowjo_biexp
#' @importFrom ggcyto scale_y_flowjo_biexp
#' @importFrom rlang .data
#' @export
#'
#' @examples 
#' 
#' data(OMIP021Samples)
#' 
#' ### 1D Examples ###
#' 
#' # simple linear scale example
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "FSC-A",
#'              xScale = "linear")
#' 
#' # with explicit linear range
#' ggplotEvents(OMIP021Samples[[1]],
#'                   xChannel = "FSC-A",
#'                   xScale = "linear",
#'                   xLinearRange = c(0, 250000))
#' 
#' # with linear scale, several flow frames
#' ggplotEvents(OMIP021Samples, xChannel = "FSC-A", xScale = "linear")
#' 
#' # simple logicle scale example
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "450/50Violet-A",
#'              xScale = "logicle")
#' 
#' # logicle scale, explicit parameters
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "450/50Violet-A",
#'              xScale = "logicle", xLogicleParams = list(
#'                  a = 1,
#'                  w = 2,
#'                  m = 7,
#'                  t = 270000))
#' 
#' # with sub-sampling
#' ggplotEvents(OMIP021Samples[[2]],
#'              xChannel = "450/50Violet-A",
#'              xScale = "logicle", nDisplayCells = 5000)
#' 
#' # tuning some plot parameters
#' ggplotEvents(OMIP021Samples[[2]],
#'              xChannel = "450/50Violet-A",
#'              xScale = "logicle", alpha = 0.5, fill = "red")
#' 
#' # examples that use a transformation list, estimated after compensation
#' compensationMatrix <- flowCore::spillover(OMIP021Samples[[1]])$SPILL
#' 
#' ffC <- runCompensation(OMIP021Samples[[1]],
#'                        spillover = compensationMatrix,
#'                        updateChannelNames = FALSE)
#' 
#' transList <- flowCore::estimateLogicle(
#'     ffC,
#'     colnames(compensationMatrix))
#' 
#' transList <-
#'     c(transList,
#'       flowCore::transformList(
#'           "FSC-A",
#'           flowCore::linearTransform(a = 0.00001)))
#' 
#' # linear example, without running the transformations on data
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "450/50Violet-A",
#'              xScale = "linear", 
#'              transList = transList,
#'              runTransforms = FALSE)
#' 
#' # linear example, now running the transformations on data
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "450/50Violet-A",
#'              xScale = "linear", 
#'              transList = transList,
#'              runTransforms = TRUE)
#' 
#' # logicle example, without running the transformations on data
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "FSC-A",
#'              xScale = "logicle", 
#'              transList = transList,
#'              runTransforms = FALSE)
#' 
#' # logicle example, now running the transformations on data
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "FSC-A",
#'              xScale = "logicle", 
#'              transList = transList,
#'              runTransforms = TRUE)
#' 
#' ### 2D examples ###
#' 
#' 
#' # simple linear example
#' ggplotEvents(OMIP021Samples[[1]],
#'                   xChannel = "FSC-A",
#'                   xScale = "linear",
#'                   yChannel = "610/20Violet-A",
#'                   yScale = "logicle")
#' 
#' # simple linear example, 2 flow frames
#' ggplotEvents(OMIP021Samples,
#'              xChannel = "FSC-A",
#'              xScale = "linear",
#'              yChannel = "SSC-A",
#'              yScale = "linear")
#' 
#' # logicle vs linear example
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "450/50Violet-A",
#'              xScale = "logicle",
#'              yChannel = "SSC-A",
#'              yScale = "linear")
#' 
#' # 2X logicle example
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "TETaGC",
#'              xScale = "logicle",
#'              yChannel = "CD27",
#'              yScale = "logicle")
#' 
#' # tuning nb of bins
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "TETaGC",
#'              xScale = "logicle",
#'              yChannel = "CD27",
#'              yScale = "logicle",
#'              bins = 128)
#' 
#' # using transformation list, not run on data
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "TETaGC",
#'              xScale = "logicle",
#'              yChannel = "CD27",
#'              yScale = "logicle",
#'              transList = transList,
#'              runTransforms = FALSE)
#' 
#' # using transformation list, run on data                  
#' ggplotEvents(OMIP021Samples[[1]],
#'              xChannel = "TETaGC",
#'              xScale = "logicle",
#'              yChannel = "CD27",
#'              yScale = "logicle",
#'              transList = transList,
#'              runTransforms = TRUE)
#' 
ggplotEvents <- function(obj,
                         xChannel,
                         yChannel = NULL,
                         nDisplayCells = Inf,
                         seed = NULL,
                         bins = 216,
                         fill = "lightblue",
                         alpha = 0.2,
                         xScale = c("linear", "logicle"),
                         yScale = c("linear", "logicle"),
                         xLogicleParams = NULL,
                         yLogicleParams = NULL,
                         xLinearRange = NULL,
                         yLinearRange = NULL,
                         transList = NULL,
                         runTransforms = FALSE) {


    #browser()
    
    xScale <- match.arg(xScale)
    yScale <- match.arg(yScale)

    isFlowSet <- FALSE
    if (inherits(obj, "flowSet")) {
        fr <- obj[[1]]
        isFlowSet <- TRUE
    } else if (inherits(obj, "flowFrame")) {
        fr <- obj
    } else {
        stop("obj type not recognized, should be a flowFrame or flowSet")
    }
    
    # selection of scale_logical function to use (temporary fix due to ggcyto)
    used_logical_scale <- "logicle"
    #used_logical_scale <- "flowjo_biexp"

    # find channel and marker names to specify axis labels

    xChMk <- flowCore::getChannelMarker(fr, xChannel)
    xChannel <- xChMk$name
    xLabel <- xChannel
    if (!is.na(xChMk$desc)) {
        xLabel <- paste0(xLabel, " : ", xChMk$desc)
    }

    if (!is.null(yChannel)) {
        yChMk <- flowCore::getChannelMarker(fr, yChannel)
        yChannel <- yChMk$name
        yLabel <- yChannel
        if (!is.na(yChMk$desc)) {
            yLabel <- paste0(yLabel, " : ", yChMk$desc)
        }
    }
    
    # reduce flow frame or flow set to the displayed channels
    selectedCols <- xChannel
    if (!is.null(yChannel)){
        selectedCols <- c(selectedCols, yChannel)
    }
    if (isFlowSet) {
        obj <- flowCore::fsApply(obj,
                                 FUN = function(ff, selectedCols){
                                     ff <- ff[, selectedCols]
                                 },
                                 selectedCols = selectedCols
        )    
    } else {
        obj <- obj[, selectedCols]
    }
    
    # build scale transforms for not linear scales,
    # when these are not directly provided by the user
    if (is.null(transList) && (xScale != "linear" || yScale != "linear")) {
        estimatedCols <- c()
        if (xScale == "logicle" && is.null(xLogicleParams)) {
            estimatedCols <- xChannel
        }
        if (!is.null(yChannel) && yScale == "logicle" && 
            is.null(yLogicleParams)) {
            estimatedCols <- c(estimatedCols, yChannel)
        }
        if (!is.null(estimatedCols)) {
            if (isFlowSet){
                transList <- flowCore::estimateLogicle(obj[[1]], estimatedCols)
            } else {
                transList <- flowCore::estimateLogicle(obj, estimatedCols)
            }
        }
        
        # add the self provided logicle parameters
        if (xScale == "logicle" && !is.null(xLogicleParams)){
            xLogicleTrans <- flowCore::logicleTransform(
                w = xLogicleParams$w,
                t = xLogicleParams$t,
                m = xLogicleParams$m,
                a = xLogicleParams$a
            )
            if (is.null(transList)) {
                transList <- flowCore::transformList(from = xChannel,
                                                     tfun = xLogicleTrans)
            } else {
                transList@transforms[[xChannel]] <- xLogicleTrans
            }
        }
        if (!is.null(yChannel) && yScale == "logicle" &&
            !is.null(yLogicleParams)) {
            yLogicleTrans <- flowCore::logicleTransform(
                w = yLogicleParams$w,
                t = yLogicleParams$t,
                m = yLogicleParams$m,
                a = yLogicleParams$a
            )
            if (is.null(transList)) {
                transList <- flowCore::transformList(from = yChannel,
                                                     tfun = yLogicleTrans)
            } else {
                transList@transforms[[yChannel]] <- yLogicleTrans
            }
        }
    }

    # perform sub-sampling if necessary

    if (nDisplayCells < Inf) {
        if (nDisplayCells < 1) {
            stop("n_display_cells should be strictly positive!")
        }
        if (isFlowSet) {
            obj <- flowCore::fsApply(obj,
                FUN = subsample,
                nEvents = nDisplayCells,
                seed = seed
            )
        } else {
            obj <- subsample(obj, nEvents = nDisplayCells, seed = seed)
        }
    }

    xTransformed <- FALSE
    yTransformed <- FALSE
    if (!is.null(transList)) {
        res <- getTransfoParams(transList, xChannel)
        if (!is.null(res)) {
            if (runTransforms) {
                xScale <- "linear"
                xTransformed <- TRUE
                if (res$type == "logicle") {
                    xLinearRange <- c(0, res$paramsList$m)
                }
            } else {
                xScale <- res$type
                if (res$type == "logicle") {
                    xLogicleParams <- res$paramsList
                }
            }
        }

        if (!is.null(yChannel)) {
            res <- getTransfoParams(transList, yChannel)
            if (!is.null(res)) {
                if (runTransforms) {
                    yScale <- "linear"
                    yTransformed <- TRUE
                    if (res$type == "logicle") {
                        yLinearRange <- c(0, res$paramsList$m)
                    }
                } else {
                    yScale <- res$type
                    if (res$type == "logicle") {
                        yLogicleParams <- res$paramsList
                    }
                }
            }
        }

        if (runTransforms) {
            # adapt scale to linear and linear_range to null
            # if trans_list is passed and is effectively run
            appliedTransList <- c()

            # remove not mentioned channels from trans_list
            appliedTransList <- transList
            transChannels <- names(transList@transforms)
            for (n in transChannels) {
                if (!(n %in% c(xChannel, yChannel))) {
                    appliedTransList@transforms[[n]] <- NULL
                }
            }
            if (isFlowSet) {
                obj <- flowCore::fsApply(obj, FUN = function(ff, transList) {
                    flowCore::transform(obj, translist = transList)
                }, transList = appliedTransList)
            } else {
                obj <- flowCore::transform(obj, translist = appliedTransList)
            }
        }
    }

    # find axis ranges depending on scales

    xScale <- match.arg(xScale)
    yScale <- match.arg(yScale)

    if (xScale == "logicle") {
        if (is.null(xLogicleParams)) {
            stop("xLogicleParams can't be NULL if xScale == logicle")
        } else {
            xAxesLimits <- c(0, xLogicleParams$m)
            myTrans <-
                do.call(flowCore::logicleTransform, args = xLogicleParams)
            myInverseTrans <- flowCore::inverseLogicleTransform(myTrans)
            myXRange <- myInverseTrans@.Data(xAxesLimits)
        }
    }
    if (yScale == "logicle") {
        if (is.null(yLogicleParams)) {
            stop("yLogicleParams can't be NULL if yScale == logicle")
        } else {
            yAxesLimits <- c(0, yLogicleParams$m)
            myTrans <-
                do.call(flowCore::logicleTransform, args = yLogicleParams)
            myInverseTrans <- flowCore::inverseLogicleTransform(myTrans)
            myYRange <- myInverseTrans@.Data(yAxesLimits)
        }
    }

    # main plot layer

    if (xTransformed) {
        xLabel <- paste0(xLabel, " (transformed)")
    }
    if (yTransformed) {
        yLabel <- paste0(yLabel, " (transformed)")
    }
    
    # browser()

    if (is.null(yChannel)) {
        p <- ggplot(
            data = obj,
            mapping = aes(x = .data[[xChannel]])
        ) +
            geom_density(fill = fill, alpha = alpha, na.rm = TRUE) +
            xlab(xLabel)
    } else {
        p <- ggplot(
            data = obj,
            mapping = aes(
                x = .data[[xChannel]],
                y = .data[[yChannel]])) +
            geom_hex(bins = bins, na.rm = TRUE) +
            scale_fill_gradientn(
            colours = grDevices::rainbow(
                n = 5, end = 0.65,
                rev = TRUE,
                s = 0.8, v = 0.8)) + 
            xlab(xLabel) +
            ylab(yLabel)
    }


    # add faceting

    p <- p + facet_wrap(~name) +
        theme(legend.position = "none")

    # add scales for x and y, as well as axis limits

    if (xScale == "logicle") {
        if (is.null(yChannel)) {
            p <- p +
                .my_scale_logicle(axis = "x", 
                                  usedScale = used_logical_scale,
                                  logicleParams = xLogicleParams,
                                  axisRange = myXRange) #+
                #coord_cartesian(xlim = myXRange)
        } else if (yScale == "logicle") {
            p <- p + 
                .my_scale_logicle(axis = "x", 
                                  usedScale = used_logical_scale,
                                  logicleParams = xLogicleParams,
                                  axisRange = myXRange) + 
                .my_scale_logicle(axis = "y", 
                                  usedScale = used_logical_scale,
                                  logicleParams = yLogicleParams,
                                  axisRange = myYRange)  #+
                # coord_cartesian(
                #     xlim = myXRange,
                #     ylim = myYRange
                # )
        } else {
            p <- p +
                .my_scale_logicle(axis = "x", 
                                  usedScale = used_logical_scale,
                                  logicleParams = xLogicleParams,
                                  axisRange = myXRange) + 
                scale_y_continuous(limits = yLinearRange)
                # coord_cartesian(
                #     xlim = myXRange,
                #     ylim = yLinearRange
                # )
        }
    } else {
        if (is.null(yChannel)) {
            p <- p +
                scale_x_continuous(limits = xLinearRange)
                #coord_cartesian(xlim = xLinearRange)
        } else if (yScale == "logicle") {
            p <- p +
                scale_x_continuous(limits = xLinearRange) + 
                .my_scale_logicle(axis = "y", 
                                  usedScale = used_logical_scale,
                                  logicleParams = yLogicleParams,
                                  axisRange = myYRange) #+
                # coord_cartesian(
                #     xlim = xLinearRange,
                #     ylim = myYRange
                # )
        } else {
            p <- p +
                scale_x_continuous(limits = xLinearRange) + 
                scale_y_continuous(limits = yLinearRange) # +
                # coord_cartesian(
                #     xlim = xLinearRange,
                #     ylim = yLinearRange
                # )
        }
    }

    p
}


#' @title plot filtered events in 2D, using ggplot
#' @description  plot events of specific channels of either :
#' flowCore::flowFrame, or
#' flowCore::flowSet
#' in 2D, showing the impact of applying a filter between :
#' - a 'pre' flowframe
#  - a 'post' flowframe
#'
#' @param ffPre a flowCore::flowFrame, before applying filter
#' @param ffPost a flowCore::flowFrame, after applying filter
#' @param xChannel channel (name or index) or marker name to be displayed
#' on x axis
#' @param yChannel channel (name or index) or marker name to be displayed
#' on y axis
#' @param nDisplayCells maximum number of events that will be plotted. If
#' the number of events exceed this number, a subsampling will be performed
#' @param seed seed used for sub-sampling (if any)
#' @param size used by geom_point()
#' @param xScale scale to be used for the x axis
#' (note "linear" corresponds to no transformation)
#' @param yScale scale to be used for the y axis
#' (note "linear" corresponds to no transformation)
#' @param xLogicleParams if (xScale == "logicle"), the parameters of the logicle
#' transformation to be used, as a list(w = ..., m = ..., a = ..., t = ...)
#' If NULL, these parameters will be estimated by flowCore::estimateLogicle()
#' @param yLogicleParams if (yScale == "logicle"), the parameters of the logicle
#' transformation to be used, as a list(w = ..., m = ..., a = ..., t = ...)
#' If NULL, these parameters will be estimated by flowCore::estimateLogicle()
#' @param xLinearRange if (xScale == "linear"), linear range to be used
#' @param yLinearRange if (yScale == "linear"), linear range to be used
#' @param transList optional list of scale transformations to be applied to each
#' channel. If it is non null, 'x/yScale', 'x/yLogicleParams' and
#' 'x/yLinear_range' will be discarded.
#' @param runTransforms (TRUE/FALSE) Will the application of non linear scale
#' result in data being effectively transformed ?
#' - If TRUE, than the data will undergo transformations prior to
#' visualization.
#' - If FALSE, the axis will be scaled but the data themselves are
#' not transformed.
#' @param interactive if TRUE, transform the scaling formats such that the
#' ggcyto::x_scale_logicle() and ggcyto::y_scale_logicle() do work with
#' plotly::ggplotly()
#'
#' @return a ggplot object
#' @import ggplot2
#' @importFrom ggcyto scale_x_logicle
#' @importFrom ggcyto scale_y_logicle
#' @importFrom ggcyto scale_x_flowjo_biexp
#' @importFrom ggcyto scale_y_flowjo_biexp
#' @importFrom rlang .data
#' @export
#'
#' @examples 
#' 
#' data(OMIP021Samples)
#' 
#' ffPre <- OMIP021Samples[[1]]
#' 
#' # creating a manual polygon gate filter based on channels L/D and FSC-A
#' 
#' LDMarker <- "L/D Aqua - Viability"
#' 
#' LDChannel <- getChannelNamesFromMarkers(ffPre, markers = LDMarker)
#' liveGateMatrix <- matrix(
#'     data = c(
#'         50000, 50000, 100000, 200000, 200000,
#'         100, 1000, 2000, 2000, 1
#'     ),
#'     ncol = 2,
#'     dimnames = list(
#'         c(),
#'         c("FSC-A", LDChannel)
#'     )
#' )
#' 
#' liveGate <- flowCore::polygonGate(
#'     filterId = "Live",
#'     .gate = liveGateMatrix
#' )
#' 
#' selectedLive <- flowCore::filter(ffPre, liveGate)
#' ffL <- flowCore::Subset(ffPre, selectedLive)
#' 
#' 
#' # show the results
#' 
#' # subsample 5000 points    
#' ggplotFilterEvents(
#'     ffPre = ffPre,
#'     ffPost = ffL,
#'     nDisplayCells = 5000,
#'     xChannel = "FSC-A", xScale = "linear",
#'     yChannel = LDMarker, yScale = "logicle") +
#'     ggplot2::ggtitle("Live gate filter - 5000 points")
#' 
#' # with all points
#' ggplotFilterEvents(
#'     ffPre = ffPre,
#'     ffPost = ffL,
#'     nDisplayCells = Inf,
#'     xChannel = "FSC-A", xScale = "linear",
#'     yChannel = LDMarker, yScale = "logicle") +
#'     ggplot2::ggtitle("Live gate filter - all points")
#'
#' 
ggplotFilterEvents <- function(ffPre, ffPost,
                               xChannel,
                               yChannel,
                               nDisplayCells = 10000,
                               seed = NULL,
                               size = 0.5,
                               xScale = c("linear", "logicle"),
                               yScale = c("linear", "logicle"),
                               xLogicleParams = NULL,
                               yLogicleParams = NULL,
                               xLinearRange = NULL,
                               yLinearRange = NULL,
                               transList = NULL,
                               runTransforms = FALSE,
                               interactive = FALSE) {
    if (!inherits(ffPre, "flowFrame")) {
        stop("ffPre type not recognized, should be a flowFrame")
    }

    if (!inherits(ffPost, "flowFrame")) {
        stop("ffPost type not recognized, should be a flowFrame")
    }
    
    xScale <- match.arg(xScale)
    yScale <- match.arg(yScale)
    
    # selection of scale_logical function to use (temporary fix due to ggcyto)
    used_logical_scale <- "logicle"
    #used_logical_scale <- "flowjo_biexp"

    # find channel and marker names to specify axis labels

    # browser()

    xChMk <- flowCore::getChannelMarker(ffPre, xChannel)
    xChannel <- xChMk$name
    xLabel <- xChannel
    if (!is.na(xChMk$desc)) {
        xLabel <- paste0(xLabel, " : ", xChMk$desc)
    }

    if (!is.null(yChannel)) {
        yChMk <- flowCore::getChannelMarker(ffPre, yChannel)
        yChannel <- yChMk$name
        yLabel <- yChannel
        if (!is.na(yChMk$desc)) {
            yLabel <- paste0(yLabel, " : ", yChMk$desc)
        }
    }

    if (!"Original_ID" %in% flowCore::colnames(ffPre)) {
        ffPre <- appendCellID(ffPre)
    }

    if (!"Original_ID" %in% flowCore::colnames(ffPost)) {
        stop(
            "ffPost should have Original_ID column in expression matrix ",
            "for the filter display to work!"
        )
    }

    # reduce flow frame or flow set to the displayed channels
    selectedCols <- xChannel
    if (!is.null(yChannel)){
        selectedCols <- c(selectedCols, yChannel)
    }
    selectedCols <- c(selectedCols, "Original_ID")

    ffPre <- ffPre[, selectedCols]
    ffPost <- ffPost[, selectedCols]


    # build scale transforms for not linear scales,
    # when these are not directly provided by the user
    if (is.null(transList) && (xScale != "linear" || yScale != "linear")) {
        estimatedCols <- c()
        if (xScale == "logicle" && is.null(xLogicleParams)) {
            estimatedCols <- xChannel
        }
        if (!is.null(yChannel) && yScale == "logicle" &&
            is.null(yLogicleParams)) {
            estimatedCols <- c(estimatedCols, yChannel)
        }
        if (!is.null(estimatedCols)) {
            transList <- flowCore::estimateLogicle(ffPre, estimatedCols)
        }

        # add the self provided logicle parameters
        if (xScale == "logicle" && !is.null(xLogicleParams)){
            xLogicleTrans <- flowCore::logicleTransform(
                w = xLogicleParams$w,
                t = xLogicleParams$t,
                m = xLogicleParams$m,
                a = xLogicleParams$a
            )
            if (is.null(transList)) {
                transList <- flowCore::transformList(from = xChannel,
                                                     tfun = xLogicleTrans)
            } else {
                transList@transforms[[xChannel]] <- xLogicleTrans
            }
        }
        if (!is.null(yChannel) && yScale == "logicle" &&
            !is.null(yLogicleParams)) {
            yLogicleTrans <- flowCore::logicleTransform(
                w = yLogicleParams$w,
                t = yLogicleParams$t,
                m = yLogicleParams$m,
                a = yLogicleParams$a
            )
            if (is.null(transList)) {
                transList <- flowCore::transformList(from = yChannel,
                                                     tfun = yLogicleTrans)
            } else {
                transList@transforms[[yChannel]] <- yLogicleTrans
            }
        }
    }

    # perform sub-sampling if necessary

    nEvents <- nrow(flowCore::exprs(ffPre))
    if (nDisplayCells < 1) stop("n_display_cells should be strictly positive!")
    i <- 0
    if (nDisplayCells < nEvents) {
        if (!is.null(seed)) {
            withr::with_seed(
                seed,
                i <- sample(nEvents, nDisplayCells)
            )
        } else {
            i <- sample(nEvents, nDisplayCells)
        }
    } else {
        i <- seq(nEvents)
    }
    
    xTransformed <- FALSE
    yTransformed <- FALSE
    if (!is.null(transList)) {
        res <- getTransfoParams(transList, xChannel)
        if (!is.null(res)) {
            if (runTransforms) {
                xScale <- "linear"
                xTransformed <- TRUE
                if (res$type == "logicle") {
                    xLinearRange <- c(0, res$paramsList$m)
                }
            } else {
                xScale <- res$type
                if (res$type == "logicle") {
                    xLogicleParams <- res$paramsList
                }
            }
        }
        
        if (!is.null(yChannel)) {
            res <- getTransfoParams(transList, yChannel)
            if (!is.null(res)) {
                if (runTransforms) {
                    yScale <- "linear"
                    yTransformed <- TRUE
                    if (res$type == "logicle") {
                        yLinearRange <- c(0, res$paramsList$m)
                    }
                } else {
                    yScale <- res$type
                    if (res$type == "logicle") {
                        yLogicleParams <- res$paramsList
                    }
                }
            }
        }
        
        if (runTransforms) {
            # adapt scale to linear and linear_range to null
            # if trans_list is passed and is effectively run
            appliedTransList <- c()
            
            # remove not mentioned channels from trans_list
            appliedTransList <- transList
            transChannels <- names(transList@transforms)
            for (n in transChannels) {
                if (!(n %in% c(xChannel, yChannel))) {
                    appliedTransList@transforms[[n]] <- NULL
                }
            }
            
            ffPre <- flowCore::transform(ffPre, translist = appliedTransList)
            ffPost <- flowCore::transform(ffPost, translist = appliedTransList)
        }
    }

    # find axis ranges depending on scales

    xScale <- match.arg(xScale)
    yScale <- match.arg(yScale)

    if (xScale == "logicle") {
        if (is.null(xLogicleParams)) {
            stop("x_logicle_params can't be NULL if xScale == logicle")
        } else {
            xAxesLimits <- c(0, xLogicleParams$m)
            myTrans <-
                do.call(flowCore::logicleTransform, args = xLogicleParams)
            myInverseTrans <- flowCore::inverseLogicleTransform(myTrans)
            myXRange <- myInverseTrans@.Data(xAxesLimits)
        }
    }
    if (yScale == "logicle") {
        if (is.null(yLogicleParams)) {
            stop("yLogicleParams can't be NULL if yScale == logicle")
        } else {
            yAxesLimits <- c(0, yLogicleParams$m)
            myTrans <-
                do.call(flowCore::logicleTransform, args = yLogicleParams)
            myInverseTrans <- flowCore::inverseLogicleTransform(myTrans)
            myYRange <- myInverseTrans@.Data(yAxesLimits)
        }
    }

    # main plot layer
    
    if (xTransformed) {
        xLabel <- paste0(xLabel, " (transformed)")
    }
    if (yTransformed) {
        yLabel <- paste0(yLabel, " (transformed)")
    }
    
    df <- data.frame(
        x = flowCore::exprs(ffPre)[, xChannel],
        y = flowCore::exprs(ffPre)[, yChannel]
    )
    
    # following 2 statements just to allow R cmd CHECK w/o note
    # x <- NULL
    # y <- NULL
    p <- ggplot(df[i, ], aes(x = .data$x, y = .data$y)) +
        geom_point(
            size = size,
            color = ifelse(flowCore::exprs(ffPre)[i, "Original_ID"] %in%
                flowCore::exprs(ffPost)[, "Original_ID"],
            "blue",
            "red"
            )
        ) +
        xlab(xLabel) +
        ylab(yLabel) +
        # theme_minimal() +
        theme(legend.position = "none")
    # ggtitle("my title")

    # add faceting
    # p <- p + facet_wrap(~name) +


    # add scales for x and y, as well as axis limits

    if (xScale == "logicle") {
        if (is.null(yChannel)) {
            p <- p +
                .my_scale_logicle(axis = "x", 
                                  usedScale = used_logical_scale,
                                  logicleParams = xLogicleParams,
                                  axisRange = myXRange,
                                  interactive = interactive) +
                coord_cartesian(xlim = myXRange)
        } else if (yScale == "logicle") {
            p <- p +
                .my_scale_logicle(axis = "x", 
                                  usedScale = used_logical_scale,
                                  logicleParams = xLogicleParams,
                                  axisRange = myXRange,
                                  interactive = interactive) +
                .my_scale_logicle(axis = "y", 
                                  usedScale = used_logical_scale,
                                  logicleParams = yLogicleParams,
                                  axisRange = myYRange,
                                  interactive = interactive) +
                coord_cartesian(
                    xlim = myXRange,
                    ylim = myYRange
                )
        } else {
            p <- p +
                .my_scale_logicle(axis = "x", 
                                  usedScale = used_logical_scale,
                                  logicleParams = xLogicleParams,
                                  axisRange = myXRange,
                                  interactive = interactive) +
                coord_cartesian(
                    xlim = myXRange,
                    ylim = yLinearRange
                )
        }
    } else {
        if (is.null(yChannel)) {
            p <- p +
                coord_cartesian(xlim = xLinearRange)
        } else if (yScale == "logicle") {
            p <- p +
                .my_scale_logicle(axis = "y", 
                                  usedScale = used_logical_scale,
                                  logicleParams = yLogicleParams,
                                  axisRange = myYRange,
                                  interactive = interactive) +
                coord_cartesian(
                    xlim = xLinearRange,
                    ylim = myYRange
                )
        } else {
            p <- p +
                coord_cartesian(
                    xlim = xLinearRange,
                    ylim = yLinearRange
                )
        }
    }
    return(p)
}
