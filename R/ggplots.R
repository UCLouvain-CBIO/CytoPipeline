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
ggplotFlowRate <- function(obj, title = "Flow Rate", timeUnit = 100)
{
  isFlowSet <- FALSE
  if (inherits(obj, "flowSet")) {
    isFlowSet <- TRUE
  } else if (inherits(obj, "flowFrame")) {
  } else {
    stop("obj type not recognized, should be a flowFrame or flowSet")
  }

  flowRateDataFrame <- function(ff, timeUnit){
    Y <- flowCore::exprs(ff)
    timeChName <- findTimeChannel(ff)
    if (is.null(timeChName)) {
      stop("No Time columns in flow frame => impossible to calculate flow rate")
    }

    h <- graphics::hist(flowCore::exprs(ff)[, timeChName],
                        breaks = seq(min(flowCore::exprs(ff)[, timeChName]),
                                     max(flowCore::exprs(ff)[, timeChName]) + timeUnit,
                                     by = timeUnit),
                        plot = FALSE)

    df <- data.frame(time = h$mids,
                     nbEvents = h$counts * 1000 / timeUnit,
                     name = flowCore::identifier(ff))

    df
  }
  
  
  if (isFlowSet) {
    df <- flowCore::fsApply(obj, FUN = flowRateDataFrame,
                               timeUnit = timeUnit)
    df <- do.call(rbind.data.frame, df)

  } else {
    df <- flowRateDataFrame(obj, timeUnit = timeUnit)
  }

  # following 2 statements just to allow R cmd CHECK w/o note
  time <- NULL
  nbEvents <- NULL
  pTime <- ggplot(df) + theme_gray() +
    geom_point(aes(x = time, y = nbEvents)) +
    geom_line(aes(x = time, y = nbEvents)) +
    ggtitle(title) +
    facet_wrap(~name) +
    xlab("Time") + ylab("Nr of events per second")
  return(pTime)

}

#' @title plot events in 1D or 2D, using ggplot2
#' @description  plot events of specific channels of either :
#' flowCore::flowFrame, or flowCore::flowSet
#' in 2D or 1D, mimicking FlowJo type of graph. \cr
#' if 1D : geom_density will be used \cr
#' if 2D : geom_hex will be used \cr
#'
#' @param obj a flowCore::flowFrame or flowCore::flowSet
#' @param xChannel channel (name or index) or marker name to be displayed on x axis
#' @param yChannel channel (name or index) or marker name to be displayed on y axis
#' @param nDisplayCells maximum number of events that will be plotted. If
#' the number of events exceed this number, a sub-sampling will be performed
#' @param seed seed used for sub-sampling (if any)
#' @param bins used in geom_hex
#' @param fill used in geom_density
#' @param alpha used in geom_density
#' @param xScale scale to be used for the x axis
#' @param yScale scale to be used for the y axis
#' @param xLogicleParams if (xScale == "logicle"), the parameters of the logicle
#' transformation to be used, as a list(w = ..., m = ..., a = ..., t = ...)
#' @param yLogicleParams if (yScale == "logicle"), the parameters of the logicle
#' transformation to be used, as a list(w = ..., m = ..., a = ..., t = ...)
#' @param xLinearRange if (xScale == "linear"), the x axis range to be used
#' @param yLinearRange if (yScale == "linear"), the y axis range to be used
#' @param transList optional list of scale transformations to be applied to each 
#' channel. If it is non null, 'x/yScale', 'x/yLogicleParams' and 
#' 'x/yLinear_range' will be discarded.
#' @param runTransforms (TRUE/FALSE) only taken into account if transList is 
#' not NULL. Will 'transList' result in data being effectively transformed ?
#' - If TRUE, than the data will undergo transformations prior to
#' visualization.
#' - If FALSE, the axis will be scaled but the data themselves are
#' not transformed.
#' @return a list of ggplot objects
#' @import ggplot2
#' @importFrom ggcyto scale_x_logicle
#' @importFrom ggcyto scale_y_logicle
#' @export
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
                         xLogicleParams = list(w = 2, m = 6.42,
                                                 a = 0, t = 262144),
                         yLogicleParams = list(w = 2, m = 6.42,
                                                 a = 0, t = 262144),
                         xLinearRange = NULL,
                         yLinearRange = NULL,
                         transList = NULL,
                         runTransforms = FALSE){


  #browser()

  isFlowSet <- FALSE
  if (inherits(obj, "flowSet")) {
    fr <- obj[[1]]
    isFlowSet <- TRUE
  } else if (inherits(obj, "flowFrame")) {
    fr <- obj
  } else {
    stop("obj type not recognized, should be a flowFrame or flowSet")
  }

  # find channel and marker names to specify axis labels

  xChMk <- flowCore::getChannelMarker(fr, xChannel)
  xChannel <- xChMk$name
  xLabel <- xChannel
  if (!is.na(xChMk$desc)) {
    xLabel <- paste0(xLabel, " : ", xChMk$desc)
  }

  if(!is.null(yChannel)){
    yChMk <- flowCore::getChannelMarker(fr, yChannel)
    yChannel <- yChMk$name
    yLabel <- yChannel
    if (!is.na(yChMk$desc)) {
      yLabel <- paste0(yLabel, " : ", yChMk$desc)
    }
  }

  # perform sub-sampling if necessary

  if (nDisplayCells < Inf) {
    if (nDisplayCells < 1) stop("n_display_cells should be strictly positive!")
    if (isFlowSet) {
      obj <- flowCore::fsApply(obj, FUN = subsample,
                     nSamples = nDisplayCells,
                     seed = seed)
    } else {
      obj <- subsample(obj, nSamples = nDisplayCells, seed = seed)
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
        if(res$type == "logicle"){
          xLogicleParams <- res$paramsList
        }
      }
    }
    
    if (!is.null(yChannel)) {
      res <- getTransfoParams(transList, yChannel)
      if(!is.null(res)){
        if (runTransforms) {
          yScale <- "linear"
          yTransformed <- TRUE
          if (res$type == "logicle") {
            yLinearRange <- c(0, res$paramsList$m)
          }
        } else {
          yScale <- res$type
          if(res$type == "logicle"){
            yLogicleParams <- res$paramsList
          }
        }
      }
    }
    
    if (runTransforms) {
      # adapt scale to linear and linear_range to null if trans_list is passed
      # and is effectively run
      appliedTransList <- c()
      
      # remove not mentioned channels from trans_list
      appliedTransList <- transList
      transChannels <- names(transList@transforms)
      for (n in transChannels) {
        if (! (n %in% c(xChannel, yChannel))) {
          appliedTransList@transforms[[n]] <- NULL
        }
      }
      if (isFlowSet) {
        obj <- flowCore::fsApply(obj, FUN = function(ff, transList){
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

  if(xScale == "logicle"){
    if(is.null(xLogicleParams)){
      stop("xLogicleParams can't be NULL if xScale == logicle")
    } else {
      xAxesLimits <- c(0, xLogicleParams$m)
      myTrans <- do.call(flowCore::logicleTransform, args = xLogicleParams)
      myInverseTrans <- flowCore::inverseLogicleTransform(myTrans)
      myXRange <- myInverseTrans@.Data(xAxesLimits)
    }
  }
  if(yScale == "logicle"){
    if(is.null(yLogicleParams)){
      stop("yLogicleParams can't be NULL if yScale == logicle")
    } else {
      yAxesLimits <- c(0, yLogicleParams$m)
      myTrans <- do.call(flowCore::logicleTransform, args = yLogicleParams)
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

  if(is.null(yChannel)){
    p <- ggplot(data = obj,
                mapping = aes_q(x = as.symbol(xChannel))) +
      geom_density(fill = fill, alpha = alpha) +
      xlab(xLabel)
  } else {
    p <- ggplot(data = obj,
                mapping = aes_q(x = as.symbol(xChannel),
                                y = as.symbol(yChannel))) +
      #geom_hex(bins = bins) + # will apply geom_hex after coord_cartesian
      xlab(xLabel) +
      ylab(yLabel)
  }


  # add faceting

  p <- p + facet_wrap(~name) +
    theme(legend.position="none")

  # add scales for x and y, as well as axis limits

  if(xScale == "logicle"){
    if(is.null(yChannel)){
      p <- p +
        do.call(scale_x_logicle, args = xLogicleParams) +
        coord_cartesian(xlim = myXRange)
    } else if(yScale == "logicle") {
      p <- p +
        do.call(scale_x_logicle, args = xLogicleParams) +
        do.call(scale_y_logicle, args = yLogicleParams) +
        coord_cartesian(xlim = myXRange,
                        ylim = myYRange)
    } else {
      p <- p +
        do.call(scale_x_logicle, args = xLogicleParams) +
        coord_cartesian(xlim = myXRange,
                        ylim = yLinearRange)
    }
  } else {
    if(is.null(yChannel)){
      p <- p +
        coord_cartesian(xlim = xLinearRange)
    } else if(yScale == "logicle") {
      p <- p +
        do.call(scale_y_logicle, args = yLogicleParams) +
        coord_cartesian(xlim = xLinearRange,
                        ylim = myYRange)
    } else {
      p <- p +
        coord_cartesian(xlim = xLinearRange,
                        ylim = yLinearRange)
    }
  }

  # if 2D, apply geom_hex at the end (after axis ranges)
  # and hex fill colour palette
  if(!is.null(yChannel)){
    p <- p + geom_hex(bins = bins) +
      scale_fill_gradientn(colours = grDevices::rainbow(n = 5, end = 0.65,
                                                        rev = TRUE,
                                                        s = 0.8, v = 0.8))
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
#' @param xChannel channel (name or index) or marker name to be displayed on x axis
#' @param yChannel channel (name or index) or marker name to be displayed on y axis
#' @param nDisplayCells maximum number of events that will be plotted. If
#' the number of events exceed this number, a subsampling will be performed
#' @param seed seed used for sub-sampling (if any)
#' @param size used by geom_point()
#' @param xScale scale to be used for the x axis
#' @param yScale scale to be used for the y axis
#' @param xLogicleParams if (xScale == "logicle"), the parameters of the logicle
#' transformation to be used, as a list(w = ..., m = ..., a = ..., t = ...)
#' @param yLogicleParams if (yScale == "logicle"), the parameters of the logicle
#' transformation to be used, as a list(w = ..., m = ..., a = ..., t = ...)
#' @param xLinearRange if (xScale == "linear"), linear range to be used
#' @param yLinearRange if (yScale == "linear"), linear range to be used
#' @param interactive if TRUE, transform the scaling formats such that the 
#' ggcyto::x_scale_logicle() and ggcyto::y_scale_logicle() do work with 
#' plotly::ggplotly()
#'
#' @return a ggplot object
#' @import ggplot2
#' @importFrom ggcyto scale_x_logicle
#' @importFrom ggcyto scale_y_logicle
#' @export
#'
ggplotFilterEvents <- function(ffPre, ffPost,
                               xChannel,
                               yChannel,
                               nDisplayCells = 10000,
                               seed = NULL,
                               size = 0.5,
                               xScale = c("linear", "logicle"),
                               yScale = c("linear", "logicle"),
                               xLogicleParams = list(w = 2, m = 6.42,
                                                       a = 0, t = 262144),
                               yLogicleParams = list(w = 2, m = 6.42,
                                                       a = 0, t = 262144),
                               xLinearRange = NULL,
                               yLinearRange = NULL,
                               interactive = FALSE){
  if (!inherits(ffPre, "flowFrame")) {
    stop("ffPre type not recognized, should be a flowFrame")
  }

  if (!inherits(ffPost, "flowFrame")) {
    stop("ffPost type not recognized, should be a flowFrame")
  }

  # find channel and marker names to specify axis labels

  #browser()

  xChMk <- flowCore::getChannelMarker(ffPre, xChannel)
  xChannel <- xChMk$name
  xLabel <- xChannel
  if (!is.na(xChMk$desc)) {
    xLabel <- paste0(xLabel, " : ", xChMk$desc)
  }
  
  if(!is.null(yChannel)){
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

  if (!"Original_ID" %in% flowCore::colnames(ffPost))
    stop("ffPost should have Original_ID column in expression matrix for the ",
         "filter display to work!")

  df <- data.frame(x = flowCore::exprs(ffPre)[,xChannel],
                   y = flowCore::exprs(ffPre)[,yChannel])

  # perform sub-sampling if necessary

  nEvents <- nrow(df)
  if(nDisplayCells < 1) stop("n_display_cells should be strictly positive!")
  i <- 0
  if(nDisplayCells < nEvents){
    if(!is.null(seed)){
      set.seed(seed)
    }
    i <- sample(nEvents, nDisplayCells)
  } else {
    i <- seq(nEvents)
  }

  # find axis ranges depending on scales

  xScale <- match.arg(xScale)
  yScale <- match.arg(yScale)

  if(xScale == "logicle"){
    if(is.null(xLogicleParams)){
      stop("x_logicle_params can't be NULL if xScale == logicle")
    } else {
      xAxesLimits <- c(0, xLogicleParams$m)
      myTrans <- do.call(flowCore::logicleTransform, args = xLogicleParams)
      myInverseTrans <- flowCore::inverseLogicleTransform(myTrans)
      myXRange <- myInverseTrans@.Data(xAxesLimits)
    }
  }
  if(yScale == "logicle"){
    if(is.null(yLogicleParams)){
      stop("yLogicleParams can't be NULL if yScale == logicle")
    } else {
      yAxesLimits <- c(0, yLogicleParams$m)
      myTrans <- do.call(flowCore::logicleTransform, args = yLogicleParams)
      myInverseTrans <- flowCore::inverseLogicleTransform(myTrans)
      myYRange <- myInverseTrans@.Data(yAxesLimits)
    }
  }

  # plot main layer
  # following 2 statements just to allow R cmd CHECK w/o note
  x <- NULL
  y <- NULL
  p <- ggplot(df[i,], aes(x = x, y = y)) +
    geom_point(size = size,
               color = ifelse(flowCore::exprs(ffPre)[i,"Original_ID"] %in%
                                flowCore::exprs(ffPost)[,"Original_ID"],
                              'blue',
                              'red')) +
    xlab(xLabel) +
    ylab(yLabel) +
    #theme_minimal() +
    theme(legend.position = "none")
    #ggtitle("my title")

  # add faceting
  #p <- p + facet_wrap(~name) +


  # add scales for x and y, as well as axis limits

  if(xScale == "logicle"){
    xArgs <- xLogicleParams
    if (interactive) {
      xArgs <- c(xArgs, label=scales::scientific_format())
    } 
    if(is.null(yChannel)){
      p <- p +
        do.call(scale_x_logicle, args = xArgs) +
        coord_cartesian(xlim = myXRange)
    } else if(yScale == "logicle") {
      yArgs <- yLogicleParams
      if (interactive) {
        yArgs <- c(yArgs, label=scales::scientific_format())
      } 
      p <- p +
        do.call(scale_x_logicle, args = xArgs) +
        do.call(scale_y_logicle, args = yArgs) +
        coord_cartesian(xlim = myXRange,
                        ylim = myYRange)
    } else {
      p <- p +
        do.call(scale_x_logicle, args = xArgs) +
        coord_cartesian(xlim = myXRange,
                        ylim = yLinearRange)
    }
  } else {
    if(is.null(yChannel)){
      p <- p +
        coord_cartesian(xlim = xLinearRange)
    } else if(yScale == "logicle") {
      yArgs <- yLogicleParams
      if (interactive) {
        yArgs <- c(yArgs, label=scales::scientific_format())
      }  
      p <- p +
        do.call(scale_y_logicle, args = yArgs) +
        coord_cartesian(xlim = xLinearRange,
                        ylim = myYRange)
    } else {
      p <- p +
        coord_cartesian(xlim = xLinearRange,
                        ylim = yLinearRange)
    }
  }
  return(p)

}





