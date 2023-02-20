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


#' @title Clean doublet events from flow cytometry data
#' @description will adjust a polygon gate aimed at cleaning doublet events
#' from the flowFrame.
#' The main idea is to use the ratio between the two indicated channel as an
#' indicator and select only the events for which this ratio is 'not too far'
#' from the median ratio. More specifically, the computed ratio is ch1/(1+ch2).
#' However, instead of looking at a constant range of this ratio, as is done in
#' PeacoQC::removeDoublets(), which leads to a semi-conic gate,
#' we apply a parallelogram shaped gate, by keeping a constant range of channel
#' 2 intensity, based on the target ratio range at the mid value of channel 1.
#' @param ff A flowCore::flowframe that contains flow cytometry data.
#' @param filterId the name for the filter that is returned
#' @param channel1 The first channel that will be used to determine the doublet
#' events. Default is "FSC-A"
#' @param channel2 The second channels that will be used to determine the
#' doublet events. Default is "FSC-H"
#' @param nmad Bandwidth above the ratio allowed (cells are kept if their ratio
#' is smaller than the median ratio + nmad times the median absolute deviation
#' of the ratios). Default is 4.
#' @param verbose If set to TRUE, the median ratio and width will be printed.
#' Default is FALSE.
#'
#' @return This function returns a flowCore::polygonGate.
#' @export
#'
#' @examples
#' 
#' data(OMIP021Samples)
#'
#' # simple example with one single singlets gate filter 
#' # FSC-A and FSC-H channels are used by default
#' 
#' mySingletsGate <- singletsGate(OMIP021Samples[[1]], nmad = 3)
#' 
#' selectedSinglets <- flowCore::filter(
#'     OMIP021Samples[[1]],
#'     mySingletsGate)
#' 
#' ff_l <- flowCore::Subset(OMIP021Samples[[1]], selectedSinglets)
#' 
#' linRange <- c(0, 250000)
#' 
#' ggplotFilterEvents(
#'     ffPre = OMIP021Samples[[1]],
#'     ffPost = ff_l,
#'     xChannel = "FSC-A", xLinearRange = linRange,
#'     yChannel = "FSC-H", yLinearRange = linRange)
#' 
#' # application of two singlets gates one after the other
#' 
#' singletsGate1 <- singletsGate(OMIP021Samples[[1]], nmad = 3)
#' singletsGate2 <- singletsGate(OMIP021Samples[[1]],
#'                               channel1 = "SSC-A",
#'                               channel2 = "SSC-H",
#'                               filterId = "Singlets2")
#' 
#' singletCombinedGate <- singletsGate1 & singletsGate2
#' 
#' selectedSinglets <- flowCore::filter(
#'     OMIP021Samples[[1]],
#'     singletCombinedGate)
#' 
#' ff_l <- flowCore::Subset(OMIP021Samples[[1]], selectedSinglets)
#' 
#' ggplotFilterEvents(
#'     ffPre = OMIP021Samples[[1]],
#'     ffPost = ff_l,
#'     xChannel = "FSC-A", xLinearRange = linRange,
#'     yChannel = "FSC-H", yLinearRange = linRange)
#' 
#' ggplotFilterEvents(
#'     ffPre = OMIP021Samples[[1]],
#'     ffPost = ff_l,
#'     xChannel = "SSC-A", xLinearRange = linRange,
#'     yChannel = "SSC-H", yLinearRange = linRange)
#' 
singletsGate <- function(ff, filterId = "Singlets",
                         channel1 = "FSC-A", channel2 = "FSC-H", nmad = 4,
                         verbose = FALSE) {
    if (!inherits(ff, "flowFrame")) {
        stop("ff type not recognized, should be a flowFrame")
    }

    ratio <- flowCore::exprs(ff)[, channel1] / (1 + flowCore::exprs(ff)[
        ,
        channel2
    ])

    ratioMedian <- stats::median(ratio)
    ratioMad <- stats::mad(ratio)
    ch1Median <- stats::median(flowCore::exprs(ff)[, channel1])
    ch1Min <- min(flowCore::exprs(ff)[, channel1])
    ch1Max <- max(flowCore::exprs(ff)[, channel1])
    ch2RefPoint <- ch1Median / ratioMedian - 1
    ch2Min <- ch1Median / (ratioMedian + nmad * ratioMad) - 1
    ch2Range <- ch2RefPoint - ch2Min
    if (verbose) {
        message(
            "Median ratio: ", ratioMedian, ", width: ", nmad * ratioMad,
            ", Ch1 median : ", ch1Median, ", Ch 2 range : ",
            ch2Range
        )
    }

    ch2TargetAtCh1Min <- ch1Min / ratioMedian - 1
    ch2TargetAtCh1Max <- ch1Max / ratioMedian - 1

    # browser()

    # parallelogram gate with vertical edges
    polygonGateMatrix <- matrix(
        data = c(
            ch1Min, ch1Min, ch1Max, ch1Max,
            ch2TargetAtCh1Min - ch2Range,
            ch2TargetAtCh1Min + ch2Range,
            ch2TargetAtCh1Max + ch2Range,
            ch2TargetAtCh1Max - ch2Range
        ),
        ncol = 2,
        dimnames = list(
            c(),
            c(channel1, channel2)
        )
    )


    singletsGate <- flowCore::polygonGate(
        filterId = filterId,
        .gate = polygonGateMatrix
    )

    return(singletsGate)
}
