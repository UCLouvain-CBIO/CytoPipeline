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

setClassUnion("DataFrameOrNull",
              members=c("data.frame", "NULL"))

#' @title CytoPipeline class
#'
#' @aliases CytoPipeline-class, CytoPipelineClass
#'
#' @name CytoPipeline-class
#'
#' @rdname CytoPipelineClass
#'
#' @description
#'
#' Class representing a flow cytometry pipeline, and composed of two processing
#' queues, i.e. lists of CytoProcessingStep objects :
#' - a list of CytoProcessingStep(s) for pre-calculation of scale
#' transformations per channel
#' - a list of CytoProcessingStep(s) for the pre-processing of flow frames
#'
#' @slot scaleTransformProcessingQueue A `list` of    
#' CytoProcessingStep objects containing the steps   
#' for obtaining the scale transformations per channel
#'
#' @slot flowFramesPreProcessingQueue A `list` of     
#' CytoProcessingStep objects containing the steps   
#' for pre-processing of the samples flow frames
#'
#' @slot experimentName A `character` containing     
#' the experiment (run) name
#'
#' @slot sampleFiles A `character` vector storing   
#' all fcs files to be run into the pipeline
#'
#' @slot pData An optional `data.frame` containing
#' additional information for each sample file. 
#' The `pData` raw names must correspond to `basename(sampleFiles)`
#' otherwise validation of the CytoPipeline object will fail!
#'
#' @exportClass CytoPipeline
#' 
#' @return nothing
#' 
#' @examples
#' 
#' ### *** EXAMPLE 1: building CytoPipeline step by step *** ###
#' 
#' rawDataDir <-
#'     system.file("extdata", package = "CytoPipeline")
#' experimentName <- "OMIP021_PeacoQC"
#' sampleFiles <- file.path(rawDataDir, list.files(rawDataDir,
#'                                              pattern = "Donor"))
#'                                              
#' outputDir <- base::tempdir()
#' 
#' # main parameters : sample files and output files
#' pipL <- CytoPipeline(experimentName = experimentName,
#'                      sampleFiles = sampleFiles)
#' 
#' ### SCALE TRANSFORMATION STEPS ###
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "scale transform",
#'                       CytoProcessingStep(
#'                           name = "flowframe_read",
#'                           FUN = "readSampleFiles",
#'                           ARGS = list(
#'                               whichSamples = "all",
#'                               truncate_max_range = FALSE,
#'                               min.limit = NULL
#'                           )
#'                       )
#'     )
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "scale transform",
#'                       CytoProcessingStep(
#'                           name = "remove_margins",
#'                           FUN = "removeMarginsPeacoQC",
#'                           ARGS = list()
#'                      )
#'     )
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "scale transform",
#'                       CytoProcessingStep(
#'                           name = "compensate",
#'                           FUN = "compensateFromMatrix",
#'                           ARGS = list(matrixSource = "fcs")
#'                       )
#'     )
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "scale transform",
#'                       CytoProcessingStep(
#'                           name = "flowframe_aggregate",
#'                           FUN = "aggregateAndSample",
#'                           ARGS = list(
#'                               nTotalEvents = 10000,
#'                               seed = 0
#'                           )
#'                       )
#'     )
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "scale transform",
#'                       CytoProcessingStep(
#'                           name = "scale_transform_estimate",
#'                           FUN = "estimateScaleTransforms",
#'                           ARGS = list(
#'                               fluoMethod = "estimateLogicle",
#'                               scatterMethod = "linear",
#'                               scatterRefMarker = "BV785 - CD3"
#'                           )
#'                       )
#'     )
#' 
#' ### PRE-PROCESSING STEPS ###
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "pre-processing",
#'                       CytoProcessingStep(
#'                           name = "flowframe_read",
#'                           FUN = "readSampleFiles",
#'                           ARGS = list(
#'                               truncate_max_range = FALSE,
#'                               min.limit = NULL
#'                           )
#'                       )
#'     )
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "pre-processing",
#'                       CytoProcessingStep(
#'                           name = "remove_margins",
#'                           FUN = "removeMarginsPeacoQC",
#'                           ARGS = list()
#'                       )
#'     )
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "pre-processing",
#'                       CytoProcessingStep(
#'                           name = "compensate",
#'                           FUN = "compensateFromMatrix",
#'                           ARGS = list(matrixSource = "fcs")
#'                       )
#'     )
#' 
#' pipL <-
#'     addProcessingStep(
#'         pipL,
#'         whichQueue = "pre-processing",
#'         CytoProcessingStep(
#'             name = "remove_debris",
#'             FUN = "removeDebrisManualGate",
#'             ARGS = list(
#'                 FSCChannel = "FSC-A",
#'                 SSCChannel = "SSC-A",
#'                 gateData =  c(73615, 110174, 213000, 201000, 126000,
#'                               47679, 260500, 260500, 113000, 35000)))
#'     )
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "pre-processing",
#'                       CytoProcessingStep(
#'                           name = "remove_dead_cells",
#'                           FUN = "removeDeadCellsManualGate",
#'                           ARGS = list(
#'                               FSCChannel = "FSC-A",
#'                               LDMarker = "L/D Aqua - Viability",
#'                               gateData = c(0, 0, 250000, 250000,
#'                                            0, 650, 650, 0)
#'                           )
#'                       )
#'     )
#' 
#' pipL <-
#'     addProcessingStep(
#'         pipL,
#'         whichQueue = "pre-processing",
#'         CytoProcessingStep(
#'             name = "perform_QC",
#'             FUN = "qualityControlPeacoQC",
#'             ARGS = list(
#'                 preTransform = TRUE,
#'                 min_cells = 150, # default
#'                 max_bins = 500, # default
#'                 step = 500, # default,
#'                 MAD = 6, # default
#'                 IT_limit = 0.55, # default
#'                 force_IT = 150, # default
#'                 peak_removal = 0.3333, # default
#'                 min_nr_bins_peakdetection = 10 # default
#'             )
#'         )
#'     )
#' 
#' pipL <-
#'     addProcessingStep(pipL,
#'                       whichQueue = "pre-processing",
#'                       CytoProcessingStep(
#'                           name = "transform",
#'                           FUN = "applyScaleTransforms",
#'                           ARGS = list()
#'                       )
#'     )
#' 
#' ### *** EXAMPLE 2: building CytoPipeline from JSON file *** ###
#' 
#' jsonDir <- system.file("extdata", package = "CytoPipeline")
#' jsonPath <- file.path(jsonDir, "pipelineParams.json")
#' 
#' pipL2 <- CytoPipeline(jsonPath,
#'                       experimentName = experimentName,
#'                       sampleFiles = sampleFiles)
#' 
setClass("CytoPipeline",
    slots = c(
        experimentName = "character",
        scaleTransformProcessingQueue = "list",
        flowFramesPreProcessingQueue = "list",
        sampleFiles = "character",
        pData = "DataFrameOrNull"
    ),
    prototype = list(
        experimentName = "default_experiment",
        scaleTransformProcessingQueue = list(),
        flowFramesPreProcessingQueue = list(),
        sampleFiles = character(),
        pData = NULL
    )
)

setValidity("CytoPipeline", function(object) {
    if (!is.null(object@pData)) {
        #browser()
        if (!inherits(object@pData, "data.frame")) {
            return("Non-null @pData slot should be a data.frame")
        }
        # if (!isTRUE(all.equal(rownames(object@pData), 
        #                       basename(object@sampleFiles))))
        if (!all(basename(object@sampleFiles) %in% rownames(object@pData))) {
            msg <- paste0("Row names of non-null @pData slot should contain ",
                          "all sample file basenames")
            return(msg)
        }
    }
    
    msg1 <- .validProcessingQueue(
        object@scaleTransformProcessingQueue,
        "scaleTransformProcessingQueue"
    )
    msg2 <- .validProcessingQueue(
        object@flowFramesPreProcessingQueue,
        "flowFramesPreProcessingQueue"
    )

    msg <- NULL
    if (length(msg1)) {
        msg <- msg1
        if (length(msg2)) {
            msg <- paste0(msg, " ; ", msg2)
        }
    } else {
        if (length(msg2)) {
            msg <- msg2
        }
    }
    if (length(msg)) {
        msg
    } else {
        TRUE
    }
})

#' @rdname CytoPipeline
#'
#' @param object a `CytoPipeline` object
#'
#' @importMethodsFrom methods show
#'
setMethod(
    "show", "CytoPipeline",
    function(object) {
        cat(
            "Pipeline object for flow cytometry experiment:",
            object@experimentName, "\n"
        )
        if (length(object@sampleFiles)) {
            cat(
                "Sample files:", length(object@sampleFiles),
                "sample file(s)\n"
            )
        } else {
            cat("No sample file\n")
        }
        if (!is.null(object@pData)) {
            cat("pheno data:\n")
            show(object@pData)
        } else {
            cat("No pheno data\n")
        }
        showProcessingSteps(object, whichQueue = "scale transform")
        showProcessingSteps(object, whichQueue = "pre-processing")
    }
)


setGeneric("CytoPipeline", function(object, ...) {
    standardGeneric("CytoPipeline")
})

#' @rdname CytoPipeline
#' @param experimentName the experiment name
#' @param sampleFiles the sample files
#' @param pData the pheno data (data.frame or NULL)
#' @export
#'
setMethod(
    "CytoPipeline", "missing",
    function(object,
             experimentName = "default_experiment",
             sampleFiles = character(),
             pData = NULL) {
        methods::new("CytoPipeline",
            experimentName = experimentName,
            scaleTransformProcessingQueue = list(),
            flowFramesPreProcessingQueue = list(),
            sampleFiles = sampleFiles,
            pData = pData
        )
    }
)

#' @rdname CytoPipeline
#'
#' @param object a `list()`
#' @param experimentName the experiment name
#' @param sampleFiles the sample files
#' @param pData the phenoData (data.frame or NULL)
#'
#' @export
#'
setMethod(
    "CytoPipeline", "list",
    function(object,
             experimentName = "default_experiment",
             sampleFiles = character(),
             pData = NULL) {
        x <- methods::new("CytoPipeline",
            experimentName = experimentName,
            scaleTransformProcessingQueue = list(),
            flowFramesPreProcessingQueue = list(),
            sampleFiles = sampleFiles,
            pData = pData
        )
        x <- .makeSlots(x, object)
        x <- .makeProcessingQueues(x, object)
        return(x)
    }
)

#' @rdname CytoPipeline
#'
#'
#' @param object a `character()` containing a JSON input
#' @param experimentName the experiment name
#' @param sampleFiles the sample files
#' @param pData the pheno Data (data.frame or NULL)
#'
#' @export
#'
setMethod(
    "CytoPipeline", "character",
    function(object,
             experimentName = "default_experiment",
             sampleFiles = character(),
             pData = NULL) {
        pipelineParams <- jsonlite::read_json(object,
            simplifyVector = TRUE,
            simplifyDataFrame = FALSE
        )
        cytoPipeline <- CytoPipeline(pipelineParams,
                                     experimentName = experimentName,
                                     sampleFiles = sampleFiles,
                                     pData = pData)
        #browser()
        return(cytoPipeline)
    }
)

#' @rdname CytoPipeline
#'
#' @param x a `CytoPipeline` object
#' @param ... additional arguments (not used here)
#'
#' @return - for `as.list.CytoPipeline`: the obtained list
#'
#' @export
#'
as.list.CytoPipeline <- function(x, ...) {
    stopifnot(inherits(x, "CytoPipeline"))

    # browser()

    slots <- methods::slotNames(x)
    slots <- slots[!slots %in% c(
        "scaleTransformProcessingQueue",
        "flowFramesPreProcessingQueue"
    )]

    outputList <- list()

    for (sl in slots) {
        outputList[[sl]] <- methods::slot(x, sl)
    }

    # now handles both processing queue
    pQueue2List <- function(pQueue) {
        pList <- list()
        for (p in seq_along(pQueue)) {
            pList[[p]] <- as.list.CytoProcessingStep(pQueue[[p]])
        }
        pList
    }

    pList1 <- pQueue2List(x@scaleTransformProcessingQueue)
    outputList[["scaleTransformProcessingSteps"]] <- pList1

    pList2 <- pQueue2List(x@flowFramesPreProcessingQueue)
    outputList[["flowFramesPreProcessingSteps"]] <- pList2

    return(outputList)
}

##' @rdname CytoPipeline
##'
##' @param x a `CytoPipeline` object
##'
##' @export
##'
experimentName <- function(x) {
    stopifnot(inherits(x, "CytoPipeline"))
    return(x@experimentName)
}

##' @rdname CytoPipeline
##' @param x a `CytoPipeline` object
##' @param value the new value to be assigned
##' @export
##'
"experimentName<-" <- function(x, value) {
    stopifnot(inherits(x, "CytoPipeline"))
    x@experimentName <- value
    return(x)
}

##' @rdname CytoPipeline
##' @param x a `CytoPipeline` object
##' @export
##'
sampleFiles <- function(x) {
    stopifnot(inherits(x, "CytoPipeline"))
    return(x@sampleFiles)
}

##' @rdname CytoPipeline
##' @param x a `CytoPipeline` object
##' @param value the new value to be assigned
##' @export
##'
"sampleFiles<-" <- function(x, value) {
    stopifnot(inherits(x, "CytoPipeline"))
    x@sampleFiles <- value
    return(x)
}

##' @rdname CytoPipeline
##' @param x a `CytoPipeline` object
##' @export
##'
pData <- function(x) {
    stopifnot(inherits(x, "CytoPipeline"))
    return(x@pData)
}

##' @rdname CytoPipeline
##' @param x a `CytoPipeline` object
##' @param value the new value to be assigned
##' @export
##'
"pData<-" <- function(x, value) {
    stopifnot(inherits(x, "CytoPipeline"))
    x@pData <- value
    if (isTRUE(methods::validObject(x))) return(x)
}

.validProcessingQueue <- function(x, queueName) {
    msg <- NULL
    if (length(x) && !all(vapply(x, inherits, "CytoProcessingStep"))) {
        msg <- paste0(
            "'",
            queueName,
            "' should only contain CytoProcessingStep objects."
        )
    }
    return(msg)
}

.makeSlots <- function(x, params) {
    stopifnot(c(
        inherits(x, "CytoPipeline"),
        is.list(params)
    ))

    # experimentName no more explicitly mandatory from the params list
    # because populated by default (default_experiment)
    #mandatory <- c("experimentName")
    mandatory <- c() 
    #optional <- c("sampleFiles")
    optional <- c("experimentName",
                  "sampleFiles",
                  "pData")

    for (m in mandatory) {
        if (is.null(params[[m]])) {
            stop("No ", m, " provided")
        }

        if (length(params[[m]]) == 0) {
            methods::slot(x, m) <-
                vector(
                    mode = mode(methods::slot(x, m)),
                    length = 0
                )
        } else {
            methods::slot(x, m) <- params[[m]]
        }
    }

    for (o in optional) {
        if (!is.null(params[[o]])) {
            if (length(params[[o]]) == 0) {
                methods::slot(x, o) <-
                    vector(
                        mode = mode(methods::slot(x, o)),
                        length = 0
                    )
            } else {
                methods::slot(x, o) <- params[[o]]
            }
        }
    }

    return(x)
}

.makeProcessingQueues <- function(x, params) {
    stopifnot(c(
        inherits(x, "CytoPipeline"),
        is.list(params)
    ))

    # browser()

    # generate scale transform processing queue
    if (!is.null(params$scaleTransformProcessingSteps)) {
        for (s in seq_along(params$scaleTransformProcessingSteps)) {
            prStep <- params$scaleTransformProcessingSteps[[s]]
            if (is.null(prStep$name) || is.null(prStep$FUN) ||
                is.null(prStep$ARGS)) {
                stop(
                    "step ", prStep,
                    " not well defined in input data, needs 'name',",
                    "'FUN' and 'ARGS'"
                )
            }
            x <- addProcessingStep(x,
                whichQueue = "scale transform",
                CytoProcessingStep(
                    name = prStep$name,
                    FUN = prStep$FUN,
                    ARGS = prStep$ARGS
                )
            )
        }
    }

    # generate flow frame pre-processing queue
    if (!is.null(params$flowFramesPreProcessingSteps)) {
        for (s in seq_along(params$flowFramesPreProcessingSteps)) {
            prStep <- params$flowFramesPreProcessingSteps[[s]]
            if (is.null(prStep$name) || is.null(prStep$FUN) ||
                is.null(prStep$ARGS)) {
                stop(
                    "step ", prStep,
                    " not well defined in input data, needs 'name',",
                    "'FUN' and 'ARGS'"
                )
            }
            x <- addProcessingStep(x,
                whichQueue = "pre-processing",
                CytoProcessingStep(
                    name = prStep$name,
                    FUN = prStep$FUN,
                    ARGS = prStep$ARGS
                )
            )
        }
    }
    return(x)
}
