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


##' @title CytoPipeline class
##'
##' @aliases CytoPipeline-class, CytoPipeline
##'
##' @name CytoPipeline-class
##'
##' @rdname CytoPipeline
##'
##' @description
##'
##' Class representing a flow cytometry pipeline, and composed of two processing
##' queues, i.e. lists of CytoProcessingStep objects :
##' - a list of CytoProcessingStep(s) for pre-calculation of scale 
##' transformations per channel
##' - a list of CytoProcessingStep(s) for the pre-processing of flow frames
##'
##' @slot scaleTransformProcessingQueue A `list()` of CytoProcessingStep objects
##' containing the steps for obtaining the scale transformations per channel
##'
##' @slot flowFramesPreProcessingQueue A `list()` of CytoProcessingStep objects
##' containing the steps for pre-processing of the samples flow frames
##' 
##' @slot experimentName A `character()` containing the experiment (run) name
##'
##' @slot sampleFiles A `character()` vector storing all fcs files to be
##' run into the pipeline
##'
# @slot savePreprocessedFiles if TRUE, will save pre-processed fcs file in
# /QC subdirectory
#
# @slot savePlotsInFiles if TRUE, will save files corresponding to default
# generated plot at each pre-processing step
#
##' @slot saveScaleTransform if TRUE, will save rds object storing scale
# transformation list generated
#
##' @slot scaleTransformFile basename of the file to use to save the scale
# transformation list (if 'saveScaleTransform' == TRUE)
#
##' @exportClass CytoPipeline
setClass("CytoPipeline",
         slots = c(experimentName = "character",
                   scaleTransformProcessingQueue = "list",
                   flowFramesPreProcessingQueue = "list",
                   sampleFiles = "character",
                   #savePreprocessedFiles = "logical",
                   #savePlotsInFiles = "logical",
                   saveScaleTransform = "logical",
                   scaleTransformFile = "character"
                   ),
         prototype = list(experimentName = "default_experiment",
                          scaleTransformProcessingQueue = list(),
                          flowFramesPreProcessingQueue = list(),
                          sampleFiles = character(),
                          #savePreprocessedFiles = FALSE,
                          #savePlotsInFiles = FALSE,
                          saveScaleTransform = FALSE,
                          scaleTransformFile = character()
                          )
         )

setValidity("CytoPipeline", function(object) {
  msg1 <- .validProcessingQueue(object@scaleTransformProcessingQueue,
                                  "scaleTransformProcessingQueue")
  msg2 <- .validProcessingQueue(object@flowFramesPreProcessingQueue,
                                  "flowFramesPreProcessingQueue")

  msg <- NULL
  if (length(msg1)){
    msg <- msg1
    if (length(msg2)){
      msg <- paste0(msg, " ; ", msg2)
    }
  } else {
    if (length(msg2)){
      msg <- msg2
    }
  }
  if (length(msg)) msg
  else TRUE
})

#' @rdname CytoPipeline
#' 
#' @param object a `CytoPipeline` object
#'
#' @importMethodsFrom methods show
#'
setMethod("show", "CytoPipeline",
          function(object) {
            cat("Pipeline object for flow cytometry experiment: ",
                object@experimentName, "\n")
            if (length(object@sampleFiles)) {
              cat("Sample files:", length(object@sampleFiles),
                  " sample file(s)\n")
            } else {
              cat("No sample file\n")
            }
            showProcessingSteps(object, whichQueue = "scale transform")
            showProcessingSteps(object, whichQueue = "pre-processing")
            cat("Other slots:\n")
            # cat("savePreprocessedFiles: ", object@savePreprocessedFiles, "\n")
            # cat("savePlotsInFiles: ", object@savePlotsInFiles, "\n")
            cat("saveScaleTransform: ", object@saveScaleTransform, "\n")
            cat("scaleTransformFile: ", object@scaleTransformFile, "\n")
          })


setGeneric("CytoPipeline", function(object, ...) 
  standardGeneric("CytoPipeline"))

#' @rdname CytoPipeline
#' @param experimentName the experiment name
#' @param sampleFiles the sample files
#' @export
#'
setMethod("CytoPipeline", "missing",
          function(object,
                   experimentName = "default_experiment",
                   sampleFiles = character()) {
            methods::new("CytoPipeline",
                         experimentName = experimentName,
                         scaleTransformProcessingQueue = list(),
                         flowFramesPreProcessingQueue = list(),
                         sampleFiles = sampleFiles)
          }
        )

#' @rdname CytoPipeline
#' 
#' @param object a `list()`
#'
#' @export
#'
setMethod("CytoPipeline", "list",
          function(object) {
            x <- methods::new("CytoPipeline",
                             experimentName = "default_experiment",
                             scaleTransformProcessingQueue = list(),
                             flowFramesPreProcessingQueue = list(),
                             sampleFiles = character())
            x <- .makeSlots(x, object)
            x <- .makeProcessingQueues(x, object)
            return(x)
          }
        )

#' @rdname CytoPipeline
#' 
#' 
#' @param object a `character()` containing a JSON input
#'
#' @export
#'
setMethod("CytoPipeline", "character",
          function(object) {
            pipelineParams <- jsonlite::read_json(object,
                                                  simplifyVector = TRUE,
                                                  simplifyDataFrame = FALSE)
            cytoPipeline <- CytoPipeline(pipelineParams)
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
as.list.CytoPipeline <- function(x,...){
  stopifnot(inherits(x, "CytoPipeline"))
  
  #browser()
  
  slots <- methods::slotNames(x)
  slots <- slots[!slots %in% c("scaleTransformProcessingQueue",
                               "flowFramesPreProcessingQueue")]
  
  outputList <- list()
  
  for (sl in slots){
    outputList[[sl]] <- methods::slot(x, sl)
  }
  
  # now handles both processing queue
  pQueue2List <- function(pQueue){
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
"experimentName<-" <- function(x, value){
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
"sampleFiles<-" <- function(x, value){
  stopifnot(inherits(x, "CytoPipeline"))
  x@sampleFiles <- value
  return(x)
}

.validProcessingQueue <- function(x, queueName) {
  msg <- NULL
  if (length(x) && !all(vapply(x, inherits, "CytoProcessingStep"))) {
    msg <- paste0("'",
                  queueName,
                  "' should only contain CytoProcessingStep objects.")
    
  }
  return(msg)
}

.makeSlots <- function(x, params) {
  stopifnot(c(inherits(x, "CytoPipeline"),
              is.list(params)))
  
  mandatory <- c("experimentName", "sampleFiles")
  optional <- c("saveScaleTransform", "scaleTransformFile",
                "", "")
  
  for (m in mandatory) {
    if (is.null(params[[m]]))
      stop("No ", m, " provided")
    
    if (length(params[[m]]) == 0) {
      methods::slot(x, m) <-
        vector(mode = mode(methods::slot(x, m)),
               length = 0)
    } else {
      methods::slot(x, m) <- params[[m]]
    }
  }
  
  for (o in optional) {
    if (!is.null(params[[o]])) {
      if (length(params[[o]]) == 0) {
        methods::slot(x, o) <-
          vector(mode = mode(methods::slot(x, o)),
                 length = 0)
      } else {
        methods::slot(x, o) <- params[[o]]
      }
      
    }
  }
  
  return(x)
}

.makeProcessingQueues <- function(x, params) {
  stopifnot(c(inherits(x, "CytoPipeline"),
              is.list(params)))
  
  #browser()
  
  # generate scale transform processing queue
  if (!is.null(params$scaleTransformProcessingSteps)) {
    for (s in seq_along(params$scaleTransformProcessingSteps)) {
      prStep <- params$scaleTransformProcessingSteps[[s]]
      if (is.null(prStep$name) || is.null(prStep$FUN) || is.null(prStep$ARGS)) {
        stop("step ", prStep, " not well defined in input data, needs 'name',",
             "'FUN' and 'ARGS'")
      }
      x <- addProcessingStep(x,
                             whichQueue = "scale transform",
                             CytoProcessingStep(name = prStep$name,
                                                FUN = prStep$FUN,
                                                ARGS = prStep$ARGS))
    }
  }
  
  # generate flow frame pre-processing queue
  if (!is.null(params$flowFramesPreProcessingSteps)) {
    for (s in seq_along(params$flowFramesPreProcessingSteps)) {
      prStep <- params$flowFramesPreProcessingSteps[[s]]
      if (is.null(prStep$name) || is.null(prStep$FUN) || is.null(prStep$ARGS)) {
        stop("step ", prStep, " not well defined in input data, needs 'name',",
             "'FUN' and 'ARGS'")
      }
      x <- addProcessingStep(x,
                             whichQueue = "pre-processing",
                             CytoProcessingStep(name = prStep$name,
                                                FUN = prStep$FUN,
                                                ARGS = prStep$ARGS))
    }
  }
  return(x)
}

