
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
##' - a list of CytoProcessingStep(s) for pre-calculation of scale transformations
##' per channel
##' - a list of CytoProcessingStep(s) for the pre-processing of flow frames
##'
##' @slot scaleTransformProcessingQueue A `list()` of CytoProcessingStep objects
##' containing the steps for obtaining the scale transformations per channel
##'
##' @slot flowFramesPreProcessingQueue A `list()` of CytoProcessingStep objects
##' containing the steps for pre-processing of the samples flow frames
##'
##' @slot sampleFiles A `character()` vector storing all fcs files to be
##' run into the pipeline
##'
##' @slot savePreprocessedFiles if TRUE, will save pre-processed fcs file in
##' /QC subdirectory
##'
##' @slot savePlotsInFiles if TRUE, will save files corresponding to default
##' generated plot at each pre-processing step
##'
##' @slot saveScaleTransform if TRUE, will save rds object storing scale
##' transformation list generated
##'
##' @slot scaleTransformFile basename of the file to use to save the scale
##' transformation list (if 'saveScaleTransform' == TRUE)
##'
##' @exportClass CytoPipeline
setClass("CytoPipeline",
         slots = c(experimentName = "character",
                   scaleTransformProcessingQueue = "list",
                   flowFramesPreProcessingQueue = "list",
                   sampleFiles = "character",
                   savePreprocessedFiles = "logical",
                   savePlotsInFiles = "logical",
                   saveScaleTransform = "logical",
                   scaleTransformFile = "character"),
         prototype = list(experimentName = "default_experiment",
                          scaleTransformProcessingQueue = list(),
                          flowFramesPreProcessingQueue = list(),
                          sampleFiles = character(),
                          savePreprocessedFiles = FALSE,
                          savePlotsInFiles = FALSE,
                          saveScaleTransform = FALSE,
                          scaleTransformFile = character())
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
#' @importMethodsFrom methods show
#'
setMethod("show", "CytoPipeline",
          function(object) {
            cat("Pipeline object for flow cytometry experiment: ",
                object@experimentName, "\n")
            if (length(object@sampleFiles)) {
              cat("Sample files:", length(object@sampleFiles),
                  " sample file(s)\n")
            }
            showProcessingSteps(object, whichQueue = "scale transform")
            showProcessingSteps(object, whichQueue = "pre-processing")
            cat("Other slots:\n")
            cat("savePreprocessedFiles: ", object@savePreprocessedFiles, "\n")
            cat("savePlotsInFiles: ", object@savePlotsInFiles, "\n")
            cat("saveScaleTransform: ", object@saveScaleTransform, "\n")
            cat("scaleTransformFile: ", object@scaleTransformFile, "\n")
          })


setGeneric("CytoPipeline", function(object, ...) standardGeneric("CytoPipeline"))

#' @rdname CytoPipeline
#'
#' @export
#'
setMethod("CytoPipeline", "missing",
          function(object,
                   experimentName = "default_experiment",
                   sampleFiles = character()) {
            new("CytoPipeline",
                experimentName = experimentName,
                scaleTransformProcessingQueue = list(),
                flowFramesPreProcessingQueue = list(),
                sampleFiles = sampleFiles)
          }
        )

#' @rdname CytoPipeline
#'
#' @export
#'
setMethod("CytoPipeline", "list",
          function(object) {
            x = new("CytoPipeline",
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
#' @export
#'
setMethod("CytoPipeline", "character",
          function(object) {
            pipelineParams <- jsonlite::read_json(object,
                                                  simplifyVector = TRUE,
                                                  simplifyDataFrame = FALSE)
            cytoPipeline = CytoPipeline(pipelineParams)
            return(cytoPipeline)
          }
)
