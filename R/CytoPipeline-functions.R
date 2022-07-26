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

##' @name handlingProcessingSteps
##' @title handling processing steps in CytoPipeline objects
##' @description functions to manipulate processing steps in processing queues 
##' of CytoPipeline objects
##' @param x a CytoPipeline object
##' @param whichQueue selects the processing queue for which we manage the processing 
##' steps
NULL

##' @title handling processing steps in CytoPipeline objects
##' @param newPS the new processing step to be added (CytoProcessingStep object)
##' @describeIn handlingProcessingSteps adds a processing step in one of the processing queues (at the
##' end), returns the modified CytoPipeline object
##'
##' @export
#'
addProcessingStep <- function(x,
                              whichQueue = c("scale transform",
                                             "pre-processing"),
                              newPS){
  stopifnot(c(inherits(x, "CytoPipeline"),
              inherits(newPS, "CytoProcessingStep")))
  whichQueue <- match.arg(whichQueue)

  if (whichQueue == "scale transform"){
    ll <- length(x@scaleTransformProcessingQueue)
    for (j in seq_along(x@scaleTransformProcessingQueue)) {
      pS <- getProcessingStep(x, whichQueue = whichQueue, index = j)
      if (pS@name == newPS@name) {
        stop("There already exist a step with name '", newPS@name,
             "' in ", whichQueue, " queue!")
      }
    }
    x@scaleTransformProcessingQueue[[ll+1]] <- newPS
  } else {
    ll <- length(x@flowFramesPreProcessingQueue)
    for (j in seq_along(x@flowFramesPreProcessingQueue)) {
      pS <- getProcessingStep(x, whichQueue = whichQueue, index = j)
      if (pS@name == newPS@name) {
        stop("There already exist a step with name '", newPS@name,
             "' in ", whichQueue, " queue!")
      }
    }
    x@flowFramesPreProcessingQueue[[ll+1]] <- newPS
  }
  return(x)
}




##' @param index index of the processing step to remove
##' @describeIn handlingProcessingSteps removes a processing step from one of the
##' processing queues, returns the modified CytoPipeline object
##' @export
##'
removeProcessingStep <- function(x,
                                 whichQueue = c("scale transform",
                                                "pre-processing"),
                                index) {
  if(!is.numeric(index)){
    stop("index should be a numeric")
  }
  stopifnot(inherits(x, "CytoPipeline"))
  whichQueue <- match.arg(whichQueue)
  if (whichQueue == "scale transform"){
    ll <- length(x@scaleTransformProcessingQueue)
    if(index<=0 || index > ll){
      stop("index out of bound")
    }
    x@scaleTransformProcessingQueue <-
      x@scaleTransformProcessingQueue[-index]
  } else {
    ll <- length(x@flowFramesPreProcessingQueue)
    if(index<=0 || index > ll){
      stop("index out of bound")
    }
    x@flowFramesPreProcessingQueue <-
      x@flowFramesPreProcessingQueue[-index]
  }
  return(x)
}

##' @describeIn handlingProcessingSteps gets the number of processing
##' steps in a processing queue
##' @export
##'
getNbProcessingSteps <- function(x,
                                 whichQueue = c("scale transform",
                                               "pre-processing")) {
  stopifnot(inherits(x, "CytoPipeline"))
  whichQueue <- match.arg(whichQueue)

  if (whichQueue == "scale transform"){
    return(length(x@scaleTransformProcessingQueue))
  } else {
    return(length(x@flowFramesPreProcessingQueue))
  }
}

##' @describeIn handlingProcessingSteps gets a processing step at a
##' specific index of a processing queue
##' @export
##'
getProcessingStep <- function(x,
                              whichQueue = c("scale transform",
                                             "pre-processing"),
                              index) {
  if(!is.numeric(index)){
    stop("index should be a numeric")
  }
  stopifnot(inherits(x, "CytoPipeline"))
  whichQueue <- match.arg(whichQueue)
  if (whichQueue == "scale transform"){
    ll <- length(x@scaleTransformProcessingQueue)
    if(index<=0 || index > ll){
      stop("index out of bound")
    }
    return(x@scaleTransformProcessingQueue[[index]])
  } else {
    ll <- length(x@flowFramesPreProcessingQueue)
    if(index<=0 || index > ll){
      stop("index out of bound")
    }
    return(x@flowFramesPreProcessingQueue[[index]])
  }
  return(x)
}

##' @describeIn handlingProcessingSteps gets a character vector of all
##' processing step names of a specific processing queue
##' @export
##'
getProcessingStepNames <- function(x, whichQueue = c("scale transform",
                                                     "pre-processing")) {
  nSteps <- getNbProcessingSteps(x, whichQueue)
  stepNames <- character()
  if(nSteps > 0) {
    stepNames <- sapply(1:nSteps,
                        FUN = function(i) {
                          pS <- getProcessingStep(x, whichQueue, i)
                          pS@name
                        })
  }

  return(stepNames)
}

##' @describeIn handlingProcessingSteps deletes all processing steps in one
##' or both processing queues, returns the modified CytoPipeline object
##' @export
##'
cleanProcessingSteps <- function(x,
                                 whichQueue = c("both", "scale transform",
                                                "pre-processing")) {
  stopifnot(inherits(x, "CytoPipeline"))
  whichQueue <- match.arg(whichQueue)
  if (whichQueue == "scale transform" || whichQueue == "both"){
    x@scaleTransformProcessingQueue <- list()
  }
  if (whichQueue == "pre-processing" || whichQueue == "both"){
    x@flowFramesPreProcessingQueue <- list()
  }
  return(x)
}

##' @describeIn handlingProcessingSteps shows all processing steps in a
##' processing queue
##' @export
##'
showProcessingSteps <- function(x,
                                whichQueue = c("scale transform",
                                               "pre-processing")) {
  stopifnot(inherits(x, "CytoPipeline"))
  whichQueue <- match.arg(whichQueue)

  if (whichQueue == "scale transform") {
    queue <- x@scaleTransformProcessingQueue
    queueName <- "Scale transformations evaluation queue"
  } else {
    queue <- x@flowFramesPreProcessingQueue
    queueName <- "Flow frames pre-processing evaluation queue"
  }

  if (length(queue)) {
    cat(queueName,": ", length(queue), " processing step(s)\n")
    for (i in seq_along(queue)) {
      ps <- queue[[i]]
      cat(i, ": ")
      show(ps)
    }
  } else {
    cat(queueName, " has no processing step\n")
  }

}


#' @title executing CytoPipeline object
#' @description this function triggers the execution of the processing queues of
#' a CytoPipeline object. 
#' First, the scale tranform processing queue is run, taking
#' the set of sample names as an implicit first input. At the end of the queue,
#' a scale transform List is assumed to be created.
#' Second, the flowFrame pre-processing queue, reapeatedly for each sample file.
#' The scale transform list generated in the previous step is taken as implicit 
#' input, together with the initial sample file. At the end of the queue run, a
#' pre-processed flowFrame is assumed to be generated.
#' No change is made on the input CytoPipeline object, all results are stored in
#' the cache.
#' @param x CytoPipeline object
#' @param path base path, a subdirectory with name equal to the experiment will
#' be created to store the output data, in particular the experiment cache
#' @param rmCache if TRUE, starts by removing the already existing cache
#' directory corresponding to the experiment
#' @returns nothing
#' @export
execute <- function(x,
                    path = ".",
                    rmCache = FALSE) {

  stopifnot(inherits(x, "CytoPipeline"))

  #browser()

  outputDir <- paste0(path, "/", x@experimentName, "/output/")

  qualityControlDir <- paste0(outputDir, "QC/")
  rdsOutputDir <- paste0(outputDir, "RDS/")

  for (newPath in c(qualityControlDir, rdsOutputDir)) {
    if (!dir.exists(newPath)) {
      dir.create(newPath, recursive = TRUE)
    }
  }

  # create or use local cache
  localCacheDir <- paste0(path, "/", x@experimentName, "/.cache")
  bfc <- BiocFileCache::BiocFileCache(localCacheDir, ask = FALSE)

  # remove cache if 'rmCache' is set to
  if (rmCache) {
    bfci <- BiocFileCache::bfcinfo(bfc)
    if(nrow(bfci) > 0) {
      warning("Found a previous cache with experiment name ",
              x@experimentName, " => removed it")
      BiocFileCache::removebfc(bfc, ask = FALSE)
      # recreate cache
      bfc <- BiocFileCache::BiocFileCache(localCacheDir, ask = FALSE)
    }
  }

  # checking consistency between to-be-run processing steps and cache
  res <- checkCytoPipelineConsistencyWithCache(x, path = path)
  if (!res$isConsistent)
    stop(res$inconsistencyMsg)

  ### first part is always to compute common transformation list ###
  message("#####################################################")
  message("### running SCALE TRANSFORMATION processing steps ###")
  message("#####################################################")

  currentTransList <- NULL

  for (s in seq_along(x@scaleTransformProcessingQueue)) {
    cacheResourceName <- getName(x@scaleTransformProcessingQueue[[s]])
    stepName <- getName(x@scaleTransformProcessingQueue[[s]])
    cacheResourceName = paste0("scaleTransform_step", s, "_",
                               stepName)

    msg <- paste0("Proceeding with step ", s, " [", stepName, "]")
    if (cacheResourceName %in% BiocFileCache::bfcinfo(bfc)$rname) {
      message(msg, ": found in cache!")
      cacheResourceFile <- BiocFileCache::bfcrpath(x = bfc, rnames = cacheResourceName)
      res <- readRDS(cacheResourceFile)
    } else {
      message(msg, " ...")
      #browser()
      if (s==1){
        # try runnning with adding sampleFiles parameters
        # if complain => run with
        res <- try(executeProcessingStep(
          x@scaleTransformProcessingQueue[[s]],
          sampleFiles = x@sampleFiles),
          silent = TRUE)
        if (methods::is(res, "try-error"))
          res <-
            executeProcessingStep(
              x@scaleTransformProcessingQueue[[s]])
      } else {
        res <-
          executeProcessingStep(
            x@scaleTransformProcessingQueue[[s]],
            res)
      }
      # store file in cache
      cacheResourceFile = BiocFileCache::bfcnew(bfc, cacheResourceName)
      saveRDS(res, unname(cacheResourceFile))

      # add one entry in cache meta data tables
      outputClass <- class(res)
      outputObjectName <- paste0(stepName, "_obj")
      stepJsonSerialize <-
        as.json.CytoProcessingStep(x@scaleTransformProcessingQueue[[s]])
      genericMeta <- data.frame(list(rid = names(cacheResourceFile),
                                  type = "scale transform",
                                  stepNb = s,
                                  stepName = stepName,
                                  stepJsonSerialize = stepJsonSerialize,
                                  outputClass = outputClass,
                                  outputObjectName = outputObjectName))

      BiocFileCache::bfcmeta(bfc, name = "generic", append = TRUE) <- genericMeta
    }
  } # end loop on steps

  # store transformList for use in flowFrames pre-processing steps
  if(inherits(res, "transformList")) {
    currentTransList <- res
    if (x@saveScaleTransform) {
      if (!length(x@scaleTransformFile))
        stop("saving scale tranformations require a 'scaleTranformFile' to be set")
      saveRDS(currentTransList,
              file = paste0(rdsOutputDir,
                            x@scaleTransformFile))
    }
  }


  for (file in x@sampleFiles) {
    #browser()
    message("#####################################################")
    message(paste0("### NOW PRE-PROCESSING FILE ", file, "..."))
    message("#####################################################")

    for (s in seq_along(x@flowFramesPreProcessingQueue)) {
      stepName <- getName(x@flowFramesPreProcessingQueue[[s]])
      cacheResourceName = paste0("preprocessing_",
                                 basename(file),
                                 "_step", s, "_",
                                 stepName)

      msg <- paste0("Proceeding with step ", s, " [", stepName, "]")

      if (cacheResourceName %in% BiocFileCache::bfcinfo(bfc)$rname) {
        message(msg, ": found in cache!")
        cacheResourceFile <- BiocFileCache::bfcrpath(x = bfc,
                                                     rnames = cacheResourceName)
        res <- readRDS(cacheResourceFile)
      } else {
        message(msg, " ...")
        #browser()
        if (s == 1) {
          res <-
            executeProcessingStep(
              x@flowFramesPreProcessingQueue[[s]],
              sampleFiles = file,
              transList = currentTransList)
        } else {
          res <-
            executeProcessingStep(
              x@flowFramesPreProcessingQueue[[s]],
              #ff = res,
              res,
              transList = currentTransList
            )
        }
         # store file in cache
        cacheResourceFile = BiocFileCache::bfcnew(bfc, cacheResourceName)
        saveRDS(res, unname(cacheResourceFile))

        # add one entry in cache meta data tables
        outputClass <- class(res)
        outputObjectName <- paste0(stepName, "_obj")
        stepJsonSerialize <-
          as.json.CytoProcessingStep(x@flowFramesPreProcessingQueue[[s]])
        genericMeta <- data.frame(list(rid = names(cacheResourceFile),
                                    type = "pre-processing",
                                    stepNb = s,
                                    stepName = stepName,
                                    stepJsonSerialize = stepJsonSerialize,
                                    outputClass = outputClass,
                                    outputObjectName = outputObjectName
                                    ))
        preprocessingMeta <- data.frame(list(rid = names(cacheResourceFile),
                                             fcsfile = basename(file)))

        BiocFileCache::bfcmeta(bfc, name = "generic", append = TRUE) <- genericMeta
        BiocFileCache::bfcmeta(bfc, name = "preprocessing", append = TRUE) <-
          preprocessingMeta

      } # if (newResource)
    } # end loop on steps
  } # end loop on sample files
}

##' @name interactingWithCytoPipelineCache
##' @title interaction between CytoPipeline object and disk cache
##' @description functions supporting the interaction between a CytoPipeline object
##' and the file cache on disk
##' @param x a CytoPipeline object
##' @param path the full path to the experiment storage on disk 
##' (without the /.cache)
##' @param experimentName the experimentName used to select the file cache on
##' disk
##' @return 
##' for `deleteCytoPipelineCache`: nothing\cr  
##' for `buildCytoPipelineCache`: the built CytoPipeline object\cr  
##' for `checkCytoPipelineConsistencyWithCache`: a list with the following values:  
##' - `isConsistent` (TRUE/FALSE)
##' - `inconsistencyMsg`: character filled in by an inconsistency message in case the 
##' cache and CytoPipeline object are not consistent with each other
##' - `scaleTransformStepStatus`: a character vector, containing, for each scale
##' transform step, a status from c("run", "not run", "inconsistent")
##' - `preProcessingStepStatus`: a character matrix, containing, for each
##' pre-processing step (rows), for each sample file (columns), a status from
##' c("run", "not run", "inconsistent")
##

NULL


##' @describeIn interactingWithCytoPipelineCache delete the whole disk cache 
##' corresponding to the experiment of a CytoPipeline object
##' @export
##'
deleteCytoPipelineCache <- function(x,  path = "."){
  stopifnot(inherits(x, "CytoPipeline"))
  localCacheDir <- paste0(path, "/", x@experimentName, "/.cache")
  bfc <- BiocFileCache::BiocFileCache(localCacheDir, ask = FALSE)
  BiocFileCache::removebfc(bfc, ask = FALSE)
}


##' @describeIn interactingWithCytoPipelineCache builds a new CytoPipeline
##' object, based on the information stored in the file cache
##' @export
##'
buildCytoPipelineFromCache <- function(experimentName, path = "."){

  #browser()

  x <- CytoPipeline(experimentName = experimentName)

  # find cache corresponding to experiment name
  experimentNames <-
    getCytoPipelineExperimentNames(path = path,
                                   pattern = x@experimentName,
                                   fixed = TRUE)
  if (length(experimentNames) == 0) {
    warning("no cache directory found for [", x@experimentName,
         "], did not build any processing step")
    return(x)
  } else {
    cacheDir <- paste0(path, "/", x@experimentName, "/.cache")

    bfc <- BiocFileCache::BiocFileCache(cacheDir, ask = FALSE)
    cacheInfo <- BiocFileCache::bfcinfo(bfc)

    # now building scale transform processing queue
    # take only steps with scale transform
    stepsInfos <-
      cacheInfo[cacheInfo$type == "scale transform",
                c("stepNb", "stepName", "stepJsonSerialize")]
    stepsInfos <- unique(stepsInfos)
    nSteps <- nrow(stepsInfos)
    if (nSteps > 0) {
      uniqueStepNbs <- unique(stepsInfos$stepNb)

      if(length(uniqueStepNbs) != nrow(stepsInfos))
        stop("more than one step having the same step nb. Cache is inconsistent,",
             " deleting it manually is advised")

      stepsInfos <- stepsInfos[order(stepsInfos$stepNb),]

      #x <- cleanProcessingSteps(x, whichQueue = "scale transform")

      for (j in 1:nSteps) {
        pS <- from.json.CytoProcessingStep(
          as.character(stepsInfos[j, "stepJsonSerialize"]))
        x <- addProcessingStep(x, whichQueue = "scale transform", pS)
      }
    }

    # now building pre-processing queue
    # take unique steps with pre-processing type

    stepsInfos <-
      cacheInfo[cacheInfo$type == "pre-processing",
                c("stepNb", "stepName", "stepJsonSerialize")]
    stepsInfos <- unique(stepsInfos)
    nSteps <- nrow(stepsInfos)
    if (nSteps > 0) {
      uniqueStepNbs <- unique(stepsInfos$stepNb)
      if(length(uniqueStepNbs) != nrow(stepsInfos))
        stop("more than one step having the same step nb. Cache is inconsistent,",
             " deleting it manually is advised")
      stepsInfos <- stepsInfos[order(stepsInfos$stepNb),]
      #x <- cleanProcessingSteps(x, whichQueue = "pre-processing")

      for (j in 1:nSteps) {
        pS <- from.json.CytoProcessingStep(
          as.character(stepsInfos[j, "stepJsonSerialize"]))
        x <- addProcessingStep(x, whichQueue = "pre-processing", pS)
      }

      # updating sample files
      x@sampleFiles <- sort(unique(stats::na.omit(cacheInfo$fcsfile)))
    } else {
      # no pre-processing step found in cache => no sample file can be updated
      x@sampleFiles <- character()
    }

    return(x)
  }

}

##' @describeIn interactingWithCytoPipelineCache check the consistency
##' between the processing steps described in a CytoPipeline object, 
##' and what is stored in the file cache
##' @export
##'
checkCytoPipelineConsistencyWithCache <- function(x, path = ".") {
  stopifnot(inherits(x, "CytoPipeline"))

  #browser()
  ret <- list(isConsistent = TRUE, inconsistencyMsg = character(0))

  nScaleTransformSteps <- length(x@scaleTransformProcessingQueue)
  ret$scaleTransformStepStatus = rep("not_run", nScaleTransformSteps)
  ret$scaleTransformStepOutputObjNames =
    rep("unknown", nScaleTransformSteps)
  ret$scaleTransformStepOutputClasses =
    rep("unknown", nScaleTransformSteps)

  if (nScaleTransformSteps >0) {
    names(ret$scaleTransformStepStatus) <-
      getProcessingStepNames(x, whichQueue = "scale transform")
  }

  nPreProcessingSteps <- length(x@flowFramesPreProcessingQueue)
  nSampleFiles <- length(x@sampleFiles)
  ret$preProcessingStepStatus =
    matrix(rep("not_run", nPreProcessingSteps * nSampleFiles),
           nrow = nPreProcessingSteps,
           ncol = nSampleFiles)
  ret$preProcessingStepOutputObjNames =
    rep("unknown", nPreProcessingSteps)
  ret$preProcessingStepOutputClasses =
    rep("unknown", nPreProcessingSteps)

  if (nPreProcessingSteps > 0 && nSampleFiles > 0) {
    rownames(ret$preProcessingStepStatus) <-
      getProcessingStepNames(x, whichQueue = "pre-processing")
    colnames(ret$preProcessingStepStatus) <- basename(x@sampleFiles)
  }

  # find cache corresponding to experiment name
  experimentNames <-
    getCytoPipelineExperimentNames(path = path,
                                   pattern = x@experimentName,
                                   fixed = TRUE)
  if (length(experimentNames) == 0)
    return(ret)

  # cache found => check consistency
  cacheDir <- paste0(path, "/", x@experimentName, "/.cache")

  bfc <- BiocFileCache::BiocFileCache(cacheDir, ask = FALSE)
  cacheInfo <- BiocFileCache::bfcinfo(bfc)

  if(nrow(cacheInfo) == 0) {
    return(ret)
  }

  # now check scale transform steps consistency
  # take only steps with scale transform type
  stepsInfos <-
      cacheInfo[cacheInfo$type == "scale transform",
                c("stepNb", "stepName", "stepJsonSerialize",
                  "outputClass", "outputObjectName")]

  stepsInfos <- unique(stepsInfos)
  nStepsInCache <- nrow(stepsInfos)
  if (nStepsInCache > 0) {
    uniqueStepNbs <- unique(stepsInfos$stepNb)

    if(length(uniqueStepNbs) != nrow(stepsInfos))
      stop("more than one step having the same step nb. Cache is inconsistent,",
           " deleting it manually is advised")

    stepsInfos <- stepsInfos[order(stepsInfos$stepNb),]
    if (nStepsInCache > nScaleTransformSteps) {
      ret$isConsistent <- FALSE
      ret$inconsistencyMsg <-
        "more scale transform steps in cache than in CytoPipeline object"
      return(ret)
    }

    for (j in 1:nStepsInCache) {
      pS <- from.json.CytoProcessingStep(
        as.character(stepsInfos[j, "stepJsonSerialize"]))
      pS2 <- getProcessingStep(x,
                               whichQueue = "scale transform",
                               index = j)
      if (identical(pS, pS2)) {
        ret$scaleTransformStepStatus[j] <- "run"
        ret$scaleTransformStepOutputClasses[j] <-
          as.character(stepsInfos[j, "outputClass"])
        ret$scaleTransformStepOutputObjNames[j] <-
          as.character(stepsInfos[j, "outputObjectName"])

        # and continue...
      } else {
        ret$scaleTransformStepStatus[j] <- "inconsistent"
        ret$isConsistent <- FALSE
        ret$isconsistencyMsg <-
          paste0("inconsistent scale transform step #", j,
                 " (different in cache)")
        return(ret)
      }
    }
  }

  #browser()

  # now check pre-processing steps consistency
  # take only steps with pre-processing type

  stepsInfos <-
    cacheInfo[cacheInfo$type == "pre-processing",
              c("stepNb", "stepName", "stepJsonSerialize",
                "outputClass", "outputObjectName")]
  stepsInfos <- unique(stepsInfos)
  nStepsInCache <- nrow(stepsInfos)
  if (nStepsInCache > 0) {
    uniqueStepNbs <- unique(stepsInfos$stepNb)

    if (length(uniqueStepNbs) != nrow(stepsInfos))
      stop("more than one step having the same step nb. Cache is inconsistent,",
           " deleting it manually is advised")

    stepsInfos <- stepsInfos[order(stepsInfos$stepNb),]
    if (nStepsInCache > nPreProcessingSteps) {
      ret$isConsistent <- FALSE
      ret$inconsistencyMsg <-
        "more pre-processing steps in cache than in CytoPipeline object"
      return(ret)
    }
  }




  if (nrow(cacheInfo[cacheInfo$type == "pre-processing",]) == 0) {
    return(ret)
  }

  for (s in seq_along(x@sampleFiles)) {
    # take only steps with the target sample file
    sampleFile <- basename(x@sampleFiles[s])
    stepsInfos <-
      cacheInfo[!is.na(cacheInfo$fcsfile) & cacheInfo$fcsfile == sampleFile,
                c("stepNb", "stepName", "stepJsonSerialize",
                  "outputClass", "outputObjectName")]
    stepsInfos <- unique(stepsInfos)
    nStepsInCache <- nrow(stepsInfos)
    if (nStepsInCache > 0) {
      uniqueStepNbs <- unique(stepsInfos$stepNb)

      if(length(uniqueStepNbs) != nrow(stepsInfos))
        stop("more than one step having the same step nb. Cache is inconsistent,",
             " deleting it manually is advised")

      stepsInfos <- stepsInfos[order(stepsInfos$stepNb),]

      for (j in 1:nStepsInCache) {
        pS <- from.json.CytoProcessingStep(
          as.character(stepsInfos[j, "stepJsonSerialize"]))
        pS2 <- getProcessingStep(x,
                                 whichQueue = "pre-processing",
                                 index = j)
        if (identical(pS, pS2)) {
          ret$preProcessingStepStatus[j, sampleFile] <- "run"
          ret$preProcessingStepOutputClasses[j] <-
            as.character(stepsInfos[j, "outputClass"])
          ret$preProcessingStepOutputObjNames[j] <-
            as.character(stepsInfos[j, "outputObjectName"])
          # and continue...
        } else {
          ret$preProcessingStepStatus[j, sampleFile] <- "inconsistent"
          ret$isConsistent <- FALSE
          ret$inconsistencyMsg <-
            paste0("inconsistent pre-processing step #", j,
                   " for sample file ", sampleFile,
                   " (different in cache)")
          return(ret)
        }
      }
    }
  } # end loop on sample files


  return(ret)

}

##' @name exportCytoPipeline
##' @title exporting CytoPipeline objects
##' @description functions to export CytoPipeline objects in various
##' formats
##' @param x a CytoPipeline object
##' @param path the full path to the name of the file to be created
NULL


##' 
##' @describeIn exportCytoPipeline exports a CytoPipeline object 
##' to a JSON file
##' @export
##'
export2JSONFile <- function(x, path){
  myList <- as.list.CytoPipeline(x)
  jsonlite::write_json(myList, path = path, pretty = TRUE, null = "null")
}


##' @title Find CytoPipeline experiments stored in a file cache
##' @description This function looks into a path for stored file caches
##' and gets the corresponding experiment names
##' @param path root path to locate the search for file caches
##' @param pattern optional pattern limiting the search for experiment
##' names
##' @param ignore.case (TRUE/FALSE) used in pattern matching (grepl)
##' @param fixed (TRUE/FALSE) used in pattern matching (grepl)
##' @return a vector of character containing found experiment names
##' 
##' @export
##'
getCytoPipelineExperimentNames <-
  function(path = ".",
           pattern = NULL,
           ignore.case = FALSE,
           fixed = FALSE) {

    subdirs <- list.dirs(path, full.names = FALSE, recursive = FALSE)
    if (!is.null(pattern)) {
      subdirs <- subdirs[grepl(pattern = pattern,
                               x = subdirs,
                               ignore.case = ignore.case,
                               fixed = fixed)]
    }
    experimentNames <- vector(mode = "character", length = 0)
    for (d in subdirs) {
      tentative <- paste0(path, "/", d, "/.cache")
      if (dir.exists(tentative)) {
        experimentNames <- c(experimentNames, d)
      }
    }
    return (experimentNames)
  }


##' @title Retrieves an object from a file cache
##' @description Given a CytoPipeline object, this function retrieves
##' a specific object in the corresponding file cache
##' @param x a CytoPipeline object
##' @param path root path to locate the search for file caches
##' @param whichQueue which queue to look into
##' @param sampleFile which sampleFile is looked for:
##' - if whichQueue == "scale transform", the sampleFile is ignored
##' - if NULL and whichQueue == "pre-processing", the sampleFile is
##' defaulted to the first one belonging to the experiment 
##' @param objectName (character) which object name to look for
##' @return the found object (or stops with an error message if the
##' target object is not found)
##' 
##' @export
##'
getCytoPipelineObjectFromCache <-
  function(x,
           path = ".",
           whichQueue = c("scale transform", "pre-processing"),
           sampleFile = NULL,
           objectName
           ) {
    #browser()
    stopifnot(inherits(x, "CytoPipeline"))
    whichQueue <- match.arg(whichQueue)
    if(whichQueue == "scale transform") {
      sampleFile <- NULL
    } else {
      if (is.numeric(sampleFile)) {
        if (sampleFile>0 && sampleFile <= length(x@sampleFiles)) {
          sampleFile <- x@sampleFiles[sampleFile]
        } else {
          stop("sampleFile out of bounds")
        }
      } else {
        sampleFileIndex <- which(basename(x@sampleFiles) == basename(sampleFile))
        if (length(sampleFileIndex) == 0) {
          stop("sampleFile not found in CytoPipeline")
        } else if (length(sampleFileIndex) > 1) {
          stop("sampleFile found multiple times in CytoPipeline ",
               "(unexpected inconsistency)")
        }
      }
    }

    # checking consistency between to-be-run processing steps and cache
    res <- checkCytoPipelineConsistencyWithCache(x, path = path)
    if (!res$isConsistent)
      stop(res$inconsistencyMsg)


    cacheDir <- paste0(path, "/", x@experimentName, "/.cache")

    bfc <- BiocFileCache::BiocFileCache(cacheDir, ask = FALSE)
    cacheInfo <- BiocFileCache::bfcinfo(bfc)
    if (!is.null(sampleFile)) {
      sampleFiles <- unique(stats::na.omit(cacheInfo$fcsfile))
      sampleFile <- basename(sampleFile)
      if (!(sampleFile %in% sampleFiles)) {
        stop(sampleFile, " not found in experiment run!")
      }
      indexInCache <- which(cacheInfo$outputObjectName == objectName &
                              cacheInfo$fcsfile == sampleFile)
      if (length(indexInCache) > 1) {
        stop("More than one line in cache corresponding to object name and sample",
             "file => unexpected inconsistency")
      }

      if(length(indexInCache) == 0) {
        stop("Combination (sample file, object name) not found in cache")
      }
    } else {
      if (!("fcsfile" %in% colnames(cacheInfo))) {
        # special case appearing when no pre-processing steps have been run
        indexInCache <- which(cacheInfo$outputObjectName == objectName)
      } else {
        indexInCache <- which(cacheInfo$outputObjectName == objectName &
                                is.na(cacheInfo$fcsfile))
      }
      
      if (length(indexInCache) > 1) {
        stop("More than one line in cache corresponding to object name and w/o",
             "sample file => unexpected inconsistency")
      }

      if(length(indexInCache) == 0) {
        stop("Object name not found in cache")
      }
    }

    cacheResourceFile <-
      BiocFileCache::bfcpath(bfc, rids = cacheInfo$rid[indexInCache])

    ret <- readRDS(cacheResourceFile)

    return(ret)
  }

##' internal function for the time being
##' @title File cache objects information
##' @description Given a CytoPipeline object, this function retrieves
##' the information related to a specific object name, i.e. 
##' object name and object class
##' @param x a CytoPipeline object
##' @param path root path to locate the search for file caches
##' @param whichQueue which queue to look into
##' @param sampleFile which sampleFile is looked for:
##' - if whichQueue == "scale transform", the sampleFile is ignored
##' - if NULL and whichQueue == "pre-processing", the sampleFile is
##' defaulted to the first one belonging to the experiment 
##' @return a dataframe with the collected information about the
##' found objects (or stops with an error message if no target object 
##' was found)
#
##' @export
# 
getCytoPipelineObjectInfos <-
  function(x,
           path = ".",
           whichQueue = c("scale transform", "pre-processing"),
           sampleFile = NULL) {
    #browser()
    stopifnot(inherits(x, "CytoPipeline"))
    whichQueue <- match.arg(whichQueue)
    if(whichQueue == "scale transform") {
      sampleFile <- NULL
    }

    # checking consistency between to-be-run processing steps and cache
    res <- checkCytoPipelineConsistencyWithCache(x, path = path)
    if (!res$isConsistent)
      stop(res$inconsistencyMsg)

    cacheDir <- paste0(path, "/", x@experimentName, "/.cache")

    bfc <- BiocFileCache::BiocFileCache(cacheDir, ask = FALSE)
    cacheInfo <- BiocFileCache::bfcinfo(bfc)
    indexesInCache <- numeric()
    if (!is.null(sampleFile)) {
      sampleFiles <- unique(stats::na.omit(cacheInfo$fcsfile))
      sampleFile <- basename(sampleFile)
      if (!(sampleFile %in% sampleFiles)) {
        stop(sampleFile, " not found in experiment run!")
      }
      indexesInCache <- which(cacheInfo$type == whichQueue &
                              cacheInfo$fcsfile == sampleFile)

      if(length(indexesInCache) == 0) {
        stop("Combination ('", whichQueue,"'", sampleFile, ") not found in cache")
      }
    } else {
      if (!("fcsfile" %in% colnames(cacheInfo))) {
        # special case appearing when no pre-processing steps have been run
        indexesInCache <- which(cacheInfo$type == whichQueue)
      } else {
        indexesInCache <- which(cacheInfo$type == whichQueue &
                                  is.na(cacheInfo$fcsfile))
      }
      
      if (length(indexesInCache) == 0) {
        stop("Combination ('", whichQueue,"', no sample file) not found in cache")
      }
    }

    objectInfos <- as.data.frame(cacheInfo[indexesInCache,
                                           c("outputObjectName",
                                             "outputClass")])
    colnames(objectInfos) <- c("ObjectName", "ObjectClass")
    return(objectInfos)
  }

##' @title Retrieves a flowFrame object from a file cache
##' @description Given a CytoPipeline object, this function retrieves
##' a specific flowCore::flowFrame object in the corresponding file cache
##' @param x a CytoPipeline object
##' @param path root path to locate the search for file caches
##' @param whichQueue which queue to look into
##' @param sampleFile which sampleFile is looked for:
##' - if whichQueue == "scale transform", the sampleFile is ignored
##' - if NULL and whichQueue == "pre-processing", the sampleFile is
##' defaulted to the first one belonging to the experiment 
##' @param objectName (character) which object name to look for
##' @return the found flowFrame (or stops with an error message if the
##' target object is not found, or if the object is no flowFrame)
##' 
##' @export
##'
getCytoPipelineFlowFrame <-
  function(x,
           path = ".",
           whichQueue = c("scale transform", "pre-processing"),
           sampleFile,
           objectName) {
    ret <-
      getCytoPipelineObjectFromCache(x = x,
                            whichQueue = whichQueue,
                            sampleFile = sampleFile,
                            objectName = objectName,
                            path = path)

    if(!inherits(ret, "flowFrame")) {
      stop("Object '", objectName,
           "' does not appear to be a flowFrame")
    }
    return(ret)
  }

##' @title Retrieves a transformList object from a file cache
##' @description Given a CytoPipeline object, this function retrieves
##' a specific flowCore::transformList object in the corresponding 
##' file cache 
##' @param x a CytoPipeline object
##' @param path root path to locate the search for file caches
##' @param whichQueue which queue to look into
##' @param sampleFile which sampleFile is looked for:
##' - if whichQueue == "scale transform", the sampleFile is ignored
##' - if NULL and whichQueue == "pre-processing", the sampleFile is
##' defaulted to the first one belonging to the experiment 
##' @param objectName (character) which object name to look for
##' @return the found flowFrame (or stops with an error message if the
##' target object is not found, or if the object is no transformList)
##' 
##' @export
##'
getCytoPipelineScaleTransform <-
  function(x,
           path = ".",
           whichQueue = c("scale transform", "pre-processing"),
           sampleFile = NULL,
           objectName) {
    ret <-
      getCytoPipelineObjectFromCache(x = x,
                            whichQueue = whichQueue,
                            sampleFile = sampleFile,
                            objectName = objectName,
                            path = path)

    if(!inherits(ret, "transformList")) {
      stop("Object '", objectName,
           "' does not appear to be a transformList")
    }

    return(ret)
  }

##' @title Plots a processing queue
##' @description This functions displays a plot of a processing
##' queue of a CytoPipeline object, using diagram::plotmat().
##' - If a step is in run state for all sample files, the 
##' corresponding box appears in green
##' - If a step is in non run state for at least one sample file, 
##' the corresponding box appears in orange
##' - If at least one step is not consistent with cache, the whole
##' set of boxes appears in red
##'
##' @param x a CytoPipeline object
##' @param path root path to locate the search for file caches
##' @param whichQueue which queue to look into
##' @param sampleFile which sampleFile is looked for:
##' - if whichQueue == "scale transform", the sampleFile is ignored
##' - if NULL and whichQueue == "pre-processing", the sampleFile is
##' defaulted to the first one belonging to the experiment 
##' @param title if TRUE, adds a title to the plot
##' @param box.type shape of label box (rect, ellipse, diamond, 
##' round, hexa, multi)
##' @param lwd default line width of arrow and box (one numeric value)
##' @param box.prop length/width ratio of label box (one numeric value)
##' @param box.cex relative size of text in boxes (one numeric value)
##' @param cex.txt relative size of arrow text (one numeric value)
##' @param box.size size of label box (one numeric value)
##' @param dtext controls the position of arrow text relative to arrowhead
##' (one numeric value)
##' @param ... other arguments passed to diagram::plotmat()
##' @return nothing
##' 
##' @export
##'
plotCytoPipelineProcessingQueue <-
  function(x,
           whichQueue = c("pre-processing", "scale transform"),
           sampleFile = NULL,
           path = ".",
           title = TRUE,
           box.type = "ellipse",
           lwd = 1,
           box.prop = 0.5,
           box.cex = 0.7,
           cex.txt = 0.7,
           box.size = 0.1,
           dtext = 0.15, 
           ...) {

    #browser()

    stopifnot(inherits(x, "CytoPipeline"))

    whichQueue <- match.arg(whichQueue)

    sampleFileIndex <- 0
    if (whichQueue == "pre-processing") {
      if (is.null(sampleFile) && length(x@sampleFiles) > 0) {
        message("no sample file passed as argument => defaulting to first ",
                "sample file")
        sampleFileIndex = 1
      } else if (is.numeric(sampleFile)) {
        if (sampleFile>0 && sampleFile <= length(x@sampleFiles)) {
          sampleFileIndex <- sampleFile
        } else {
          stop("sampleFile out of bounds")
        }
      } else {
        sampleFileIndex <- which(basename(x@sampleFiles) == basename(sampleFile))
        if (length(sampleFileIndex) == 0) {
          stop("sampleFile not found in CytoPipeline")
        } else if (length(sampleFileIndex) > 1) {
          stop("sampleFile found multiple times in CytoPipeline ",
               "(unexpected inconsistency)")
        }
      }
    } else {
      if (!is.null(sampleFile)) {
        message("sample file passed but not taken into account as the ",
                "processing queue is 'scale transform'")
      }
    }

    nSteps <- getNbProcessingSteps(x, whichQueue = whichQueue)
    steps <- getProcessingStepNames(x, whichQueue = whichQueue)
    names <- c("raw data", rep(" ", nSteps))
    box.lcol <- c("black", rep("orange", nSteps))

    res <- checkCytoPipelineConsistencyWithCache(x, path = path)
    if (!res$isConsistent) {
      warning("CytoPipeline object not consistent with cache: ",
              res$inconsistencyMsg)
      box.lcol <- c("black", rep("red", nSteps))
    } else if (sampleFileIndex == 0 && whichQueue == "pre-processing") {
      warning("No sample file in CytoPipeline object")
      box.lcol <- c("black", rep("red", nSteps))
    } else if (nSteps > 0) {
      for (i in 1:nSteps) {
        stepStatus <- ""
        if (whichQueue == "pre-processing") {
          if (unname(res$preProcessingStepStatus[i, sampleFileIndex]) == "run") {
            box.lcol[i+1] <- "green"
          }
          if (res$preProcessingStepOutputObjNames[i] != "unknown") {
            names[i+1] <- res$preProcessingStepOutputObjNames[i]
            if (res$preProcessingStepOutputClasses[i] != "unknown") {
              names[i+1] <- paste0(names[i+1], "\n(",
                                   res$preProcessingStepOutputClasses[i], ")")
            }
          }
        } else {
          if (unname(res$scaleTransformStepStatus[i]) == "run") {
            box.lcol[i+1] <- "green"
          }
          if (res$scaleTransformStepOutputObjNames[i] != "unknown") {
            names[i+1] <- res$scaleTransformStepOutputObjNames[i]
            if (res$scaleTransformStepOutputClasses[i] != "unknown") {
              names[i+1] <- paste0(names[i+1], "\n(",
                                   res$scaleTransformStepOutputClasses[i], ")")
            }
          }
        }
      }
    }

    #browser()

    M <- matrix(nrow = nSteps+1, ncol = nSteps+1, byrow = TRUE, data = 0)

    if(nSteps >=1) {
      for (j in seq_along(steps)) {
        M[j+1 ,j] <- steps[j]
      }
    }

    if (nSteps %% 2 == 0) {
      pos <- c(rep(2, nSteps/2), 1)
    } else {
      pos <- rep(2, (nSteps+1)/2)
    }

    diagram::plotmat(M, pos = pos, curve = 0, name = names,
                     box.type = box.type,
                     box.lcol = box.lcol,
                     lwd = lwd,
                     box.prop = box.prop,
                     box.cex = box.cex,
                     cex.txt = cex.txt,
                     box.size = box.size,
                     dtext = dtext,
                     ...)


    if (title) {
      theTitle <- paste0("Experiment: ", x@experimentName,
                         "\nProcessing queue: ", whichQueue)
      if (sampleFileIndex != 0) {
        theTitle <- paste0(theTitle,
                           "\nSample:",
                           basename(sampleFiles(x)[sampleFileIndex]))
      }
      graphics::title(main = theTitle)
    }

  }











