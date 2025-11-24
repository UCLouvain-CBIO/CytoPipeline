# executing CytoPipeline object

this function triggers the execution of the processing queues of a
CytoPipeline object. First, the scale tranform processing queue is run,
taking the set of sample names as an implicit first input. At the end of
the queue, a scale transform List is assumed to be created. Second, the
flowFrame pre-processing queue, reapeatedly for each sample file. The
scale transform list generated in the previous step is taken as implicit
input, together with the initial sample file. At the end of the queue
run, a pre-processed flowFrame is assumed to be generated. No change is
made on the input CytoPipeline object, all results are stored in the
cache.

## Usage

``` r
execute(
  x,
  path = ".",
  rmCache = FALSE,
  useBiocParallel = FALSE,
  BPPARAM = BiocParallel::bpparam(),
  BPOPTIONS = BiocParallel::bpoptions(packages = c("flowCore")),
  saveLastStepFF = TRUE,
  saveFFSuffix = "_preprocessed",
  saveFFFormat = c("fcs", "csv"),
  saveFFCsvUseChannelMarker = TRUE,
  saveScaleTransforms = FALSE
)
```

## Arguments

- x:

  CytoPipeline object

- path:

  base path, a subdirectory with name equal to the experiment will be
  created to store the output data, in particular the experiment cache

- rmCache:

  if TRUE, starts by removing the already existing cache directory
  corresponding to the experiment

- useBiocParallel:

  if TRUE, use BiocParallel for computation of the sample file
  pre-processing in parallel (one file per worker at a time). Note the
  BiocParallel function used is `bplapply()`

- BPPARAM:

  if `useBiocParallel` is TRUE, sets the BPPARAM back-end to be used for
  the computation. If not provided, will use the top back-end on the
  [`BiocParallel::registered()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  stack.

- BPOPTIONS:

  if `useBiocParallel` is TRUE, sets the BPOPTIONS to be passed to
  `bplapply()` function. Note that if you use a `SnowParams` back-end,
  you need to specify all the packages that need to be loaded for the
  different CytoProcessingStep to work properly (visibility of
  functions). As a minimum, the `flowCore` package needs to be loaded.
  (hence the default `BPOPTIONS = bpoptions(packages = c("flowCore"))` )

- saveLastStepFF:

  if TRUE, save the final result of the pre-processing, for each file.
  By convention, these output files are stored in
  `path`/`x@experimentName`/output/, the file names used are the same as
  the initial fcs file basenames, concatenated with `saveFFSuffix`, and
  with file extension corresponding to `saveFFFormat`.

- saveFFSuffix:

  FF file name suffix

- saveFFFormat:

  either `fcs` or `csv`

- saveFFCsvUseChannelMarker:

  if TRUE (default), converts the channels to the corresponding marker
  names (where the Marker is not NA). This setting is only applicable to
  export in csv format.

- saveScaleTransforms:

  if TRUE (default FALSE), save on disk (in RDS format) the
  [`flowCore::transformList`](https://rdrr.io/pkg/flowCore/man/transformList-class.html)
  object obtained after running the scaleTransform processing queue. The
  file name is hardcoded to
  `path`/`experimentName`/`RDS`/`scaleTransformList.rds`

## Value

nothing

## Examples

``` r

### *** EXAMPLE 1: building CytoPipeline step by step *** ###

rawDataDir <-
    system.file("extdata", package = "CytoPipeline")
experimentName <- "OMIP021_PeacoQC"
sampleFiles <- file.path(rawDataDir, list.files(rawDataDir,
                                             pattern = "Donor"))
                                             
outputDir <- base::tempdir()

# main parameters : sample files and output files
pipelineParams <- list()
pipelineParams$experimentName <- experimentName
pipelineParams$sampleFiles <- sampleFiles
pipL <- CytoPipeline(pipelineParams)

### SCALE TRANSFORMATION STEPS ###

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "scale transform",
                      CytoProcessingStep(
                          name = "flowframe_read",
                          FUN = "readSampleFiles",
                          ARGS = list(
                              whichSamples = "all",
                              truncate_max_range = FALSE,
                              min.limit = NULL
                          )
                      )
    )

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "scale transform",
                      CytoProcessingStep(
                          name = "remove_margins",
                          FUN = "removeMarginsPeacoQC",
                          ARGS = list()
                     )
    )

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "scale transform",
                      CytoProcessingStep(
                          name = "compensate",
                          FUN = "compensateFromMatrix",
                          ARGS = list(matrixSource = "fcs")
                      )
    )

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "scale transform",
                      CytoProcessingStep(
                          name = "flowframe_aggregate",
                          FUN = "aggregateAndSample",
                          ARGS = list(
                              nTotalEvents = 10000,
                              seed = 0
                          )
                      )
    )

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "scale transform",
                      CytoProcessingStep(
                          name = "scale_transform_estimate",
                          FUN = "estimateScaleTransforms",
                          ARGS = list(
                              fluoMethod = "estimateLogicle",
                              scatterMethod = "linear",
                              scatterRefMarker = "BV785 - CD3"
                          )
                      )
    )

### PRE-PROCESSING STEPS ###

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "pre-processing",
                      CytoProcessingStep(
                          name = "flowframe_read",
                          FUN = "readSampleFiles",
                          ARGS = list(
                              truncate_max_range = FALSE,
                              min.limit = NULL
                          )
                      )
    )

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "pre-processing",
                      CytoProcessingStep(
                          name = "remove_margins",
                          FUN = "removeMarginsPeacoQC",
                          ARGS = list()
                      )
    )

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "pre-processing",
                      CytoProcessingStep(
                          name = "compensate",
                          FUN = "compensateFromMatrix",
                          ARGS = list(matrixSource = "fcs")
                      )
    )

pipL <-
addProcessingStep(
    pipL,
    whichQueue = "pre-processing",
    CytoProcessingStep(
        name = "remove_debris",
        FUN = "removeDebrisManualGate",
        ARGS = list(
            FSCChannel = "FSC-A",
            SSCChannel = "SSC-A",
            gateData =  c(73615, 110174, 213000, 201000, 126000,
                          47679, 260500, 260500, 113000, 35000)
                   )
   )
)

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "pre-processing",
                      CytoProcessingStep(
                          name = "remove_dead_cells",
                          FUN = "removeDeadCellsManualGate",
                          ARGS = list(
                              FSCChannel = "FSC-A",
                              LDMarker = "L/D Aqua - Viability",
                              gateData = c(0, 0, 250000, 250000,
                                           0, 650, 650, 0)
                          )
                      )
    )

pipL <-
    addProcessingStep(
        pipL,
        whichQueue = "pre-processing",
        CytoProcessingStep(
            name = "perform_QC",
            FUN = "qualityControlPeacoQC",
            ARGS = list(
                preTransform = TRUE,
                min_cells = 150, # default
                max_bins = 500, # default
                step = 500, # default,
                MAD = 6, # default
                IT_limit = 0.55, # default
                force_IT = 150, # default
                peak_removal = 0.3333, # default
                min_nr_bins_peakdetection = 10 # default
            )
        )
    )

pipL <-
    addProcessingStep(pipL,
                      whichQueue = "pre-processing",
                      CytoProcessingStep(
                          name = "transform",
                          FUN = "applyScaleTransforms",
                          ARGS = list()
                      )
    )

# execute pipeline, remove cache if existing with the same experiment name
suppressWarnings(execute(pipL, rmCache = TRUE, path = outputDir))
#> #####################################################
#> ### running SCALE TRANSFORMATION processing steps ###
#> #####################################################
#> Proceeding with step 1 [flowframe_read] ...
#> Proceeding with step 2 [remove_margins] ...
#> Removing margins from file : Donor1.fcs
#> Removing margins from file : Donor2.fcs
#> Proceeding with step 3 [compensate] ...
#> Proceeding with step 4 [flowframe_aggregate] ...
#> Proceeding with step 5 [scale_transform_estimate] ...
#> #####################################################
#> ### NOW PRE-PROCESSING FILE /__w/_temp/Library/CytoPipeline/extdata/Donor1.fcs...
#> #####################################################
#> Proceeding with step 1 [flowframe_read] ...
#> Proceeding with step 2 [remove_margins] ...
#> Removing margins from file : Donor1.fcs
#> Proceeding with step 3 [compensate] ...
#> Proceeding with step 4 [remove_debris] ...
#> Proceeding with step 5 [remove_dead_cells] ...
#> Proceeding with step 6 [perform_QC] ...
#> Applying PeacoQC method...
#> Starting quality control analysis for Donor1.fcs
#> Calculating peaks
#> MAD analysis removed 16.54% of the measurements
#> The algorithm removed 16.54% of the measurements
#> Proceeding with step 7 [transform] ...
#> #####################################################
#> ### NOW PRE-PROCESSING FILE /__w/_temp/Library/CytoPipeline/extdata/Donor2.fcs...
#> #####################################################
#> Proceeding with step 1 [flowframe_read] ...
#> Proceeding with step 2 [remove_margins] ...
#> Removing margins from file : Donor2.fcs
#> Proceeding with step 3 [compensate] ...
#> Proceeding with step 4 [remove_debris] ...
#> Proceeding with step 5 [remove_dead_cells] ...
#> Proceeding with step 6 [perform_QC] ...
#> Applying PeacoQC method...
#> Starting quality control analysis for Donor2.fcs
#> Calculating peaks
#> MAD analysis removed 5.4% of the measurements
#> The algorithm removed 5.4% of the measurements
#> Proceeding with step 7 [transform] ...

# re-execute as is without removing cache => all results found in cache!
suppressWarnings(execute(pipL, rmCache = FALSE, path = outputDir))
#> #####################################################
#> ### running SCALE TRANSFORMATION processing steps ###
#> #####################################################
#> Proceeding with step 1 [flowframe_read]: found in cache!
#> Proceeding with step 2 [remove_margins]: found in cache!
#> Proceeding with step 3 [compensate]: found in cache!
#> Proceeding with step 4 [flowframe_aggregate]: found in cache!
#> Proceeding with step 5 [scale_transform_estimate]: found in cache!
#> #####################################################
#> ### NOW PRE-PROCESSING FILE /__w/_temp/Library/CytoPipeline/extdata/Donor1.fcs...
#> #####################################################
#> Proceeding with step 1 [flowframe_read]: found in cache!
#> Proceeding with step 2 [remove_margins]: found in cache!
#> Proceeding with step 3 [compensate]: found in cache!
#> Proceeding with step 4 [remove_debris]: found in cache!
#> Proceeding with step 5 [remove_dead_cells]: found in cache!
#> Proceeding with step 6 [perform_QC]: found in cache!
#> Proceeding with step 7 [transform]: found in cache!
#> #####################################################
#> ### NOW PRE-PROCESSING FILE /__w/_temp/Library/CytoPipeline/extdata/Donor2.fcs...
#> #####################################################
#> Proceeding with step 1 [flowframe_read]: found in cache!
#> Proceeding with step 2 [remove_margins]: found in cache!
#> Proceeding with step 3 [compensate]: found in cache!
#> Proceeding with step 4 [remove_debris]: found in cache!
#> Proceeding with step 5 [remove_dead_cells]: found in cache!
#> Proceeding with step 6 [perform_QC]: found in cache!
#> Proceeding with step 7 [transform]: found in cache!

### *** EXAMPLE 2: building CytoPipeline from JSON file *** ###

jsonDir <- system.file("extdata", package = "CytoPipeline")
jsonPath <- file.path(jsonDir, "pipelineParams.json")

pipL2 <- CytoPipeline(jsonPath, 
                      experimentName = experimentName,
                      sampleFiles = sampleFiles)

# note we temporarily set working directory into package root directory
# needed as json path mentions "./" path for sample files
suppressWarnings(execute(pipL2, rmCache = TRUE, path = outputDir))
#> #####################################################
#> ### running SCALE TRANSFORMATION processing steps ###
#> #####################################################
#> Proceeding with step 1 [flowframe_read] ...
#> Proceeding with step 2 [remove_margins] ...
#> Removing margins from file : Donor1.fcs
#> Removing margins from file : Donor2.fcs
#> Proceeding with step 3 [compensate] ...
#> Proceeding with step 4 [flowframe_aggregate] ...
#> Proceeding with step 5 [scale_transform_estimate] ...
#> #####################################################
#> ### NOW PRE-PROCESSING FILE /__w/_temp/Library/CytoPipeline/extdata/Donor1.fcs...
#> #####################################################
#> Proceeding with step 1 [flowframe_read] ...
#> Proceeding with step 2 [remove_margins] ...
#> Removing margins from file : Donor1.fcs
#> Proceeding with step 3 [compensate] ...
#> Proceeding with step 4 [remove_doublets] ...
#> Proceeding with step 5 [remove_debris] ...
#> Proceeding with step 6 [remove_dead_cells] ...
#> Proceeding with step 7 [perform_QC] ...
#> Applying PeacoQC method...
#> Starting quality control analysis for Donor1.fcs
#> Calculating peaks
#> MAD analysis removed 30.75% of the measurements
#> The algorithm removed 30.75% of the measurements
#> Proceeding with step 8 [transform] ...
#> #####################################################
#> ### NOW PRE-PROCESSING FILE /__w/_temp/Library/CytoPipeline/extdata/Donor2.fcs...
#> #####################################################
#> Proceeding with step 1 [flowframe_read] ...
#> Proceeding with step 2 [remove_margins] ...
#> Removing margins from file : Donor2.fcs
#> Proceeding with step 3 [compensate] ...
#> Proceeding with step 4 [remove_doublets] ...
#> Proceeding with step 5 [remove_debris] ...
#> Proceeding with step 6 [remove_dead_cells] ...
#> Proceeding with step 7 [perform_QC] ...
#> Applying PeacoQC method...
#> Starting quality control analysis for Donor2.fcs
#> Calculating peaks
#> MAD analysis removed 24.38% of the measurements
#> The algorithm removed 24.38% of the measurements
#> Proceeding with step 8 [transform] ...

### *** EXAMPLE 3: building CytoPipeline from cache (previously run) *** ###

experimentName <- "OMIP021_PeacoQC"
pipL3 <- buildCytoPipelineFromCache(
    experimentName = experimentName,
    path = outputDir)

suppressWarnings(execute(pipL3,
        rmCache = FALSE,
        path = outputDir))
#> #####################################################
#> ### running SCALE TRANSFORMATION processing steps ###
#> #####################################################
#> Proceeding with step 1 [flowframe_read]: found in cache!
#> Proceeding with step 2 [remove_margins]: found in cache!
#> Proceeding with step 3 [compensate]: found in cache!
#> Proceeding with step 4 [flowframe_aggregate]: found in cache!
#> Proceeding with step 5 [scale_transform_estimate]: found in cache!
#> #####################################################
#> ### NOW PRE-PROCESSING FILE Donor1.fcs...
#> #####################################################
#> Proceeding with step 1 [flowframe_read]: found in cache!
#> Proceeding with step 2 [remove_margins]: found in cache!
#> Proceeding with step 3 [compensate]: found in cache!
#> Proceeding with step 4 [remove_doublets]: found in cache!
#> Proceeding with step 5 [remove_debris]: found in cache!
#> Proceeding with step 6 [remove_dead_cells]: found in cache!
#> Proceeding with step 7 [perform_QC]: found in cache!
#> Proceeding with step 8 [transform]: found in cache!
#> #####################################################
#> ### NOW PRE-PROCESSING FILE Donor2.fcs...
#> #####################################################
#> Proceeding with step 1 [flowframe_read]: found in cache!
#> Proceeding with step 2 [remove_margins]: found in cache!
#> Proceeding with step 3 [compensate]: found in cache!
#> Proceeding with step 4 [remove_doublets]: found in cache!
#> Proceeding with step 5 [remove_debris]: found in cache!
#> Proceeding with step 6 [remove_dead_cells]: found in cache!
#> Proceeding with step 7 [perform_QC]: found in cache!
#> Proceeding with step 8 [transform]: found in cache!
```
