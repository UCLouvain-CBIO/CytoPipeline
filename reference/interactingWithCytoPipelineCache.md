# interaction between CytoPipeline object and disk cache

functions supporting the interaction between a CytoPipeline object and
the file cache on disk

## Usage

``` r
deleteCytoPipelineCache(x, path = ".")

buildCytoPipelineFromCache(experimentName, path = ".")

checkCytoPipelineConsistencyWithCache(
  x,
  path = ".",
  whichQueue = c("both", "scale transform", "pre-processing"),
  sampleFile = NULL
)
```

## Arguments

- x:

  a CytoPipeline object

- path:

  the full path to the experiment storage on disk (without the /.cache)

- experimentName:

  the experimentName used to select the file cache on disk

- whichQueue:

  which processing queue to check the consistency of

- sampleFile:

  if whichQueue == "pre-processing" or "both": which sample file(s) to
  check on the disk cache

## Value

for `deleteCytoPipelineCache`: TRUE if successfully removed  
for `buildCytoPipelineFromCache`: the built CytoPipeline object  
for `checkCytoPipelineConsistencyWithCache`: a list with the following
values:

- `isConsistent` (TRUE/FALSE)

- `inconsistencyMsg`: character filled in by an inconsistency message in
  case the cache and CytoPipeline object are not consistent with each
  other

- `scaleTransformStepStatus`: a character vector, containing, for each
  scale transform step, a status from c("run", "not run",
  "inconsistent")

- `preProcessingStepStatus`: a character matrix, containing, for each
  pre-processing step (rows), for each sample file (columns), a status
  from c("run", "not run", "inconsistent")

## Functions

- `deleteCytoPipelineCache()`: delete the whole disk cache corresponding
  to the experiment of a CytoPipeline object

- `buildCytoPipelineFromCache()`: builds a new CytoPipeline object,
  based on the information stored in the file cache

- `checkCytoPipelineConsistencyWithCache()`: check the consistency
  between the processing steps described in a CytoPipeline object, and
  what is stored in the file cache

## Examples

``` r

# preliminary run:
# build CytoPipeline object using json input, run and store results in cache
rawDataDir <-
    system.file("extdata", package = "CytoPipeline")
experimentName <- "OMIP021_PeacoQC"
sampleFiles <- file.path(rawDataDir, list.files(rawDataDir,
                                             pattern = "Donor"))
                                             
jsonDir <- system.file("extdata", package = "CytoPipeline")
jsonPath <- file.path(jsonDir, "pipelineParams.json")
outputDir <- base::tempdir()
pipL <- CytoPipeline(jsonPath,
                     experimentName = experimentName,
                     sampleFiles = sampleFiles)

# note we temporarily set working directory into package root directory
# needed as json path mentions "./" path for sample files
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
     

# rebuild CytoPipeline from stored results in cache, for a specific 
# experiment

experimentName <- "OMIP021_PeacoQC"
pipL2 <- buildCytoPipelineFromCache(
    experimentName = experimentName,
    path = outputDir)


# checking consistency between CytoPipeline object and cache
res <- checkCytoPipelineConsistencyWithCache(pipL2)
#res

suppressWarnings(execute(pipL2, rmCache = FALSE, path = outputDir))
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
# (everything is already stored in cache)

# deleting cache related to a specific experiment
pipL3 <- CytoPipeline(experimentName = experimentName)
deleteCytoPipelineCache(pipL3, path = outputDir)
#> [1] TRUE
```
