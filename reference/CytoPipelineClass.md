# CytoPipeline class

Class representing a flow cytometry pipeline, and composed of two
processing queues, i.e. lists of CytoProcessingStep objects :

- a list of CytoProcessingStep(s) for pre-calculation of scale
  transformations per channel

- a list of CytoProcessingStep(s) for the pre-processing of flow frames

## Usage

``` r
# S4 method for class 'CytoPipeline'
show(object)

# S4 method for class 'missing'
CytoPipeline(
  object,
  experimentName = "default_experiment",
  sampleFiles = character(),
  pData = NULL
)

# S4 method for class 'list'
CytoPipeline(
  object,
  experimentName = "default_experiment",
  sampleFiles = character(),
  pData = NULL
)

# S4 method for class 'character'
CytoPipeline(
  object,
  experimentName = "default_experiment",
  sampleFiles = character(),
  pData = NULL
)

# S3 method for class 'CytoPipeline'
as.list(x, ...)

experimentName(x)

experimentName(x) <- value

sampleFiles(x)

sampleFiles(x) <- value

pData(x)

pData(x) <- value

sampleDisplayNames(x, sampleFiles = NULL)

sampleNameFromDisplayName(x, displayName)
```

## Arguments

- object:

  a [`character()`](https://rdrr.io/r/base/character.html) containing a
  JSON input

- experimentName:

  the experiment name

- sampleFiles:

  a character (e.g. sampleFileNames) or a numeric vector (e.g. indices
  of sample files). If NULL, all samples will be displayed.

- pData:

  the pheno Data (data.frame or NULL)

- x:

  a `CytoPipeline` object

- ...:

  additional arguments (not used here)

- value:

  the new value to be assigned. the `pData<-` setter is a bit more
  liberal than it used to be:

  1.  It can accept new pData containing more rows than existing sample
      names (the corresponding subset of pData is taken).

  2.  It can accept pData with row names pointing to either sample file
      full paths or base file names

  3.  It can accept pData with no row names provided the number of rows
      correspond to the number of sample files. Row names are then set
      by default to sample file base names (if unique), or sample file
      full paths.

- displayName:

  a character

## Value

nothing

- for `as.list.CytoPipeline`: the obtained list

&nbsp;

- for `sampleDisplayNames`: a character vector of sample display names

&nbsp;

- for `sampleNameFromDisplayName`: the sample name corresponding to the
  specified display name. of sample display names

## Slots

- `scaleTransformProcessingQueue`:

  A `list` of CytoProcessingStep objects containing the steps for
  obtaining the scale transformations per channel

- `flowFramesPreProcessingQueue`:

  A `list` of CytoProcessingStep objects containing the steps for
  pre-processing of the samples flow frames

- `experimentName`:

  A `character` containing the experiment (run) name

- `sampleFiles`:

  A `character` vector storing all fcs files to be run into the pipeline

- `pData`:

  An optional `data.frame` containing additional information for each
  sample file. The `pData` raw names should correspond to the sample
  files (using full paths or base paths). If the `pData` contains a
  columns with name 'displayName', this will have an impact in the
  `sampleDisplayNames()` function, i.e. sample display names will be the
  one mentioned in `pData`, instead of typically base file names (or
  larger paths if base file names are not unique)

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
pipL <- CytoPipeline(experimentName = experimentName,
                     sampleFiles = sampleFiles)

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
                              47679, 260500, 260500, 113000, 35000)))
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

### *** EXAMPLE 2: building CytoPipeline from JSON file *** ###

jsonDir <- system.file("extdata", package = "CytoPipeline")
jsonPath <- file.path(jsonDir, "pipelineParams.json")

pipL2 <- CytoPipeline(jsonPath,
                      experimentName = experimentName,
                      sampleFiles = sampleFiles)
```
