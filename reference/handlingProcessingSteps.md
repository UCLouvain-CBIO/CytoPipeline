# handling processing steps in CytoPipeline objects

functions to manipulate processing steps in processing queues of
CytoPipeline objects

## Usage

``` r
addProcessingStep(
  x,
  whichQueue = c("scale transform", "pre-processing"),
  newPS
)

removeProcessingStep(
  x,
  whichQueue = c("scale transform", "pre-processing"),
  index
)

getNbProcessingSteps(x, whichQueue = c("scale transform", "pre-processing"))

getProcessingStep(
  x,
  whichQueue = c("scale transform", "pre-processing"),
  index
)

getProcessingStepNames(x, whichQueue = c("scale transform", "pre-processing"))

cleanProcessingSteps(
  x,
  whichQueue = c("both", "scale transform", "pre-processing")
)

showProcessingSteps(x, whichQueue = c("scale transform", "pre-processing"))
```

## Arguments

- x:

  a CytoPipeline object

- whichQueue:

  selects the processing queue for which we manage the processing steps

- newPS:

  the new processing step to be added (CytoProcessingStep object)

- index:

  index of the processing step to remove

## Value

- for `addProcessingStep`: the updated CytoPipeline object

&nbsp;

- for `removeProcessingStep`: the updated CytoPipeline object

&nbsp;

- for `getNbProcessingSteps`: the number of processing steps present in
  the target queue

&nbsp;

- for `getProcessingStep`: the obtained CytoProcessingStep object

&nbsp;

- for `getProcessingStepNames`: the vector of step names

&nbsp;

- for `cleanProcessingSteps`: the updated CytoPipeline object

&nbsp;

- for `showProcessingSteps`: nothing (only console display side effect
  is required)

## Functions

- `addProcessingStep()`: adds a processing step in one of the processing
  queues (at the end), returns the modified CytoPipeline object

- `removeProcessingStep()`: removes a processing step from one of the
  processing queues, returns the modified CytoPipeline object

- `getNbProcessingSteps()`: gets the number of processing steps in a
  processing queue

- `getProcessingStep()`: gets a processing step at a specific index of a
  processing queue

- `getProcessingStepNames()`: gets a character vector of all processing
  step names of a specific processing queue

- `cleanProcessingSteps()`: deletes all processing steps in one or both
  processing queues, returns the modified CytoPipeline object

- `showProcessingSteps()`: shows all processing steps in a processing
  queue

## Examples

``` r

rawDataDir <-
    system.file("extdata", package = "CytoPipeline")
experimentName <- "OMIP021_PeacoQC"
sampleFiles <- file.path(rawDataDir, list.files(rawDataDir,
                                             pattern = "Donor"))
transListPath <- 
    file.path(system.file("extdata", package = "CytoPipeline"), 
              "OMIP021_TransList.rds")

# main parameters : sample files and experiment name
pipelineParams <- list()
pipelineParams$experimentName <- experimentName
pipelineParams$sampleFiles <- sampleFiles

# create CytoPipeline object (no step defined yet)
pipL <- CytoPipeline(pipelineParams)

# add a processing step in scale tranformation queue
pipL <- addProcessingStep(pipL,
                          whichQueue = "scale transform",
                          CytoProcessingStep(
                              name = "scale_transform_read",
                              FUN = "readRDS",
                              ARGS = list(file = transListPath)
                          ))

getNbProcessingSteps(pipL, "scale transform") # returns 1
#> [1] 1

# add another processing step in scale transformation queue
pipL <- addProcessingStep(pipL,
                          whichQueue = "scale transform",
                          CytoProcessingStep(
                              name = "scale_transform_sum",
                              FUN = "sum",
                              ARGS = list()
                          )
)

getNbProcessingSteps(pipL, "scale transform") # returns 2
#> [1] 2

getProcessingStepNames(pipL, whichQueue = "scale transform")
#> [1] "scale_transform_read" "scale_transform_sum" 

# removes second processing step in scale transformation queue
pipL <- removeProcessingStep(pipL,
                             whichQueue = "scale transform",
                             index = 2)

# get processing step object
pS <- getProcessingStep(pipL, whichQueue = "scale transform", index = 1)
getCPSName(pS) #"scale_transform_read"
#> [1] "scale_transform_read"

# add a processing step in pre-processing queue
pipL <- addProcessingStep(pipL,
                          whichQueue = "pre-processing",
                          CytoProcessingStep(
                              name = "pre-processing_sum",
                              FUN = "sum",
                              ARGS = list()
                          ))
getNbProcessingSteps(pipL, "scale transform") # returns 1
#> [1] 1
getNbProcessingSteps(pipL, "pre-processing") # returns also 1
#> [1] 1

showProcessingSteps(pipL, whichQueue = "scale transform")
#> Scale transformations evaluation queue : 1 processing step(s)
#> 1 :Object of class "CytoProcessingStep"
#>  Name: scale_transform_read
#>  Function: readRDS
#>  Arguments:
#>   o file = /__w/_temp/Library/CytoPipeline/extdata/OMIP021_TransList.rds
showProcessingSteps(pipL, whichQueue = "pre-processing")
#> Flow frames pre-processing evaluation queue : 1 processing step(s)
#> 1 :Object of class "CytoProcessingStep"
#>  Name: pre-processing_sum
#>  Function: sum

# cleans both processing queues
pipL <- cleanProcessingSteps(pipL)
pipL
#> Pipeline object for flow cytometry experiment: OMIP021_PeacoQC 
#> Sample files: 2 sample file(s)
#> head(samples):
#>   displayName                                         sampleFile
#> 1  Donor1.fcs /__w/_temp/Library/CytoPipeline/extdata/Donor1.fcs
#> 2  Donor2.fcs /__w/_temp/Library/CytoPipeline/extdata/Donor2.fcs
#> No pheno data
#> Scale transformations evaluation queue has no processing step
#> Flow frames pre-processing evaluation queue has no processing step
```
