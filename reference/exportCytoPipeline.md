# exporting CytoPipeline objects

functions to export CytoPipeline objects in various formats

## Usage

``` r
export2JSONFile(x, path)
```

## Arguments

- x:

  a CytoPipeline object

- path:

  the full path to the name of the file to be created

## Value

- for `export2JSONFile`: nothing

## Functions

- `export2JSONFile()`: exports a CytoPipeline object to a JSON file
  (writing the file = side effect)

## Examples

``` r

outputDir <- base::tempdir()

rawDataDir <-
    system.file("extdata", package = "CytoPipeline")
experimentName <- "OMIP021_PeacoQC"
sampleFiles <- file.path(rawDataDir, list.files(rawDataDir,
                                             pattern = "Donor"))

# build CytoPipeline object using json input
jsonPath <- file.path(system.file("extdata", package = "CytoPipeline"), 
                      "pipelineParams.json")
  
pipL <- CytoPipeline(jsonPath,
                     experimentName = experimentName,
                     sampleFiles = sampleFiles)

# remove the last pre-processing step
nPreProcessing <- getNbProcessingSteps(pipL, whichQueue = "pre-processing")
pipL <- removeProcessingStep(pipL, whichQueue = "pre-processing", 
                                   index = nPreProcessing)

# export back to json file    
export2JSONFile(pipL, path = file.path(outputDir, "newFile.json")) 
```
