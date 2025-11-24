# remove debris from a flowFrame using manual gating

remove debris from a flowFrame, using manual gating in the FSC-A, SSC-A
2D representation. The function internally uses flowCore::polygonGate()

## Usage

``` r
removeDebrisManualGate(ff, FSCChannel, SSCChannel, gateData, ...)
```

## Arguments

- ff:

  a flowCore::flowFrame

- FSCChannel:

  a character containing the exact name of the forward scatter channel

- SSCChannel:

  a character containing the exact name of the side scatter channel

- gateData:

  a numerical vector containing the polygon gate coordinates first the
  `FSCChannel` channel coordinates of each points of the polygon gate,
  then the `SSCChannel` channel coordinates of each points.

- ...:

  additional parameters passed to flowCore::polygonGate()

## Value

a flowCore::flowFrame with removed debris events from the input

## Examples

``` r

rawDataDir <-
    system.file("extdata", package = "CytoPipeline")
sampleFiles <-
    file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))

truncateMaxRange <- FALSE
minLimit <- NULL

# create flowCore::flowSet with all samples of a dataset
fsRaw <- readSampleFiles(
    sampleFiles = sampleFiles,
    whichSamples = "all",
    truncate_max_range = truncateMaxRange,
    min.limit = minLimit)

suppressWarnings(ff_m <- removeMarginsPeacoQC(x = fsRaw[[2]]))
#> Removing margins from file : Donor2.fcs
    
ff_c <-
    compensateFromMatrix(ff_m,
                         matrixSource = "fcs")        


remDebrisGateData <- c(73615, 110174, 213000, 201000, 126000,
                       47679, 260500, 260500, 113000, 35000)

ff_cells <-
    removeDebrisManualGate(ff_c,
                           FSCChannel = "FSC-A",
                           SSCChannel = "SSC-A",
                           gateData = remDebrisGateData)

```
