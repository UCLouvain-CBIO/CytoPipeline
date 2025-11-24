# remove dead cells from a flowFrame using manual gating

remove dead cells from a flowFrame, using manual gating in the FSC-A,
'(a)Live/Dead' 2D representation. The function uses
flowCore::polygonGate()

## Usage

``` r
removeDeadCellsManualGate(
  ff,
  preTransform = FALSE,
  transList = NULL,
  FSCChannel,
  LDMarker,
  gateData,
  ...
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- preTransform:

  boolean, if TRUE: the transList list of scale transforms will be
  applied first on the LD channel.

- transList:

  applied in conjunction with preTransform == TRUE

- FSCChannel:

  a character containing the exact name of the forward scatter channel

- LDMarker:

  a character containing the exact name of the marker corresponding to
  (a)Live/Dead channel, or the Live/Dead channel name itself

- gateData:

  a numerical vector containing the polygon gate coordinates first the
  `FSCChannel` channel coordinates of each points of the polygon gate,
  then the LD channel coordinates of each points (prior to scale
  transform)

- ...:

  additional parameters passed to flowCore::polygonGate()

## Value

a flowCore::flowFrame with removed dead cells from the input

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
                         
remDeadCellsGateData <- c(0, 0, 250000, 250000,
                          0, 650, 650, 0)  

ff_lcells <-
    removeDeadCellsManualGate(ff_c,
                              FSCChannel = "FSC-A",
                              LDMarker = "L/D Aqua - Viability",
                              gateData = remDeadCellsGateData)
   
```
