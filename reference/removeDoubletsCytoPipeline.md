# remove doublets from a flowFrame, using CytoPipeline custom algorithm

Wrapper around CytoPipeline::singletGate(). Can apply the flowStats
function subsequently on several channel pairs, e.g. (FSC-A, FSC-H) and
(SSC-A, SSC-H)

## Usage

``` r
removeDoubletsCytoPipeline(ff, areaChannels, heightChannels, nmads, ...)
```

## Arguments

- ff:

  a flowCore::flowFrame

- areaChannels:

  a character vector containing the name of the 'area type' channels one
  wants to use

- heightChannels:

  a character vector containing the name of the 'height type' channels
  one wants to use

- nmads:

  a numeric vector with the bandwidth above the ratio allowed, per
  channels pair (cells are kept if the ratio between -A channel\[i\] and
  -H channel\[i\] is smaller than the median ratio + nmad\[i\] times the
  median absolute deviation of the ratios). Default is 4, for all
  channel pairs.

- ...:

  additional parameters passed to CytoPipeline::singletGate()

## Value

a flowCore::flowFrame with removed doublets events from the input

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

ff_s <-
    removeDoubletsCytoPipeline(ff_c,
                               areaChannels = c("FSC-A", "SSC-A"),
                               heightChannels = c("FSC-H", "SSC-H"),
                               nmads = c(3, 5))
                            
```
