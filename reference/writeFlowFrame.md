# write flowFrame to disk

wrapper around flowCore::write.FCS() or utils::write.csv that discards
any additional parameter passed in (...)

## Usage

``` r
writeFlowFrame(
  ff,
  dir = ".",
  useFCSFileName = TRUE,
  prefix = "",
  suffix = "",
  format = c("fcs", "csv"),
  csvUseChannelMarker = TRUE,
  ...
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- dir:

  an existing directory to store the flowFrame,

- useFCSFileName:

  if TRUE filename used will be based on original fcs filename

- prefix:

  file name prefix

- suffix:

  file name suffix

- format:

  either fcs or csv

- csvUseChannelMarker:

  if TRUE (default), converts the channels to the corresponding marker
  names (where the Marker is not NA). This setting is only applicable to
  export in csv format.

- ...:

  other arguments (not used)

## Value

nothing

## Examples

``` r
rawDataDir <-
    system.file("extdata", package = "CytoPipeline")
sampleFiles <-
    file.path(rawDataDir, list.files(rawDataDir, pattern = "Donor"))

truncateMaxRange <- FALSE
minLimit <- NULL

# create flowCore::flowSet with all samples of a dataset
res <- readSampleFiles(
    sampleFiles = sampleFiles,
    whichSamples = "all",
    truncate_max_range = truncateMaxRange,
    min.limit = minLimit)
    
ff_c <- compensateFromMatrix(res[[2]], matrixSource = "fcs") 
outputDir <- base::tempdir()
writeFlowFrame(ff_c, 
               dir = outputDir,
               suffix = "_fcs_export",
               format = "csv")
```
