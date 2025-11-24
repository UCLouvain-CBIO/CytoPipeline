# Read fcs sample files

Wrapper around flowCore::read.fcs() or flowCore::read.flowSet(). Also
adds a "Cell_ID" additional column, used in flowFrames comparison

## Usage

``` r
readSampleFiles(
  sampleFiles,
  whichSamples = "all",
  nSamples = NULL,
  seed = NULL,
  channelMarkerFile = NULL,
  ...
)
```

## Arguments

- sampleFiles:

  a vector of character path to sample files

- whichSamples:

  one of:

  - 'all' if all sample files need to be read

  - 'random' if some samples need to be chosen randomly (in that case,
    using `nSamples` and `seed`)

  - a vector of indexes pointing to the sampleFiles vector

- nSamples:

  number of samples to randomly select (if `whichSamples == "random"`).
  If `nSamples` is higher than nb of available samples, the output will
  be all samples

- seed:

  an optional seed parameters (provided to ease reproducibility).

- channelMarkerFile:

  an optional path to a csv file which provides the mapping between
  channels and markers. If provided, this csv file should contain a
  `Channel` column, and a `Marker` column. Optionally a 'Used' column
  can be provided as well (TRUE/FALSE). Channels for which the 'Used'
  column is set to FALSE will not be incorporated in the created
  flowFrame.

- ...:

  additional parameters passed to flowCore file reading functions.

## Value

either a flowCore::flowSet or a flowCore::flowFrame if
length(sampleFiles) == 1

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

#res

# create a flowCore::flowFrame with one single sample
res2 <- readSampleFiles(
    sampleFiles = sampleFiles,
    whichSamples = 2,
    truncate_max_range = truncateMaxRange,
    min.limit = minLimit)

#res2
```
