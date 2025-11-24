# extract compensation matrix from a flowCore::flowFrame

helper function retrieving the compensation matrix stored in fcs file
(if any). It scans the following keywords: \$SPILL, \$spillover and
\$SPILLOVER

## Usage

``` r
getAcquiredCompensationMatrix(ff)
```

## Arguments

- ff:

  a flowCore::flowFrame

## Value

the found compensation matrix

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
compensationMatrix <- getAcquiredCompensationMatrix(fsRaw[[2]])
```
