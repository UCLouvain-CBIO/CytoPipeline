# compensation of fcs file(s) from matrix

executes the classical compensation function on a flowSet or flowFrame,
given a compensation matrix. The matrix can be either retrieved in the
fcs files themselves or provided as a csv file.

## Usage

``` r
compensateFromMatrix(
  x,
  matrixSource = c("fcs", "import"),
  matrixPath = NULL,
  updateChannelNames = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  a
  [`flowCore::flowFrame`](https://rdrr.io/pkg/flowCore/man/flowFrame-class.html)
  or
  [`flowCore::flowSet`](https://rdrr.io/pkg/flowCore/man/flowSet-class.html)

- matrixSource:

  if "fcs", the compensation matrix will be fetched from the fcs files
  (different compensation matrices can then be applied by fcs file) if
  "import", uses `matrixPath` to read the matrix (should be a csv file)

- matrixPath:

  if matrixSource == "import", will be used as the input csv file path

- updateChannelNames:

  if TRUE, updates the fluo channel names by prefixing them with "comp-"

- verbose:

  if TRUE, displays information messages

- ...:

  additional arguments (not used)

## Value

the compensated flowSet or flowFrame

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
```
