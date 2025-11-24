# compensate with additional options

: this is a simple wrapper around the flowCore::compensate() utility,
allowing to trigger an update of the fluo channel names with a prefix
'comp-' (as in FlowJo)

## Usage

``` r
runCompensation(obj, spillover, updateChannelNames = TRUE)
```

## Arguments

- obj:

  a flowCore::flowFrame or flowCore::flowSet

- spillover:

  compensation object or spillover matrix or a list of compensation
  objects

- updateChannelNames:

  if TRUE, add a 'comp-' prefix to all fluorochrome channels (hence does
  not impact the columns related to FSC, SSC, or other specific keyword
  like TIME, Original_ID, File,...) Default TRUE.

## Value

a new object with compensated data, and possibly updated column names

## Examples

``` r

data(OMIP021Samples)

ff <- OMIP021Samples[[1]]
compMatrix <- flowCore::spillover(ff)$SPILL
ff <- runCompensation(ff, 
                      spillover = compMatrix, 
                      updateChannelNames = TRUE)
```
