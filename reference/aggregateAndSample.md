# Aggregate and sample multiple flow frames of a flow set together

Aggregate multiple flow frames in order to analyze them simultaneously.
A new FF, which contains about nTotalEvents cells, nTotalEvents/nFiles
cells from each file. Two new columns are added: a column indicating the
original file by index, and a noisy version of this, for better plotting
opportunities, This function is based on PeacoQC::AggregateFlowframes()
where file names inputs have been replaced by a flowSet input.

## Usage

``` r
aggregateAndSample(
  fs,
  nTotalEvents,
  setup = c("forceNEvent", "forceBalance"),
  seed = NULL,
  channels = NULL,
  writeOutput = FALSE,
  outputFile = "aggregate.fcs",
  keepOrder = FALSE
)
```

## Arguments

- fs:

  a flowCore::flowset

- nTotalEvents:

  Total number of cells to select from the input flow frames

- setup:

  How to proceed when nTotalEvents/nFiles is too high for some of the
  flow frames:

  - forceBalance (default): compute the minimum nb of events per flow
    frame, and keep that amount of events from each flow frame.

  - forceNEvents: try to be as balanced as possible, but force a total
    of nTotalEvents if possible, i.e. takes all events from the flow
    frame with too low nb of events, and then fill in the total with
    events from the bigger flow frames in a balanced way. However, if
    nTotalEvents is greater than the sum of all events, take all events
    only once.

- seed:

  seed to be set before sampling for reproducibility. Default NULL does
  not set any seed.

- channels:

  Channels/markers to keep in the aggregate. Default NULL takes all
  channels of the first file.

- writeOutput:

  Whether to write the resulting flowframe to a file. Default FALSE

- outputFile:

  Full path to output file. Default "aggregate.fcs"

- keepOrder:

  If TRUE, the random subsample will be ordered in the same way as they
  were originally ordered in the file. Default = FALSE.

## Value

returns a new flowCore::flowFrame

## Examples

``` r

data(OMIP021Samples)

nCells <- 1000
agg <- aggregateAndSample(
    fs = OMIP021Samples,
    nTotalEvents = nCells)
```
