# sub-sampling of a flowFrame

: sub-samples a flowFrame with the specified number of samples, without
replacement. adds also a column 'Original_ID' if not already present in
flowFrame.

## Usage

``` r
subsample(ff, nEvents, seed = NULL, keepOriginalCellIDs = TRUE, ...)
```

## Arguments

- ff:

  a flowCore::flowFrame

- nEvents:

  number of events to be obtained using sub-sampling

- seed:

  can be set for reproducibility of event sub-sampling

- keepOriginalCellIDs:

  if TRUE, adds (if not already present) a 'OriginalID' column
  containing the initial IDs of the cell (from 1 to nrow prior to
  subsampling). if FALSE, does the same, but takes as IDs (1 to nrow
  after subsampling)

- ...:

  additional parameters (currently not used)

## Value

new flowCore::flowFrame with the obtained subset of samples

## Examples

``` r

data(OMIP021Samples)

# take first sample of dataset, subsample 100 events and create new flowFrame
ff <- subsample(OMIP021Samples[[1]], nEvents = 100)

```
