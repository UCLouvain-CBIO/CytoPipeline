# reset 'Original_ID' column in a flowframe

: on a flowCore::flowFrame, reset 'Original_ID' column. This column can
be used in plots comparing the events pre and post gating. If the
'Original_ID' column already exists, the function replaces the existing
IDs by the user provided ones. If not, an
[`appendCellID()`](https://uclouvain-cbio.github.io/CytoPipeline/reference/appendCellID.md)
is called.

## Usage

``` r
resetCellIDs(ff, eventIDs = seq_len(flowCore::nrow(ff)))
```

## Arguments

- ff:

  a flowCore::flowFrame

- eventIDs:

  an integer vector containing the values to be set in expression
  matrix, as Original ID's.

## Value

new flowCore::flowFrame containing the amended (or added) 'Original_ID'
column

## Examples

``` r

data(OMIP021Samples)

ff <- appendCellID(OMIP021Samples[[1]])

subsample_ff <- subsample(ff, 100, keepOriginalCellIDs = TRUE)

# re-create a sequence of IDs, ignoring the ones before subsampling
reset_ff <- resetCellIDs(subsample_ff)
```
