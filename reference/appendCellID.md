# append 'Original_ID' column to a flowframe

: on a flowCore::flowFrame, append a 'Original_ID' column. This column
can be used in plots comparing the events pre and post gating. If the
'Original_ID' column already exists, the function does nothing

## Usage

``` r
appendCellID(ff, eventIDs = seq_len(flowCore::nrow(ff)))
```

## Arguments

- ff:

  a flowCore::flowFrame

- eventIDs:

  an integer vector containing the values to be added in expression
  matrix, as Original ID's.

## Value

new flowCore::flowFrame containing the added 'Original_ID' column

## Examples

``` r

data(OMIP021Samples)

retFF <- appendCellID(OMIP021Samples[[1]])
```
