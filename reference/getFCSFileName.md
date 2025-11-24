# get fcs file name

get basename of \$FILENAME keyword if exists

## Usage

``` r
getFCSFileName(ff)
```

## Arguments

- ff:

  a flowCore::flowFrame

## Value

the basename of \$FILENAME keyword

## Examples

``` r

data(OMIP021Samples)

fName <- getFCSFileName(OMIP021Samples[[1]])
```
