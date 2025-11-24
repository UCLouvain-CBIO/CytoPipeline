# find flow frame columns that represent true signal

: find flow frame columns that represent true signal

## Usage

``` r
areSignalCols(
  x,
  toRemovePatterns = c("Time", "Original_ID", "File", "SampleID")
)
```

## Arguments

- x:

  a flowCore::flowFrame or a flowCore::flowSet

- toRemovePatterns:

  a vector of string patterns that are to be considered as non signal

## Value

a vector of booleans of which the dimension is equal to the number of
columns in ff

## Examples

``` r

data(OMIP021Samples)

areSignalCols(OMIP021Samples)
#>          FSC-A          FSC-H          SSC-A          SSC-H 450/50Violet-A 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> 525/50Violet-A 540/30Violet-A 585/15Violet-A 610/20Violet-A 670/30Violet-A 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#>    670/14Red-A   730//45Red-A    780/60Red-A   530/30Blue-A   710/50Blue-A 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#> 582/15Yellow-A 610/20Yellow-A 670/30Yellow-A 710/50Yellow-A 780/60Yellow-A 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#>           Time    Original_ID 
#>          FALSE          FALSE 
```
