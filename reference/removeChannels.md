# remove channels from a flowFrame

: in a flowCore::flowFrame, remove the channels of the given names.

## Usage

``` r
removeChannels(ff, channels)
```

## Arguments

- ff:

  a flowCore::flowFrame

- channels:

  the channel names to be removed

## Value

a new flowCore::flowFrame with the removed channels

## Examples

``` r

data(OMIP021Samples)

retFF <- removeChannels(OMIP021Samples[[1]],
                        channel = "FSC-A")
```
