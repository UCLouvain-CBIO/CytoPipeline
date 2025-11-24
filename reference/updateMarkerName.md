# update marker name of a given flowFrame channel

: in a flowCore::flowFrame, update the marker name (stored in 'desc' of
parameters data) of a given channel. Also update the corresponding
keyword in the flowFrame.

## Usage

``` r
updateMarkerName(ff, channel, newMarkerName)
```

## Arguments

- ff:

  a flowCore::flowFrame

- channel:

  the channel for which to update the marker name

- newMarkerName:

  the new marker name to be given to the selected channel

## Value

a new flowCore::flowFrame with the updated marker name

## Examples

``` r

data(OMIP021Samples)

retFF <- updateMarkerName(OMIP021Samples[[1]],
                          channel = "FSC-A",
                          newMarkerName = "Fwd Scatter-A")
```
