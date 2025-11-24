# estimates scale tranformations

this function estimates the scale transformations to be applied on a
flowFrame to obtain 'good behaving' distributions, i.e. the best
possible separation between + population and - population. It
distinguishes between scatter channels, where either linear, or no
transform is applied, and fluo channels, where either logicle transform

- using flowCore::estimateLogicle - is estimated, or no transform is
  applied.

The idea of linear transform of scatter channels is as follows: a
reference channel (not a scatter one) is selected and a linear transform
(Y = AX + B) is applied to all scatter channel, as to align their 5 and
95 percentiles to those of the reference channel For the estimateLogicle
function, see flowCore documentation.

## Usage

``` r
estimateScaleTransforms(
  ff,
  fluoMethod = c("estimateLogicle", "none"),
  scatterMethod = c("none", "linearQuantile"),
  scatterRefMarker = NULL,
  specificScatterChannels = NULL,
  verbose = FALSE
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- fluoMethod:

  method to be applied to all fluo channels

- scatterMethod:

  method to be applied to all scatter channels

- scatterRefMarker:

  the reference channel that is used to align the

- specificScatterChannels:

  vector of scatter channels for which we still want to apply the fluo
  method (and not the scatter Method)

- verbose:

  if TRUE, send messages to the user at each step

## Value

a flowCore::flowFrame with removed low quality events from the input

## Examples

``` r

data(OMIP021Samples)

compMatrix <- flowCore::spillover(OMIP021Samples[[1]])$SPILL
ff_c <- runCompensation(OMIP021Samples[[1]], spillover = compMatrix)

transList <- 
    estimateScaleTransforms(        
        ff = ff_c,
        fluoMethod = "estimateLogicle",
        scatterMethod = "linear",
        scatterRefMarker = "BV785 - CD3")
```
