# compute linear transformation of scatter channels found in ff, based on 5% and 95% of referenceChannel, set as target. If there is a transformation defined in transList for referenceChannel, it is applied first, before computing quantiles. Then the computed linear transformations (or each scatter channel) are added into the transfo_list. -A channels are computed, and same linear transformation is then applied to corresponding -W and -H channels (if they exist in ff).

based on a referenceChannel

## Usage

``` r
computeScatterChannelsLinearScale(
  ff,
  transList = NULL,
  referenceChannel,
  silent = TRUE
)
```

## Arguments

- ff:

  a flowCore::flowFrame

- transList:

  an initial flowCore::transformList

- referenceChannel:

  the reference channel to take target quantile values from. Can be
  defined as marker or channel name.

- silent:

  if FALSE, will output some information on the computed linear
  transformations

## Value

the transList with added linear scale transformations

## Examples

``` r

data(OMIP021Samples)

ff <- OMIP021Samples[[1]]
refMarker <- "APCCy7 - CD4"
refChannel <- "780/60Red-A"
transList <- flowCore::estimateLogicle(ff,
                                       channels = refChannel)
retTransList <-
    computeScatterChannelsLinearScale(ff,
                                      transList = transList,
                                      referenceChannel = refMarker,
                                      silent = TRUE
    )
```
