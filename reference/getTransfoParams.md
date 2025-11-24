# get tranformation parameters for a specific channel

investigates a flowCore::tranformList object to get the type and
parameters of the transformation applying to a specific channel

## Usage

``` r
getTransfoParams(transList, channel)
```

## Arguments

- transList:

  a flowCore::transformList

- channel:

  channel name

## Value

If the transformation exists for the specified channel, and is either
recognized as a logicle transfo or a linear transfo, a list with two
slots:

- \$type a character containing the transfo type ('logicle' or 'linear')

- \$params_list a list of named numeric, according to transfo type

Otherwise, NULL is returned.

## Examples

``` r

data(OMIP021Samples)

# set-up a hybrid transformation list :
# - two channels are logicle-ly transformed with automatic param estimates
# - one channel has explicit logicle transfo with default parameters
# - one channel has linear transformation
# - other channels have no transformation
translist <- flowCore::estimateLogicle(
    OMIP021Samples[[1]],
    c("450/50Violet-A", "525/50Violet-A")
)
translist <- c(
    translist,
    flowCore::transformList(
        "FSC-A",
        flowCore::linearTransform(
            a = 0.1,
            b = 0
       )
    ),
    flowCore::transformList(
        "540/30Violet-A",
        flowCore::logicleTransform()
    )
)

ret1 <- getTransfoParams(translist, channel = "FSC-A")
ret1$type # "linear"
#> [1] "linear"
ret1$paramsList # a = 0.1, b = 0.
#> $a
#> [1] 0.1
#> 
#> $b
#> [1] 0
#> 

ret2 <- getTransfoParams(translist, channel = "525/50Violet-A")
ret2$type # "logicle"
#> [1] "logicle"
ret2$paramsList # a = 0., w = 0.2834, m = 4.5, t = 262143
#> $a
#> [1] 0
#> 
#> $w
#> [1] 0.2628157
#> 
#> $m
#> [1] 4.5
#> 
#> $t
#> [1] 262143
#> 

ret3 <- getTransfoParams(translist, channel = "540/30Violet-A")
ret3$type # "logicle
#> [1] "logicle"
ret3$paramsList # a = 0., w = 0.5, m = 4.5, t = 262144
#> $a
#> [1] 0
#> 
#> $w
#> [1] 0.5
#> 
#> $m
#> [1] 4.5
#> 
#> $t
#> [1] 262144
#> 
```
