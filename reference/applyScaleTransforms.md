# apply scale transforms

wrapper around flowCore::transform() that discards any additional
parameter passed in (...) Additionally, some checks regarding channels
correspondance is done: if `transList` contains transformations for
channels that are not present in `x`, then these transformations are
first removed.

## Usage

``` r
applyScaleTransforms(x, transList, verbose = FALSE, ...)
```

## Arguments

- x:

  a flowCore::flowSet or a flowCore::flowFrame

- transList:

  a flowCore::transformList

- verbose:

  if TRUE, send a message per flowFrame transformed

- ...:

  other arguments (not used)

## Value

the transformed flowFrame

## Examples

``` r

data(OMIP021Samples)

transListPath <- file.path(system.file("extdata", 
                                       package = "CytoPipeline"),
                           "OMIP021_TransList.rds") 

transList <- readRDSObject(transListPath)

ff_c <- compensateFromMatrix(OMIP021Samples[[1]],
                             matrixSource = "fcs")  

ff_t <- applyScaleTransforms(ff_c, transList = transList)
```
