# read RDS object

wrapper around readRDS, which discards any additional parameters passed
in (...)

## Usage

``` r
readRDSObject(RDSFile, ...)
```

## Arguments

- RDSFile:

  a RDS file containing a R object object

- ...:

  other arguments (not used)

## Value

the read R object

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
