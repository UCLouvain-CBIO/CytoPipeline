# Cyto Processing step

Class containing the function and arguments to be applied in a
lazy-execution framework.

Objects of this class are created using the `CytoProcessingStep()`
function. The processing step is executed with the
`executeProcessingStep()` function.

## Usage

``` r
CytoProcessingStep(name = character(), FUN = character(), ARGS = list())

# S4 method for class 'CytoProcessingStep'
show(object)

executeProcessingStep(x, ...)

getCPSName(x)

getCPSFUN(x)

getCPSARGS(x)

# S3 method for class 'CytoProcessingStep'
as.list(x, ...)

as.json.CytoProcessingStep(x, pretty = FALSE)

from.json.CytoProcessingStep(jsonString)
```

## Arguments

- name:

  `character` denoting a name to the step, which can be different from
  the function name

- FUN:

  `function` or `character` representing a function name.

- ARGS:

  `list` of arguments to be passed along to `FUN`.

- object:

  a `CytoProcessingStep` object.

- x:

  a `CytoProcessingStep` object.

- ...:

  other arguments (not used)

- pretty:

  formatting set-up (see jsonlite::toJSON doc)

- jsonString:

  a [`character()`](https://rdrr.io/r/base/character.html) containing a
  JSON string.

## Value

The `CytoProcessingStep` function returns and object of type
`CytoProcessingStep`.

## Details

This object contains all relevant information of a data analysis
processing step, i.e. the function and all of its arguments to be
applied to the data.

## Examples

``` r

## Create a simple processing step object
ps1 <- CytoProcessingStep("summing step", sum)

getCPSName(ps1)
#> [1] "summing step"

getCPSFUN(ps1)
#> function (..., na.rm = FALSE)  .Primitive("sum")

getCPSARGS(ps1)
#> list()

executeProcessingStep(ps1, 1:10)
#> [1] 55

as.list(ps1)
#> $name
#> [1] "summing step"
#> 
#> $FUN
#> [1] "sum"
#> 
#> $ARGS
#> list()
#> 

js_str <- as.json.CytoProcessingStep(ps1)

ps2 <- from.json.CytoProcessingStep(js_str)

identical(ps1, ps2)
#> [1] FALSE
```
