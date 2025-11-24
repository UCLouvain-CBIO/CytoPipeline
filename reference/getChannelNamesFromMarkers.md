# get channel names from markers

finds name of channels corresponding to user provided markers

## Usage

``` r
getChannelNamesFromMarkers(ff, markers)
```

## Arguments

- ff:

  a flowCore::flowFrame

- markers:

  a vector of markers, either provided as :

  - an array of booleans (referring to flowFrame columns)

  - an array of integers (indices in flowFrame columns)

  - an array of characters (exact markers or channel patterns)

## Value

a character vector, containing the names of the corresponding channels

## Examples

``` r

data(OMIP021Samples)

# with existing markers
ret <- getChannelNamesFromMarkers(
    OMIP021Samples[[1]],
    c(
        "FSC-A",
        "L/D Aqua - Viability",
        "FITC - gdTCR",
        "PECy5 - CD28"
    ))
    
ret # c("FSC-A", "525/50Violet-A", "530/30Blue-A", "670/30Yellow-A")
#> [1] "FSC-A"          "525/50Violet-A" "530/30Blue-A"   "670/30Yellow-A"

# with boolean vector
indices <- c(1, 6, 14, 18)
boolInput <- rep(FALSE, 21)
boolInput[indices] <- TRUE
ret2 <- getChannelNamesFromMarkers(
    OMIP021Samples[[1]],
    boolInput)
    
ret2 # c("FSC-A", "525/50Violet-A", "530/30Blue-A", "670/30Yellow-A")
#> [1] "FSC-A"          "525/50Violet-A" "530/30Blue-A"   "670/30Yellow-A"

# with indices vector
ret3 <- getChannelNamesFromMarkers(
    OMIP021Samples[[1]],
    indices
)
ret3 # c("FSC-A", "525/50Violet-A", "530/30Blue-A", "670/30Yellow-A")
#> [1] "FSC-A"          "525/50Violet-A" "530/30Blue-A"   "670/30Yellow-A"

```
