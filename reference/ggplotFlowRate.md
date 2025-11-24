# plot flow rate as a function of time, using ggplot2

plot flow rate as a function of time, using ggplot2

## Usage

``` r
ggplotFlowRate(obj, title = "Flow Rate", timeUnit = 100)
```

## Arguments

- obj:

  a flowCore::flowFrame or flowCore::flowSet

- title:

  a title for the graph

- timeUnit:

  which time interval is used to calculate "instant" flow rate (default
  = 100 ms)

## Value

a ggplot graph

## Examples

``` r

data(OMIP021Samples)

# single flowFrame plot
ggplotFlowRate(OMIP021Samples[[1]])


# two flowFrames plot 
ggplotFlowRate(OMIP021Samples)


# single plot with title
ggplotFlowRate(OMIP021Samples[[1]], title = "Test Flow Rate plot")


# explicit time unit
ggplotFlowRate(OMIP021Samples[[1]], timeUnit = 50)

```
