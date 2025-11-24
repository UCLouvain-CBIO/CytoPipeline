# find time channel in flowSet/flowFrame

tries to find a channel in a flowSet/flowFrame that could be the time
channel. First tries to identify a channel name containing the 'time'
string, then tries to identify a single monotonically increasing
channel.

## Usage

``` r
findTimeChannel(obj, excludeChannels = c())
```

## Arguments

- obj:

  a flowCore::flowFrame or flowCore::flowSet

- excludeChannels:

  vector of column names to exclude in the search

## Value

a character, name of the found channel that should be representing time.
If not found, returns NULL.

## Examples

``` r

data(OMIP021Samples)

ret <- findTimeChannel(OMIP021Samples[[1]])
ret # "Time"
#> [1] "Time"
```
