# CytoPipeline - Copyright (C) <2022>
# <Université catholique de Louvain (UCLouvain), Belgique>
#
#   Description and complete License: see LICENSE file.
#
# This program (CytoPipeline) is free software:
#   you can redistribute it and/or modify it under the terms of the GNU General
#   Public License as published by the Free Software Foundation, either
#   version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (<http://www.gnu.org/licenses/>).


## code to prepare `OMIP021UT` dataset (light version of `OMIP021` used in
## unit testing (testthat))

library(CytoPipeline)

data(OMIP021Samples)

# sub-sample equal nb of events in each fcs
sampleSize <- 100
OMIP021UTSamples <- fsApply(
    x = OMIP021Samples,
    FUN = function(ff) {
        subsample(ff, nSamples = sampleSize, seed = 1)
    }
)
