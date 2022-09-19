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


## code to prepare `OMIP021` dataset

# flow cytometry core librairies
library(flowCore)
library(CytoPipeline)

# in the following statement, specify the folder containing the dataset 
# manually downloaded from FlowRepository.org
# note that automatic download API for FlowRepository does not work anymore 
# for the time being (May 2022)
datasetPath <- "/Datasets/OMIP-021/FlowRepository_FR-FCM-ZZ9H_files"
datasetPath <- "C:/CBIO/2021-phd-philippe-hauchamps/Datasets/OMIP-021/FlowRepository_FR-FCM-ZZ9H_files"


files <- list.files(datasetPath, pattern = "Donor", recursive = TRUE)

nDonors <- length(files)

# init Flow Set

OMIP021Samples <- read.flowSet(
    paste0(datasetPath, "/", files),
    transformation = "linearize",
    alter.names = FALSE,
    name = "OMIP-21",
    truncate_max_range = TRUE,
    min.limit = NULL
)




pData(OMIP021Samples)$Donor <- paste0("Donor_", seq(nDonors))

# sub-sample equal nb of events in each fcs
sampleSize <- 5000
OMIP021Samples <- fsApply(
    x = OMIP021Samples,
    FUN = function(ff) {
        subsample(ff, nSamples = sampleSize, seed = 1)
    }
)



usethis::use_data(OMIP021Samples, overwrite = TRUE)

fsApply(
    x = OMIP021Samples,
    FUN = function(ff) {
        basefilename <- paste0(
            "sample_",
            flowCore::identifier(ff)
        )
        path <- system.file("extdata",
            package = "CytoPipeline"
        )

        flowCore::keyword(ff)[["$FIL"]] <- basefilename
        flowCore::keyword(ff)[["FILENAME"]] <- basefilename
        flowCore::write.FCS(ff,
            filename = paste0(path, "/", basefilename)
        )
    }
)
