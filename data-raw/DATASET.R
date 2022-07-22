## code to prepare `DATASET` dataset goes here

# flow cytometry core librairies
library(flowCore)
library(flowWorkspace)
library(cytophil)

dirDataset <- paste0("C:/CBIO/2021-phd-philippe-hauchamps/Datasets/OMIP-021/",
                     "FlowRepository_FR-FCM-ZZ9H_files")

files <- list.files(dirDataset, pattern = "Donor", recursive = TRUE)

nDonors <- length(files)

# init Flow Set

OMIP021Samples <- read.flowSet(
  paste0(dirDataset, "/", files), transformation = "linearize",
  alter.names = FALSE,
  name = "OMIP-21",
  truncate_max_range = TRUE,
  min.limit = NULL)

pData(OMIP021Samples)$Donor = paste0("Donor_", seq(nDonors))

# sub-sample equal nb of events in each fcs
sampleSize <- 100000
OMIP021Samples <- fsApply(x = OMIP021Samples,
                          FUN = function(ff){
                            FF_subsample(ff, nSamples = sampleSize, seed = 1)
                          })



usethis::use_data(OMIP021Samples, overwrite = TRUE)

fsApply(x = OMIP021Samples,
        FUN = function(ff){
          basefilename <- paste0("sample_",
                             flowCore::identifier(ff))
          path <- system.file("extdata",
                              package = "cytophil")

          flowCore::keyword(ff)[["$FIL"]] <- basefilename
          flowCore::keyword(ff)[["FILENAME"]] <- basefilename
          flowCore::write.FCS(ff,
                              filename = paste0(path,"/", basefilename))
        })

