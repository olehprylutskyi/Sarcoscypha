rm(list = ls()) # Reset R`s brain

# Optional: set working directory (change PATH to yours)
# No need if RStudio's Rproj is used
# setwd("~/R/sarcoscypha-urnula_legacy/Maxent_Sarco_Ukraine_2023")

datadir <- "./data"
outputdir <- "./outputs"
# set minimal value for coordinate uncertainty in meters
minCoordPrecision <- 200

# Loading libraries
library(dplyr)
library(tidyr)

# Read Ukrainian citizen science data
ukraine_cs <- read.csv(
  paste0(datadir, "/raw_occ_data/Citizen_science_Sarco-2025-05-19.csv")
) %>%
  dplyr::select(
    coordinateUncertaintyInMeters, 
    Longitude, 
    Latitude, 
    Source
  ) %>%
  as.data.frame() %>%
  filter(
    coordinateUncertaintyInMeters <= minCoordPrecision
  )

colnames(ukraine_cs) <- c("coordUncert", "x", "y", "source")

# Read iNaturalist data from local download
inat_ua <- read.csv(
  paste0(datadir, "/raw_occ_data/inat_observations-2025-05-19.csv")
) %>%
  dplyr::select(
    positional_accuracy, longitude, latitude
  ) %>%
  as.data.frame() %>%
  filter(
    is.na(positional_accuracy) | 
      positional_accuracy <= minCoordPrecision
  ) %>%
  mutate(Source = "iNaturalist")

# Rename columns
colnames(inat_ua) <- c("coordUncert", "x", "y", "source")

# Read GBIF data from local download
gbif_ua <- read.csv(
  paste0(datadir, "/raw_occ_data/gbif_observations-2025-05-19.csv")
) %>% 
  dplyr::select(
    coordinateUncertaintyInMeters, 
    decimalLongitude, 
    decimalLatitude
  ) %>%
  as.data.frame() %>%
  filter(
    is.na(coordinateUncertaintyInMeters) | 
      coordinateUncertaintyInMeters <= minCoordPrecision
  ) %>%
  mutate(Source = "GBIF")

colnames(gbif_ua) <- c("coordUncert", "x", "y", "source")

# Read PlutoF data from local download
plutof_ua <- read.csv(
  paste0(datadir, "/raw_occ_data/plutof_data_2025-05-19.csv")
) %>%
  dplyr::select(
    Sampling.area.Accuracy, 
    Sampling.area.Longitude, 
    Sampling.area.Latitude) %>%
  as.data.frame() %>%
  filter(
    is.na(Sampling.area.Accuracy) | 
      Sampling.area.Accuracy <= minCoordPrecision
  ) %>%
  mutate(Source = "PlutoF")

colnames(plutof_ua) <- c("coordUncert", "x", "y", "source")

# Merge all data and calculate number of occurrences by source
occ <- rbind(ukraine_cs, inat_ua, gbif_ua, plutof_ua)

input_datasources <- as.data.frame.table(
  table(occ$source)
) %>%
  rename(
    Source = Var1,
    NumObs = Freq
  )

input_datasources

# Export results to csv file
write.csv(
  input_datasources, 
  file = paste0(outputdir, "/input_data_sources.csv")
)

# Prepare data for modeling
# How many data points has precise coordinates?
hist(occ$coordUncert)


occ %>%
  filter(!is.na(coordUncert)) %>% 
  mutate(scientific_name = "Sarcoscypha sp") %>%
  relocate(scientific_name) %>% 
  write.csv(
    file = paste0(outputdir, "/sarcoscypha_occ.csv"),
    row.names = FALSE
  )

# Clean the environment
rm(list = ls())

# End of the script ####
