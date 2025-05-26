# Prepare the environment ####
rm(list = ls()) # Reset R`s brain

# Optional: set working directory (change PATH to yours)
# No need if RStudio's Rproj is used
# setwd("~/R/sarcoscypha-urnula_legacy/Maxent_Sarco_Ukraine_2023")

datadir <- "./data"
outputdir <- "./outputs"
modeldir <- "./models/"
evaldir <- paste0(modeldir, "evaluations/")

# Loading libraries
library(dismo)
# library(raster)
# library(sp)
library(rJava)
library(gbm)

# set a stable temporary directory where rasters won’t be deleted mid-session
rasterOptions()
# raster::rasterOptions(tmpdir = "/Users/oleh/R/sarcoscypha-urnula_legacy/Maxent_Sarco_Ukraine_2023/tempdir")
# raster::rasterOptions(default = TRUE)



# Maxent modeling (dismo) ####
# Load pre-processed covariates
load(file = paste0(outputdir, "/maxent_ready/covariates.Rdata"))
# Load pre-processed presence data
load(file = paste0(outputdir, "/maxent_ready/presence.Rdata"))
# Load pre-processed background data
load(file = paste0(outputdir, "/maxent_ready/background.Rdata"))


for (i in 1:100) {
  # Splitting data to train/test
  # Dividing into train and test subsamples for both data and background
  group <- kfold(presvals, 5) # randomly divide the data into 5 equal parts
  pres_train <- presvals[group != 1, ] # 80% (parts 1-5) as training data
  pres_test <- presvals[group == 1, ] # 20% (part 1) as testing data
  group1 <- kfold(backgvals, 5)
  backg_train <- backgvals[group1 != 1, ]
  backg_test <- backgvals[group1 == 1, ]
  
  # Build a maxent model
  m_maxent <- maxent(x = covariates,
                     p = pres_train,
                     a = backg_train)
  # Evaluate the model
  e_maxent <- evaluate(pres_test, backg_test, m_maxent, covariates)
  # Save the model
  save(m_maxent, file = paste0(modeldir, paste0("m", i), ".Rdata"))
  save(e_maxent, file = paste0(evaldir, paste0("e_m", i), ".Rdata"))
  gc() # force garbage collection
}

# End of the script ####