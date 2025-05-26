# Prepare the environment ####
rm(list = ls()) # Reset R`s brain

# Optional: set working directory (change PATH to yours)
# No need if RStudio's Rproj is used
# setwd("~/R/sarcoscypha-urnula_legacy/Maxent_Sarco_Ukraine_2023")

datadir <- "./data"
outputdir <- "./outputs"
# modeldir <- "./models/"
# evaldir <- paste0(modeldir, "evaluations/")
# Covariates (climate + land cover)
rasterdir <- "/CHELSA_1981-2010_V.2.1"
# # Averaged climate scenaria
# rasterdir_126 <- "/CHELSA_2041-2070_average_ssp126_V.2.1"
# rasterdir_370 <- "/CHELSA_2041-2070_average_ssp370_V.2.1"
# rasterdir_585 <- "/CHELSA_2041-2070_average_ssp585_V.2.1"

# the minimal distance that you want records to be separated by
thinDist <- 5 # in kilometers

# Loading libraries
# library(dismo)
library(sp)
library(sf)
# library(rJava)
# library(gbm)
library(raster)
library(usdm) # Variance Inflation Factor
library(corrplot)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load state boundary of Ukraine
load(file = paste0(datadir, "/Ukraine.Rdata"))

# Mask water bodies
waterbodies <- sf::st_read(
  paste0(datadir, "/waterbodies_ua.shp")
) %>% 
  filter(!st_is_empty(.))

waterbodies <- as(waterbodies, "Spatial")

# Covariates ####
# Load input rasters
files <- list.files(
  path =  paste0(datadir, rasterdir), 
  pattern = "tif"
)

rasPath <- file.path(
  paste0(datadir, rasterdir),
  files
)

covariates <- raster::stack(rasPath)

layer_names <- sub(".tif", "", files) # drop ".tif" part from the filenames

names(covariates) <- layer_names # rename layers (get band1 names by default)

# Mask covariates by Ukraine boundary polygon
covariates <- mask(covariates, Ukraine)

# Mask covariates by waterbodies spatial polygon
covariates <- mask(covariates, waterbodies, inverse = TRUE)


# Explore covariates ####
# Plot all covariates to explore
for (i in 1:nlayers(covariates)) {
  png(
    paste0(
      outputdir, 
      "/covariates_plots/",
      names(covariates[[i]]),
      ".png"
    ),
    width = 16, height = 12, units = "cm" , res = 150
  )
  plot(covariates[[i]], main = names(covariates[[i]]))
  dev.off()
}

# Drop bioclim layers containing visible issues
covariates <- dropLayer(covariates, c("bio08", "bio09"))

names(covariates)



# ------------------

# ## Variable selection
# included_vars <- c(
#   "bio05", # Max Temperature of Warmest Month
#   "bio06", # Min Temperature of Coldest Month
#   "bio15", # Precipitation Seasonality (Coefficient of Variation)
#   "bio18", # Precipitation of Warmest Quarter
#   "bio19", # Precipitation of Coldest Quarter
#   "TREES_BD", # Broadleaved deciduous tree cover percentage (ESA landcover)
#   "TREES_NE" # Needleleaved evergreen tree cover percentage (ESA landcover)
# )

# covs <- subset(covariates, included_vars)

# names(covs)

# Occurrence data ####
# Load text data
occ <- read.csv(paste0(outputdir, "/sarcoscypha_occ.csv"))

# delete service 1st var with counter and last var with sources
occ <- occ[c("scientific_name", "x", "y")]

# subset to just records with latitude and longitude
occ <- occ[!is.na(occ$x) & !is.na(occ$y),]

# round longitude and latitude with 5 digits
occ$x <- round(occ$x, 5)    
occ$y <- round(occ$y, 5)

# make it spatial
coordinates(occ) <- ~ x + y
occ@proj4string <- CRS(SRS_string = "EPSG:4326")
plot(occ)

# proj4string(occ) <- CRS("+init=epsg:4326") # this is WGS84

# Remove duplicates
occ <- remove.duplicates(occ)

# Filter occurrences by bounding box
occ <- crop(occ, Ukraine) # crop occs by polygon extent

## Remove duplicates and NAs
# Extract predictors` values from raster stack
occs_vals_Ss <- as.data.frame(
  raster::extract(covariates, occ@coords, cellnumbers = TRUE)
)

# Remove duplicated same cell values
occ <- occ[!duplicated(occs_vals_Ss[, 1]), ] # from occ data
occs_vals_Ss <- occs_vals_Ss[!duplicated(occs_vals_Ss[, 1]), -1] # from extracted predictors

# remove occurrence records with NA environmental values
occ <- occ[!(rowSums(is.na(occs_vals_Ss)) >= 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Ss <- na.omit(occs_vals_Ss)
# add columns for env variable values for each occurrence record
occs_Ss <- cbind(occ, occs_vals_Ss)

## Spatial thinning

occs <- cbind(occs_Ss@data, occs_Ss@coords)
output <- spThin::thin(
  loc.data = occs, lat.col = "y",
  long.col = "x", spec.col = 'scientific_name',
  thin.par = thinDist, reps = 100,
  locs.thinned.list.return = TRUE, write.files = FALSE,
  write.log.file = FALSE, verbose = FALSE
)

# pull thinned dataset with max records, not just the first in the list
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
# if more than one max, pick first
maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]
presvals <- occs[as.numeric(rownames(maxThin)),]
# make it spatial back
coordinates(presvals) <- ~ x + y
# proj4string(presvals) <- CRS("+init=epsg:4326") # this is WGS84
presvals@proj4string <- CRS(SRS_string = "EPSG:4326")
# drop "scientific_name" variable
presvals <- presvals[-1]

nrow(presvals)
View(presvals@data)

# Set number of background points simulated. as 5 times num. of presence points
bgNumPoints <- nrow(presvals@data) * 100
backgr <- spsample (
  Ukraine,
  bgNumPoints,
  type = "random"
) 

## Remove duplicates and NAs from the background
# Extract predictors` values from raster stack
backgr_vals_Ss <- as.data.frame(
  raster::extract(covariates, backgr@coords, cellnumbers = TRUE)
)

# Remove duplicated same cell values
# from background spatial points
backgr <- backgr[!duplicated(backgr_vals_Ss[, 1]), ]
# from extracted predictors
backgr_vals_Ss <- backgr_vals_Ss[!duplicated(backgr_vals_Ss[, 1]), -1]

# remove points with NA environmental values
backgr <- backgr[!(rowSums(is.na(backgr_vals_Ss)) >= 1), ]
# also remove variable value rows with NA environmental values
backgr_vals_Ss <- na.omit(backgr_vals_Ss)

# add columns for env variable values for each background record
# Since no slot of name "data" for backgr {sp} object, we need additional steps
backgvals <- cbind(backgr_vals_Ss, backgr@coords)
# make it spatial back
coordinates(backgvals) <- ~ x + y
# proj4string(backgvals) <- CRS("+init=epsg:4326") # this is WGS84
backgvals@proj4string <- CRS(SRS_string = "EPSG:4326")

nrow(backgvals)

View(backgvals@data)
# -----------------

# Collinearity control ####
# Define and exclude highly correlated variables
# Calculate the correlation matrix for the numeric columns
cor_matrix <- cor(
  backgvals@data,
  use = "complete.obs",
  method = "pearson"
)

# Plot correlation among all predictor variables
png(
  paste0(
    outputdir,
    "/covariate_selection/corrplot_all.png"
  ),
  width = 36, height = 32, units = "cm", res = 150
)
corrplot(
  cor_matrix,
  method = "color",            # Use colored squares for correlation
  type = "upper",              # Show upper triangle only
  order = "hclust",            # Reorder variables hierarchically
  addCoef.col = "black",       # Show correlation coefficients in black
  number.cex = 0.5,            # Reduce the size of correlation labels
  tl.col = "black",            # Text label color
  tl.srt = 30,                 # Rotate labels slightly for readability
  tl.cex = 0.5,                # Reduce text size of variable labels (set smaller valu)
  cl.cex = 0.8,                # Reduce text size of color legend
  diag = FALSE,                # Hide diagonal
  col = colorRampPalette(c("#11aa96", "#61c6fa", "#f6aa70"))(200),
  sig.level = 0.01, insig = "blank"
)
dev.off()

# Variance Inflation Factor
v <- vifstep(backgvals@data, th = 10)
v

# ---------- VIFs of the remained variables -------- 
#     Variables      VIF
# 1       bio02 2.231996
# 2       bio15 4.807178
# 3       bio19 2.824250
# 4        EVIs 6.563164
# 5        EVIw 2.975911
# 6         fcf 2.552801
# 7         gsl 3.264493
# 8   hurs_mean 4.163897
# 9  hurs_range 2.831591
# 10      NDWIs 4.905632
# 11      NDWIw 1.609681
# 12        scd 4.258895
# 13   TREES_BD 1.954181
# 14   TREES_NE 2.307035

# List remained variables
remained_vars <- v@results[1] %>% pull()
# [1] "bio02"      "bio15"      "bio19"      "EVIs"       "EVIw"       "fcf"
# [7] "gsl"        "hurs_mean"  "hurs_range" "NDWIs"      "NDWIw"      "scd"
# [13] "TREES_BD"   "TREES_NE"




# library(PerformanceAnalytics)
# png(paste0(figuredir, "/correlation_covariates.png"), width = 20, height = 20, units = "cm", res = 300)
# chart.Correlation(occ.pr, histogram = TRUE, method = "pearson")
# dev.off()




# Plot correlation among remained predictor variables
png(
  paste0(
    outputdir,
    "/covariate_selection/corrplot_vif.png"
  ),
  width = 16, height = 12, units = "cm", res = 150
)
corrplot(
  cor(
    backgvals@data[, remained_vars], # exclude vars based on VIF
    use = "complete.obs", 
    method = "pearson"
  ),
  method = "color",            # Use colored squares for correlation
  type = "upper",              # Show upper triangle only
  order = "hclust",            # Reorder variables hierarchically
  addCoef.col = "black",       # Show correlation coefficients in black
  number.cex = 0.5,            # Reduce the size of correlation labels
  tl.col = "black",            # Text label color
  tl.srt = 30,                 # Rotate labels slightly for readability
  tl.cex = 0.5,                # Reduce text size of variable labels (set smaller valu)
  cl.cex = 0.8,                # Reduce text size of color legend
  diag = FALSE,                # Hide diagonal
  col = colorRampPalette(c("#11aa96", "#61c6fa", "#f6aa70"))(200),
  sig.level = 0.01, insig = "blank"
)
dev.off()

# High Pearson correalation among remaining variables
# EVIs / NDWIs = 0.85


# Drop bioclim layers containing visible issues
covariates <- dropLayer(covariates, c(v@excluded, "NDWIs"))

names(covariates)

# Plot correlation among remained predictor variables
png(
  paste0(
    outputdir,
    "/covariate_selection/corrplot_vif_cor.png"
  ),
  width = 16, height = 12, units = "cm", res = 150
)
corrplot(
  cor(
    backgvals@data[, names(covariates)], # exclude vars based on VIF
    use = "complete.obs", 
    method = "pearson"
  ),
  method = "color",            # Use colored squares for correlation
  type = "upper",              # Show upper triangle only
  order = "hclust",            # Reorder variables hierarchically
  addCoef.col = "black",       # Show correlation coefficients in black
  number.cex = 0.5,            # Reduce the size of correlation labels
  tl.col = "black",            # Text label color
  tl.srt = 30,                 # Rotate labels slightly for readability
  tl.cex = 0.5,                # Reduce text size of variable labels (set smaller valu)
  cl.cex = 0.8,                # Reduce text size of color legend
  diag = FALSE,                # Hide diagonal
  col = colorRampPalette(c("#11aa96", "#61c6fa", "#f6aa70"))(200),
  sig.level = 0.01, insig = "blank"
)
dev.off()

selected_covs <- names(covariates)

print(names(covariates))
# [1] "bio02"      "bio15"      "bio19"      "EVIs"       "EVIw"       "fcf"       
# [7] "gsl"        "hurs_mean"  "hurs_range" "NDWIw"      "scd"        "TREES_BD"  
# [13] "TREES_NE" 

presvals@data <- presvals@data[, names(covariates)]

# Leave only selected covariates
presvals@data <- presvals@data[, names(presvals@data) %in% names(covariates)]
backgvals@data <- backgvals@data[, names(backgvals@data) %in% names(covariates)]

# Check it
View(presvals@data)
View(backgvals@data)

# Export fully prepared data
save(covariates, file = paste0(outputdir, "/maxent_ready/covariates.Rdata"))
save(presvals, file = paste0(outputdir, "/maxent_ready/presence.Rdata"))
save(backgvals, file = paste0(outputdir, "/maxent_ready/background.Rdata"))

# End of the script ####
