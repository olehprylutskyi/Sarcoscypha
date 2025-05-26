# This script adjusts rasters other than CHELSA to the same
# coordinate reference system and spatial resolution as CHELSA

# Vegetation indices (NDVI, NDWI, EVA) should be prepared in
# Google Earth Engine (https://earthengine.google.com/) code editor,
# then downloaded from Google Drive associated with the GEE account.

# ESA tree cover is supposed to be manually acquaired from the 
# https://catalogue.ceda.ac.uk/uuid/26a0f46c95ee4c29b5c650b129aab788/

# All those files must be stored in the ./rescaling_external_rasters
# directory before executing this script.

# Prepare the environment ####
rm(list = ls()) # Reset R`s brain

# Optional: set working directory (change PATH to yours)
# No need if RStudio's Rproj is used
# setwd("~/R/sarcoscypha-urnula_legacy/Maxent_Sarco_Ukraine_2023")

datadir <- "./data"

library(raster)
library(terra)

# Load state boundary of Ukraine
load(file = paste0(datadir, "/Ukraine.Rdata"))
# Convert to SpatVector for using as a mask in terra::mask
Ukraine <- terra::vect(Ukraine)

# rasters to be rescaled
ndvi_summer <- rast("./rescaling_external_rasters/MODIS_NDVI_Summer.tif")
plot(ndvi_summer)

# ndvi_autumn <- rast("./rescaling_external_rasters/MODIS_NDVI_Autumn_500m.tif")
# plot(ndvi_autumn)

ndvi_winter <- rast("./rescaling_external_rasters/MODIS_NDVI_Winter.tif")
plot(ndvi_winter)

ndwi_summer <- rast("./rescaling_external_rasters/MODIS_NDWI_Summer.tif")
plot(ndwi_summer)

# ndwi_autumn <- rast("./rescaling_external_rasters/MODIS_NDWI_autumn.tif")
# plot(ndwi_autumn)

ndwi_winter <- rast("./rescaling_external_rasters/MODIS_NDWI_Winter.tif")
plot(ndwi_winter)

evi_summer <- rast("./rescaling_external_rasters/MODIS_EVI_Summer.tif")
plot(evi_summer)

# evi_autumn <- rast("./rescaling_external_rasters/MODIS_EVI_Autumn.tif")
# plot(evi_autumn)

evi_winter <- rast("./rescaling_external_rasters/MODIS_EVI_Winter.tif")
plot(evi_winter)

trees_bd <- rast("./rescaling_external_rasters/ESA_LC_TREES-BD.tif")
plot(trees_bd)

trees_ne <- rast("./rescaling_external_rasters/ESA_LC_TREES-NE.tif")
plot(trees_ne)

# raster by which we are going to rescale first raster
chelsa <- rast("./data/CHELSA_1981-2010_V.2.1/bio15.tif")
plot(chelsa)

# NDVI Summer
# using another raster as CRS, extent, and resolution source
y <- terra::project(ndvi_summer, chelsa)
plot(y)

# Mask raster by another raster
x <- terra::mask(y, chelsa)
plot(x)

# Mask raster by Ukraine boundary polygon
z <- mask(x, Ukraine)
plot(z)

# write raster as tiff file
writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/NDVIs.tif",
            # format = "GTiff",
            overwrite = TRUE)

# # NDVI Autumn
# # using another raster as CRS, extent, and resolution source
# y <- terra::project(ndvi_autumn, chelsa)
# plot(y)

# # Mask raster by another raster
# x <- terra::mask(y, chelsa)
# plot(x)

# # Mask raster by Ukraine boundary polygon
# z <- mask(x, Ukraine)
# plot(z)

# # write raster as tiff file
# writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/NDVIa.tif",
#             # format = "GTiff",
#             overwrite = TRUE)

# NDVI Winter
# using another raster as CRS, extent, and resolution source
y <- terra::project(ndvi_winter, chelsa)
plot(y)

# Mask raster by another raster
x <- terra::mask(y, chelsa)
plot(x)

# Mask raster by Ukraine boundary polygon
z <- mask(x, Ukraine)
plot(z)

# write raster as tiff file
writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/NDVIw.tif",
            # format = "GTiff",
            overwrite = TRUE)


# NDWI Summer
# using another raster as CRS, extent, and resolution source
y <- terra::project(ndwi_summer, chelsa)
plot(y)

# Mask raster by another raster
x <- terra::mask(y, chelsa)
plot(x)

# Mask raster by Ukraine boundary polygon
z <- mask(x, Ukraine)
plot(z)

# write raster as tiff file
writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/NDWIs.tif",
            # format = "GTiff",
            overwrite = TRUE)

# # NDWI Autumn
# # using another raster as CRS, extent, and resolution source
# y <- terra::project(ndwi_autumn, chelsa)
# plot(y)

# # Mask raster by another raster
# x <- terra::mask(y, chelsa)
# plot(x)

# # Mask raster by Ukraine boundary polygon
# z <- mask(x, Ukraine)
# plot(z)

# # write raster as tiff file
# writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/NDWIa.tif",
#             # format = "GTiff",
#             overwrite = TRUE)

# NDWI Winter
# using another raster as CRS, extent, and resolution source
y <- terra::project(ndwi_winter, chelsa)
plot(y)

# Mask raster by another raster
x <- terra::mask(y, chelsa)
plot(x)

# Mask raster by Ukraine boundary polygon
z <- mask(x, Ukraine)
plot(z)

# write raster as tiff file
writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/NDWIw.tif",
            # format = "GTiff",
            overwrite = TRUE)


# EVI Summer
# using another raster as CRS, extent, and resolution source
y <- terra::project(evi_summer, chelsa)
plot(y)

# Mask raster by another raster
x <- terra::mask(y, chelsa)
plot(x)

# Mask raster by Ukraine boundary polygon
z <- mask(x, Ukraine)
plot(z)

# write raster as tiff file
writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/EVIs.tif",
            # format = "GTiff",
            overwrite = TRUE)

# # EVI Autumn
# # using another raster as CRS, extent, and resolution source
# y <- terra::project(evi_autumn, chelsa)
# plot(y)

# # Mask raster by another raster
# x <- terra::mask(y, chelsa)
# plot(x)

# # Mask raster by Ukraine boundary polygon
# z <- mask(x, Ukraine)
# plot(z)

# # write raster as tiff file
# writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/EVIa.tif",
#             # format = "GTiff",
#             overwrite = TRUE)

# EVI Winter
# using another raster as CRS, extent, and resolution source
y <- terra::project(evi_winter, chelsa)
plot(y)

# Mask raster by another raster
x <- terra::mask(y, chelsa)
plot(x)

# Mask raster by Ukraine boundary polygon
z <- mask(x, Ukraine)
plot(z)

# write raster as tiff file
writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/EVIw.tif",
            # format = "GTiff",
            overwrite = TRUE)

# ESA Land Cover
# tree cover broadleaved deciduous
# using another raster as CRS, extent, and resolution source
y <- terra::project(trees_bd, chelsa)
plot(y)

# Mask raster by another raster
x <- terra::mask(y, chelsa)
plot(x)

# Mask raster by Ukraine boundary polygon
z <- mask(x, Ukraine)
plot(z)

# write raster as tiff file
writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/TREES_BD.tif",
            # format = "GTiff",
            overwrite = TRUE)

# tree cover needleleaved evergreen
# using another raster as CRS, extent, and resolution source
y <- terra::project(trees_ne, chelsa)
plot(y)

# Mask raster by another raster
x <- terra::mask(y, chelsa)
plot(x)

# Mask raster by Ukraine boundary polygon
z <- mask(x, Ukraine)
plot(z)

# write raster as tiff file
writeRaster(z, filename = "./data/CHELSA_1981-2010_V.2.1/TREES_NE.tif",
            # format = "GTiff",
            overwrite = TRUE)

# End of the script ####
