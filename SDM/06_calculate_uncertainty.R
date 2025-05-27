# Prepare the environment ####
rm(list = ls()) # Reset R`s brain

# Optional: set working directory (change PATH to yours)
# No need if RStudio's Rproj is used
# setwd("~/R/sarcoscypha-urnula_legacy/Maxent_Sarco_Ukraine_2023")

datadir <- "./data"
outputdir <- "./outputs"
modeldir <- "./models"
evaldir <- paste0(modeldir, "/evaluations")

# Loading libraries
library(dismo)     # for MaxEnt routine
library(sp)        # Spatial Data
library(rJava)     # support MaxEnt
library(gbm)       # Calculation of AUC
library(raster)    # raster manipulations
library(sf)        # spatial data for ggplot2
library(patchwork) # plates
library(dplyr)     # data manipulations
library(ggplot2)   # visualizations
library(tidyr)     # data manipulations
library(doBy)      # statsummary for ggplot2

# Load state boundary of Ukraine
load(file = paste0(datadir, "/Ukraine.Rdata"))
# Convert to {sf}
Ukraine <- sf::st_as_sf(Ukraine)

# Load pre-processed covariates
load(file = paste0(outputdir, "/maxent_ready/covariates.Rdata"))
names(covariates)

# [1] "bio02"      "bio15"      "bio19"      "EVIs"       "EVIw"       "fcf"       
# [7] "gsl"        "hurs_mean"  "hurs_range" "NDWIw"      "scd"        "TREES_BD"  
# [13] "TREES_NE" 


# Get models
# Get a list of all .Rdata files in the folder
model_files <- list.files(
  path = modeldir, 
  pattern = "Rdata"
)
modelPath <- file.path(modeldir, model_files)

# Create an empty list to store the loaded data
models <- list()

# Load models as a list
for (i in 1:length(model_files)) {
  load(file = modelPath[i])
  models[i] <- m_maxent
}

# Get models` evaluations
# Get a list of all .Rdata files in the folder
eval_files <- list.files(
  path = evaldir, 
  pattern = "Rdata"
)
evalPath <- file.path(evaldir, eval_files)

# Create an empty list to store loaded data in
evals <- list()
# Store model evaluations as a list
for (i in 1:length(eval_files)) {
  load(file = evalPath[i])
  evals[i] <- e_maxent
}

# Area under receiver operator curve (AUC)
# Extract AUCs
AUCs <- list()
for (i in 1:length(evals)) {
  auc <- round(evals[[i]]@auc, 3)
  AUCs[i] <- auc
}

auc_df <- as.data.frame(AUCs)
auc_df <- t(auc_df) # vector of AUC for 100 models
# auc_df <- as.data.frame(auc_df)

# Calculate summary statistics for AUCs
auc_mean <- round(mean(auc_df), 2)
auc_max <- round(max(auc_df), 2)
auc_min <- round(min(auc_df), 2)
auc_sd <- round(sd(auc_df), 2)


# Variable contributions ####
# Preview results output to define necessry variables' indices
models[[1]]@results

vc_list <- list()
for (i in 1:length(models)) {
  vc <- models[[i]]@results[7:(6+length(names(covariates))),]
  vc_list[[i]] <- vc
}

vc_df <- as.data.frame(vc_list)
vc_df <- t(vc_df)
vc_df <- as.data.frame(vc_df)

str(vc_df)

# Rename variables to get covariates' names back
names(vc_df) <- sub("\\.contribution$", "", names(vc_df))

# Transform wide dataframe into long
long_vc <- pivot_longer(
  vc_df, 
  cols = 1:ncol(vc_df), # columns to pivot
  names_to = "variable", # new column name for variable names
  values_to = "var_contr" # new column name for variable values
)

# # Boxplot
# ggplot(long_vc, aes(x = variable, y = var_contr)) +
#   geom_boxplot() +
#   coord_flip()

# Min values
MIN <- summaryBy(var_contr ~ variable, FUN = min, data = long_vc)
# Max values
MAX <- summaryBy(var_contr ~ variable, FUN = max, data = long_vc)
# Mean values
MEAN <- summaryBy(var_contr ~ variable, FUN = mean, data = long_vc)
# Standard errors
SE <- summaryBy(
  var_contr ~ variable, 
  FUN = function(x) sd(x)/sqrt(length(x)),
  data = long_vc
)
# Pivot table with statistical summaries
dat <- data.frame(MIN[1], MIN[2], MAX[2], MEAN[2], SE[2])
names(dat) <- c("variable", "Min", "Max", "Mean", "SE")
dat$variable

dat$variable <- c(
  "Mean diurnal air temperature range",
  "Precipitation seasonality",
  "Mean monthly precipitation amount of the coldest quarter",
  "Enhanced Vegetation Index (summer)",
  "Enhanced Vegetation Index (winter)",
  "Frost change frequency",
  "Growing season length",
  "Mean monthly near-surface relative humidity",
  "Annual range of monthly near-surface relative humidity",
  "Normalized Difference Water Index (winter)",
  "Snow cover days",
  "Broadleaved deciduous tree cover",
  "Needleleaved evergreen tree cover"
)

# Export variable contributions to csv
write.csv(
  dat,
  file = paste0(outputdir, "/variable_importance.csv")
)

# Make variable importance plot
varimp_p <- ggplot(
  data = dat,
  aes(x = reorder(variable, Min), ymin = Min, ymax = Max)
) +
  geom_pointrange(aes(y = Mean), colour = "#CC101F") +
  coord_flip() +
  labs(title = "",
       x = "",
       y = "Variable contribution, %",
       caption = "100 model replication\n* Mean (Min – Max)") +
  annotate("text",  x = -Inf, y = 50, vjust = -5,
           label = paste0("AUC = ", auc_mean, " (", auc_min, " – ", auc_max, ")*"),
           size = 4) +
  theme_bw()

varimp_p

# Save variable importance plot
png(paste0(outputdir, "/variable_importance.png"), width = 160, height = 120,
    units = "mm", res = 300)
varimp_p
dev.off()

# Prediction ####
predictions <- list()

for (i in 1:length(models)) {
  maxentpred <- dismo::predict(covariates, models[[i]])
  predictions[[i]] <- maxentpred
}

predictions
# Export prediction rasters
save(predictions, file = paste0(outputdir, "/predictions.Rdata"))

# Merge predictions into one rasterBrick 
r_brick <- brick(predictions)

# Delete unused objects and free-up memory
rm(predictions)
gc()

# Calculate median
r_median <- calc(r_brick, median)
plot(r_median)

# Calculate Standard Deviation
# as prediction uncertainty metric
r_sd <- calc(r_brick, sd)
plot(r_sd)


# Visualisation ####
# Median prediction
median_pred.df =  as.data.frame(r_median, xy = TRUE) %>% drop_na()

# Plotting
p_median_preds <- ggplot() +
  geom_raster(
    aes(x = x, y = y, fill = layer),
    data = median_pred.df
  ) +
  scale_fill_viridis_c("Probability\nof presence") +
  labs(
    x = "", y = "",
    title = "Pedicted distribution of Sarcoscypha spp.",
    caption = "Maxent, 100 replications"
  ) +
  theme_minimal()

p_median_preds

# Concatenate the vector elements into a comma-separated string
my_covariates_4caption <- paste(names(covariates), collapse = "-")

png(
  paste0(outputdir, "/publication_ready/Maxent_2km_chelsa_prediction_", my_covariates_4caption, ".png"),
  width = 160, height = 120, units = "mm", res = 300
)
p_median_preds
dev.off()

# Uncertainty
sd_pred.df =  as.data.frame(r_sd, xy = TRUE) %>% drop_na()

# Plotting
p_sd_preds <- ggplot() +
  geom_raster(
    aes(x = x, y = y, fill = layer),
    data = sd_pred.df
  ) +
  # scale_fill_viridis_c("Standard\ndeviation") +
  scale_fill_viridis_c("Coefficient\nof variation") +
  labs(x = "", y = "",
       title = "Prediction uncertainty for Sarcoscypha spp.",
       caption = "Maxent, 100 replications") +
  theme_minimal()

p_sd_preds

# # Concatenate the vector elements into a comma-separated string
# my_covariates_4caption <- paste(names(covariates), collapse = "-")

png(
  paste0(outputdir, "/publication_ready/Maxent_2km_chelsa_uncertainty_", my_covariates_4caption, ".png"),
  width = 160, height = 120, units = "mm", res = 300
)
p_sd_preds
dev.off()


# Save rasters for future ####
# write raster as tiff file
writeRaster(
  r_median, 
  filename = paste0(outputdir, "/output_rasters/pred_median_100rep.tiff"),
  # format = "GTiff",
  overwrite = TRUE
)

writeRaster(
  r_sd,
  filename = paste0(outputdir, "/output_rasters/pred_sd_100rep.tiff"),
  # format = "GTiff",
  overwrite = TRUE
)


# SDM Plate ####
# # Read saved rasters
# r_median <- raster(paste0(outputdir, "/pred_median_100rep.tiff"))
# r_sd <- raster(paste0(outputdir, "/pred_sd_100rep.tiff"))

# Make ggplot-ready dataframes
median_pred.df <- as.data.frame(r_median, xy = TRUE) %>% drop_na()
sd_pred.df <- as.data.frame(r_sd, xy = TRUE) %>% drop_na()

# Plot ensemble prediction
p_median_preds <- ggplot() +
  geom_sf(data = st_as_sf(Ukraine)) +
  geom_raster(
    aes(x = x, y = y, fill = layer),
    data = median_pred.df
  ) +
  scale_fill_viridis_c(
    "Probability\nof presence\n ",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1)
  ) +
  labs(
    x = "",
    y = "",
    title = "(b)"
  ) +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank()
  )

# Show ensemble prediction
p_median_preds

# Plot prediction uncertainty
p_sd_preds <- ggplot() +
  geom_sf(data = st_as_sf(Ukraine)) +
  geom_raster(
    aes(x = x, y = y, fill = layer),
    data = sd_pred.df
  ) +
  scale_fill_viridis_c(
    "Standard\ndeviation\n ",
    breaks = c(0.1, 0.2),
    option = "B"
  ) +
  labs(
    x = "", y = "", title = "(c)"
  ) +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank()
  )

# Show prediction uncertainty
p_sd_preds

# Make variable importance plot
varimp_p <- ggplot(
  data = dat, 
  aes(x = reorder(variable, Min), ymin = Min, ymax = Max)
) +
  geom_pointrange(
    aes(y = Mean),
    colour = "#CC101F"
  ) +
  coord_flip() +
  scale_y_continuous(
    breaks = c(0, 50, 100)
  ) +
  labs(
    title = "(a)",
    x = "",
    y = "Variable contribution, %",
    caption = "* Mean (Min – Max)"
  ) +
  annotate(
    "text",  x = -Inf, y = 35, vjust = -5,
    label = paste0("AUC = ", auc_mean, " (", auc_min, "–", auc_max, ")*"),
    size = 4
  ) +
  theme(
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   panel.background = element_blank(),
   panel.border = element_rect(fill = NA),
   axis.title.y = element_text(size = 10),
   axis.title.x = element_text(size = 9)
  )


varimp_p

# Export the final plate
png(
  paste0(outputdir, "/publication_ready/Maxent_chelsa_plate.png"),
  width = 230, height = 140, units = "mm", res = 300
)
varimp_p + p_median_preds/p_sd_preds
dev.off()

# End of the script ####
