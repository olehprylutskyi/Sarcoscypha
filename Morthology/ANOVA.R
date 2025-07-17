# Load required libraries
library(ggplot2)
library(dplyr)
library(ggExtra)
library(scales)  # for color palettes

# 1. Load and clean the dataset
# Replace with your actual file name if needed
dat <- read_csv("measurements_morphology.csv")

# Keep only relevant columns and rename for clarity
dat_clean <- dat %>%
  select(length, width, species, Specimen_ID) %>%
  filter(!is.na(length), !is.na(width), !is.na(species), !is.na(Specimen_ID)) %>%
  rename(
    spore_length = length,
    spore_width = width,
    specimen = Specimen_ID
  )

# 2 faktor ANOVA
mod <- aov(spore_length ~ species * specimen, data = dat_clean)
summary(mod)

# can be also done like car::Anova for type II or III:
# car::Anova(mod, type = 2)

# 3. Create custom color palette for each specimen, grouped by species

# Get unique specimen IDs for each species
austriaca_ids <- dat_clean %>% filter(species == "S. austriaca") %>% pull(specimen) %>% unique()
coccinea_ids  <- dat_clean %>% filter(species == "S. coccinea") %>% pull(specimen) %>% unique()

# Generate distinguishable shades of red for S. austriaca
austriaca_cols <- setNames(
  sequential_hcl(length(austriaca_ids), h = 0, c = 80, l = seq(35, 85, length.out = length(austriaca_ids))),
  austriaca_ids
)

# Generate distinguishable shades of ochre/yellow for S. coccinea
coccinea_cols <- setNames(
  sequential_hcl(length(coccinea_ids), h = 50, c = 80, l = seq(40, 85, length.out = length(coccinea_ids))),
  coccinea_ids
)

# Create color palettes for each specimen within each species

# Get unique specimen IDs for each species
austriaca_ids <- dat_clean %>% filter(species == "S. austriaca") %>% pull(specimen) %>% unique()
coccinea_ids <- dat_clean %>% filter(species == "S. coccinea") %>% pull(specimen) %>% unique()

# Create gradients of red shades for S. austriaca
austriaca_cols <- setNames(
  colorRampPalette(c("firebrick3", "tomato", "indianred"))(length(austriaca_ids)),
  austriaca_ids
)

# Create gradients of ochre shades for S. coccinea
coccinea_cols <- setNames(
  colorRampPalette(c("goldenrod3", "khaki4", "darkgoldenrod1"))(length(coccinea_ids)),
  coccinea_ids
)

# Combine into one color palette by specimen
specimen_colors <- c(austriaca_cols, coccinea_cols)

# Base plot
p <- ggplot(dat_clean, aes(x = spore_length, y = spore_width, shape = species, color = specimen)) +
  geom_point(alpha = 0.9, size = 2.2) +
  scale_color_manual(values = specimen_colors) +
  scale_shape_manual(
    values = c("S. austriaca" = 17, "S. coccinea" = 16),
    labels = c("S. austriaca" = expression(italic("S. austriaca")),
               "S. coccinea" = expression(italic("S. coccinea")))
  ) +
  labs(
    x = "Spore Length",
    y = "Spore Width",
    color = "Specimen",
    shape = "Species"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# Add marginal density plots
ggMarginal(p, type = "density", groupColour = FALSE, groupFill = TRUE)

# 3.1. Combine into one color palette by specimen
specimen_colors <- c(austriaca_cols, coccinea_cols)

# Base plot with facet wrap
ggplot(dat_clean, aes(x = spore_length, y = spore_width, shape = species, color = specimen)) +
  geom_point(alpha = 0.9, size = 2.2) +
  scale_color_manual(values = specimen_colors) +
  scale_shape_manual(
    values = c("S. austriaca" = 17, "S. coccinea" = 16),
    labels = c("S. austriaca" = expression(italic("S. austriaca")),
               "S. coccinea" = expression(italic("S. coccinea")))
  ) +
  facet_wrap(~ specimen) +
  labs(
    x = "Spore Length",
    y = "Spore Width",
    color = "Specimen",
    shape = "Species"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
