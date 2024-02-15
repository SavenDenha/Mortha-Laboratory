library(tidyverse)
library(phenoptr)
library(ggplot2)
library(dplyr)
library(ggforce)

# Read the single-cell (csd) data for both conditions
csd_control <- Measurements_C
csd_diseased <- Measurements_D

# Count occurrences of each phenotype in the diseased condition
csd_diseased %>% count(Phenotype)

# Filter for the desired columns in both conditions
csd_control <- csd_control %>%
    filter(Phenotype %in% c("B220", "CD3", "KLRG1"))

csd_diseased <- csd_diseased %>%
    filter(Phenotype %in% c("B220", "CD3", "KLRG1"))

# Measure distance to the nearest cell from each phenotype in both conditions
distances_control <- find_nearest_distance(csd_control)
distances_diseased <- find_nearest_distance(csd_diseased)

# Join the distances to the left of the original csd for both conditions
csd_with_distance_control <- bind_cols(csd_control, distances_control)
csd_with_distance_diseased <- bind_cols(csd_diseased, distances_diseased)

# Create a new variable indicating the condition
csd_with_distance_control$Condition <- "Control"
csd_with_distance_diseased$Condition <- "Diseased"

# Combine the data for both conditions
combined_data <- rbind(csd_with_distance_control, csd_with_distance_diseased)

# Filter for CD3 and B220 phenotypes
filtered_data <- combined_data %>%
    filter(Phenotype %in% c("KLRG1"))

# Generate a summary of the average distances to the specific phenotype for both conditions
summary_distances_combined <- filtered_data %>%
    group_by(Condition) %>%
    summarize(
        Average_Distance = mean(`Distance to CD3`, na.rm = TRUE),
        Std_Distance = sd(`Distance to CD3`, na.rm = TRUE)
    )
print(summary_distances_combined)

# Generate a summary of the average distances to the specific phenotype for both conditions
summary_distances_combined <- filtered_data %>%
    group_by(Condition) %>%
    summarize(
        Average_Distance = mean(`Distance to B220`, na.rm = TRUE),
        Std_Distance = sd(`Distance to B220`, na.rm = TRUE)
    )
print(summary_distances_combined)

# Create a single density plot for Distance to CD3 and Distance to B220
combined_density_plot <- ggplot(filtered_data, aes(x = `Distance to CD3`)) +
    geom_density(aes(color = Condition, linetype = "CD3"), alpha = 0.5) +
    geom_density(aes(x = `Distance to B220`, color = Condition, linetype = "B220"), alpha = 0.5) +
    scale_color_manual(values = c("Control" = "black", "Diseased" = "red")) +
    scale_linetype_manual(values = c("CD3" = "solid", "B220" = "dashed")) +
    labs(
        x = "Distance",
        color = "Condition",
        linetype = "Phenotype"
    ) +
    xlab("Distance from KLRG1 to:") +  # Add x-axis label here
      theme(
      panel.background = element_blank(),  # Remove background
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black"),  # Set axis lines to black
      legend.position = "bottom"  # Move legend to the bottom
    )

# Display the combined density plot
print(combined_density_plot)
