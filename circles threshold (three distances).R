library(tidyverse)
library(phenoptr)
library(ggplot2)
library(dplyr)
library(ggforce)

# read the single cell *csd) data
csd <- Measurements_D
csd %>% count(Phenotype)

# filter for the desired columns
csd <- csd %>%
  filter(Phenotype %in% c("B220", "CD3", "KLRG1", "CD11b"))

# measure distance to the nearest cell from each phenotype
distances <- find_nearest_distance(csd)

# join the distances to the left of the origional csd 
csd_with_distance <- bind_cols(csd, distances)

# plot the distribution of the distances ****to a specific phenotype****
ggplot(csd_with_distance, aes(`Distance to CD3`, color=Phenotype)) +
  geom_density(size=1)

# generate a summary of teh average distances ****to the specific phenotype****
summary_distances <- csd_with_distance %>%
    group_by(Phenotype) %>%
    summarize(
        Average_Distance = mean(`Distance to KLRG1`, na.rm = TRUE),
        Std_Distance = sd(`Distance to B220`, na.rm = TRUE)
    )

# Print the summary_distances
print(summary_distances)

# Plot the ggplot with average distances and standard deviations
ggplot(csd_with_distance, aes(`Distance to B220`, color = Phenotype)) +
    geom_density(size = 1) +
    labs(title = 'Density Plot of Distance to B220 by Phenotype') +
    geom_vline(
        data = summary_distances,
        aes(xintercept = Average_Distance, color = Phenotype),
        linetype = 'dashed',
        size = 1
    ) 

# create a table for each phenotype
B220_cells <- csd_with_distance %>%
  filter(Phenotype == 'B220')
CD3_cells <- csd_with_distance %>%
  filter(Phenotype == 'CD3')
KLRG1_cells <- csd_with_distance %>%
  filter(Phenotype == 'KLRG1')

# create a table for the triplet interactions **from the viewpoint of each phenotype**
B220_to_CD3 = B220_cells %>% left_join(CD3_cells, by=c('Cell ID CD3'='Cell ID'),
                                    suffix=c('', '.CD3'))
B220_to_KLRG1_CD3 = B220_to_CD3 %>% left_join(KLRG1_cells, by=c('Cell ID KLRG1'='Cell ID'),
                                    suffix=c('', '.KLRG1'))

KLRG1_to_CD3 = KLRG1_cells %>% left_join(CD3_cells, by=c('Cell ID CD3'='Cell ID'),
                                    suffix=c('', '.CD3'))
KLRG1_to_B220_CD3 = KLRG1_to_CD3 %>% left_join(B220_cells, by=c('Cell ID B220'='Cell ID'),
                                    suffix=c('', '.B220'))

CD3_to_B220 = CD3_cells %>% left_join(B220_cells, by=c('Cell ID B220'='Cell ID'),
                                    suffix=c('', '.B220'))
CD3_to_KLRG1_B220 = CD3_to_B220 %>% left_join(KLRG1_cells, by=c('Cell ID KLRG1'='Cell ID'),
                                    suffix=c('', '.KLRG1'))

# find the triplet interactions (the same exact three cells) in common between the three view points
common_in_B220 <- inner_join(CD3_to_KLRG1_B220, KLRG1_to_B220_CD3, by = c("Name.B220"))
common_in_all <- inner_join(common_in_B220, B220_to_KLRG1_CD3, by = c("Name.CD3"), "Name.KLRG1")

# leave only one copy of the common triplet interactions
CD3_to_KLRG1_B220 <- anti_join(CD3_to_KLRG1_B220, common_in_all, by = c("Name" = "Name.x"))
KLRG1_to_B220_CD3 <- anti_join(KLRG1_to_B220_CD3, common_in_all, by = c("Name" = "Name.y"))

# Calculate the Euclidean distance between teh other two cells (ex. if we already have the distance in the pointview of B220 to CD3 nad KLRG1, 
# then we need the distance between CD3 and KLRG1) and add it as a new column
# Calculate the Euclidean distance and add it as a new column
B220_to_KLRG1_CD3 <- B220_to_KLRG1_CD3 %>%
  mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.KLRG1", "Cell Y Position.KLRG1"), as.numeric) %>%
  mutate(
    `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.KLRG1`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.KLRG1`)^2)
  )

KLRG1_to_B220_CD3 <- KLRG1_to_B220_CD3 %>%
  mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.B220", "Cell Y Position.B220"), as.numeric) %>%
  mutate(
    `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.B220`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.B220`)^2)
  )

CD3_to_KLRG1_B220 <- CD3_to_KLRG1_B220 %>%
  mutate_at(vars("Cell X Position.B220", "Cell Y Position.B220", "Cell X Position.KLRG1", "Cell Y Position.KLRG1"), as.numeric) %>%
  mutate(
    `3rd distance` = sqrt((`Cell X Position.KLRG1` - `Cell X Position.B220`)^2 + (`Cell Y Position.KLRG1` - `Cell Y Position.B220`)^2)
  )

# get the triplet distance
B220_to_KLRG1_CD3 <- B220_to_KLRG1_CD3 %>%
  mutate(sum = `Distance to KLRG1` + `Distance to CD3` + `3rd distance`)
KLRG1_to_B220_CD3 <- KLRG1_to_B220_CD3 %>%
  mutate(sum = `Distance to CD3` + `Distance to B220` + `3rd distance`)
CD3_to_KLRG1_B220 <- CD3_to_KLRG1_B220 %>%
  mutate(sum = `Distance to KLRG1` + `Distance to B220` + `3rd distance`)

# introduce a threshold value
threshold_distance <- 10 

# get the triplet distance
B220_to_KLRG1_CD3 <- B220_to_KLRG1_CD3 %>%
  mutate(sum = `Distance to KLRG1` + `Distance to CD3` + `3rd distance`)
KLRG1_to_B220_CD3 <- KLRG1_to_B220_CD3 %>%
  mutate(sum = `Distance to CD3` + `Distance to B220` + `3rd distance`)
CD3_to_KLRG1_B220 <- CD3_to_KLRG1_B220 %>%
  mutate(sum = `Distance to KLRG1` + `Distance to B220` + `3rd distance`)

# filter the triplet interaction tables based on the threshold (the triplet distance is under a threshold)
B220_to_KLRG1_CD3 <- B220_to_KLRG1_CD3 %>%
  filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to KLRG1` < threshold_distance)

KLRG1_to_B220_CD3 <- KLRG1_to_B220_CD3 %>%
  filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to B220` < threshold_distance)

CD3_to_KLRG1_B220 <- CD3_to_KLRG1_B220 %>%
  filter(`3rd distance` < threshold_distance & `Distance to B220` < threshold_distance & `Distance to KLRG1` < threshold_distance)

# add all tables together to get the total unique triplet interactions
column_names <- names(B220_to_KLRG1_CD3)
KLRG1_to_B220_CD3 <- KLRG1_to_B220_CD3 %>%
  setNames(column_names)
CD3_to_KLRG1_B220 <- CD3_to_KLRG1_B220 %>%
  setNames(column_names)
triplet <- rbind(B220_to_KLRG1_CD3, KLRG1_to_B220_CD3, CD3_to_KLRG1_B220)

cat("Number of triplet interactions with an average distance ", threshold_distance, " microns: ", nrow(triplet), "\n")


background_path <- "C:/Users/Saven Denha/Desktop/Kyle/February 2, 2024/Diseased/diseased_processed.jpg"
background <- jpeg::readJPEG(background_path) %>% as.raster()


base_plot <- ggplot(mapping=aes(`Cell X Position`, `Cell Y Position`)) %>% 
  phenoptr:::add_scales_and_background(background, ncol(background)*0.2071607, nrow(background)*0.2071607, scale_color='white') +
  labs(x='Cell X Position', y='Cell Y Position') +
  scale_color_manual('Phenotype', values=c('CD3'='green', 'B220'='red', 'KLRG1'='cyan2')) +
  theme_void() +
  theme(legend.position = "none")

final_plot = base_plot + 
    geom_circle(data = triplet,
                aes(x0 = `Cell X Position`, y0 = `Cell Y Position`, r = threshold_distance-1),
                color = 'white', fill = NA, size = 0.1)

ggsave("triplet_interaction.png", final_plot, width = 6, height = 8, units = "in", dpi = 2000)

radius <- threshold_distance + 7

# Number of points on the circumference
num_points <- 100

for (i in 1:nrow(triplet)) {
    triplet_row <- triplet[i, ]
    
    # Set the center and radius of the circle
    center_x <- triplet_row$`Cell X Position`
    center_y <- triplet_row$`Cell Y Position`
    
    # Calculate points on the circumference
    theta <- seq(0, 2 * pi, length.out = num_points)
    x_points <- (center_x + radius * cos(theta))/0.2071607
    y_points <- (center_y + radius * sin(theta))/0.2071607
    
    # Create a data frame for the circle
    circle_coordinates <- data.frame(X = x_points, Y = y_points)
    
    # Save the data frame to a text file
    write.table(circle_coordinates, paste0("circle_", i, "_coordinates.txt"), row.names = FALSE, col.names = TRUE, sep = "\t")
}