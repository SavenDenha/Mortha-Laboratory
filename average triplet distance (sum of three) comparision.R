library(tidyverse)
library(phenoptr)
library(ggplot2)
library(dplyr)
library(ggforce)
library(raster)
library(cowplot)
library(magick)

# read the single-cell *csd* data for the control group
csd_control <- Measurements_C

# read the single-cell *csd* data for the diseased group
csd_diseased <- Measurements_D

# filter for the desired columns in both control and diseased groups
csd_control <- csd_control %>%
    filter(Phenotype %in% c("B220", "CD3", "CD11b"))

csd_diseased <- csd_diseased %>%
    filter(Phenotype %in% c("B220", "CD3", "CD11b"))

# measure distance to the nearest cell from each phenotype for both groups
distances_control <- find_nearest_distance(csd_control)
distances_diseased <- find_nearest_distance(csd_diseased)

# join the distances to the left of the original csd for both groups
csd_with_distance_control <- bind_cols(csd_control, distances_control)
csd_with_distance_diseased <- bind_cols(csd_diseased, distances_diseased)

# create a table for each phenotype for both groups
B220_cells_control <- csd_with_distance_control %>%
    filter(Phenotype == 'B220')
CD3_cells_control <- csd_with_distance_control %>%
    filter(Phenotype == 'CD3')
CD11b_cells_control <- csd_with_distance_control %>%
    filter(Phenotype == 'CD11b')

B220_cells_diseased <- csd_with_distance_diseased %>%
    filter(Phenotype == 'B220')
CD3_cells_diseased <- csd_with_distance_diseased %>%
    filter(Phenotype == 'CD3')
CD11b_cells_diseased <- csd_with_distance_diseased %>%
    filter(Phenotype == 'CD11b')

# create a table for the triplet interactions **from the viewpoint of each phenotype** for both groups
B220_to_CD3_control = B220_cells_control %>% left_join(CD3_cells_control, by=c('Cell ID CD3'='Cell ID'),
                                                       suffix=c('', '.CD3'))
B220_to_CD11b_CD3_control = B220_to_CD3_control %>% left_join(CD11b_cells_control, by=c('Cell ID CD11b'='Cell ID'),
                                                              suffix=c('', '.CD11b'))

CD11b_to_CD3_control = CD11b_cells_control %>% left_join(CD3_cells_control, by=c('Cell ID CD3'='Cell ID'),
                                                         suffix=c('', '.CD3'))
CD11b_to_B220_CD3_control = CD11b_to_CD3_control %>% left_join(B220_cells_control, by=c('Cell ID B220'='Cell ID'),
                                                               suffix=c('', '.B220'))

CD3_to_B220_control = CD3_cells_control %>% left_join(B220_cells_control, by=c('Cell ID B220'='Cell ID'),
                                                      suffix=c('', '.B220'))
CD3_to_CD11b_B220_control = CD3_to_B220_control %>% left_join(CD11b_cells_control, by=c('Cell ID CD11b'='Cell ID'),
                                                              suffix=c('', '.CD11b'))

# create a table for the triplet interactions **from the viewpoint of each phenotype** for both groups
B220_to_CD3_diseased = B220_cells_diseased %>% left_join(CD3_cells_diseased, by=c('Cell ID CD3'='Cell ID'),
                                                         suffix=c('', '.CD3'))
B220_to_CD11b_CD3_diseased = B220_to_CD3_diseased %>% left_join(CD11b_cells_diseased, by=c('Cell ID CD11b'='Cell ID'),
                                                                suffix=c('', '.CD11b'))

CD11b_to_CD3_diseased = CD11b_cells_diseased %>% left_join(CD3_cells_diseased, by=c('Cell ID CD3'='Cell ID'),
                                                           suffix=c('', '.CD3'))
CD11b_to_B220_CD3_diseased = CD11b_to_CD3_diseased %>% left_join(B220_cells_diseased, by=c('Cell ID B220'='Cell ID'),
                                                                 suffix=c('', '.B220'))

CD3_to_B220_diseased = CD3_cells_diseased %>% left_join(B220_cells_diseased, by=c('Cell ID B220'='Cell ID'),
                                                        suffix=c('', '.B220'))
CD3_to_CD11b_B220_diseased = CD3_to_B220_diseased %>% left_join(CD11b_cells_diseased, by=c('Cell ID CD11b'='Cell ID'),
                                                                suffix=c('', '.CD11b'))

# find the triplet interactions (the same exact three cells) in common between the three viewpoints for the control group
common_in_B220_control <- inner_join(CD3_to_CD11b_B220_control, CD11b_to_B220_CD3_control, by = c("Name.B220"))
common_in_all_control <- inner_join(common_in_B220_control, B220_to_CD11b_CD3_control, by = c("Name.CD3"), "Name.CD11b")

# find the triplet interactions (the same exact three cells) in common between the three viewpoints for the diseased group
common_in_B220_diseased <- inner_join(CD3_to_CD11b_B220_diseased, CD11b_to_B220_CD3_diseased, by = c("Name.B220"))
common_in_all_diseased <- inner_join(common_in_B220_diseased, B220_to_CD11b_CD3_diseased, by = c("Name.CD3"), "Name.CD11b")

# leave only one copy of the common triplet interactions for the control group
CD3_to_CD11b_B220_control <- anti_join(CD3_to_CD11b_B220_control, common_in_all_control, by = c("Name" = "Name.x"))
CD11b_to_B220_CD3_control <- anti_join(CD11b_to_B220_CD3_control, common_in_all_control, by = c("Name" = "Name.y"))

# leave only one copy of the common triplet interactions for the diseased group
CD3_to_CD11b_B220_diseased <- anti_join(CD3_to_CD11b_B220_diseased, common_in_all_diseased, by = c("Name" = "Name.x"))
CD11b_to_B220_CD3_diseased <- anti_join(CD11b_to_B220_CD3_diseased, common_in_all_diseased, by = c("Name" = "Name.y"))

# Calculate the Euclidean distance between the other two cells for the control group
B220_to_CD11b_CD3_control <- B220_to_CD11b_CD3_control %>%
    mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.CD11b", "Cell Y Position.CD11b"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.CD11b`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.CD11b`)^2)
    )

CD11b_to_B220_CD3_control <- CD11b_to_B220_CD3_control %>%
    mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.B220", "Cell Y Position.B220"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.B220`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.B220`)^2)
    )

CD3_to_CD11b_B220_control <- CD3_to_CD11b_B220_control %>%
    mutate_at(vars("Cell X Position.B220", "Cell Y Position.B220", "Cell X Position.CD11b", "Cell Y Position.CD11b"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD11b` - `Cell X Position.B220`)^2 + (`Cell Y Position.CD11b` - `Cell Y Position.B220`)^2)
    )

# Calculate the Euclidean distance between the other two cells for the diseased group
B220_to_CD11b_CD3_diseased <- B220_to_CD11b_CD3_diseased %>%
    mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.CD11b", "Cell Y Position.CD11b"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.CD11b`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.CD11b`)^2)
    )

CD11b_to_B220_CD3_diseased <- CD11b_to_B220_CD3_diseased %>%
    mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.B220", "Cell Y Position.B220"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.B220`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.B220`)^2)
    )

CD3_to_CD11b_B220_diseased <- CD3_to_CD11b_B220_diseased %>%
    mutate_at(vars("Cell X Position.B220", "Cell Y Position.B220", "Cell X Position.CD11b", "Cell Y Position.CD11b"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD11b` - `Cell X Position.B220`)^2 + (`Cell Y Position.CD11b` - `Cell Y Position.B220`)^2)
    )

# get the triplet distance for the control group
B220_to_CD11b_CD3_control <- B220_to_CD11b_CD3_control %>%
    mutate(sum = `Distance to CD11b` + `Distance to CD3` + `3rd distance`)
CD11b_to_B220_CD3_control <- CD11b_to_B220_CD3_control %>%
    mutate(sum = `Distance to CD3` + `Distance to B220` + `3rd distance`)
CD3_to_CD11b_B220_control <- CD3_to_CD11b_B220_control %>%
    mutate(sum = `Distance to CD11b` + `Distance to B220` + `3rd distance`)

# get the triplet distance for the diseased group
B220_to_CD11b_CD3_diseased <- B220_to_CD11b_CD3_diseased %>%
    mutate(sum = `Distance to CD11b` + `Distance to CD3` + `3rd distance`)
CD11b_to_B220_CD3_diseased <- CD11b_to_B220_CD3_diseased %>%
    mutate(sum = `Distance to CD3` + `Distance to B220` + `3rd distance`)
CD3_to_CD11b_B220_diseased <- CD3_to_CD11b_B220_diseased %>%
    mutate(sum = `Distance to CD11b` + `Distance to B220` + `3rd distance`)

# add all tables together to get the total unique triplet interactions
column_names <- names(B220_to_CD11b_CD3_control)
CD11b_to_B220_CD3_control <- CD11b_to_B220_CD3_control %>%
  setNames(column_names)
CD3_to_CD11b_B220_control <- CD3_to_CD11b_B220_control %>%
  setNames(column_names)
triplet_control_CD11b <- rbind(B220_to_CD11b_CD3_control, CD11b_to_B220_CD3_control, CD3_to_CD11b_B220_control)

# add all tables together to get the total unique triplet interactions
column_names <- names(B220_to_CD11b_CD3_diseased)
CD11b_to_B220_CD3_diseased <- CD11b_to_B220_CD3_diseased %>%
  setNames(column_names)
CD3_to_CD11b_B220_diseased <- CD3_to_CD11b_B220_diseased %>%
  setNames(column_names)
triplet_diseased_CD11b <- rbind(B220_to_CD11b_CD3_diseased, CD11b_to_B220_CD3_diseased, CD3_to_CD11b_B220_diseased)

# Calculate the average and standard deviation of the "sum" column
sum_summary_control_CD11b <- summary(triplet_control_CD11b$sum)
average_sum_control_CD11b <- mean(triplet_control_CD11b$sum, na.rm = TRUE)
std_dev_sum_control_CD11b <- sd(triplet_control_CD11b$sum, na.rm = TRUE)

# Calculate the average and standard deviation of the "sum" column
sum_summary_diseased_CD11b <- summary(triplet_diseased_CD11b$sum)
average_sum_diseased_CD11b <- mean(triplet_diseased_CD11b$sum, na.rm = TRUE)
std_dev_sum_diseased_CD11b <- sd(triplet_diseased_CD11b$sum, na.rm = TRUE)

############################################################################################################


# read the single-cell *csd* data for the control group
csd_control <- Measurements_C

# read the single-cell *csd* data for the diseased group
csd_diseased <- Measurements_D

# filter for the desired columns in both control and diseased groups
csd_control <- csd_control %>%
    filter(Phenotype %in% c("B220", "CD3", "KLRG1"))

csd_diseased <- csd_diseased %>%
    filter(Phenotype %in% c("B220", "CD3", "KLRG1"))

# measure distance to the nearest cell from each phenotype for both groups
distances_control <- find_nearest_distance(csd_control)
distances_diseased <- find_nearest_distance(csd_diseased)

# join the distances to the left of the original csd for both groups
csd_with_distance_control <- bind_cols(csd_control, distances_control)
csd_with_distance_diseased <- bind_cols(csd_diseased, distances_diseased)

# create a table for each phenotype for both groups
B220_cells_control <- csd_with_distance_control %>%
    filter(Phenotype == 'B220')
CD3_cells_control <- csd_with_distance_control %>%
    filter(Phenotype == 'CD3')
KLRG1_cells_control <- csd_with_distance_control %>%
    filter(Phenotype == 'KLRG1')

B220_cells_diseased <- csd_with_distance_diseased %>%
    filter(Phenotype == 'B220')
CD3_cells_diseased <- csd_with_distance_diseased %>%
    filter(Phenotype == 'CD3')
KLRG1_cells_diseased <- csd_with_distance_diseased %>%
    filter(Phenotype == 'KLRG1')

# create a table for the triplet interactions **from the viewpoint of each phenotype** for both groups
B220_to_CD3_control = B220_cells_control %>% left_join(CD3_cells_control, by=c('Cell ID CD3'='Cell ID'),
                                                       suffix=c('', '.CD3'))
B220_to_KLRG1_CD3_control = B220_to_CD3_control %>% left_join(KLRG1_cells_control, by=c('Cell ID KLRG1'='Cell ID'),
                                                              suffix=c('', '.KLRG1'))

KLRG1_to_CD3_control = KLRG1_cells_control %>% left_join(CD3_cells_control, by=c('Cell ID CD3'='Cell ID'),
                                                         suffix=c('', '.CD3'))
KLRG1_to_B220_CD3_control = KLRG1_to_CD3_control %>% left_join(B220_cells_control, by=c('Cell ID B220'='Cell ID'),
                                                               suffix=c('', '.B220'))

CD3_to_B220_control = CD3_cells_control %>% left_join(B220_cells_control, by=c('Cell ID B220'='Cell ID'),
                                                      suffix=c('', '.B220'))
CD3_to_KLRG1_B220_control = CD3_to_B220_control %>% left_join(KLRG1_cells_control, by=c('Cell ID KLRG1'='Cell ID'),
                                                              suffix=c('', '.KLRG1'))

# create a table for the triplet interactions **from the viewpoint of each phenotype** for both groups
B220_to_CD3_diseased = B220_cells_diseased %>% left_join(CD3_cells_diseased, by=c('Cell ID CD3'='Cell ID'),
                                                         suffix=c('', '.CD3'))
B220_to_KLRG1_CD3_diseased = B220_to_CD3_diseased %>% left_join(KLRG1_cells_diseased, by=c('Cell ID KLRG1'='Cell ID'),
                                                                suffix=c('', '.KLRG1'))

KLRG1_to_CD3_diseased = KLRG1_cells_diseased %>% left_join(CD3_cells_diseased, by=c('Cell ID CD3'='Cell ID'),
                                                           suffix=c('', '.CD3'))
KLRG1_to_B220_CD3_diseased = KLRG1_to_CD3_diseased %>% left_join(B220_cells_diseased, by=c('Cell ID B220'='Cell ID'),
                                                                 suffix=c('', '.B220'))

CD3_to_B220_diseased = CD3_cells_diseased %>% left_join(B220_cells_diseased, by=c('Cell ID B220'='Cell ID'),
                                                        suffix=c('', '.B220'))
CD3_to_KLRG1_B220_diseased = CD3_to_B220_diseased %>% left_join(KLRG1_cells_diseased, by=c('Cell ID KLRG1'='Cell ID'),
                                                                suffix=c('', '.KLRG1'))

# find the triplet interactions (the same exact three cells) in common between the three viewpoints for the control group
common_in_B220_control <- inner_join(CD3_to_KLRG1_B220_control, KLRG1_to_B220_CD3_control, by = c("Name.B220"))
common_in_all_control <- inner_join(common_in_B220_control, B220_to_KLRG1_CD3_control, by = c("Name.CD3"), "Name.KLRG1")

# find the triplet interactions (the same exact three cells) in common between the three viewpoints for the diseased group
common_in_B220_diseased <- inner_join(CD3_to_KLRG1_B220_diseased, KLRG1_to_B220_CD3_diseased, by = c("Name.B220"))
common_in_all_diseased <- inner_join(common_in_B220_diseased, B220_to_KLRG1_CD3_diseased, by = c("Name.CD3"), "Name.KLRG1")

# leave only one copy of the common triplet interactions for the control group
CD3_to_KLRG1_B220_control <- anti_join(CD3_to_KLRG1_B220_control, common_in_all_control, by = c("Name" = "Name.x"))
KLRG1_to_B220_CD3_control <- anti_join(KLRG1_to_B220_CD3_control, common_in_all_control, by = c("Name" = "Name.y"))

# leave only one copy of the common triplet interactions for the diseased group
CD3_to_KLRG1_B220_diseased <- anti_join(CD3_to_KLRG1_B220_diseased, common_in_all_diseased, by = c("Name" = "Name.x"))
KLRG1_to_B220_CD3_diseased <- anti_join(KLRG1_to_B220_CD3_diseased, common_in_all_diseased, by = c("Name" = "Name.y"))

# Calculate the Euclidean distance between the other two cells for the control group
B220_to_KLRG1_CD3_control <- B220_to_KLRG1_CD3_control %>%
    mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.KLRG1", "Cell Y Position.KLRG1"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.KLRG1`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.KLRG1`)^2)
    )

KLRG1_to_B220_CD3_control <- KLRG1_to_B220_CD3_control %>%
    mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.B220", "Cell Y Position.B220"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.B220`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.B220`)^2)
    )

CD3_to_KLRG1_B220_control <- CD3_to_KLRG1_B220_control %>%
    mutate_at(vars("Cell X Position.B220", "Cell Y Position.B220", "Cell X Position.KLRG1", "Cell Y Position.KLRG1"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.KLRG1` - `Cell X Position.B220`)^2 + (`Cell Y Position.KLRG1` - `Cell Y Position.B220`)^2)
    )

# Calculate the Euclidean distance between the other two cells for the diseased group
B220_to_KLRG1_CD3_diseased <- B220_to_KLRG1_CD3_diseased %>%
    mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.KLRG1", "Cell Y Position.KLRG1"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.KLRG1`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.KLRG1`)^2)
    )

KLRG1_to_B220_CD3_diseased <- KLRG1_to_B220_CD3_diseased %>%
    mutate_at(vars("Cell X Position.CD3", "Cell Y Position.CD3", "Cell X Position.B220", "Cell Y Position.B220"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.CD3` - `Cell X Position.B220`)^2 + (`Cell Y Position.CD3` - `Cell Y Position.B220`)^2)
    )

CD3_to_KLRG1_B220_diseased <- CD3_to_KLRG1_B220_diseased %>%
    mutate_at(vars("Cell X Position.B220", "Cell Y Position.B220", "Cell X Position.KLRG1", "Cell Y Position.KLRG1"), as.numeric) %>%
    mutate(
        `3rd distance` = sqrt((`Cell X Position.KLRG1` - `Cell X Position.B220`)^2 + (`Cell Y Position.KLRG1` - `Cell Y Position.B220`)^2)
    )

# get the triplet distance for the control group
B220_to_KLRG1_CD3_control <- B220_to_KLRG1_CD3_control %>%
    mutate(sum = `Distance to KLRG1` + `Distance to CD3` + `3rd distance`)
KLRG1_to_B220_CD3_control <- KLRG1_to_B220_CD3_control %>%
    mutate(sum = `Distance to CD3` + `Distance to B220` + `3rd distance`)
CD3_to_KLRG1_B220_control <- CD3_to_KLRG1_B220_control %>%
    mutate(sum = `Distance to KLRG1` + `Distance to B220` + `3rd distance`)

# get the triplet distance for the diseased group
B220_to_KLRG1_CD3_diseased <- B220_to_KLRG1_CD3_diseased %>%
    mutate(sum = `Distance to KLRG1` + `Distance to CD3` + `3rd distance`)
KLRG1_to_B220_CD3_diseased <- KLRG1_to_B220_CD3_diseased %>%
    mutate(sum = `Distance to CD3` + `Distance to B220` + `3rd distance`)
CD3_to_KLRG1_B220_diseased <- CD3_to_KLRG1_B220_diseased %>%
    mutate(sum = `Distance to KLRG1` + `Distance to B220` + `3rd distance`)

# add all tables together to get the total unique triplet interactions
column_names <- names(B220_to_KLRG1_CD3_control)
KLRG1_to_B220_CD3_control <- KLRG1_to_B220_CD3_control %>%
  setNames(column_names)
CD3_to_KLRG1_B220_control <- CD3_to_KLRG1_B220_control %>%
  setNames(column_names)
triplet_control_KLRG1 <- rbind(B220_to_KLRG1_CD3_control, KLRG1_to_B220_CD3_control, CD3_to_KLRG1_B220_control)

# add all tables together to get the total unique triplet interactions
column_names <- names(B220_to_KLRG1_CD3_diseased)
KLRG1_to_B220_CD3_diseased <- KLRG1_to_B220_CD3_diseased %>%
  setNames(column_names)
CD3_to_KLRG1_B220_diseased <- CD3_to_KLRG1_B220_diseased %>%
  setNames(column_names)
triplet_diseased_KLRG1 <- rbind(B220_to_KLRG1_CD3_diseased, KLRG1_to_B220_CD3_diseased, CD3_to_KLRG1_B220_diseased)

# Calculate the average and standard deviation of the "sum" column
sum_summary_control_KLRG1 <- summary(triplet_control_KLRG1$sum)
average_sum_control_KLRG1 <- mean(triplet_control_KLRG1$sum, na.rm = TRUE)
std_dev_sum_control_KLRG1 <- sd(triplet_control_KLRG1$sum, na.rm = TRUE)

# Calculate the average and standard deviation of the "sum" column
sum_summary_diseased_KLRG1 <- summary(triplet_diseased_KLRG1$sum)
average_sum_diseased_KLRG1 <- mean(triplet_diseased_KLRG1$sum, na.rm = TRUE)
std_dev_sum_diseased_KLRG1 <- sd(triplet_diseased_KLRG1$sum, na.rm = TRUE)

#print
cat("Average of sum column control KLRG1:", average_sum_control_KLRG1, "\n")
cat("Standard deviation of sum column control KLRG1:", std_dev_sum_control_KLRG1, "\n")

cat("Average of sum column diseased KLRG1:", average_sum_diseased_KLRG1, "\n")
cat("Standard deviation of sum column diseased KLRG1:", std_dev_sum_diseased_KLRG1, "\n")


cat("Average of sum column control CD11b:", average_sum_control_CD11b, "\n")
cat("Standard deviation of sum column control CD11b:", std_dev_sum_control_CD11b, "\n")

cat("Average of sum column diseased CD11b:", average_sum_diseased_CD11b, "\n")
cat("Standard deviation of sum column diseased CD11b:", std_dev_sum_diseased_CD11b, "\n")


# Create a data frame with the average and standard deviation values for KLRG1 and CD11b in the control and diseased groups
data <- data.frame(
  Group = rep(c("CD11b", "KLRG1"), each = 2),
  Phenotype = rep(c("Control", "Diseased"), 2),
  Average = c(average_sum_control_CD11b, average_sum_diseased_CD11b, average_sum_control_KLRG1, average_sum_diseased_KLRG1),
  StdDev = c(std_dev_sum_control_CD11b, std_dev_sum_diseased_CD11b, std_dev_sum_control_KLRG1, std_dev_sum_diseased_KLRG1)
)

ggplot(data, aes(x = Group, y = Average, fill = Phenotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Average - StdDev, ymax = Average + StdDev), width = 0.25, position = position_dodge(0.9)) +
  labs(title = "Average and Standard Deviation of KLRG1 and CD11b",
       x = "Group",
       y = "Average Value") +
  scale_fill_manual(values = c("Control" = "blue", "Diseased" = "red")) +
  theme_minimal()

