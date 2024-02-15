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


# Initialize empty vectors to store threshold distances and corresponding Tn values for the diseased group
threshold_values <- seq(0, 550, by = 1)
Tn_values_control <- c()
Tn_values_diseased <- c()

# Loop through each threshold distance for the diseased group
for (threshold_distance in threshold_values) {
    
    # filter the triplet interaction tables based on the threshold for the control group
    B220_to_CD11b_CD3_filtered_control <- B220_to_CD11b_CD3_control %>%
        filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to CD11b` < threshold_distance)
    
    CD11b_to_B220_CD3_filtered_control <- CD11b_to_B220_CD3_control %>%
        filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to B220` < threshold_distance)
    
    CD3_to_CD11b_B220_filtered_control <- CD3_to_CD11b_B220_control %>%
        filter(`3rd distance` < threshold_distance & `Distance to B220` < threshold_distance & `Distance to CD11b` < threshold_distance)
    
    # filter the triplet interaction tables based on the threshold for the diseased group
    B220_to_CD11b_CD3_filtered_diseased <- B220_to_CD11b_CD3_diseased %>%
        filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to CD11b` < threshold_distance)
    
    CD11b_to_B220_CD3_filtered_diseased <- CD11b_to_B220_CD3_diseased %>%
        filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to B220` < threshold_distance)
    
    CD3_to_CD11b_B220_filtered_diseased <- CD3_to_CD11b_B220_diseased %>%
        filter(`3rd distance` < threshold_distance & `Distance to B220` < threshold_distance & `Distance to CD11b` < threshold_distance)
    
    # add all filtered tables together to get the total unique triplet interactions for the control group
    column_names_filtered_control <- names(B220_to_CD11b_CD3_filtered_control)
    CD11b_to_B220_CD3_filtered_control <- CD11b_to_B220_CD3_filtered_control %>%
        setNames(column_names_filtered_control)
    CD3_to_CD11b_B220_filtered_control <- CD3_to_CD11b_B220_filtered_control %>%
        setNames(column_names_filtered_control)
    triplet_filtered_control <- rbind(B220_to_CD11b_CD3_filtered_control, CD11b_to_B220_CD3_filtered_control, CD3_to_CD11b_B220_filtered_control)
    
    # add all filtered tables together to get the total unique triplet interactions for the diseased group
    column_names_filtered_diseased <- names(B220_to_CD11b_CD3_filtered_diseased)
    CD11b_to_B220_CD3_filtered_diseased <- CD11b_to_B220_CD3_filtered_diseased %>%
        setNames(column_names_filtered_diseased)
    CD3_to_CD11b_B220_filtered_diseased <- CD3_to_CD11b_B220_filtered_diseased %>%
        setNames(column_names_filtered_diseased)
    triplet_filtered_diseased <- rbind(B220_to_CD11b_CD3_filtered_diseased, CD11b_to_B220_CD3_filtered_diseased, CD3_to_CD11b_B220_filtered_diseased)
    
    # get the current Tn value for the diseased group
    Tn_control <- nrow(triplet_filtered_control)
    
    # Store the current Tn value in the vector for the diseased group
    Tn_values_control <- c(Tn_values_control, Tn_control)
    
    # get the current Tn value for the diseased group
    Tn_diseased <- nrow(triplet_filtered_diseased)
    
    # Store the current Tn value in the vector for the diseased group
    Tn_values_diseased <- c(Tn_values_diseased, Tn_diseased)
}

# Create data frames with threshold values and corresponding Tn values for both groups
result_df_control_CD11b <- data.frame(threshold_distance = threshold_values, Tn = Tn_values_control)
result_df_diseased_CD11b <- data.frame(threshold_distance = threshold_values, Tn = Tn_values_diseased)

result_df_control_CD11b$Group <- 'Control'
result_df_diseased_CD11b$Group <- 'Diseased'


###################################################################################################

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


# Initialize empty vectors to store threshold distances and corresponding Tn values for the diseased group
threshold_values <- seq(0, 550, by = 1)
Tn_values_control <- c()
Tn_values_diseased <- c()

# Loop through each threshold distance for the diseased group
for (threshold_distance in threshold_values) {
    
    # filter the triplet interaction tables based on the threshold for the control group
    B220_to_KLRG1_CD3_filtered_control <- B220_to_KLRG1_CD3_control %>%
        filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to KLRG1` < threshold_distance)
    
    KLRG1_to_B220_CD3_filtered_control <- KLRG1_to_B220_CD3_control %>%
        filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to B220` < threshold_distance)
    
    CD3_to_KLRG1_B220_filtered_control <- CD3_to_KLRG1_B220_control %>%
        filter(`3rd distance` < threshold_distance & `Distance to B220` < threshold_distance & `Distance to KLRG1` < threshold_distance)
    
    # filter the triplet interaction tables based on the threshold for the diseased group
    B220_to_KLRG1_CD3_filtered_diseased <- B220_to_KLRG1_CD3_diseased %>%
        filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to KLRG1` < threshold_distance)
    
    KLRG1_to_B220_CD3_filtered_diseased <- KLRG1_to_B220_CD3_diseased %>%
        filter(`3rd distance` < threshold_distance & `Distance to CD3` < threshold_distance & `Distance to B220` < threshold_distance)
    
    CD3_to_KLRG1_B220_filtered_diseased <- CD3_to_KLRG1_B220_diseased %>%
        filter(`3rd distance` < threshold_distance & `Distance to B220` < threshold_distance & `Distance to KLRG1` < threshold_distance)
    
    # add all filtered tables together to get the total unique triplet interactions for the control group
    column_names_filtered_control <- names(B220_to_KLRG1_CD3_filtered_control)
    KLRG1_to_B220_CD3_filtered_control <- KLRG1_to_B220_CD3_filtered_control %>%
        setNames(column_names_filtered_control)
    CD3_to_KLRG1_B220_filtered_control <- CD3_to_KLRG1_B220_filtered_control %>%
        setNames(column_names_filtered_control)
    triplet_filtered_control <- rbind(B220_to_KLRG1_CD3_filtered_control, KLRG1_to_B220_CD3_filtered_control, CD3_to_KLRG1_B220_filtered_control)
    
    # add all filtered tables together to get the total unique triplet interactions for the diseased group
    column_names_filtered_diseased <- names(B220_to_KLRG1_CD3_filtered_diseased)
    KLRG1_to_B220_CD3_filtered_diseased <- KLRG1_to_B220_CD3_filtered_diseased %>%
        setNames(column_names_filtered_diseased)
    CD3_to_KLRG1_B220_filtered_diseased <- CD3_to_KLRG1_B220_filtered_diseased %>%
        setNames(column_names_filtered_diseased)
    triplet_filtered_diseased <- rbind(B220_to_KLRG1_CD3_filtered_diseased, KLRG1_to_B220_CD3_filtered_diseased, CD3_to_KLRG1_B220_filtered_diseased)
    
    # get the current Tn value for the diseased group
    Tn_control <- nrow(triplet_filtered_control)
    
    # Store the current Tn value in the vector for the diseased group
    Tn_values_control <- c(Tn_values_control, Tn_control)
    
    # get the current Tn value for the diseased group
    Tn_diseased <- nrow(triplet_filtered_diseased)
    
    # Store the current Tn value in the vector for the diseased group
    Tn_values_diseased <- c(Tn_values_diseased, Tn_diseased)
}

# Create data frames with threshold values and corresponding Tn values for both groups
result_df_control_KLRG1 <- data.frame(threshold_distance = threshold_values, Tn = Tn_values_control)
result_df_diseased_KLRG1 <- data.frame(threshold_distance = threshold_values, Tn = Tn_values_diseased)

result_df_control_KLRG1$Group <- 'Control'
result_df_diseased_KLRG1$Group <- 'Diseased' 

# Create a combined plot for CD11b and KLRG1 with merged color and linetype legends
combined_plot <- ggplot() +
    geom_line(data = result_df_control_CD11b, aes(x = threshold_distance, y = Tn, color = "Control - CD11b", linetype = "Control - CD11b")) +
    geom_line(data = result_df_diseased_CD11b, aes(x = threshold_distance, y = Tn, color = "Diseased - CD11b", linetype = "Diseased - CD11b")) +
    geom_line(data = result_df_control_KLRG1, aes(x = threshold_distance, y = Tn, color = "Control - KLRG1", linetype = "Control - KLRG1")) +
    geom_line(data = result_df_diseased_KLRG1, aes(x = threshold_distance, y = Tn, color = "Diseased - KLRG1", linetype = "Diseased - KLRG1")) +
    scale_color_manual(values = c("Control - CD11b" = "black", "Diseased - CD11b" = "red", "Control - KLRG1" = "black", "Diseased - KLRG1" = "red")) +
    scale_linetype_manual(values = c("Control - CD11b" = "solid", "Diseased - CD11b" = "solid", "Control - KLRG1" = "dashed", "Diseased - KLRG1" = "dashed")) +
    labs(
        title = "CD11b and KLRG1: Control and Diseased Groups - Tn vs Td",
        x = "Threshold Distance (Td)",
        y = "Triplet Interactions (Tn)",
        color = "Group",
        linetype = "Group"  # Use the same label for color and linetype legends
    ) +
    theme(
        panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Set axis lines to black
        legend.position = "bottom"  # Move legend to the bottom
    )

# Display the combined plot
print(combined_plot)