# This code was written by TC and is used to produce plots in the paper 
# Interpretable Clustering: Balancing Representative Groupings with Interpretability.

# The code contained in this script is protected by a V3 creative commons licence
# If you have any questions regarding this code, please contact Tommy Cheale at Tommycheale@gmail.com

# Load in required packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(fda)
library(cluster)     
library(factoextra)  
library(dendextend)  
library(geometry)    
library(plotly)      
library(patchwork)
library(stringr)     
library(purrr)
library(cluster)
library(NbClust)
library(grid)
library(gridBase)
library(ggrepel)
library(reshape2)
library(FactoMineR)
library(mclust)
library(ggdendro)
library(MASS)
library(gridExtra)
library(readr)

###############################################################################
################################ Data Processing ##############################
###############################################################################

# Load data from 2019 SOCBR
SOCB_2019 <- read_csv("All_SOCB_2019_population_indices.csv") %>%
  select(Species = English_Name, year, index) %>%
  arrange(Species, year)

# Load group data
SOCB_2019_habitat <- read_csv("SOCB_2019_habitat_groups.csv")

# Filter to after 1974
SOCB_2019 <- SOCB_2019 %>%
  filter(year >= 1975)

# Identify species with incomplete data
SOCB_2019 <- SOCB_2019 %>%
  group_by(Species) %>%
  mutate(Complete = all(1975:2016 %in% year)) %>%
  ungroup()

incomplete_species <- SOCB_2019 %>% 
  filter(!Complete) |> 
  pull(Species) |> 
  unique()
  
# Remove species with incomplete data
SOCB_2019 <- SOCB_2019 %>%
  filter(Complete) %>%
  select(-Complete)

# Compute max-rescaled indices
SOCB_2019 <- SOCB_2019 %>%
  group_by(Species) %>%
  mutate(mr_index = index / max(index)) |> 
  ungroup()

# Create wide data
SOCB_2019_wide <- SOCB_2019 |> 
  select(Species, year, mr_index) |> 
  pivot_wider(names_from = "Species", values_from = "mr_index")

# Create data matrix required by funcitions in the fda package
SOCB_2019_matrix <- SOCB_2019_wide |> 
  select(-year) |> 
  as.matrix()

rownames(SOCB_2019_matrix) <- SOCB_2019_wide |> pull(year)
            
################################################################################
################################ Data Explorattion #############################
################################################################################

# Compute mean trends by group
SOCB_2019_mean_guild <- SOCB_2019 |> 
  left_join(SOCB_2019_habitat, by = "Species") |> 
  group_by(Group, year) |> 
  summarize(mr_mean = mean(mr_index))

# Plot rescaled curves for select groups
SOCB_2019 |> 
  left_join(SOCB_2019_habitat, by = "Species") |> 
  filter(Group %in% (c("Forest Birds", "Aerial Insectivores", "Shorebirds"))) |> 
  ggplot(aes(x = year, y = mr_index, group = Species, color = Group)) +
  geom_line(alpha = 0.4) +
  facet_wrap(~ Group) +  # individual species lines
  geom_line(
    data = SOCB_2019_mean_guild |> filter(Group %in% (c("Forest Birds", "Aerial Insectivores", "Shorebirds"))),
    aes(x = year, y = mr_mean, color = Group),
    inherit.aes = FALSE,
    size = 1.5
  ) +
  scale_color_manual(values = c(
    "Shorebirds" = "#ff7f0e",
    "Forest Birds" = "#8c564b",
    "Aerial Insectivores" = "#7f7f7f"
  )) +
  labs(
    title = "Rescaled Abundance Over Time for Selected Guilds",
    y = "Rescaled Annual Abundance",
    x = "Year",
    color = "Guild"
  ) + 
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 25),
        plot.title = element_text(size = 35, face = "bold"),        # main title
        axis.title = element_text(size = 25, face = "bold"),        # x and y axis titles
        axis.text = element_text(size = 18)  
  )

ggsave("plots/Examples.pdf", width = 12, height =8)


###############################################################################
################################## Resmoothing ################################
###############################################################################

# Fit B-splines with between 12 and 18 basis functions

time_points <- 1:nrow(SOCB_2019_matrix)

bspline_fits <- 4:42 |> 
  purrr::map(\(k){
    basis_obj <- create.bspline.basis(rangeval = range(time_points), nbasis = k)
    smooth_result <- smooth.basis(argvals = time_points, 
                                  y = SOCB_2019_matrix, 
                                  fdParobj = basis_obj)
    fd_obj <- smooth_result$fd
    fit <- eval.fd(time_points,fd_obj)
    mae <- max(abs(SOCB_2019_matrix - fit))
    mse <- mean((SOCB_2019_matrix - fit)^2)
    list(df = k, fd = fd_obj, fit = fit, mae = mae, mse = mse)
  })

# Compare mean absolute error and mean squared error
bspline_fits |> 
  purrr::map(\(fit) tibble(df = fit$df, mae = fit$mae, mse = fit$mse)) |> 
  bind_rows() |> 
  print(n = Inf)

###############################################################################
################################ Stage 1: FPCA ################################
###############################################################################

# Conduct FPCA using selected number of basis functions
k <- 15
fpca_result <- pca.fd(bspline_fits[[k - 3]]$fd, nharm = k) # Perform FPCA

# Scree plot
scree_df <- tibble(eigenvalue = fpca_result$values,
                   sd = sqrt(eigenvalue)) |> 
  mutate(pc = 1:n(),
         pve = eigenvalue / sum(eigenvalue),
         cum_pve = cumsum(pve))

scree_df |> 
  filter(pc <= 8) |> 
  ggplot(aes(x = pc, y = pve)) +
  geom_bar(stat = "identity", fill = "blue", color = "blue", alpha = 0.7) +
  geom_line(aes(y = cum_pve, group = 1), color = "darkgreen", linewidth = 1) +
  geom_point(aes(y = cum_pve), color = "darkgreen", size = 3, shape = 17) +
  labs(title = "Proportion of Variance Captured by Each Functional Principal Component",
       x = "Functional Principal Component",
       y = "Proportion of Variance Explained") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 22),  # Increase x-axis text size
    axis.text.y = element_text(size = 22),  # Increase y-axis text size
    axis.title.x = element_text(size = 22),  # Increase x-axis title size
    axis.title.y = element_text(size = 22),  # Increase y-axis title size
    plot.title = element_text(size = 25, hjust = 0.5)  # Increase title size and center it
  ) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Cumulative Variance"))

# Compute Cos2 scores (proportion of variance explained by first 2 PCs)
Total_var<- fpca_result$scores^2 %*% fpca_result$values
Partial_var<- fpca_result$scores[,1:2]^2 %*% fpca_result$values[1:2]

Cos2_scores <- Partial_var/Total_var

# Construct output data frame
fpca_output <- tibble(Species = colnames(SOCB_2019_matrix),
                      PC1 = fpca_result$scores[,1],
                      PC2 = fpca_result$scores[,2],
                      Cos2 = Cos2_scores) |> 
  left_join(SOCB_2019_habitat, by = "Species")

################################################################################
#Filter species with low cos 2 scores (specifically those that don't reach 0.75)
################################################################################

low_cos2_species <- DF_obj1_1975[DF_obj1_1975$Cos_2_scores_12 < 0.75, ]
low_cos2_species$Name <- rownames(low_cos2_species)

############################################################################
########################### Plot original curves ###########################
############################################################################

#Use of the re-scale maximum function here
PCA_Max_1975<-rescale_Max(PCA_Form_123_1975)

#plot all the data
Max_pop_tibble_1975 <- as.data.frame(PCA_Max_1975)

# Generate years based on the number of rows in the matrix
years <- 1975:2016

# Add years as a column
Max_pop_tibble_1975$year <- years

# Melt the tibble into long format
Max_pop_tibble_1975 <- Max_pop_tibble_1975 %>%
  pivot_longer(cols = -year, names_to = "species", values_to = "Max_rescaled_population")

cleaned_data$species<-cleaned_data$Name

test<-cleaned_data[cleaned_data$species == "Gray Jay",]

# Join the two datasets by species
joined_data <- Max_pop_tibble_1975 %>%
  left_join(cleaned_data, by = "species")

joined_data
ars<-joined_data[joined_data$species == "Alder Flycatcher",]


guilds <- c("Waterfowl", "Birds of Prey", "Wetland Birds", "Forest Birds",
            "Aerial Insectivores", "Grassland Birds", "Shorebirds",
            "Seabirds", "Duck Species", "Native",
            "Goose Species", "Forest Crop Specilists",
            "Tolerant to Farming",
            "Dependent on Grassland", 
            "Long distance shorebirds",
            "Short distance shorebirds", "Forest birds wintering in SA", 
            "Forest birds wintering in CA")

guild_vectors_PC1<-list()


# Loop over each guild and calculate the mean vector for that guild
for (guild in guilds) {
  # Filter species that belong to the current guild
  guild_data <- joined_data %>%
    filter(get(guild) == 1)
  
  # Ensure you're referencing the correct population column
  guild_vector <- guild_data %>%
    group_by(year) %>%
    summarise(mean_PC1 = mean(PC1, na.rm = TRUE))  
  # Store the result in the list
  guild_vectors_PC1[[guild]] <- guild_vector
}


guild_vectors_PC2<-list()


# Loop over each guild and calculate the mean vector for that guild
for (guild in guilds) {
  # Filter species that belong to the current guild
  guild_data <- joined_data %>%
    filter(get(guild) == 1)
  
  # Ensure you're referencing the correct population column
  guild_vector <- guild_data %>%
    group_by(year) %>%
    summarise(mean_PC2 = mean(PC2, na.rm = TRUE))  
  # Store the result in the list
  guild_vectors_PC2[[guild]] <- guild_vector
}


guild_vectors_PC3<-list()

# Loop over each guild and calculate the mean vector for that guild
for (guild in guilds) {
  # Filter species that belong to the current guild
  guild_data <- joined_data %>%
    filter(get(guild) == 1)
  
  # Ensure you're referencing the correct population column
  guild_vector <- guild_data %>%
    group_by(year) %>%
    summarise(mean_PC3 = mean(PC3, na.rm = TRUE))  
  # Store the result in the list
  guild_vectors_PC3[[guild]] <- guild_vector
}


##############################################################################
##################################Find mean curves for the groups ############
##############################################################################

# Extract the mean_PC1 values and create a dataframe for each guild
mean_PC1_df <- do.call(rbind, lapply(names(guild_vectors_PC1), function(guild) {
  data.frame(
    guild = guild,
    mean_PC1 = guild_vectors_PC1[[guild]]$mean_PC1  # Assuming mean_population represents mean_PC1
  )
}))

# Extract the mean_PC2 values and create a dataframe for each guild
mean_PC2_df <- do.call(rbind, lapply(names(guild_vectors_PC2), function(guild) {
  data.frame(
    guild = guild,
    mean_PC2 = guild_vectors_PC2[[guild]]$mean_PC2  # Assuming mean_population represents mean_PC2
  )
}))


# Extract the mean_PC2 values and create a dataframe for each guild
mean_PC3_df <- do.call(rbind, lapply(names(guild_vectors_PC3), function(guild) {
  data.frame(
    guild = guild,
    mean_PC3 = guild_vectors_PC3[[guild]]$mean_PC3  # Assuming mean_population represents mean_PC2
  )
}))

# Merge the two dataframes on the 'guild' column
final_df_it <- merge(mean_PC1_df, mean_PC2_df, by = "guild")
final_df <- merge(final_df_it, mean_PC3_df,  by="guild")

df_unique_guilds <- final_df[!duplicated(final_df), ]


rotations<-pca_L_Scale_max_1975$rotation

rotations<-data.frame(rotations)

rot1<-rotations$PC1
rot2<-rotations$PC2
rot3<-rotations$PC3


cleaned_data$Rotation_12 <- mapply(
  function(pc1, pc2) pc1 * rot1 + pc2 * rot2,
  cleaned_data$PC1, cleaned_data$PC2, SIMPLIFY = FALSE
)


cleaned_data$Rotation_13 <- mapply(
  function(pc1, pc3) pc1 * rot1 + pc3 * rot3,
  cleaned_data$PC1, cleaned_data$PC3, SIMPLIFY = FALSE
)


cleaned_data$Rotation_23 <- mapply(
  function(pc2, pc3) pc2 * rot2 + pc3 * rot3,
  cleaned_data$PC2, cleaned_data$PC3, SIMPLIFY = FALSE
)


cleaned_data$Rotation_123 <- mapply(
  function(pc1, pc2, pc3) pc1 * rot1 + pc2 * rot2 + pc3 * rot3,
  cleaned_data$PC1, cleaned_data$PC2, cleaned_data$PC3, SIMPLIFY = FALSE
)


PC1<-cleaned_data$PC1
PC2<-cleaned_data$PC2
PC3<-cleaned_data$PC3

cleaned_data

Forest_Nam<-cleaned_data[cleaned_data$`Forest Birds` == 1,]$species

Shorebirds_Nam<-cleaned_data[cleaned_data$Shorebirds == 1,]$species

AI_Nam<-cleaned_data[cleaned_data$`Aerial Insectivores` == 1,]$species

all_species <- c(Forest_Nam, Shorebirds_Nam, AI_Nam)

no_filter_rescaled_df <- no_filter_filtered_df_1975 %>%
  rowwise() %>%
  mutate(across(-English_Name, ~ . / max(c_across(-English_Name), na.rm = TRUE))) %>%
  ungroup()


no_filter_rescaled_df <- no_filter_rescaled_df %>%
  drop_na()

# Filter dataset for these species
species_df <- no_filter_rescaled_df %>%
  filter(English_Name %in% all_species)

# Convert from wide to long format
species_long <- species_df %>%
  pivot_longer(
    cols = -English_Name,
    names_to = "Year",
    values_to = "Abundance"
  ) %>%
  mutate(Year = as.numeric(Year))

# Add a Guild column for coloring/faceting
species_long <- species_long %>%
  mutate(Guild = case_when(
    English_Name %in% Forest_Nam ~ "Forest birds",
    English_Name %in% Shorebirds_Nam ~ "Shorebirds",
    English_Name %in% AI_Nam ~ "Aerial insectivores",
    TRUE ~ "Other"
  ))

species_long <- species_long %>%
  mutate(
    Guild_Color = case_when(
      Guild == "Shorebirds" ~ "#ff7f0e",
      Guild == "Forest birds" ~ "#8c564b",
      Guild == "Aerial insectivores" ~ "#7f7f7f",
      TRUE ~ "#000000"
    )
  )

###############################################################################################
#########################################Plot Aerial insectivores, Forest Birds, and Shorebirds
###############################################################################################
# ggplot(species_long, aes(x = Year, y = Abundance, group = English_Name, color = Guild)) +
#   geom_line(alpha = 0.6) +  # lighter lines for individual species
#   facet_wrap(~ Guild, scales = "free_y") +
#   scale_color_manual(values = c(
#     "Shorebirds" = "#ff7f0e",
#     "Forest birds" = "#8c564b",
#     "Aerial insectivores" = "#7f7f7f",
#     "Other" = "#000000"
#   )) +
#   labs(
#     title = "Abundance Over Time for Selected Species",
#     y = "Abundance",
#     color = "Guild"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "none")


SOCB_2019_mean_guild <- species_long %>%
  group_by(Guild, Year) %>%
  summarize(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

################################################################################
##############################Plot with mean curves ############################
################################################################################
ggplot(species_long, aes(x = Year, y = Abundance, group = English_Name, color = Guild)) +
  geom_line(alpha = 0.4) +  # individual species lines
  geom_line(
    data = SOCB_2019_mean_guild,
    aes(x = Year, y = Mean_Abundance, color = Guild),
    inherit.aes = FALSE,
    size = 1.5
  ) +
  facet_wrap(~ Guild, scales = "free_y") +
  scale_color_manual(values = c(
    "Shorebirds" = "#ff7f0e",
    "Forest birds" = "#8c564b",
    "Aerial insectivores" = "#7f7f7f",
    "Other" = "#000000"
  )) +
  labs(
    title = "Rescaled Abundance Over Time for Selected Guilds",
    y = "Rescaled Annual Abundance",
    color = "Guild"
  ) + 
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 25),
        plot.title = element_text(size = 35, face = "bold"),        # main title
        axis.title = element_text(size = 25, face = "bold"),        # x and y axis titles
        axis.text = element_text(size = 18)  
  )


#Final tag analysis
uniquely_tagged_species <- cleaned_data

# Define category columns
categories <- c("Grassland Birds", "Birds of Prey", "Shorebirds", "Aerial Insectivores", 
                "Forest Birds", "Seabirds", "Wetland Birds")

# Function to reset all categories for a species
reset_categories <- function(df, species_name) {
  df[df$Name == species_name, categories] <- 0
  return(df)
}

#############################################################################################
######################################Code to manually assign species to tags ###############
#############################################################################################

# Assign species in both Birds of Prey and Grassland Birds to Birds of Prey
bp_gb <- uniquely_tagged_species$`Birds of Prey` == 1 & uniquely_tagged_species$`Grassland Birds` == 1
uniquely_tagged_species$`Grassland Birds`[bp_gb] <- 0  # Remove from Grassland Birds

# Assign specific species in Grassland Birds and Shorebirds to their respective groups
species_gb <- c("Mountain Plover", "Upland Sandpiper", "Eskimo Curlew")
species_sb <- c("Willet", "Long-billed Curlew")

for (sp in species_gb) {
  uniquely_tagged_species <- reset_categories(uniquely_tagged_species, sp)
  uniquely_tagged_species[uniquely_tagged_species$Name == sp, "Grassland Birds"] <- 1
}

for (sp in species_sb) {
  uniquely_tagged_species <- reset_categories(uniquely_tagged_species, sp)
  uniquely_tagged_species[uniquely_tagged_species$Name == sp, "Shorebirds"] <- 1
}

# Assign species in both Aerial Insectivores and Forest Birds to Aerial Insectivores
ai_fb <- uniquely_tagged_species$`Aerial Insectivores` == 1 & uniquely_tagged_species$`Forest Birds` == 1
uniquely_tagged_species$`Forest Birds`[ai_fb] <- 0  # Remove from Forest Birds

# Assign species in both Shorebirds and Forest Birds to Shorebirds
sb_fb <- uniquely_tagged_species$`Shorebirds` == 1 & uniquely_tagged_species$`Forest Birds` == 1
uniquely_tagged_species$`Forest Birds`[sb_fb] <- 0  # Remove from Forest Birds

# Assign species in both Seabirds and Wetland Birds to Seabirds

species_web <- c("Arctic Tern")
species_sea <- c("Pomarine Jaeger","Parasitic Jaeger","Long-tailed Jaeger",
                 "Sabine's Gull","Iceland Gull","Common Tern")

for (sp in species_web) {
  uniquely_tagged_species <- reset_categories(uniquely_tagged_species, sp)
  uniquely_tagged_species[uniquely_tagged_species$Name == sp, "Wetland Birds"] <- 1
}

for (sp in species_sea) {
  uniquely_tagged_species <- reset_categories(uniquely_tagged_species, sp)
  uniquely_tagged_species[uniquely_tagged_species$Name == sp, "Seabirds"] <- 1
}
################################################################################
################################################################################

################################################################################
################ Final Tags Analysis ###########################################
################################################################################

group_columns <- c( "Waterfowl", "Birds of Prey", "Wetland Birds", 
                    "Seabirds", "Forest Birds", "Other",
                    "Shorebirds", "Aerial Insectivores", "Grassland Birds"
)

# Convert each group column to numeric
uniquely_tagged_species[, group_columns] <- lapply(uniquely_tagged_species[, group_columns], function(x) as.numeric(as.character(x)))

# Now calculate the sums
species_per_group <- colSums(uniquely_tagged_species[, group_columns], na.rm = TRUE)

# Identify species that are in more than one group:
uniquely_tagged_species$group_count <- rowSums(uniquely_tagged_species[, group_columns])
overlap_species <- uniquely_tagged_species[uniquely_tagged_species$group_count > 1, ]

# View the specific overlaps, for example, by selecting the relevant columns:
head(overlap_species[, c("Name", group_columns)])


n <- length(group_columns)
overlap_matrix <- matrix(0, nrow = n, ncol = n)
rownames(overlap_matrix) <- group_columns
colnames(overlap_matrix) <- group_columns

# Calculate counts for each pair of groups
for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      # Count species that are exclusively in group i:
      overlap_matrix[i, j] <- sum(uniquely_tagged_species[, group_columns[i]] == 1 & rowSums(uniquely_tagged_species[, group_columns]) == 1, na.rm = TRUE)
    } else {
      # Count species that are in both group i and group j:
      overlap_matrix[i, j] <- sum(uniquely_tagged_species[, group_columns[i]] == 1 & uniquely_tagged_species[, group_columns[j]] == 1, na.rm = TRUE)
    }
  }
}

# View the matrix
print(overlap_matrix)


# Melt the matrix into long format for ggplot2
overlap_df <- melt(overlap_matrix)
names(overlap_df) <- c("Group1", "Group2", "Count")

# Heatmap to check that there are no overlaps left. This should be just a bunch of zeros.

# ggplot(overlap_df, aes(x = Group1, y = Group2, fill = Count)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = Count), color = "black", size = 3) +
#   scale_fill_gradient(low = "lightblue", high = "darkblue") +
#   theme_minimal() +
#   labs(title = "Species Overlap Between Groups",
#        x = "Group",
#        y = "Group") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Select other
main_Func<-uniquely_tagged_species[uniquely_tagged_species$Other == "1",]


################################################################################
###############################Forming initial groups #### #####################
################################################################################

cleaned_data<-uniquely_tagged_species

tag_cols <- c("Duck Species", "Goose Species", 
              "Forest Crop Specilists", "Forest birds wintering in SA", 
              "Forest birds wintering in CA", "Short distance shorebirds", 
              "Long distance shorebirds", "Tolerant to Farming", "Dependent on Grassland")

Habitiat_category_columns <- c("Waterfowl", "Birds of Prey", "Wetland Birds", 
                               "Seabirds", "Forest Birds",
                               "Shorebirds", "Aerial Insectivores", "Grassland Birds")

#Convert to numeric
cleaned_data[Habitiat_category_columns] <- cleaned_data[Habitiat_category_columns] %>%
  mutate(across(everything(), ~ as.numeric(as.character(.)))) %>%
  replace(is.na(.), 0)

# Find species that belong to NO category
no_category_species <- cleaned_data %>%
  filter(rowSums(dplyr::select(., all_of(Habitiat_category_columns))) == 0)

print(no_category_species[, c("Name", "PC1", "PC2")])


# Filter out species that belong to no category at all
cleaned_data <- cleaned_data %>%
  filter(rowSums(dplyr::select(., all_of(Habitiat_category_columns))) > 0)

cleaned_data <- cleaned_data %>%
  dplyr::mutate(num_categories = rowSums(dplyr::select(., all_of(Habitiat_category_columns))))

single_category_species <- cleaned_data %>%
  filter(num_categories == 1) %>%
  dplyr::mutate(best_group = apply(dplyr::select(., all_of(Habitiat_category_columns)), 1, function(row) {
    Habitiat_category_columns[which(row == 1)]
  }))



cleaned_data <- cleaned_data %>%
  dplyr::mutate(num_categories1 = rowSums(dplyr::select(., all_of(Habitiat_category_columns))))

multi_tagged_species <- cleaned_data %>%
  filter(num_categories1 > 1)

cleaned_data <- cleaned_data %>%
  dplyr::mutate(Original_Group = apply(dplyr::select(., all_of(Habitiat_category_columns)), 1, function(row) {
    Habitiat_category_columns[which(row == 1)]
  }))




################################################################################################
##################################### Breaking down clusters ###################################
################################################################################################


clustering_results <- list()
plots_list <- list()
species_table_list <- list()
tags_table_list <- list()
centers_table_list <- list()
silhouette_plots_list <- list()
elbow_plots_list <- list()

# Loop over each habitat (Original_Group)
for (habitat in Habitiat_category_columns) {
  # Subset data for the current habitat
  habitat_data <- cleaned_data %>% filter(Original_Group == habitat)
  
  # Convert tag columns to numeric
  habitat_data <- habitat_data %>%
    mutate(across(all_of(tag_cols), ~ as.numeric(as.character(.))))
  
  
  #habitat_data <- habitat_data %>% filter(!Outlier)  # Remove KMM outliers
  #habitat_data <- habitat_data %>% filter(!Outlier_LOF)  # Remove LOF outliers 
  
  
  # Initialize a data frame to store average silhouette scores for each candidate k
  sil_summary <- data.frame(k = integer(), avg_sil = double())
  
  # If there is only one species in this habitat, assign it to cluster 1.
  if(nrow(habitat_data) < 2) {
    habitat_data$Cluster <- 1
    optimal_k <- 1
  } else {
    # Consider k from 2 to up to 10 (or less if not enough rows)
    k_range <- 2:min(10, nrow(habitat_data))
    silhouette_avgs <- numeric(length(k_range))
    
    # Pre-compute the distance matrix on PC1 and PC2
    d <- dist(habitat_data[, c("PC1", "PC2")])
    
    # Loop over each candidate k to calculate average silhouette width
    for (i in seq_along(k_range)) {
      k <- k_range[i]
      set.seed(123)  # for reproducibility
      km_res <- kmeans(habitat_data[, c("PC1", "PC2")], centers = k, nstart = 25)
      sil <- silhouette(km_res$cluster, d)
      silhouette_avgs[i] <- mean(sil[, 3])
      sil_summary <- rbind(sil_summary, data.frame(k = k, avg_sil = mean(sil[, 3])))
    }
    
    # Create an elbow plot: average silhouette score vs. k
    elbow_plot <- ggplot(sil_summary, aes(x = k, y = avg_sil)) +
      geom_line() +
      geom_point(size = 3) +
      labs(title = paste(habitat),
           x = "Number of Clusters (k)",
           y = "Average Silhouette Score") +
      theme_minimal()+
      theme(
        axis.title.x = element_text(size = 20),  # Increase x-axis title size
        axis.title.y = element_text(size = 20),  # Increase y-axis title size
        plot.title = element_text(size = 35, hjust = 0.5)  # Increase title size and center it
      )
    
    elbow_plots_list[[habitat]] <- elbow_plot
    
    # Select the k with the maximum average silhouette width
    optimal_index <- which.max(silhouette_avgs)
    optimal_k <- k_range[optimal_index]
    
    # Run final KMeans clustering using the optimal k
    set.seed(123)
    km_final <- kmeans(habitat_data[, c("PC1", "PC2")], centers = optimal_k, nstart = 25)
    habitat_data$Cluster <- km_final$cluster
    
    # Compute the silhouette object for the optimal clustering (for the ordered silhouette plot)
    sil_obj <- silhouette(km_final$cluster, d)
    sil_df <- as.data.frame(sil_obj)
    sil_df <- sil_df %>% arrange(cluster, sil_width) %>%
      mutate(order = row_number())
    
    # Create the ordered silhouette plot using ggplot2
    sil_plot <- ggplot(sil_df, aes(x = order, y = sil_width, fill = factor(cluster))) +
      geom_bar(stat = "identity", width = 1) +
      labs(title = paste("Ordered Silhouette Plot for", habitat, "- k =", optimal_k),
           x = "Ordered Observations", y = "Silhouette Width", fill = "Cluster") +
      theme_minimal()+
      theme(
        axis.title.x = element_text(size = 15),  # Increase x-axis title size
        axis.title.y = element_text(size = 15),  # Increase y-axis title size
        plot.title = element_text(size = 25, hjust = 0.5)  # Increase title size and center it
      )
    
    # Save the silhouette plot for the current habitat
    silhouette_plots_list[[habitat]] <- sil_plot
  }
  
  # Compute cluster centers and a dispersion measure (std dev of distances from center)
  centers_df <- habitat_data %>%
    group_by(Cluster) %>%
    summarize(
      center_PC1 = mean(PC1),
      center_PC2 = mean(PC2),
      Count = n(),
      dispersion = sd(sqrt((PC1 - mean(PC1))^2 + (PC2 - mean(PC2))^2))
    ) %>%
    ungroup() %>%
    mutate(Original_Group = habitat)
  
  # Create a table listing species per cluster
  species_table <- habitat_data %>%
    group_by(Cluster) %>%
    summarize(
      Species = paste(Name, collapse = ", "),
      Count = n()
    ) %>%
    ungroup() %>%
    mutate(Original_Group = habitat)
  
  # Create a table with the breakdown of tag counts per cluster
  tags_table <- habitat_data %>%
    group_by(Cluster) %>%
    summarize(across(all_of(tag_cols), sum)) %>%
    ungroup() %>%
    mutate(Original_Group = habitat)
  
  # Generate a scatter plot for this habitat with clusters and centers
  p <- ggplot(habitat_data, aes(x = PC1, y = PC2, color = factor(Cluster))) +
    geom_point(size = 2, alpha = 0.7) +
    geom_point(data = centers_df, aes(x = center_PC1, y = center_PC2),
               color = "black", size = 4, shape = 8) +
    labs(title = paste("KMeans Clustering for", habitat, "- Optimal k =", optimal_k),
         color = "Cluster") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 20),  # Increase x-axis title size
      axis.title.y = element_text(size = 20),  # Increase y-axis title size
      plot.title = element_text(size = 30, hjust = 0.5)  # Increase title size and center it
    )
  
  # Save the results in the respective lists
  plots_list[[habitat]] <- p
  clustering_results[[habitat]] <- habitat_data
  species_table_list[[habitat]] <- species_table
  tags_table_list[[habitat]] <- tags_table
  centers_table_list[[habitat]] <- centers_df
}

# Combine results from all habitats into one data frame per result type
all_clusters <- do.call(rbind, clustering_results)
all_species_table <- do.call(rbind, species_table_list)
all_tags_table <- do.call(rbind, tags_table_list)
all_centers_table <- do.call(rbind, centers_table_list)

# Display the species per cluster table, the tag breakdown table, and the cluster centers
print(all_species_table)
print(all_tags_table)
print(all_centers_table)

# Display the plots for each habitat: scatter plot, elbow plot, and ordered silhouette plot (if available)
for (habitat in names(plots_list)) {
  print(plots_list[[habitat]])
  if (!is.null(elbow_plots_list[[habitat]])) {
    print(elbow_plots_list[[habitat]])
  }
  if (!is.null(silhouette_plots_list[[habitat]])) {
    print(silhouette_plots_list[[habitat]])
  }
}

# Print the elbow plots for the initial break
elbow_plots_list


#############################################################################
###################################Build patchwork for plots ################
#############################################################################

# Ensure all plots share the same x and y limits
for (i in seq_along(elbow_plots_list)) {
  elbow_plots_list[[i]] <- elbow_plots_list[[i]] +
    scale_x_continuous(limits = c(1, 10), breaks = 1:10) +
    scale_y_continuous(limits = c(0.2, 0.6)) +  # Adjust as needed
    theme_minimal(base_size = 10) +
    theme(
      axis.title = element_blank(),
      plot.title = element_text(size = 25, hjust = 0.5),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
}

# Determine which plots are on the left and bottom for labels
left_indices <- c(1, 4, 7)
bottom_indices <- c(6, 7, 8)

for (i in seq_along(elbow_plots_list)) {
  if (i %in% left_indices) {
    elbow_plots_list[[i]] <- elbow_plots_list[[i]] +
      theme(axis.text.y = element_text(size = 18, hjust = 0.5)) 
  }
  if (i %in% bottom_indices) {
    elbow_plots_list[[i]] <- elbow_plots_list[[i]] +
      theme(axis.text.x = element_text(size = 18, hjust = 0.5)) 
  }
}

left_tit <-  c(4)
bot_tit <- c(8)

for (i in seq_along(elbow_plots_list)) {
  if (i %in% left_tit) {
    elbow_plots_list[[i]] <- elbow_plots_list[[i]] +
      theme( axis.title.y = element_text(size = 40, hjust = 0.5)) +
      ylab("Average Silhouette Score")
  }
  if (i %in% bot_tit) {
    elbow_plots_list[[i]] <- elbow_plots_list[[i]] +
      theme(axis.title.x = element_text(size = 40, hjust = 0.5)) +
      xlab("Number of Clusters (k)")
  }
}

# Add an empty plot for the 9th space
empty_plot <- ggplot() + theme_void()

# Combine using patchwork
combined_plot <- wrap_plots(c(elbow_plots_list, list(empty_plot)),
                            nrow = 3, ncol = 3) +
  plot_annotation(
    title = "Average Silhouette Score Computed for k = {2,...,10} Clusters Over Each Guild",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold"))
  )

################################################################################
###################################### 3 by 3 plot for elbow plots #############
################################################################################
print(combined_plot)
################################################################################

print(all_species_table)

# Save this as a file
#write.table(all_species_table, "Cluster_Membership.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Table of clusters and species
print(all_tags_table)

# Table of centres
print(all_centers_table)

# Save these as a file if you want
#write.table(all_tags_table, "Tags_of_Breaks.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(all_centers_table, "Summary_of_Breaks.txt", sep = "\t", row.names = FALSE, quote = FALSE)




#####################################################################################
###################################Removing single clusters##########################
#####################################################################################

# Reshape the all_species_table so that each species is on its own row.
all_species_long <- all_species_table %>%
  separate_rows(Species, sep = ",") %>%
  mutate(Species = trimws(Species))

# Merge with cleaned_data by matching species names.
merged_data <- cleaned_data %>%
  inner_join(all_species_long, by = c("Name" = "Species"))

# Remove species in clusters that have only one species.
merged_data <- merged_data %>%
  group_by(Original_Group.y, Cluster) %>%
  mutate(cluster_size = n()) %>%
  ungroup()

# Identify species in singleton clusters
singleton_species <- merged_data %>%
  filter(cluster_size == 1) %>%
  pull(Name) %>%
  unique()

# Keep only non-singleton clusters
non_singletons <- merged_data %>%
  filter(cluster_size > 1)

###############################################################################
# Identify outlier species within each cluster using Mahalanobis distance.
# Mahalanobis requires the covariance matrix, and vector of means per cluster.
###############################################################################

# Compute Mahalanobis distance safely within clusters
mahal_cluster_outliers <- non_singletons %>%
  group_by(Original_Group.y, Cluster) %>%
  group_modify(~ {
    cluster_data <- .x
    mat <- cluster_data[, c("PC1", "PC2")]
    
    if (nrow(mat) < 3) {
      cluster_data <- cluster_data %>%
        mutate(
          mahal_dist = NA,
          dist_threshold = NA,
          outlier_detection_skipped = TRUE,
          skip_reason = "Too few points"
        )
      return(cluster_data)
    }
    
    cov_matrix <- cov(mat)
    if (det(cov_matrix) == 0) {
      cluster_data <- cluster_data %>%
        mutate(
          mahal_dist = NA,
          dist_threshold = NA,
          outlier_detection_skipped = TRUE,
          skip_reason = "Singular covariance"
        )
      return(cluster_data)
    }
    
    cluster_data <- cluster_data %>%
      mutate(
        mahal_dist = mahalanobis(mat, colMeans(mat), cov_matrix),
        dist_threshold = qchisq(0.95, df = 2),
        outlier_detection_skipped = FALSE,
        skip_reason = NA_character_
      )
    return(cluster_data)
  }) %>%
  ungroup()


skipped_clusters_summary <- mahal_cluster_outliers %>%
  filter(outlier_detection_skipped) %>%
  distinct(Original_Group.y, Cluster, skip_reason)

print(skipped_clusters_summary)

# Filter out outliers
outlier_species <- mahal_cluster_outliers %>%
  filter(!is.na(mahal_dist) & mahal_dist > dist_threshold) %>%
  pull(Name) %>%
  unique()

# Combine all species to be removed.
removed_species <- union(singleton_species, outlier_species)

# Create the new cleaned_data_no_outlier by removing these species.
cleaned_data_no_outlier <- cleaned_data %>%
  filter(!Name %in% removed_species)

# Optionally, save the list of removed species:
removed_species_list <- removed_species

# Display results of outlier species
print(cleaned_data_no_outlier)
print(removed_species_list)

# Plotting and output remain unchanged
all_species_long <- all_species_table %>%
  separate_rows(Species, sep = ",") %>%
  mutate(Species = trimws(Species))

joined_data <- cleaned_data %>%
  inner_join(all_species_long, by = c("Name" = "Species"))


###############################################################################
######### 2-D Projection coefficents after first split with outliers ##########
###############################################################################

ggplot(joined_data, aes(x = PC1, y = PC2, color = factor(Cluster))) +
  geom_point(size=2.5, alpha = 0.6) +
  geom_point(data = joined_data %>% filter(Name %in% removed_species_list),
             color = "red", shape = 17, size = 4, show.legend = FALSE) +
  geom_text(data = joined_data %>% filter(Name %in% removed_species_list),
            aes(label = Name), vjust = 1.05, hjust = 0.58, size = 4.75, , show.legend = FALSE) +
  facet_wrap(~ Original_Group.y) +
  theme_minimal() +
  labs(title = "Species Projection Coefficients by Guild Sub-Cluster for FPC1 and FPC2",
       x = "FPC1",
       y = "FPC2",
       color = "Sub-Cluster")+
  theme(
    axis.title.x = element_text(size = 30),  # Increase x-axis title size
    axis.title.y = element_text(size = 30),  # Increase y-axis title size
    plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),  # Increase title size and center it
    # Legend title
    legend.title = element_text(size = 30),
    legend.text = element_text(size =18),
    # Subplot (facet) titles
    legend.position = "none",  # This removes the legend
    strip.text = element_text(size = 25)) + 
    ylim(c(-3.5,3.5)) +
    xlim(c(-3.5,3.5)) +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0)


removed_data <- joined_data %>% filter(Name %in% removed_species_list)
removed_data %>% dplyr::select(Name, PC1, PC2)

############################################################################
############################ Plot the outliers #############################
############################################################################
# Create a data frame for plotting
plot_data <- tibble(
  Year = years,
  `Evening Grosbeak` = removed_data$PC1[1] * rot1 + removed_data$PC2[1] * rot2,
  `Rusty Blackbird`  = removed_data$PC1[2] * rot1 + removed_data$PC2[2] * rot2,
  `Trumpeter Swan`   = removed_data$PC1[3] * rot1 + removed_data$PC2[3] * rot2,
  `Wandering Tattler`= removed_data$PC1[4] * rot1 + removed_data$PC2[4] * rot2
) %>%
  pivot_longer(
    cols = -Year,
    names_to = "Species",
    values_to = "Abundance"
  )
# Plot
ggplot(plot_data, aes(x = Year, y = Abundance, color = Species)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~ Species, ncol = 2, scales = "fixed") +
  scale_color_manual(values = c(
    "Evening Grosbeak"   = "steelblue",
    "Rusty Blackbird"    = "darkgreen",
    "Trumpeter Swan"     = "firebrick",
    "Wandering Tattler"  = "purple"
  )) +
  scale_y_continuous(limits = c(-0.75, 0.45)) +
  labs(
    x = "Year",
    y = "Reconstructed Abundance Function"
  ) +
  theme_minimal(base_size = 22) +  # makes everything bigger by default
  theme(
    legend.position = "none",
    strip.text = element_text(size = 26, face = "bold"),   # facet titles
    axis.title = element_text(size = 24, face = "bold"),   # axis labels
    axis.text = element_text(size = 20),                   # tick labels
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5)
  )

##########################################################################################
##########################################################################################
#Second Run of Analysis (after Identification of Outliers)################################
##########################################################################################
##########################################################################################
cleaned_data_no_outlier
clustering_results <- list()
plots_list <- list()
species_table_list <- list()
tags_table_list <- list()
centers_table_list <- list()
silhouette_plots_list <- list()
elbow_plots_list <- list()

# Loop over each habitat (Original_Group)
for (habitat in Habitiat_category_columns) {
  # Subset data for the current habitat
  habitat_data <- cleaned_data_no_outlier %>% filter(Original_Group == habitat)
  
  # Convert tag columns to numeric
  habitat_data <- habitat_data %>%
    mutate(across(all_of(tag_cols), ~ as.numeric(as.character(.))))
  
  
  #habitat_data <- habitat_data %>% filter(!Outlier)  # Remove KMM outliers
  #habitat_data <- habitat_data %>% filter(!Outlier_LOF)  # Remove LOF outliers 
  
  
  # Initialize a data frame to store average silhouette scores for each candidate k
  sil_summary <- data.frame(k = integer(), avg_sil = double())
  
  # If there is only one species in this habitat, assign it to cluster 1.
  if(nrow(habitat_data) < 2) {
    habitat_data$Cluster <- 1
    optimal_k <- 1
  } else {
    # Consider k from 2 to up to 10 (or less if not enough rows)
    k_range <- 2:min(10, nrow(habitat_data))
    silhouette_avgs <- numeric(length(k_range))
    
    # Pre-compute the distance matrix on PC1 and PC2
    d <- dist(habitat_data[, c("PC1", "PC2")])
    
    # Loop over each candidate k to calculate average silhouette width
    for (i in seq_along(k_range)) {
      k <- k_range[i]
      set.seed(123)  # for reproducibility
      km_res <- kmeans(habitat_data[, c("PC1", "PC2")], centers = k, nstart = 25)
      sil <- silhouette(km_res$cluster, d)
      silhouette_avgs[i] <- mean(sil[, 3])
      sil_summary <- rbind(sil_summary, data.frame(k = k, avg_sil = mean(sil[, 3])))
    }
    
    # Create an elbow plot: average silhouette score vs. k
    elbow_plot <- ggplot(sil_summary, aes(x = k, y = avg_sil)) +
      geom_line() +
      geom_point(size = 3) +
      labs(title = paste("Elbow Plot (Silhouette Score) for", habitat),
           x = "Number of Clusters (k)",
           y = "Average Silhouette Score") +
      theme_minimal()
    elbow_plots_list[[habitat]] <- elbow_plot
    
    # Select the k with the maximum average silhouette width
    optimal_index <- which.max(silhouette_avgs)
    optimal_k <- k_range[optimal_index]
    
    # Run final KMeans clustering using the optimal k
    set.seed(123)
    km_final <- kmeans(habitat_data[, c("PC1", "PC2")], centers = optimal_k, nstart = 25)
    habitat_data$Cluster <- km_final$cluster
    
    # Compute the silhouette object for the optimal clustering (for the ordered silhouette plot)
    sil_obj <- silhouette(km_final$cluster, d)
    sil_df <- as.data.frame(sil_obj)
    sil_df <- sil_df %>% arrange(cluster, sil_width) %>%
      mutate(order = row_number())
    
    # Create the ordered silhouette plot using ggplot2
    sil_plot <- ggplot(sil_df, aes(x = order, y = sil_width, fill = factor(cluster))) +
      geom_bar(stat = "identity", width = 1) +
      labs(title = paste("Ordered Silhouette Plot for", habitat, "- k =", optimal_k),
           x = "Ordered Observations", y = "Silhouette Width", fill = "Cluster") +
      theme_minimal()
    
    # Save the silhouette plot for the current habitat
    silhouette_plots_list[[habitat]] <- sil_plot
  }
  
  # Compute cluster centers and a dispersion measure (std dev of distances from center)
  centers_df <- habitat_data %>%
    group_by(Cluster) %>%
    summarize(
      center_PC1 = mean(PC1),
      center_PC2 = mean(PC2),
      Count = n(),
      dispersion = sd(sqrt((PC1 - mean(PC1))^2 + (PC2 - mean(PC2))^2))
    ) %>%
    ungroup() %>%
    mutate(Original_Group = habitat)
  
  # Create a table listing species per cluster
  species_table <- habitat_data %>%
    group_by(Cluster) %>%
    summarize(
      Species = paste(Name, collapse = ", "),
      Count = n()
    ) %>%
    ungroup() %>%
    mutate(Original_Group = habitat)
  
  # Create a table with the breakdown of tag counts per cluster
  tags_table <- habitat_data %>%
    group_by(Cluster) %>%
    summarize(across(all_of(tag_cols), sum)) %>%
    ungroup() %>%
    mutate(Original_Group = habitat)
  
  # Generate a scatter plot for this habitat with clusters and centers
  p <- ggplot(habitat_data, aes(x = PC1, y = PC2, color = factor(Cluster))) +
    geom_point(size = 2, alpha = 0.7) +
    geom_point(data = centers_df, aes(x = center_PC1, y = center_PC2),
               color = "black", size = 4, shape = 8) +
    labs(title = paste("KMeans Clustering for", habitat, "- Optimal k =", optimal_k),
         color = "Cluster") +
    theme_minimal()
  
  # Save the results in the respective lists
  plots_list[[habitat]] <- p
  clustering_results[[habitat]] <- habitat_data
  species_table_list[[habitat]] <- species_table
  tags_table_list[[habitat]] <- tags_table
  centers_table_list[[habitat]] <- centers_df
}

# Combine results from all habitats into one data frame per result type
all_clusters <- do.call(rbind, clustering_results)
all_species_table <- do.call(rbind, species_table_list)
all_tags_table <- do.call(rbind, tags_table_list)
all_centers_table <- do.call(rbind, centers_table_list)

# Display the species per cluster table, the tag breakdown table, and the cluster centers
print(all_species_table)
print(unique(all_species_table$Original_Group))
print(all_tags_table)
print(all_centers_table)

# Display the plots for each habitat: scatter plot, elbow plot, and ordered silhouette plot (if available)
for (habitat in names(plots_list)) {
  print(plots_list[[habitat]])
  if (!is.null(elbow_plots_list[[habitat]])) {
    print(elbow_plots_list[[habitat]])
  }
  if (!is.null(silhouette_plots_list[[habitat]])) {
    print(silhouette_plots_list[[habitat]])
  }
}



print(all_species_table)
print(all_tags_table)
print(all_centers_table)

# write.csv(all_tags_table, "all_tags_table.csv", row.names = FALSE, fileEncoding = "UTF-8")
# write.table(all_centers_table, "Summary_of_Breaks.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Reorder the columns
export_table <- all_species_table %>%
  dplyr::select(Cluster, Original_Group, Count, Species)

# Save to CSV
#write.csv(export_table, "all_species_table.csv", row.names = FALSE, fileEncoding = "UTF-8")

######################################################################################################################
######################################################################################################################

######################################################################################
#############################Rejoining of the Clusters################################
######################################################################################

# --- PREPARATION: Build the initial clusters as leaves ---
all_clusters

# Add a unique label for each initial cluster (e.g., "Forest Birds 1")
all_clusters$Subcluster <- paste(all_clusters$Original_Group, all_clusters$Cluster, sep = " ")

unique(all_clusters$Subcluster)

initial_clusters <- split(all_clusters, all_clusters$Subcluster)
initial_clusters

tag_cols <- c("Duck Species", "Goose Species",
              "Forest Crop Specilists", "Forest birds wintering in SA",
              "Forest birds wintering in CA", "Short distance shorebirds",
              "Long distance shorebirds", "Tolerant to Farming", "Dependent on Grassland")

clusters_list <- list()
for (cl in names(initial_clusters)) {
  df <- initial_clusters[[cl]]
  clusters_list[[cl]] <- list(
    species = df,
    lineage = cl,
    tag_counts = colSums(df[, tag_cols, drop = FALSE])
  )
}

length(clusters_list)

m <- length(clusters_list)
cluster_indices <- setNames(-seq_len(m), names(clusters_list))
next_index <- 1

merge_matrix <- NULL
merge_heights <- c()
merge_steps <- list()
scatter_plots <- list()
PC_points<-list()
tag_counts_history <- list()  # to track tag counts per iteration

iteration <- 0


# Main hierarchical merging loop
while (length(clusters_list) > 1) {
  iteration <- iteration + 1
  
  # Extract names and prepare an empty distance matrix
  cluster_names <- names(clusters_list)
  n <- length(cluster_names)
  dist_mat <- matrix(Inf, n, n, dimnames = list(cluster_names, cluster_names))
  
  # Compute pairwise distances between clusters
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      sp1 <- clusters_list[[ cluster_names[i] ]]$species[, c("PC1", "PC2"), drop = FALSE]
      sp2 <- clusters_list[[ cluster_names[j] ]]$species[, c("PC1", "PC2"), drop = FALSE]
      
      if (nrow(sp1) > 0 & nrow(sp2) > 0) {
        # Max pairwise distance between clusters (complete linkage)
        d_ij <- max(as.matrix(dist(rbind(sp1, sp2)))[1:nrow(sp1), (nrow(sp1) + 1):nrow(rbind(sp1, sp2))])
        dist_mat[i, j] <- d_ij
        dist_mat[j, i] <- d_ij
      }
    }
  }
  
  # Identify closest pair of clusters to merge
  min_idx <- which(dist_mat == min(dist_mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
  name1 <- cluster_names[min_idx[1]]
  name2 <- cluster_names[min_idx[2]]
  merge_distance <- dist_mat[name1, name2]
  
  new_species <- rbind(clusters_list[[name1]]$species, clusters_list[[name2]]$species)
  new_lineage <- paste(clusters_list[[name1]]$lineage, clusters_list[[name2]]$lineage, sep = " + ")
  new_tag_counts <- colSums(new_species[, tag_cols, drop = FALSE])
  
  new_cluster_name <- new_lineage
  new_cluster <- list(
    species = new_species,
    lineage = new_lineage,
    tag_counts = new_tag_counts
  )
  
  tag_counts_history[[iteration]] <- data.frame(
    Iteration = iteration,
    Cluster = new_cluster_name,
    t(as.data.frame(new_tag_counts))  # transpose to make tags as columns
  )
  
  
  # Register merge indices for dendrogram reconstruction
  idx1 <- cluster_indices[name1]
  idx2 <- cluster_indices[name2]
  pair_indices <- c(idx1, idx2)
  new_index <- next_index
  next_index <- next_index + 1
  
  merge_matrix <- rbind(merge_matrix, c(pair_indices[1], pair_indices[2]))
  merge_heights <- c(merge_heights, merge_distance)
  
  merge_steps[[iteration]] <- list(
    merged = c(name1, name2),
    new_lineage = new_lineage,
    distance = merge_distance,
    tag_counts = new_tag_counts
  )
  
  # Remove old clusters and add the new merged one
  clusters_list[[name1]] <- NULL
  clusters_list[[name2]] <- NULL
  clusters_list[[new_cluster_name]] <- new_cluster
  
  cluster_indices <- cluster_indices[!names(cluster_indices) %in% c(name1, name2)]
  cluster_indices[new_cluster_name] <- new_index
  
  
  
  current_species <- do.call(rbind, lapply(clusters_list, function(cl) cl$species))
  current_species$CurrentCluster <- rep(names(clusters_list), sapply(clusters_list, function(cl) nrow(cl$species)))
  
  PC_points[[iteration]]<-current_species
  
  
  centroids <- do.call(rbind, lapply(names(clusters_list), function(nm) {
    cl <- clusters_list[[nm]]
    data.frame(
      PC1 = mean(cl$species$PC1),
      PC2 = mean(cl$species$PC2),
      CurrentCluster = nm
    )
  }))
  
  
  
  # Scatter plot with correct coloring
  p_scatter <- ggplot() +
    # Plot individual points
    geom_point(data = current_species, aes(x = PC1, y = PC2, color = factor(CurrentCluster)), size = 2, alpha = 0.7) +
    # Plot centroids as squares in the same color
    geom_point(data = centroids, aes(x = PC1, y = PC2, color = factor(CurrentCluster)), shape = 15, size = 5) +
    # Ensure consistent colors
    
    scale_color_manual(values = rainbow(length(unique(current_species$CurrentCluster)))) +
    
    # Wrap legend labels to prevent long labels from cutting off
    guides(color = guide_legend(override.aes = list(size = 4), title = "Cluster (Lineage)")) +
    theme_minimal() +
    labs(
      title = paste("Merge", iteration),
      subtitle = paste("Merge distance:", round(merge_distance, 2)),
      color = "Cluster (Lineage)",
    ) +
    theme(
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      axis.title.x = element_text(size=18, face="bold"),
      axis.title.y = element_text(size=18, face="bold"),
      legend.box = "vertical"     # spacing between items
      
      
      
    )+
    scale_color_discrete(labels = function(x) str_wrap(x, width = 20))
  # Wrap legend labels at ~20 characters
  
  scatter_plots[[iteration]] <- p_scatter
  
}

n_orig <- length(initial_clusters)
hclust_obj <- list(
  merge = merge_matrix,
  height = merge_heights,
  order = 1:n_orig,
  labels = names(initial_clusters),
  method = "complete-linkage",
  call = match.call()
)
class(hclust_obj) <- "hclust"

dendro_data <- dendro_data(hclust_obj)

merge_yvals <- unique(dendro_data$segments$yend)
merge_yvals <- merge_yvals[merge_yvals > 0]  # Filter out tips (0 height)

# Sort merge heights from lowest to highest
merge_order <- sort(merge_yvals)

# Match each height with an order number
merge_labels_df <- data.frame(
  y = merge_order,
  label = seq_along(merge_order)
)

# For x positions, estimate midpoint of each merge from the dendro segments
midpoints <- dendro_data$segments %>%
  filter(y == yend & y > 0) %>% # horizontal lines at merge points
  group_by(y) %>%
  summarize(x = mean(c(x, xend)), .groups = "drop")  # midpoint x of each horizontal merge line

# Merge with label data
merge_labels_df <- left_join(merge_labels_df, midpoints, by = "y")

merge_labels_df

label_offset <- max(dendro_data$segments$yend) * 0.11

dendro_plot <- ggplot(segment(dendro_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(dendro_data),
            aes(x = x, y = -label_offset, label = label),
            hjust = -0, vjust=-0.8 , angle = 90, size = 4.5) +
  labs(title = "Dendrogram of Merged Clusters", x = "Original Cluster", y = "Height") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(9, 40, 10, 10),
  ) +
  geom_text(data = merge_labels_df,
            aes(x = x, y = y - (label_offset*0.2), label = label),
            color = "blue", size = 4, fontface = "bold")+
  ylim(-label_offset * 1.05, max(dendro_data$segments$yend))+
geom_hline(yintercept = 2.555, colour = "blue", linetype = "dashed", linewidth = 1) +
geom_hline(yintercept = 3.14, colour = "orange", linetype = "dashed", linewidth = 1)


print(dendro_plot)
tag_counts_history <- do.call(rbind, tag_counts_history)

#Write to file if you want
#write.csv(tag_counts_history, "Merge_Tags.csv", row.names = FALSE, fileEncoding = "UTF-8")

###################################################################################
#########################################Species Level Plots#######################
###################################################################################
##

cluster_membership_to_be_used_later<-initial_clusters

# Define the years vector (1975 to 2016)
years <- 1975:2016  # 42 years

# List to store the initial time series plots for each leaf cluster.
initial_time_series_plots <- list()

# Assume "initial_clusters" is a list where each element is a data frame
# of species for that initial cluster (as defined earlier using your Subcluster labels).
for (nm in names(initial_clusters)) {
  df <- initial_clusters[[nm]]
  
  # For each species in the cluster, compute its time series:
  # Use the PC1 and PC2 values from the species, and compute: PC1*rot1 + PC2*rot2.
  # Compute the time series for each species in the cluster.
  # Each row of species_ts will be a species' time series (length 42).
  species_ts <- t(apply(df[, c("PC1", "PC2")], 1, function(x) {
   ( x[1] * rot1 + x[2] * rot2 )
  }))
  
  # Set the column names to the years (1975-2016)
  colnames(species_ts) <- as.character(1975:2016)
  
  # Convert the matrix to a data frame.
  species_ts_df <- as.data.frame(species_ts)
  # Add a column for species names (assuming rownames are set or you can assign a unique identifier)
  species_ts_df$Species <- rownames(species_ts_df)
  
  # Melt the data frame so each row is a time point for a species.
  species_ts_long <- melt(species_ts_df, id.vars = "Species", variable.name = "year", value.name = "Abundance")
  # Convert year to numeric if needed:
  species_ts_long$year <- as.numeric(as.character(species_ts_long$year))
  
  # Compute the cluster centre time series (using the mean PC1 and PC2 of the cluster)
  centre <- c(mean(df$PC1), mean(df$PC2))
  centre_ts <- ( centre[1] * rot1 + centre[2] * rot2)
  centre_df <- data.frame(year = years, Abundance = centre_ts)
  
  # Create the plot: individual species in faint orange, centre as thick dark-blue line.
  p <- ggplot() +
    geom_line(data = species_ts_long, aes(x = year, y = Abundance, group = Species),
              color = "orange", size = 0.4) +
    ylim(-0.75, 0.45)+
    geom_line(data = centre_df, aes(x = year, y = Abundance),
              color = "darkblue", size = 1.5) +
    labs(title = paste(nm),
         x = "Year", y = "Reconstructed Abundance Function") +
    scale_x_continuous(breaks = seq(1975, 2016, by = 10)) +
    theme_minimal()+
    theme(
      plot.title = element_text(hjust = 0.5, size=18),
      axis.title.y = element_text(size=25 ),
      axis.text = element_text(size = 11, face="bold"),
      axis.title.x = element_text(size=30 ),)
  
  initial_time_series_plots[[nm]] <- p
}

# For example, view the plot for one of the initial clusters:


group_colors <- c("Wetland Birds" = "#1f77b4",
                  "Shorebirds" = "#ff7f0e",
                  "Birds of Prey" = "#2ca02c",
                  "Seabirds" = "#d62728",
                  "Waterfowl" = "#9467bd",
                  "Forest Birds" = "#8c564b",
                  "Grassland Birds" = "#e377c2",
                  "Aerial Insectivores" = "#7f7f7f")

# Print individual layers if you want
# print(initial_time_series_plots[[1]])


# Blank plot for padding
blank_plot <- ggplot() + theme_void()

# Helper function to prep each group
prepare_group <- function(plots, total = 4, show_x = FALSE, show_y_first_only = TRUE) {
  n <- length(plots)
  # Pad with blank plots
  if (n < total) {
    plots <- c(plots, rep(list(blank_plot), total - n))
  }
  
  # Style each plot
  plots <- lapply(seq_along(plots), function(i) {
    p <- plots[[i]] +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank()) 
          
    if (!inherits(p, "gg")) return(p)  # Skip blank
    if (!show_x) {
      p <- p + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
    }
    if (show_y_first_only && i != 1) {
      p <- p + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
    }
    return(p)
  })
  
  return(wrap_plots(plots, ncol = total))
}

# Define your groups (replace `initial_time_series_plots` with your actual list)
group1 <- prepare_group(initial_time_series_plots[1:3])
group2 <- prepare_group(initial_time_series_plots[4:5])
group3 <- prepare_group(initial_time_series_plots[6:8])
group4 <- prepare_group(initial_time_series_plots[9:10])
group5 <- prepare_group(initial_time_series_plots[11:12])
group6 <- prepare_group(initial_time_series_plots[13:15])
group8 <- prepare_group(initial_time_series_plots[16:19], show_x = TRUE)
group7 <- prepare_group(initial_time_series_plots[20:23])

# Combine all groups into one plot (without shared axis titles yet)
final_plot <- group1 / group2 / group3 / group4 / group5 / group6 / group7 / group8 +
  plot_layout(ncol = 1)

final_plot <- final_plot +
  theme(
      plot.margin = margin(5, 5, 5, 5),  # Add space for axis titles
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
  )

# --- SHARED AXIS TITLES SECTION ---
# Convert to grob for layout manipulation
g <- patchwork::patchworkGrob(final_plot)

# Create shared axis titles
x_label <- textGrob("Years", gp = gpar(fontsize = 25, fontface = "bold"))
y_label <- textGrob("Reconstructed Centered Trajectory Function", rot = 90, gp = gpar(fontsize = 25, fontface = "bold"))


main_title <- textGrob(
  "Reconstructed Centered Trajectory Function for Each Guild Sub-Cluster",
  gp = gpar(fontsize = 25, fontface = "bold")
)

# Arrange the full layout
grid.arrange(
  main_title,
  arrangeGrob(
    arrangeGrob(
      y_label, g, ncol = 2, widths = unit.c(unit(1.5, "lines"), unit(1, "npc") - unit(1.5, "lines"))
    ),
    x_label,
    nrow = 2,
    heights = unit.c(unit(1, "npc") - unit(1.5, "lines"), unit(1.5, "lines"))
  ),
  nrow = 2,
  heights = unit.c(unit(1.5, "lines"), unit(1, "npc") - unit(1.5, "lines"))
)

######################################################
######################################################
#############Save this plot###########################
######################################################

output_dir <- "plots"
# Create the folder if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir)

for (i in seq_along(initial_time_series_plots)) {
  plot_name <- paste0(output_dir, "/time_series_plot_", i, ".png")
  p <- initial_time_series_plots[[i]] + 
    theme_minimal() + 
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(filename = plot_name, plot = p, width = 8, height = 6, dpi = 300)
}



PC_points[[18]]$CurrentCluster
PC_points[[18]]$Subcluster
PC_points[[18]]$Original_Group
PC_points[[18]]$Name

PC_points





main_time_ser <- function(Main_data, Iter) {
  Iter_data <- list()
  
  last_cluster <- Main_data[[Iter]]$CurrentCluster[nrow(Main_data[[Iter]])]
  
  # Now, filter the data based on the last value of CurrentCluster
  Main_spec <- Main_data[[Iter]][Main_data[[Iter]]$CurrentCluster == last_cluster, ]
  
  print(Main_spec)
  
  for (i in 1:nrow(Main_spec)) {
    vecz <- NULL
    vecz <- Main_spec[i,]$PC1 * rot1 + Main_spec[i,]$PC2 * rot2
    Iter_data[[i]] <- list(Values = vecz, Original_Group = Main_spec[i,]$Original_Group)
  }
  return(Iter_data)
}

N <- 22  # You can change this as needed
plotlist <- list()

for (I in 1:N) {
  
  # Dynamically generate Iter_1_dat for each I
  Iter_1_dat <- main_time_ser(PC_points, I)
  
  # Create the dataframe iter_df_I, ensuring we include Original_Group
  iter_df_I <- do.call(rbind, lapply(1:length(Iter_1_dat), function(i) {
    data.frame(Year = years, Value = Iter_1_dat[[i]]$Values, 
               Group = as.factor(i), 
               Original_Group = as.factor(Iter_1_dat[[i]]$Original_Group))  # Ensure Original_Group is retained
  }))
  
  # Calculate the mean time series
  mean_df <- iter_df_I %>%
    group_by(Year) %>%
    summarize(Mean_Value = mean(Value))
  
  # Create the plot for this iteration I
  plot_title <- paste("Time Series Plots of Cluster Formed After Merge", I)
  
  Plot12 <- ggplot() +
    # Plot individual time series, now colored by Original_Group
    geom_line(data = iter_df_I, aes(x = Year, y = Value, group = Group, color = Original_Group), alpha = 0.6) +
    # Plot the mean time series with blue thick lines
    geom_line(data = mean_df, aes(x = Year, y = Mean_Value), color = "darkblue", size = 1.5) + 
    labs(title = plot_title, 
         x = "Year", y = "Rescaled Abundance Index", color = "Original Group") +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size=25, face="bold"),
          axis.title.y = element_text(size=20,face = "bold"),
          axis.title.x = element_text(size=20,face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
    )
  
  plotlist[[I]] <- Plot12
}

plotlist
#Same Colour Scales

# Define fixed levels and colors
group_levels <- c("Wetland Birds", "Shorebirds", "Birds of Prey", "Seabirds", 
                  "Waterfowl", "Forest Birds", "Grassland Birds", "Aerial Insectivores")

group_colors <- c("Wetland Birds" = "#1f77b4",
                  "Shorebirds" = "#ff7f0e",
                  "Birds of Prey" = "#2ca02c",
                  "Seabirds" = "#d62728",
                  "Waterfowl" = "#9467bd",
                  "Forest Birds" = "#8c564b",
                  "Grassland Birds" = "#e377c2",
                  "Aerial Insectivores" = "#7f7f7f")

N <- 22
plotlist <- list()

for (I in 1:N) {
  
  Iter_1_dat <- main_time_ser(PC_points, I)
  
  iter_df_I <- do.call(rbind, lapply(1:length(Iter_1_dat), function(i) {
    data.frame(Year = years,
               Value = Iter_1_dat[[i]]$Values,
               Group = as.factor(i),
               Original_Group = factor(Iter_1_dat[[i]]$Original_Group, levels = group_levels))
  }))
  
  mean_df <- iter_df_I %>%
    group_by(Year) %>%
    summarize(Mean_Value = mean(Value), .groups = "drop")
  
  plot_title <- paste("Cluster at Merge", I)
  
  Plot12 <- ggplot() +
    geom_line(data = iter_df_I, aes(x = Year, y = Value, group = Group, color = Original_Group), alpha = 0.6) +
    geom_line(data = mean_df, aes(x = Year, y = Mean_Value), color = "darkblue", size = 1.8) + 
    scale_color_manual(values = group_colors, drop = FALSE) +  # <- Lock color mapping
    labs(title = plot_title, 
         x = "Year", y = "Rescaled Abundance Function", color = "Original Group") +
    ylim(-0.75,0.45)+
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.title.y = element_text(size = 30, face = "bold"),
          axis.title.x = element_text(size = 30, face = "bold"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22))
  
  plotlist[[I]] <- Plot12
}



#########################################################################
########################### Recover the three and eight cluster solution#
#########################################################################

#Add franklin gull to the plots too

Franklin<-cleaned_data[cleaned_data$Name=="Franklin's Gull",]

pc1<-Franklin$PC1
pc2<-Franklin$PC2

Rot_fra<- (pc1 * rot1 + pc2 * rot2)

plotlist[16]
plot(y=Rot_fra, x=years)

new_line_df <- data.frame(years = years, Rot_fra = Rot_fra)


p1 <- plotlist[[16]] +
  geom_line(data = new_line_df, aes(x = years, y = Rot_fra), 
            color = "red", size = 1, linetype = "dashed", show.legend = FALSE)+
  labs(title = "Red Cluster (Merge 16)", 
       x = "Year", y = "Reconstructed Centered Trajectory Function", color = "Original Group") + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 20))

p2 <- plotlist[[19]] +
  labs(title = "Blue Cluster (Merge 19)", 
       x = "Year", y = "Rescaled Abundance Function", color = "Original Group") + 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text = element_text(size = 20))

p3 <- plotlist[[20]] +
  labs(title = "Green Cluster (Merge 20)", 
       x = "Year", y = "Rescaled Abundance Function", color = "Original Group") + 
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text = element_text(size = 20))

# Combine plots horizontally
final_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)

# Show the final plot
print(final_plot)


#Save the plot if you want

# Set output directory
output_dir <- "merge_plots"
# Create the folder if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir)


for (i in seq_along(plotlist)) {
  filename <- paste0(output_dir, "/merge_", i, ".png")
  plot_clean <- plotlist[[i]] + 
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(filename = filename,
         plot = plot_clean,
         width = 8, height = 6, dpi = 300)
}




# Sanity check to see if reuslts look okay
# rub<-All_SOCB_2019[All_SOCB_2019$English_Name=="Rusty Blackbird",]
# rub$index
# plot(x=1970:2016,y=rub$index, type="l", main="Rusty Blackbird Raw Index")
# write.csv(rub, "Rusty_Blackbird.csv", row.names = FALSE, fileEncoding = "UTF-8")


#Colour by Groups

years <- 1975:2016  # 42 years

# List to store the initial time series plots for each leaf cluster.
initial_time_series_plots <- list()

# Assume "initial_clusters" is a list where each element is a data frame
# of species for that initial cluster (as defined earlier using your Subcluster labels).
for (nm in names(initial_clusters)) {
  df <- initial_clusters[[nm]]
  
  # For each species in the cluster, compute its time series:
  # Use the PC1 and PC2 values from the species, and compute: PC1*rot1 + PC2*rot2.

  species_ts <- t(apply(df[, c("PC1", "PC2")], 1, function(x) {
    x[1] * rot1 + x[2] * rot2
  }))
  
  # Set the column names to the years (1975-2016)
  colnames(species_ts) <- as.character(1975:2016)
  
  # Convert the matrix to a data frame.
  species_ts_df <- as.data.frame(species_ts)
  # Add a column for species names (assuming rownames are set or you can assign a unique identifier)
  species_ts_df$Species <- rownames(species_ts_df)
  
  # Melt the data frame so each row is a time point for a species.
  species_ts_long <- melt(species_ts_df, id.vars = "Species", variable.name = "year", value.name = "Abundance")
  # Convert year to numeric if needed:
  species_ts_long$year <- as.numeric(as.character(species_ts_long$year))
  
  # Compute the cluster centre time series (using the mean PC1 and PC2 of the cluster)
  centre <- c(mean(df$PC1), mean(df$PC2))
  centre_ts <- centre[1] * rot1 + centre[2] * rot2
  centre_df <- data.frame(year = years, Abundance = centre_ts)
  
  # Create the plot: individual species in faint orange, centre as thick dark-blue line.
  p <- ggplot() +
    geom_line(data = species_ts_long, aes(x = year, y = Abundance, group = Species),
              color = "#8c564b", size = 0.4) +
    ylim(-0.75, 0.45)+
    geom_line(data = centre_df, aes(x = year, y = Abundance),
              color = "darkblue", size = 1.5) +
    labs(title = paste(nm),
         x = "Year", y = "Reconstructed Centered Trajectory Function") +
    scale_x_continuous(breaks = seq(1975, 2016, by = 10)) +
    theme_minimal()+
    theme(
      plot.title = element_text(hjust = 0.5, size=25),
      axis.title.y = element_text(size=25 ),
      axis.text = element_text(size = 11, face="bold"),
      axis.title.x = element_text(size=30 ),)
  
  initial_time_series_plots[[nm]] <- p
}

# For example, view the plot for one of the initial clusters:


group_colors <- c("Wetland Birds" = "#1f77b4",
                  "Shorebirds" = "#ff7f0e",
                  "Birds of Prey" = "#2ca02c",
                  "Seabirds" = "#d62728",
                  "Waterfowl" = "#9467bd",
                  "Forest Birds" = "#8c564b",
                  "Grassland Birds" = "#e377c2",
                  "Aerial Insectivores" = "#7f7f7f")


Wetland_3<-initial_time_series_plots[[22]]

Seabirds_2<-print(initial_time_series_plots[[12]])

Forest_Birds_2<-initial_time_series_plots[[7]]

plotlist[[9]]=plotlist[[9]]+theme(legend.position = "none", 
                                  plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                                  axis.title.y = element_text(size = 30, face = "bold"),
                                  axis.title.x = element_text(size = 30, face = "bold"),
                                  legend.text = element_text(size = 22),
                                  axis.text = element_text(size = 15, face ="bold"),
                                  legend.title = element_text(size = 22))

plotlist[[11]]=plotlist[[11]]+theme(legend.position = "none",
                                    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                                    axis.title.y = element_text(size = 30, face = "bold"),
                                    axis.title.x = element_text(size = 30, face = "bold"),
                                    legend.text = element_text(size = 22),
                                    axis.text = element_text(size = 15, face ="bold"),
                                    legend.title = element_text(size = 22))

plotlist[[10]]=plotlist[[10]]+theme(legend.position = "none",
                                    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                                    axis.title.y = element_text(size = 30, face = "bold"),
                                    axis.title.x = element_text(size = 30, face = "bold"),
                                    legend.text = element_text(size = 22),
                                    axis.text = element_text(size = 15, face ="bold"),
                                    legend.title = element_text(size = 22))

plotlist[[14]]=plotlist[[14]]+labs( x = "Year", y = "Reconstructed Centered Trajectory Function")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.title.y = element_text(size = 30, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 22),
        axis.text = element_text(size = 15, face ="bold"),
        legend.title = element_text(size = 22))

plotlist[[15]]=plotlist[[15]]+theme(legend.position = "none",
                                    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                                    axis.title.y = element_text(size = 30, face = "bold"),
                                    axis.title.x = element_text(size = 30, face = "bold"),
                                    legend.text = element_text(size = 22),
                                    axis.text = element_text(size = 15, face ="bold"),
                                    legend.title = element_text(size = 22))

Wetland_3=Wetland_3+theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                          axis.title.y = element_text(size = 30, face = "bold"),
                          axis.title.x = element_text(size = 30, face = "bold"),
                          legend.text = element_text(size = 22),
                          axis.text = element_text(size = 15, face ="bold"),
                          legend.title = element_text(size = 22))

Seabirds_2=Seabirds_2+theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                            axis.title.y = element_text(size = 30, face = "bold"),
                            axis.title.x = element_text(size = 30, face = "bold"),
                            legend.text = element_text(size = 22),
                            axis.text = element_text(size = 15, face ="bold"),
                            legend.title = element_text(size = 22))

Forest_Birds_2=Forest_Birds_2+theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                                    axis.title.y = element_text(size = 30, face = "bold"),
                                    axis.title.x = element_text(size = 30, face = "bold"),
                                    legend.text = element_text(size = 22),
                                    axis.text = element_text(size = 15, face ="bold"),
                                    legend.title = element_text(size = 22))


# Prepare the plots
p1 <- plotlist[[9]]+theme(legend.position = "none",
                          axis.title.y = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x = element_blank() )
p1

p2 <- plotlist[[11]]+theme(legend.position = "none",
                           axis.title.y = element_blank(),
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank())
p2

p3 <- plotlist[[10]]+theme(legend.position = "none",
                           axis.title.y = element_blank(),
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank())
p3

p4 <- plotlist[[14]]+theme(legend.position = "none",
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank())


p5 <- plotlist[[15]]+theme(legend.position = "none",
                           axis.title.y = element_blank(),
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank())

p6 <- Wetland_3 +theme(legend.position = "none",
                       axis.title.y = element_blank(),
                       axis.title.x = element_blank(),
                       axis.text.y = element_blank() )


p7 <- Seabirds_2 + theme(legend.position = "none",
                         axis.title.y = element_blank(),
                         axis.title.x = element_blank()
)


p8 <- Forest_Birds_2 + theme(legend.position = "none",
                             axis.title.y = element_blank(),
                             axis.text.y = element_blank() )
p9 <- plot_spacer()     # Empty placeholder for the 9th slot


# Arrange in 3x3 grid
final_grid <- (p1 | p2 | p3) /
  (p4 | p5 | p6) /
  (p7 | p8 | p9)

# Display it
print(final_grid)




###################################################################
###################################################################
#################### Extra plot at Simon's request ################
###################################################################
###################################################################


#Need

# Raw rescaled indices

Raw_rescaled_for_later

Raw_rescaled_for_later$year <- years

# categories of clustering

cluster_membership_to_be_used_later

merged_df <- do.call(rbind, cluster_membership_to_be_used_later)

filtered_cluster_categories<-merged_df %>%
  dplyr::select("Subcluster","Original_Group","Name")

# Merge data

# Make sure your species columns include non-syntactic names (spaces, etc.) — pivot_longer handles them well.
long_joined_df <- 
  tidyr::pivot_longer(
    data   = Raw_rescaled_for_later,
    cols   = -c(year, SampleID),           # all species columns
    names_to  = "Species",
    values_to = "Value"
  ) |>
  dplyr::mutate(
    # optional hygiene: trim spaces just in case
    Species = stringr::str_trim(Species)
  ) |>
  # 2) inner join on species name, keeping only matches
  dplyr::inner_join(
    y  = filtered_cluster_categories |>
      dplyr::mutate(Name = stringr::str_trim(Name)),
    by = c("Species" = "Name")
  )



# Result columns:
# SampleID, Year, Species, Value, Subcluster, Original_Group

# Plot as before

#-- Assumes long_joined_df contains:
#   SampleID, Year, Species, Value,
#   Subcluster, Original_Group

plots_by_subcluster <- long_joined_df |>
  dplyr::group_by(Subcluster) |>
  dplyr::group_split() |>
  purrr::map(function(df_sub) {
    
    # Subcluster name
    sc <- unique(df_sub$Subcluster)
    og <- unique(df_sub$Original_Group)
    
    # ---- Compute centre (mean trajectory) ----
    centre_df <-
      df_sub |>
      dplyr::group_by(year) |>
      dplyr::summarise(
        centre_value = mean(Value, na.rm = TRUE),
        .groups = "drop"
      )
    
    # ---- Make the plot ----
    ggplot() +
      geom_line(
        data = df_sub,
        aes(x = year, y = Value, group = Species),
        color = "orange",
        size = 0.5,
        alpha = 0.8
      ) +
      geom_line(
        data = centre_df,
        aes(x = year, y = centre_value),
        color = "darkblue",
        size = 1.5
      ) +
      labs(
        title = paste0(sc),
        x = "Years",
        y = "Max-Rescaled Trajectories"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          size  = 16,
          face  = "bold"
        ),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text    = element_text(size = 10)
      )
  })

names(plots_by_subcluster) <- long_joined_df |>
  dplyr::distinct(Subcluster, Original_Group) |>
  dplyr::transmute(name = paste0(Subcluster)) |>
  dplyr::pull(name)


plots_by_subcluster

# 
# desired_order <- names(initial_time_series_plots)
# 
# # keep only the ones that exist, in the desired order
# plots_by_subcluster_ordered <- plots_by_subcluster[desired_order]
# 
# # If any desired names are missing, they'll be NA; drop them safely:
# plots_by_subcluster_ordered <- plots_by_subcluster_ordered[!vapply(plots_by_subcluster_ordered, is.null, logical(1))]


# Blank plot for padding
blank_plot <- ggplot() + theme_void()

# Helper function to prep each group
prepare_group <- function(plots, total = 4, show_x = FALSE, show_y_first_only = TRUE) {
  n <- length(plots)
  # Pad with blank plots
  if (n < total) {
    plots <- c(plots, rep(list(blank_plot), total - n))
  }
  
  # Style each plot
  plots <- lapply(seq_along(plots), function(i) {
    p <- plots[[i]] +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    if (!inherits(p, "gg")) return(p)  # Skip blank
    if (!show_x) {
      p <- p + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
    }
    if (show_y_first_only && i != 1) {
      p <- p + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
    }
    return(p)
  })
  
  return(wrap_plots(plots, ncol = total))
}

# Define your groups (replace `initial_time_series_plots` with your actual list)
group1 <- prepare_group(plots_by_subcluster[1:3])
group2 <- prepare_group(plots_by_subcluster[4:5])
group3 <- prepare_group(plots_by_subcluster[6:8])
group4 <- prepare_group(plots_by_subcluster[9:10])
group5 <- prepare_group(plots_by_subcluster[11:12])
group6 <- prepare_group(plots_by_subcluster[13:15])
group8 <- prepare_group(plots_by_subcluster[16:19], show_x = TRUE)
group7 <- prepare_group(plots_by_subcluster[20:23])

# Combine all groups into one plot (without shared axis titles yet)
final_plot <- group1 / group2 / group3 / group4 / group5 / group6 / group7 / group8 +
  plot_layout(ncol = 1)

final_plot <- final_plot +
  theme(
      plot.margin = margin(5, 5, 5, 5),  # Add space for axis titles
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) 

# --- SHARED AXIS TITLES SECTION ---
# Convert to grob for layout manipulation
g <- patchwork::patchworkGrob(final_plot)

# Create shared axis titles
x_label <- textGrob("Years", gp = gpar(fontsize = 25, fontface = "bold"))
y_label <- textGrob("Max-Rescaled Trajectories", rot = 90, gp = gpar(fontsize = 25, fontface = "bold"))


main_title <- textGrob(
  "Raw Max-Rescaled Trajectories for Each Guild Sub-Cluster",
  gp = gpar(fontsize = 25, fontface = "bold")
)

# Arrange the full layout
final_plot<-grid.arrange(
  main_title,
  arrangeGrob(
    arrangeGrob(
      y_label, g, ncol = 2, widths = unit.c(unit(1.5, "lines"), unit(1, "npc") - unit(1.5, "lines"))
    ),
    x_label,
    nrow = 2,
    heights = unit.c(unit(1, "npc") - unit(1.5, "lines"), unit(1.5, "lines"))
  ),
  nrow = 2,
  heights = unit.c(unit(1.5, "lines"), unit(1, "npc") - unit(1.5, "lines"))
)


# Recommended: cairo_pdf for better font rendering + vector output
grDevices::cairo_pdf(
  filename = "plots/raw_max_rescaled_trajectories_by_guild_subcluster.pdf",
  width = 11.69,  # A4 landscape width in inches
  height = 16.54  # A4 landscape height in inches (tall for many rows)
  # For A3 portrait you might try: width = 11.69, height = 16.54
  # For A3 landscape: width = 16.54, height = 11.69
)

# Draw the grob created earlier
grid::grid.newpage()
grid::grid.draw(final_plot)

dev.off()


################################################################################
############################### Save some files and end code####################
################################################################################

# sain<-main_Func$Name
# 
# filtered_earliest_latest_1975
# 
# filt1<-data.frame(filtered_earliest_latest_1975)
# 
# write.csv(filt1, "Earliest_Latest.csv", row.names = FALSE, fileEncoding = "UTF-8")
# 
# Nata<-data.frame(low_cos2_species$Name, low_cos2_species$Cos_2_scores_12)
# 
# write.csv(Nata, "Los_cos_2.csv", row.names = FALSE, fileEncoding = "UTF-8")
