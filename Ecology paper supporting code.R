# This code was written by TC and is used to produce plots in the paper 
# Interpretable Clustering: Balancing Representative Groupings with Interpretability.

# The code contained in this script is protected by a V3 creative commons licence
# If you have any questions regarding this code, please contact Tommy Cheale at Tommycheale@gmail.com

# Load in required packages
library(cluster)
library(fda)
library(ggdendro)
library(ggh4x)
library(patchwork)
library(tidyverse)

# Set control parameters
n_pc <- 4 # Number of components used in approximation
cos2_tol <- -1 # Cos^2 score tolerance
figure_export <- FALSE # If TRUE, then save figures to disk
figure_format <- "jpg" # File format for saved figures: jpg or pdf.

###############################################################################
################################ Data Processing ##############################
###############################################################################

# Load data from 2019 SOCBR
SOCB_2019 <- read_csv("All_SOCB_2019_population_indices.csv") |>
  select(Species = English_Name, year, index) |>
  arrange(Species, year)

# Load guild data
SOCB_2019_habitat <- read_csv("SOCB_2019_habitat_guilds.csv")

# Merge data
SOCB_2019 <- SOCB_2019 |> 
  full_join(SOCB_2019_habitat, by = "Species")

# Filter to after 1974
SOCB_2019 <- SOCB_2019 |>
  filter(year >= 1975)

# Identify species with incomplete data
SOCB_2019 <- SOCB_2019 |>
  group_by(Species) |>
  mutate(Complete = all(1975:2016 %in% year)) |>
  ungroup()

removed_species <- SOCB_2019 |> 
  filter(!Complete) |> 
  group_by(Species) |> 
  summarize(Start = min(year)) |> 
  mutate(Reason = "Incomplete",
         Colour = "BrickRed",
         Description = paste("No data before",Start)) |> 
  select(-Start)    
  
# Remove species with incomplete data
SOCB_2019 <- SOCB_2019 |>
  filter(Complete) |>
  select(-Complete)

SOCB_2019_habitat <- SOCB_2019_habitat |> 
  filter(! (Species %in% pull(removed_species, "Species")))

# Identify species with only Other habitat tag or No habitat tags
removed_species <- SOCB_2019_habitat |> 
  filter(Guild %in% c("None")) |> 
  mutate(Reason = Guild,
         Colour = ifelse(Guild == "Other", "ForestGreen", "BurntOrange"),
         Description = ifelse(Guild == "Other", "Tagged as other", "No guild tag")) |> 
  select(-Guild) |> 
  bind_rows(removed_species)

# Remove species with only Other habitat tag or No habitat tags
SOCB_2019  <- SOCB_2019 |> 
  filter(! (Species %in% pull(removed_species, "Species")))

SOCB_2019_habitat <- SOCB_2019_habitat |> 
  filter(! (Species %in% pull(removed_species, "Species")))

###############################################################################
################################## Resmoothing ################################
###############################################################################

# Create data matrix
SOCB_2019_matrix <- SOCB_2019 |> 
  select(Species, year, index) |> 
  pivot_wider(names_from = "Species", values_from = "index") |> 
  select(-year) |> 
  as.matrix()

# Fit B-splines 
K <- nrow(SOCB_2019_matrix)
time_points <- 1:K
basis_obj <- create.bspline.basis(rangeval = range(time_points), nbasis = K)
smooth_result <- smooth.basis(argvals = time_points, 
                              y = SOCB_2019_matrix, 
                              fdParobj = basis_obj)

# Extract resmoothed fits
SOCB_2019_resmoothed <- eval.fd(1:K,smooth_result$fd)

# Apply max-rescaling
SOCB_2019_max_rescaled <- SOCB_2019_resmoothed |>
  apply(2, \(y) y/max(y))

# Append to original data
SOCB_2019 <- SOCB_2019_max_rescaled |> 
  as_tibble() |> 
  add_column(year = 1975:2016) |> 
  pivot_longer(-year,names_to = "Species", values_to = "mr_index") |> 
  left_join(SOCB_2019, by = c("Species","year"))

################################################################################
################################ Data Exploration ##############################
################################################################################

# Compute mean trends by guild
SOCB_2019_mean_guild <- SOCB_2019 |> 
  group_by(Guild, year) |> 
  summarize(mr_mean = mean(mr_index),
            .groups = "drop")

# Plot rescaled curves for select guilds
SOCB_2019 |> 
  filter(Guild %in% (c("Forest Birds", "Aerial Insectivores", "Shorebirds"))) |> 
  ggplot(aes(x = year, y = mr_index, group = Species, color = Guild)) +
  geom_line(alpha = 0.4) +
  facet_wrap(~ Guild) +  # individual species lines
  geom_line(
    data = SOCB_2019_mean_guild |> filter(Guild %in% (c("Forest Birds", "Aerial Insectivores", "Shorebirds"))),
    aes(x = year, y = mr_mean, color = Guild),
    inherit.aes = FALSE,
    linewidth = 1.5
  ) +
  scale_color_manual(values = c(
    "Shorebirds" = "#ff7f0e",
    "Forest Birds" = "#8c564b",
    "Aerial Insectivores" = "#7f7f7f"
  )) +
  labs(
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

if(figure_export){
  "plots/Examples" |> 
    paste0(".",figure_format) |> 
    ggsave(width = 12, height =8)
}

###############################################################################
################################ Stage 1: FPCA ################################
###############################################################################

# Reit B-splines to max-rescaled trajectories
smooth_result <- smooth.basis(argvals = time_points, 
                              y = SOCB_2019_max_rescaled, 
                              fdParobj = basis_obj)

# Conduct FPCA 
fpca_result <- pca.fd(smooth_result$fd, nharm = 42) # Perform FPCA

# Plot harmonics
harmonics <- eval.fd(time_points, fpca_result$harmonics) 

harmonics[,1:4] |> 
  as_tibble() |> 
  add_column(Year = 1975:2016, .before = 1) |> 
  pivot_longer(-Year, names_to = "PC", values_to = "Value") |> 
  mutate(PC = str_sub(PC, start = 3)) |> 
  ggplot(aes(x = Year, y = Value)) + 
  geom_line() + 
  geom_hline(yintercept = 0, lty = 2) + 
  facet_wrap(vars(PC)) +
  ylab("") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 20),
    axis.text.x = element_text(hjust = 1, size = 22),  # Increase x-axis text size
    axis.text.y = element_text(size = 22),  # Increase y-axis text size
    axis.title.x = element_text(size = 22),  # Increase x-axis title size
    axis.title.y = element_text(size = 22),  # Increase y-axis title size
    plot.title = element_text(size = 25, hjust = 0.5)  # Increase title size and center it
  ) 




# Scree plot
scree_df <- tibble(pc = 1:42,
                   pve = fpca_result$varprop,
                   cum_pve = cumsum(pve))

scree_df |> 
  filter(pc <= 8) |> 
  ggplot(aes(x = pc, y = pve)) +
  geom_bar(stat = "identity", fill = "blue", color = "blue", alpha = 0.7) +
  geom_line(aes(y = cum_pve, group = 1), color = "darkgreen", linewidth = 1) +
  geom_point(aes(y = cum_pve), color = "darkgreen", size = 3, shape = 17) +
  labs(x = "Functional Principal Component",
       y = "Proportion of Variance Explained") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1, size = 22),  # Increase x-axis text size
    axis.text.y = element_text(size = 22),  # Increase y-axis text size
    axis.title.x = element_text(size = 22),  # Increase x-axis title size
    axis.title.y = element_text(size = 22),  # Increase y-axis title size
    plot.title = element_text(size = 25, hjust = 0.5)  # Increase title size and center it
  ) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Cumulative Variance"))

if(figure_export){
  "plots/Scree" |> 
    paste0(".",figure_format) |> 
    ggsave(width = 12, height = 8)
}

# Compute Cos2 scores 
Explained_var <- fpca_result$scores[,1:5]^2 |> 
  apply(1,cumsum)

Total_var <- fpca_result$scores^2 |> 
  apply(1,sum)

Cos2_scores <- (t(Explained_var)/Total_var) |> 
  as_tibble() |> 
  set_names(paste0("Cos2_",1:5))

# Construct output data frame
fpca_output <- tibble(Species = colnames(SOCB_2019_max_rescaled),
                      PC = fpca_result$scores[,1:n_pc]) |> 
  bind_cols(Cos2 = Cos2_scores) |> 
  left_join(SOCB_2019_habitat, by = "Species")

# Boxplot of cos2 scores vs number of components
fpca_output |> 
  select(Species, starts_with("Cos2")) |> 
  pivot_longer(-Species, names_to = "Order", values_to = "Value") |> 
  separate(Order, c("Dummy","Order"), "_") |> 
  select(-Dummy) |> 
  filter(Order > 1) |> 
  ggplot(aes(x = Order, y = Value)) +
  geom_boxplot() +
  geom_hline(yintercept = .75, colour = "orange", lty = 2, linewidth = 1.5) +
  ylab("Cos2 Score") + 
  xlab("Order of Approximation") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 30, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title=element_blank())

if(figure_export){
  "plots/Cos2_boxplot" |> 
    paste0(".",figure_format) |> 
    ggsave(width = 12, height =8)
}

# Identify species with low cos2 scores
removed_species <- fpca_output |> 
  filter(Cos2_2 < cos2_tol) |> 
  select(Species, Cos2_2) |> 
  mutate(Reason = "Cos2",
         Colour = "RoyalPurple",
         Description = paste("Cos$^2$ score is",round(Cos2_2,2))) |> 
  select(-starts_with("Cos2")) |> 
  bind_rows(removed_species)

# Create LaTeX table of removed species
removed_species |> 
  arrange(Species) |> 
  mutate(Text = paste0("\\textcolor{",Colour,"}{",Species,"}", "&",
                       "\\textcolor{",Colour,"}{",Description,"}")) |> 
  select(Text) |>
  rowid_to_column("Row") |> 
  mutate(Row = ceiling(Row/2)) |> 
  group_by(Row) |> 
  reframe(Text = paste(paste(Text,collapse = "&"),"\\\\")) |> 
  select(Text) |> 
  write_csv("output/removed_species.tex", col_names = FALSE)
  

# Compute reconstructed trajectories
fpca_mean <- eval.fd(time_points, fpca_result$meanfd) |> 
  as.vector()

SOCB_2019 <- (fpca_result$scores[,1:n_pc] %*% t(harmonics[,1:n_pc])) |> 
  as_tibble(.name_repair = "minimal") |> 
  set_names(1975:2016) |> 
  add_column(Species = SOCB_2019_habitat |> pull("Species"), .before = 1) |> 
  pivot_longer(-Species, names_to = "year", values_to = "fit") |> 
  mutate(year = as.numeric(year), fit = fit + fpca_mean) |> 
  left_join(SOCB_2019, by = c("Species", "year"))
  
################################################################################################
##################################### Stage 2: Subdividing #####################################
################################################################################################

mahal_dist <- function(M, pc){
  if(nrow(M) <= pc){
    return(rep(0, nrow(M)))
  }
  else{
    m <- colMeans(M)
    V <- var(M)
    mahalanobis(M, m, V)
  }
}

subdivide1 <- fpca_output |> 
  filter(Cos2_2 > cos2_tol) |> 
  group_by(Guild) |> 
  group_map(\(guild_data,keys){
    # Extract guild name
    guild <- keys[[1]]

    # Compute distance matrix
    d <- dist(guild_data[, "PC"])
    
    # Apply k-means clustering for each possible k
    max_k <- min(10, nrow(guild_data) - 1)
    km_res <- lapply(2:max_k, \(k) kmeans(guild_data[, "PC"], centers = k, nstart = 25))
    
    # Compute silhouette scores
    sil <- lapply(km_res, \(km) silhouette(km$cluster, d))
    
    # Compute average silhouette scores
    sil_summary <- tibble(Guild = guild,
                          k = 2:max_k,
                          avg_sil = sapply(sil, \(s) mean(s[,"sil_width"])))
    
    # Select the k with the maximum average silhouette width
    optimal_k <- sil_summary |> 
      filter(avg_sil == max(avg_sil)) |> 
      pull("k")
    
    # Check for outliers/singletons
    output <- guild_data |> 
      add_column(Guild = guild,
                 Cluster = km_res[[optimal_k - 1]]$cluster) |> 
      group_by(Cluster) |>
      group_map(\(clust,keys){
        clust |> 
          mutate(Cluster = keys[[1]],
                 Singleton = (n() == 1),
                 Mahal = mahal_dist(clust |> select(PC), n_pc))
      }) |> 
      bind_rows() |> 
      mutate(Outlier = Mahal > qchisq(0.95, df = n_pc))
    
    
    list(output = output,
         kmeans = km_res,
         k = optimal_k,
         silhouette = sil,
         silhouette_summary = sil_summary)
  })

# Elbow plots
subdivide1 |> 
  map(\(res) res$silhouette_summary) |> 
  bind_rows() |> 
  ggplot(aes(x = k, y = avg_sil)) +
  geom_line() +
  geom_point(size = 3) +
  labs(x = "Number of Clusters (k)",
       y = "Average Silhouette Score") +
  theme_minimal()+
  facet_wrap(vars(Guild)) +
  scale_x_continuous(limits = c(1, 10), breaks = 1:10) +
  theme(
    axis.title.x = element_text(size = 20),  # Increase x-axis title size
    axis.title.y = element_text(size = 20),  # Increase y-axis title size
    plot.title = element_text(size = 35, hjust = 0.5, face = "bold"),  # Increase title size and center it
    strip.text = element_text(size = 25))

if(figure_export){
  "plots/Elbow" |> 
    paste0(".",figure_format) |> 
    ggsave(width = 12, height =8)
}

# Extract data on clusters for all guilds
initial_clusters <- subdivide1 |> 
  map(\(res) res$output) |>
  bind_rows()

# Plot trajectories for outlying species
initial_clusters |> 
  filter(Singleton | Outlier) |> 
  left_join(SOCB_2019, by = c("Species")) |> 
  ggplot(aes(x = year, y = mr_index, color = Species)) +
  geom_line(linewidth = 1.5) +
  geom_line(mapping = aes(x = year, y = fit), lty = 2,linewidth = 1.5) +
  facet_wrap(~ Species, ncol = 2, scales = "fixed") +
  scale_color_manual(values = c(
    "Evening Grosbeak"   = "steelblue",
    "Rusty Blackbird"    = "darkgreen",
    "Trumpeter Swan"     = "firebrick",
    "Eared Grebe"  = "purple"
  )) +
  labs(
    x = "Year",
    y = "Max-Rescaled Abundance"
  ) +
  theme_minimal(base_size = 22) +  # makes everything bigger by default
  theme(
    legend.position = "none",
    strip.text = element_text(size = 25, face = "bold"),   # facet titles
    axis.title = element_text(size = 24, face = "bold"),   # axis labels
    axis.text = element_text(size = 20),                   # tick labels
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5)
  )

if(figure_export){
  "plots/Removed_Out" |> 
    paste0(".",figure_format) |> 
    ggsave(width = 12, height = 6)
}

# Repeat clustering excluding outliers and singletons
subdivide2 <- fpca_output |> 
  left_join(initial_clusters |> select(Species, Singleton, Outlier), by = c("Species")) |> 
  filter(!is.na(Guild),
         !Singleton,
         !Outlier,
         Cos2_2 > cos2_tol) |> 
  group_by(Guild) |> 
  group_map(\(guild_data,keys){
    # Extract guild name
    guild <- keys[[1]]
    
    # Compute distance matrix
    d <- dist(guild_data[, "PC"])
    
    # Apply k-means clustering for each possible k
    max_k <- min(10, nrow(guild_data) - 1)
    km_res <- lapply(2:max_k, \(k) kmeans(guild_data[, "PC"], centers = k, nstart = 25))
    
    # Compute silhouette scores
    sil <- lapply(km_res, \(km) silhouette(km$cluster, d))
    
    # Compute average silhouette scores
    sil_summary <- tibble(Guild = guild,
                          k = 2:max_k,
                          avg_sil = sapply(sil, \(s) mean(s[,"sil_width"])))
    
    # Select the k with the maximum average silhouette width
    optimal_k <- sil_summary |> 
      filter(avg_sil == max(avg_sil)) |> 
      pull("k")
    
    # Check for outliers/singletons
    output <- guild_data |> 
      add_column(Guild = guild,
                 Cluster = km_res[[optimal_k - 1]]$cluster) 
    
    list(output = output,
         kmeans = km_res,
         k = optimal_k,
         silhouette = sil,
         silhouette_summary = sil_summary)
  })

final_clusters <- subdivide2 |> 
  map(\(res) res$output) |>
  bind_rows()

# Compute fits at cluster means
fpca_means <- final_clusters |> 
  group_by(Guild, Cluster)  |> 
  summarize(PC = t(apply(PC, 2, mean)),
            .groups = "keep") |> 
  reframe(year = 1975:2016,
          Mean = (harmonics[,1:n_pc] %*% t(PC) + fpca_mean)[,1])

# Reorder the clusters
cluster_order <- fpca_means |> 
  filter(year %in% c(1975, 2016)) |> 
  pivot_wider(names_from = year, values_from = Mean) |> 
  mutate(Diff = `2016` - `1975`) |> 
  arrange(Guild, Diff)|> 
  group_by(Guild)  |> 
  mutate(New_Cluster = 1:n()) |> 
  select(Guild, Cluster, New_Cluster)

final_clusters <- final_clusters |> 
  full_join(cluster_order, by = c("Guild","Cluster")) |> 
  mutate(Cluster = New_Cluster) |> 
  select(-New_Cluster)

fpca_means <- fpca_means |> 
  full_join(cluster_order, by = c("Guild","Cluster")) |> 
  mutate(Cluster = New_Cluster) |> 
  select(-New_Cluster)


# Add a unique label for each cluster (e.g., "Forest Birds 1")
final_clusters <- final_clusters |> 
  mutate(Subcluster = paste(Guild, Cluster)) 

fpca_means <- fpca_means |> 
  mutate(Subcluster = paste(Guild, Cluster)) 

# Plot clusters labelling outliers and singletons
final_clusters |> 
  ggplot(aes(x = PC[,1], y = PC[,2], color = factor(Cluster))) +
  geom_point(size=2.5, alpha = 0.6) +
  geom_point(data = initial_clusters |> filter(Singleton | Outlier),
             color = "red", shape = 17, size = 4, show.legend = FALSE) +
  geom_text(data = initial_clusters |> filter(Singleton | Outlier),
            aes(label = Species), color = "red", vjust = -1.05, hjust = 0.58, size = 4.75, , show.legend = FALSE) +
  facet_wrap(~ Guild) +
  theme_minimal() +
  labs(x = "FPC1",
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

if(figure_export){
  "plots/Cluster_Out" |> 
    paste0(".",figure_format) |> 
    ggsave(width = 12, height = 8)
}

# Plot original trajectories by subcluster
layout <- list(c(1, 2, 3, NA, NA, NA), # Aerial Insectivores
               c(NA, NA, 4, 5, NA, NA), # Birds of Prey
               c(6, 7, 8, 9, 10, 11), # Forest Birds
               c(12, NA, 13, NA, NA, NA), # Grassland Birds
               c(14, 15, 16, 17, 18, NA), # Other
               c(NA, NA, 19, NA, 20, NA), # Seabirds
               c(21, 22, 23, NA, NA, NA), # Shorebirds
               c(NA, NA, 24, 25, 26, 27), # Waterfowl
               c(28, 29, NA, 30, 31, NA)) |> # Wetland Birds
  do.call(what = "rbind")

final_clusters |> 
  filter(!(Singleton | Outlier)) |> 
  select(Species, Guild, Cluster, Subcluster) |> 
  left_join(SOCB_2019, by = "Species") |> 
  ggplot(aes(x = year, y = mr_index)) +
  geom_line(aes(group = Species), color = "orange", linewidth = 0.4) +
  geom_line(data = fpca_means, aes(x = year, y = Mean),
            color = "darkblue", linewidth = 1.5) +
  #facet_manual(~Subcluster, design = layout) +
  #facet_grid(Guild ~ Subcluster) + 
  facet_wrap(vars(Subcluster))+
  labs(x = "Year", y = "Relative Abundance Function") +
  scale_x_continuous(breaks = seq(1975, 2016, by = 10)) +
  theme_minimal()+
  theme(
    strip.text = element_text(size=10, face = "bold"),
    axis.title.y = element_text(size=25 ),
    axis.text = element_text(size = 11, face="bold"),
    axis.title.x = element_text(size=30 ))

if(figure_export){
  "plots/All_Subgroups" |> 
    paste0(".",figure_format) |> 
    ggsave(width = 12, height = 8)
}

# Generate LaTeX output
final_clusters |> 
  group_by(Guild) |> 
  group_map(\(df,keys){
    tibble(Text = paste("\\multicolumn{2}{l}{\\textbf{",keys[1],"}}\\\\")) |> 
      bind_rows(df |> 
                  group_by(Cluster) |> 
                  summarize(Text = paste(Cluster[1]," & ", paste(Species, collapse = ", "), "\\\\")) |> 
                  select(Text))
  }) |> 
  bind_rows() |> 
  write_csv("output/subcluster_membership.tex", col_names = FALSE, quote = "none")

######################################################################################
############################# Stage 3: Rejoining #####################################
######################################################################################

# Create cluster list structure
subclusters <- final_clusters |> 
  pull("Subcluster") |> 
  unique() |> 
  sort()

clusters_list <- subclusters |> 
  map(\(sc){
    list(species = final_clusters |> filter(Subcluster == sc) |> add_column(Merge = 0),
         lineage = sc)
  }) |> 
  structure(names = subclusters)

cluster_indices <- setNames(-seq_len(length(clusters_list)), names(clusters_list))
next_index <- 1

PC_points<-list()
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
      sp1 <- clusters_list[[ cluster_names[i] ]]$species[, "PC", drop = FALSE]
      sp2 <- clusters_list[[ cluster_names[j] ]]$species[, "PC", drop = FALSE]
      
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

  new_cluster_name <- new_lineage
  new_cluster <- list(
    species = new_species,
    lineage = new_lineage
  )
  
  # Register merge indices for dendrogram reconstruction
  idx1 <- cluster_indices[name1]
  idx2 <- cluster_indices[name2]
  pair_indices <- c(idx1, idx2)
  new_index <- next_index
  next_index <- next_index + 1
  
  if(iteration == 1)
    merge_df <- tibble(Merge = 1,
                       Cluster1 = pair_indices[1],
                       Cluster2 = pair_indices[2],
                       Height = merge_distance)
  else
    merge_df <- tibble(Merge = iteration,
                       Cluster1 = pair_indices[1],
                       Cluster2 = pair_indices[2],
                       Height = merge_distance) |> 
    bind_rows(merge_df) |> 
    arrange(Merge)
  
  # Remove old clusters and add the new merged one
  clusters_list[[name1]] <- NULL
  clusters_list[[name2]] <- NULL
  clusters_list[[new_cluster_name]] <- new_cluster
  
  clusters_list[[new_cluster_name]]$species$Merge <- iteration
  
  cluster_indices <- cluster_indices[!names(cluster_indices) %in% c(name1, name2)]
  cluster_indices[new_cluster_name] <- new_index
  
  current_species <- do.call(rbind, lapply(clusters_list, function(cl) cl$species))
  current_species$CurrentCluster <- rep(names(clusters_list), sapply(clusters_list, function(cl) nrow(cl$species)))
  
  PC_points[[iteration]]<-current_species
}

# Construct hclust obj to make use of existing functions for visualizing clusters
n_orig <- length(subclusters)

hclust_obj <- list(
  merge = merge_df |> select(Cluster1, Cluster2) |> as.matrix(),
  height = merge_df |> pull(Height),
  order = 1:n_orig,
  labels = subclusters,
  method = "complete-linkage",
  call = match.call()
) |> 
  structure(class = "hclust")

# Construct dendrogram
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
midpoints <- dendro_data$segments |>
  filter(y == yend & y > 0) |> # horizontal lines at merge points
  group_by(y) |>
  summarize(x = mean(c(x, xend)), .groups = "drop")  # midpoint x of each horizontal merge line

# Merge with label data
merge_labels_df <- left_join(merge_labels_df, midpoints, by = "y")

label_offset <- max(dendro_data$segments$yend) * 0.11

# Identify solutions to focus on
merge_df |> 
  mutate(Height1 = lead(Height),
         Diff = Height1 - Height) |> 
  arrange(desc(Diff))

merges <-  c(26,18)
clusters <- nrow(merge_df) - merges +1 

# Construct dendrogram

line_heights <- (merge_df$Height[merges] + merge_df$Height[merges + 1])/2

ggplot(segment(dendro_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(dendro_data),
            aes(x = x, y = -label_offset, label = label),
            hjust = -0, vjust=-0.8 , angle = 90, size = 4.5) +
  labs(x = "Sub-cluster", y = "Height") +
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
geom_hline(yintercept = line_heights, colour = c("blue","orange"), linetype = "dashed", linewidth = 1)

if(figure_export){
  "plots/Dend" |> 
    paste0(".",figure_format) |> 
    ggsave(width = 12, height = 8)
}

# Plot results for the selected solutions

for(i in 1:length(clusters)){

  cluster <- clusters[i]
  merge <- merges[i]
  
  centroids <- PC_points[[merge]] |> 
    group_by(CurrentCluster, Merge) |> 
    summarise(PC = t(apply(PC, 2, mean))) |> 
    mutate(label = ifelse(Merge == 0, CurrentCluster, paste0("Merge ",Merge, ": ",CurrentCluster)))
  
  PC_points[[merge]] |> 
    mutate(label = ifelse(Merge == 0, CurrentCluster, paste0("Merge ",Merge, ": ",CurrentCluster))) |> 
    ggplot(aes(x = PC[,1], y = PC[,2], color = factor(label))) +
    # Plot individual points
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) +
    geom_point(size = 2, alpha = 0.7) +
    # Plot centroids as squares in the same color
    geom_point(data = centroids, shape = 15, size = 5) +
    # Ensure consistent colors
    scale_color_manual(values = rainbow(length(unique(current_species$CurrentCluster)))) +
    # Wrap legend labels to prevent long labels from cutting off
    guides(color = guide_legend(override.aes = list(size = 4), title = "Cluster (Lineage)")) +
    xlim(c(-3,3)) + 
    ylim(c(-3,3)) + 
    theme_minimal() +
    labs(
      title = paste("Merge", merge),
      subtitle = paste("Height:", merge_df |> filter(Merge == merge) |> pull("Height") |> round(2)),
      color = "Cluster (Lineage)",
    ) +
    theme(
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      axis.title.x = element_text(size=18, face="bold"),
      axis.title.y = element_text(size=18, face="bold"),
      legend.box = "vertical"     # spacing between items
    )+
    scale_color_discrete(labels = function(x) str_wrap(x, width = 67))

  if(figure_export){  
    paste0("plots/",cluster,"_cluster_scatter",".",figure_format) |> 
      ggsave(width = 12, height = 6)
  }
  
  rejoined_means <- centroids |> 
    group_by(CurrentCluster, Merge) |> 
    reframe(year = 1975:2016,
            Mean = (harmonics[,1:n_pc] %*% t(PC) + fpca_mean)[,1]) |> 
    mutate(label = ifelse(Merge == 0, CurrentCluster, paste0("Merge ",Merge)))
  
  group_colors <- c("Wetland Birds" = "#1f77b4",
                    "Shorebirds" = "#ff7f0e",
                    "Birds of Prey" = "#2ca02c",
                    "Seabirds" = "#d62728",
                    "Waterfowl" = "#9467bd",
                    "Forest Birds" = "#8c564b",
                    "Grassland Birds" = "#e377c2",
                    "Aerial Insectivores" = "#7f7f7f")
  
  join_data <- PC_points[[merge]] |> 
    select(Species, Cluster, Subcluster, Merge, CurrentCluster) |> 
    right_join(SOCB_2019, by = "Species") |> 
    filter(!is.na(CurrentCluster)) |> 
    mutate(label = ifelse(Merge == 0, CurrentCluster, paste0("Merge ", Merge)))
  
  join_data |> 
    ggplot(aes(x = year, y = mr_index)) +
    geom_line(aes(group = Species, colour = Guild), alpha = 0.6) + 
    geom_line(data = rejoined_means, aes(y = Mean), color = "darkblue", linewidth = 1.8) +
    scale_color_manual(values = group_colors, drop = FALSE) +  # <- Lock color mapping
    facet_wrap(vars(label), nrow = (cluster - 1) %/% 4 + 1) +
    labs(x = "Year", y = "Rescaled Abundance", color = "Original Group") +
    theme_minimal() +
    theme(strip.text = element_text(size=16, face = "bold"),
          axis.title.y = element_text(size = 30, face = "bold"),
          axis.title.x = element_text(size = 30, face = "bold"),
          legend.text = element_text(size = 14),
          legend.position = "bottom",
          legend.title=element_blank())
  
  if(figure_export){
    paste0("plots/",cluster,"_cluster_line", ".",figure_format) |> 
      ggsave(width = 12, height = 6 + 3 * (clusters == 8))
  }
}