# This code was written by TC and is used to produce plots in the paper 
# Interpretable Clustering: Balancing Representative Groupings with Interpretability.

# The code contained in this script is protected by a V3 creative commons licence
# If you have any questions regarding this code, please contact Tommy Cheale at Tommycheale@gmail.com

# Load in required packages
library(tidyverse)
library(fda)
library(cluster)     
library(factoextra)  
library(dendextend)  
library(geometry)    
library(plotly)      
library(patchwork)
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
library(ggh4x)

###############################################################################
################################ Data Processing ##############################
###############################################################################

# Load data from 2019 SOCBR
SOCB_2019 <- read_csv("All_SOCB_2019_population_indices.csv") %>%
  select(Species = English_Name, year, index) %>%
  arrange(Species, year)

# Load Guild data
SOCB_2019_habitat <- read_csv("SOCB_2019_habitat_guilds.csv")

# Filter to after 1974
SOCB_2019 <- SOCB_2019 %>%
  filter(year >= 1975)

# Identify species with incomplete data
SOCB_2019 <- SOCB_2019 %>%
  group_by(Species) %>%
  mutate(Complete = all(1975:2016 %in% year)) %>%
  ungroup()

removed_species <- SOCB_2019 %>% 
  filter(!Complete) |> 
  group_by(Species) |> 
  summarize(Start = min(year)) |> 
  mutate(Reason = "Incomplete",
         Description = paste("No data before",Start)) |> 
  select(-Start)    
  
# Remove species with incomplete data
SOCB_2019 <- SOCB_2019 %>%
  filter(Complete) %>%
  select(-Complete)

SOCB_2019_habitat <- SOCB_2019_habitat |> 
  filter(! (Species %in% pull(removed_species, "Species")))

# Identify species with only Other habitat tag or No habitat tags
removed_species <- SOCB_2019_habitat |> 
  filter(Guild %in% c("Other","None")) |> 
  mutate(Reason = Guild,
         Description = ifelse(Guild == "Other", "Tagged as other", "No guild tag")) |> 
  select(-Guild) |> 
  bind_rows(removed_species)

# Remove species with only Other habitat tag or No habitat tags
SOCB_2019  <- SOCB_2019 |> 
  filter(! (Species %in% pull(removed_species, "Species")))

SOCB_2019_habitat <- SOCB_2019_habitat |> 
  filter(! (Species %in% pull(removed_species, "Species")))

# Compute max-rescaled indices
SOCB_2019 <- SOCB_2019 %>%
  group_by(Species) %>%
  mutate(mr_index = index / max(index)) |> 
  ungroup()

# Create wide data
SOCB_2019_wide <- SOCB_2019 |> 
  select(Species, year, mr_index) |> 
  pivot_wider(names_from = "Species", values_from = "mr_index")

# Create data matrix required by functions in the fda package
SOCB_2019_matrix <- SOCB_2019_wide |> 
  select(-year) |> 
  as.matrix()

rownames(SOCB_2019_matrix) <- SOCB_2019_wide |> pull(year)
            
################################################################################
################################ Data Explorattion #############################
################################################################################

# Compute mean trends by guild
SOCB_2019_mean_guild <- SOCB_2019 |> 
  left_join(SOCB_2019_habitat, by = "Species") |> 
  group_by(Guild, year) |> 
  summarize(mr_mean = mean(mr_index))

# Plot rescaled curves for select guilds
SOCB_2019 |> 
  left_join(SOCB_2019_habitat, by = "Species") |> 
  filter(Guild %in% (c("Forest Birds", "Aerial Insectivores", "Shorebirds"))) |> 
  ggplot(aes(x = year, y = mr_index, group = Species, color = Guild)) +
  geom_line(alpha = 0.4) +
  facet_wrap(~ Guild) +  # individual species lines
  geom_line(
    data = SOCB_2019_mean_guild |> filter(Guild %in% (c("Forest Birds", "Aerial Insectivores", "Shorebirds"))),
    aes(x = year, y = mr_mean, color = Guild),
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
    axis.text.x = element_text(hjust = 1, size = 22),  # Increase x-axis text size
    axis.text.y = element_text(size = 22),  # Increase y-axis text size
    axis.title.x = element_text(size = 22),  # Increase x-axis title size
    axis.title.y = element_text(size = 22),  # Increase y-axis title size
    plot.title = element_text(size = 25, hjust = 0.5)  # Increase title size and center it
  ) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Cumulative Variance"))

ggsave("plots/Scree.pdf", width = 12, height = 8)

# Compute Cos2 scores (proportion of variance explained by first 2 PCs)
Total_var<- fpca_result$scores^2 %*% fpca_result$values |> 
  as.vector()

Explained_var<- fpca_result$scores[,1:2]^2 %*% fpca_result$values[1:2] |> 
  as.vector()

Cos2_scores <- Explained_var/Total_var

# Construct output data frame
fpca_output <- tibble(Species = colnames(SOCB_2019_matrix),
                      PC1 = fpca_result$scores[,1],
                      PC2 = fpca_result$scores[,2],
                      Cos2 = Cos2_scores) |> 
  left_join(SOCB_2019_habitat, by = "Species")

# Identify species with low cos2 scores
removed_species <- fpca_output |> 
  filter(Cos2 < .80) |> 
  select(Species, Cos2) |> 
  mutate(Reason = "Cos2",
         Description = paste("Cos2 score is",Cos2)) |> 
  select(-Cos2) |> 
  bind_rows(removed_species)

# Compute reconstructed trajectories based on first 2 PCs
harmonics <- fda::eval.fd(time_points, fpca_result$harmonics)

fpca_mean <- fda::eval.fd(time_points, fpca_result$meanfd) |> 
  as.vector()

fpca_fit <- (fpca_result$scores[,1:2] %*% t(harmonics[,1:2]) + fpca_mean) |> 
  as_tibble() |> 
  set_names(1975:2016) |> 
  add_column(Species = SOCB_2019_habitat |> pull("Species"), .before = 1) |> 
  pivot_longer(-Species, names_to = "year", values_to = "fit") |> 
  mutate(year = as.numeric(year))
  
################################################################################################
##################################### Stage 2: Subdividing #####################################
################################################################################################

mahal_dist <- function(M){
  if(nrow(M) == 1){
    return(NA)
  }
  if(nrow(M) == 2){
    return(c(0,0))
  }
  else{
    m <- colMeans(M)
    V <- var(M)
    mahalanobis(M, m, V)
  }
}

subdivide1 <- fpca_output |> 
  filter(! (Guild %in% c("Other","None")),
         Cos2 > .80) |> 
  group_by(Guild) |> 
  group_map(\(guild_data,keys){
    # Extract guild name
    guild <- keys[[1]]
 
    # Compute distance matrix
    d <- dist(guild_data[, c("PC1", "PC2")])
    
    # Apply k-means clustering for each possible k
    km_res <- lapply(2:10, \(k) kmeans(guild_data[, c("PC1", "PC2")], centers = k, nstart = 25))
    
    # Compute silhouette scores
    sil <- lapply(km_res, \(km) silhouette(km$cluster, d))
    
    # Compute average silhouette scores
    sil_summary <- tibble(Guild = guild,
                          k = 2:10,
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
                 Mahal = mahal_dist(clust |> select(PC1,PC2)))
      }) |> 
      bind_rows() |> 
      mutate(Outlier = Mahal > qchisq(0.95, df = 2))
    
    
    list(output = output,
         kmeans = km_res,
         k = optimal_k,
         silhouette = sil,
         silhouette_summary = sil_summary)
  })

# Elbow plots
subdivide1 |> 
  purrr::map(\(res) res$silhouette_summary) |> 
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
    plot.title = element_text(size = 35, hjust = 0.5, face = "bold")  # Increase title size and center it
  ) + 
  plot_annotation(title = "Average Silhouette Score Computed for k = {2,...,10} Clusters Over Each Guild"
)

ggsave("plots/Elbow.pdf", width = 12, height =8)

# Extract data on clusters for all guilds
initial_clusters <- subdivide1 |> 
  purrr::map(\(res) res$output) |>
  bind_rows()


# Plot trajectories for outlying species
initial_clusters |> 
  filter(Singleton | Outlier) |> 
  left_join(SOCB_2019, by = c("Species")) |> 
  mutate(mr_index = mr_index) |> 
  left_join(fpca_fit, by = c("Species","year")) |> 
  ggplot(aes(x = year, y = mr_index, color = Species)) +
  geom_line(linewidth = 1.5) +
  geom_line(mapping = aes(x = year, y = fit), lty = 2,linewidth = 1.5) +
  facet_wrap(~ Species, ncol = 2, scales = "fixed") +
  scale_color_manual(values = c(
    "Evening Grosbeak"   = "steelblue",
    "Rusty Blackbird"    = "darkgreen",
    "Trumpeter Swan"     = "firebrick",
    "Wandering Tattler"  = "purple"
  )) +
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

ggsave("plots/Removed_Out.pdf", width = 12, height = 8)

# Repeat clustering excluding outliers and singletons
subdivide2 <- fpca_output |> 
  left_join(initial_clusters |> select(Species, Singleton, Outlier), by = c("Species")) |> 
  filter(!is.na(Guild),
         !Singleton,
         !Outlier,
         Cos2 > .80) |> 
  group_by(Guild) |> 
  group_map(\(guild_data,keys){
    # Extract guild name
    guild <- keys[[1]]
    
    # Compute distance matrix
    d <- dist(guild_data[, c("PC1", "PC2")])
    
    # Apply k-means clustering for each possible k
    km_res <- lapply(2:10, \(k) kmeans(guild_data[, c("PC1", "PC2")], centers = k, nstart = 25))
    
    # Compute silhouette scores
    sil <- lapply(km_res, \(km) silhouette(km$cluster, d))
    
    # Compute average silhouette scores
    sil_summary <- tibble(Guild = guild,
                          k = 2:10,
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
  purrr::map(\(res) res$output) |>
  bind_rows()

# Plot clusters labelling outliers and singletons
final_clusters |> 
  ggplot(aes(x = PC1, y = PC2, color = factor(Cluster))) +
  geom_point(size=2.5, alpha = 0.6) +
  geom_point(data = initial_clusters %>% filter(Singleton | Outlier),
             color = "red", shape = 17, size = 4, show.legend = FALSE) +
  geom_text(data = initial_clusters %>% filter(Singleton | Outlier),
            aes(label = Species), color = "red", vjust = -1.05, hjust = 0.58, size = 4.75, , show.legend = FALSE) +
  facet_wrap(~ Guild) +
  theme_minimal() +
  labs(title = "Species Projection Coefficients by Guild Sub-Cluster",
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

ggsave("plots/Cluster_Out.pdf", width = 12, height = 8)

# List species by guild and cluster
final_clusters |> 
  group_by(Guild, Cluster) |> 
  summarize(Species = str_c(Species, collapse = ", "))

# Add a unique label for each cluster (e.g., "Forest Birds 1")
final_clusters <- final_clusters |> 
  mutate(Subcluster = paste(Guild, Cluster))

# Compute fits at cluster means
fpca_means <- final_clusters |> 
  group_by(Guild, Cluster, Subcluster) |> 
  summarize(PC1 = mean(PC1),
            PC2 = mean(PC2),
            .groups = "keep") |> 
  reframe(year = 1975:2016,
          Mean = PC1 * harmonics[,1] + PC2 * harmonics[,2] + fpca_mean)

# Plot original trajectories by subcluster
layout <- list(c(1, 2, 3, NA),
               c(4, 5, NA, NA),
               c(6, 7, 8, NA),
               c(9, 10, NA, NA),
               c(11, 12, NA, NA),
               13:16,
               17:20,
               21:24) |> 
  do.call(what = "rbind")

final_clusters |> 
  filter(!(Singleton | Outlier)) |> 
  select(Species, Guild, Cluster, Subcluster) |> 
  left_join(SOCB_2019, by = "Species") |> 
  ggplot(aes(x = year, y = mr_index)) +
  geom_line(aes(group = Species), color = "orange", linewidth = 0.4) +
  geom_line(data = fpca_means, aes(x = year, y = Mean),
            color = "darkblue", linewidth = 1.5) +
  facet_manual(~Subcluster, design = layout) +
  labs(x = "Year", y = "Relative Abundance Function") +
  scale_x_continuous(breaks = seq(1975, 2016, by = 10)) +
  theme_minimal()+
  theme(
    strip.text = element_text(size=10, face = "bold"),
    axis.title.y = element_text(size=25 ),
    axis.text = element_text(size = 11, face="bold"),
    axis.title.x = element_text(size=30 ))

ggsave("plots/All_Subgroups.pdf", width = 12, height = 8)


######################################################################################
############################# Stage 3: Rejoining #####################################
######################################################################################

# Create cluster list structure
subclusters <- final_clusters |> 
  pull("Subcluster") |> 
  unique() |> 
  sort()

clusters_list <- subclusters |> 
  purrr::map(\(sc){
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
midpoints <- dendro_data$segments %>%
  filter(y == yend & y > 0) %>% # horizontal lines at merge points
  group_by(y) %>%
  summarize(x = mean(c(x, xend)), .groups = "drop")  # midpoint x of each horizontal merge line

# Merge with label data
merge_labels_df <- left_join(merge_labels_df, midpoints, by = "y")

label_offset <- max(dendro_data$segments$yend) * 0.11

# Construct dendrogram
ggplot(segment(dendro_data)) +
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

ggsave("plots/Dend.pdf", width = 12, height = 8)

# Plot results for the 3-cluster or 8-cluster solution

# clusters <- 3 # Select one
clusters <- 8

merge <- 24 - clusters

centroids <- PC_points[[merge]] |> 
  group_by(CurrentCluster, Merge) |> 
  summarize(PC1 = mean(PC1),PC2 = mean(PC2)) |> 
  mutate(label = ifelse(Merge == 0, CurrentCluster, paste0("Merge ",Merge, ": ",CurrentCluster)))
            
PC_points[[merge]] |> 
  mutate(label = ifelse(Merge == 0, CurrentCluster, paste0("Merge ",Merge, ": ",CurrentCluster))) |> 
  ggplot(aes(x = PC1, y = PC2, color = factor(label))) +
  # Plot individual points
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_point(size = 2, alpha = 0.7) +
  geom_point(data = PC_points[[merge]] |> filter(Species == "Franklin's Gull"), color = "red", shape = 17, size = 4) + 
  geom_text(data = PC_points[[merge]] |> filter(Species == "Franklin's Gull"),
            aes(label = Species), color = "red", vjust = 2.05, hjust = 0.5, size = 4.75, , show.legend = FALSE) +
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
    subtitle = paste("Merge distance:", merge_df |> filter(Merge == merge) |> pull("Height") |> round(2)),
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

ggsave(paste0("plots/",clusters,"_cluster_scatter.pdf"), width = 12, height = 6)

rejoined_means <- centroids |> 
  group_by(CurrentCluster, Merge) |> 
  reframe(year = 1975:2016,
          Mean = PC1 * harmonics[,1] + PC2 * harmonics[,2] + fpca_mean) |> 
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
  right_join(SOCB_2019, by = "Species") |> 
  filter(!is.na(CurrentCluster)) |> 
  mutate(label = ifelse(Merge == 0, CurrentCluster, paste0("Merge ", Merge)))

join_data |> 
  ggplot(aes(x = year, y = mr_index)) +
  geom_line(aes(group = Species, colour = Guild), alpha = 0.6) + 
  geom_line(data = rejoined_means, aes(y = Mean), color = "darkblue", linewidth = 1.8) +
  scale_color_manual(values = group_colors, drop = FALSE) +  # <- Lock color mapping
  facet_wrap(vars(label), ncol = 3) +
  labs(x = "Year", y = "Rescaled Abundance", color = "Original Group") +
  theme_minimal() +
  theme(strip.text = element_text(size=25, face = "bold"),
        axis.title.y = element_text(size = 30, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "bottom",
        legend.title=element_blank()) +
  geom_line(data = join_data |> filter(Species == "Franklin's Gull"), colour ="red", lty = 2, linewidth = 1.5)

ggsave(paste0("plots/",clusters,"_cluster_line.pdf"), width = 12, height = 6 + 6 * (clusters == 8))

# SJB 200426



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

# Assume "final_clusters" is a list where each element is a data frame
# of species for that initial cluster (as defined earlier using your Subcluster labels).
for (nm in names(final_clusters)) {
  df <- final_clusters[[nm]]
  
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

