library("readxl")
library("Ckmeans.1d.dp")
library("factoextra")
library(tidyverse)
library(dplyr)
library(ggplot2)
library(dtwclust)
library(ggdendro)
library(dendextend)

#depends on the base_clustering.R

#ACTUAL RUN
group_of_interest <- "KO_t."
scaled_data <- get_scaled_count_data_for_group_of_interest(group_of_interest)

#get a list of significantly changed genes
fold_changes <- read.csv(file = 'clustering/data/fold_changes.csv', row.names = 1)
  
#Now subset all genes by the top 100 significant ones
genes_of_interest <- get_top100_genes(fold_changes,'t24_KO_vs_t0_KO')
data_to_cluster <- scaled_data[rownames(scaled_data) %in% genes_of_interest, ]

#Cluster the samples themselves
png(file=paste0(results_dir_sample_clusters, "sample_dendogram_KO.png" ), width=1000, height=600)
p <- cluster_samples(collapsed_replicates)
dev.off()

#KMEANS
cluster_number_kmeans <- get_optimal_number_of_clusters(data_to_cluster, kmeans)
km_clusters <-  run_kmeans(data_to_cluster, cluster_number_kmeans, folder_CC)
genes_per_cluster <- extract_clusters(data_to_cluster, km_clusters, cluster_number_kmeans)
save_clusters(genes_per_cluster, paste0(results_dir_kmeans, folder_top100))

#HCLUST
cluster_number_hcut <- get_optimal_number_of_clusters(data_to_cluster, hcut)
cluster_number_hcut <- 10
hc_clusters <-  run_hclust(data_to_cluster, cluster_number_hcut, folder_CC)
genes_per_cluster <- extract_clusters(data_to_cluster, hc_clusters, cluster_number_hcut)
save_clusters(genes_per_cluster, paste0(results_dir_hclust, folder_top100, cluster_number_hcut, '/'))



