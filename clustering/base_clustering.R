library("readxl")
library("Ckmeans.1d.dp")
library("factoextra")
library(tidyverse)
library(dplyr)
library(ggplot2)
library(dtwclust)
library(ggdendro)
library(dendextend)

#Depends on utilities.R

results_dir_sample_clusters <- 'clustering/results/Sample_Clusters/'
results_dir_kmeans <- 'clustering/results/Kmeans/'
results_dir_hclust <- 'clustering/results/HClust/'
folder_CC <- "KO_time_CC/"
folder_top100 <- 'KO_time_top100/'
#KO_time, WT_time

cluster_samples <- function(samples){
  #Sample clustering, just to see if there are outliers and if the order kind of makes sense
  hc <- hclust(as.dist(1-cor(samples, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
  TreeC = as.dendrogram(hc, method="average")
  return(plot(TreeC,
       main = "Sample Clustering",
       ylab = "Height"
  ))
}

get_top100_genes <- function(fold_changes, comparison){
  fold_change_column <- paste0('log2FoldChange.', comparison)
  padj_column <- paste0('padj.', comparison)
  #We want to only observe the significant genes between t24 vs t0, order them by significance and from those we also only take a subset
  row.names(fold_changes[, c(fold_change_column, padj_column)][order(fold_changes[[padj_column]]), ][1:100, ])
}


#type can be kmeans or hcut
get_optimal_number_of_clusters <- function(clustering_data, clustering_type){
  #Get optimal number of clusters
  n_clust <- fviz_nbclust(clustering_data, clustering_type, method = "silhouette", k.max = 8) +
    theme_minimal()
  n_clust<-n_clust$data
  max_cluster<-as.numeric(n_clust$clusters[which.max(n_clust$y)])
  return(max_cluster)
}

scale_data <- function(data){
  #Scaling
  data <- data %>% 
    # transpose the matrix so genes are as columns
    t() %>% 
    # apply scalling to each column of the matrix (genes)
    scale() %>% 
    # transpose back so genes are as rows again
    t() %>% 
    na.omit()
}

extract_clusters <- function(cluster_data, clusters, cluster_number){
  
  genes_per_cluster <- list()
  #Get the clusters with the actual genes in them
  clustered_genes <- as.data.frame(cbind(gene = row.names(cluster_data), cluster = clusters))
  
  for (i in 1:cluster_number){
    cluster <- data.frame(genes = clustered_genes[clustered_genes$cluster ==i,]$gene)
    genes_per_cluster <- append(genes_per_cluster,cluster)
  }
  
  return(genes_per_cluster)
  
}


run_kmeans <- function(sign_genes_matrix, cluster_number, folder_for_plot){
  #Actual Clustering K-means
  set.seed(20)
  km_res <- kmeans(sign_genes_matrix, centers=cluster_number, nstart = 1000, iter.max = 20)
  km_clusters <- km_res$cluster
  
  #Plotting
  fviz_cluster(km_res, data = sign_genes_matrix,
               geom = "point",
               ellipse.type = "convex", 
               ggtheme = theme_bw()
  )
  ggsave(paste0(results_dir_kmeans, folder_for_plot, "clusterplot.png" ))
  
  return(km_clusters)
}

run_hclust <- function(sign_genes_matrix, cluster_number, folder_for_plot){
  #Clustering hclust
  
  #Distance Matrix
  dist_mat <- dist(sign_genes_matrix, method="euclidean")
  
  #Clustering for time series
  hc_model <- hclust(dist_mat, method = "ward.D")
  
  #Extract Clusters
  hc_clusters <- cutree(hc_model, k = cluster_number)
  
  #Plotting
  
  png(file=paste0(results_dir_hclust, folder_for_plot, "clusters_dendogram.png" ), width=1000, height=600)
  par(mar = c(10,1,1,1))
  dend <- hc_model %>% as.dendrogram %>%
    set("branches_k_color", k=4) %>%
    set("labels_cex", 0.8)
  dend %>% 
    plot(main="Dendrogram KO over time gene clustering")
  dev.off()
  
  return(hc_clusters)
}

save_clusters <- function(clusters, dir){
  for (i in seq_along(clusters)) {
    write.csv(clusters[i],paste0(dir, 'Cluster_', i,".csv"))
  }
}

count_number_of_clusters_in_folder <- function(path){
  all_files <- list.files(path)
  
  # Use grepl to get files with csv extension
  csv_files <- all_files[grepl("\\.csv$", all_files)]
  
  # Count the number of csv files
  return(length(csv_files))
}
