#depends on the base_clustering.R
group_of_interest <- "KO_t."
scaled_data <- get_scaled_count_data_for_group_of_interest(group_of_interest)

all_genes <- read.csv(file = paste0('CC_checkpoint_genes.csv'), sep =';')
all_stages <- unique(all_genes$Phase)

#Do it for all the phases
for (stage in all_stages){
  genes_of_interest <-subset(all_genes, Phase==stage)
  
  data_to_cluster <- scaled_data[rownames(scaled_data) %in% genes_of_interest$genes, ]
  
  #HCLUST
  cluster_number_hcut <- get_optimal_number_of_clusters(data_to_cluster, hcut)
  hc_clusters <-  run_hclust(data_to_cluster, cluster_number_hcut, paste0(folder_CC, stage, '/'))
  genes_per_cluster <- extract_clusters(data_to_cluster, hc_clusters, cluster_number_hcut)
  save_clusters(genes_per_cluster, paste0(results_dir_hclust, folder_CC, '/', stage, '/'))
}


