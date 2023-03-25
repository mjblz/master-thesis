group_of_interest <- "KO_t."
scaled_data <- get_scaled_count_data_for_group_of_interest(group_of_interest)


#Top100
#For Clusters
cluster_number <- 10

for(i in 1:cluster_number){
  genes <- read.csv(file = paste0('clustering/results/HClust/', folder_top100, cluster_number,'/Cluster_',i,'.csv'), row.names = 1)
  plot_specific_genes_on_timeline(
    genes_to_plot = extract_genes_and_prep_for_plotting(genes,scaled_data), 
    y = 'expression',
    path = paste0('timeline_genes/Clusters/top100/',cluster_number, '/Cluster_',i,'.png'), 
    legend = TRUE, 
    plot_title = '',
    group = 'gene',
    color = 'gene')
}

#Stages
all_genes <- read.csv(file = paste0('CC_checkpoint_genes.csv'), sep =';')
all_stages <- unique(all_genes$Phase)

for(stage in all_stages){
  stage_folder <- paste0('clustering/results/HClust/', folder_CC,stage, '/')
  cluster_number <- count_number_of_clusters_in_folder(stage_folder)
  for(i in 1:cluster_number){
    genes <- read.csv(file = paste0('clustering/results/HClust/', folder_CC, '/',stage, '/Cluster_',i,'.csv'), row.names = 1)
    plot_specific_genes_on_timeline(
      genes_to_plot = extract_genes_and_prep_for_plotting(genes,scaled_data), 
      y = 'expression',
      path = paste0('Clusters/CC/',stage, '/Cluster_',i,'.png'), 
      legend = TRUE, 
      plot_title = '',
      group = 'gene',
      color = 'gene')
  }
}
