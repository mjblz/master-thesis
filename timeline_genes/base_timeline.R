library(ggplot2)
library(tidyr)


# fold_changes <- read.csv(file = 'clustering/data/fold_changes.csv', row.names = 1)
# 
# #Select out relevant genes
# 
# 
# fold_changes_subset <- subset(fold_changes_subset, 
#                     select = 
#                       grepl("log2FoldChange", names(fold_changes_subset)) & 
#                       grepl("KO", names(fold_changes_subset)) & 
#                       !grepl("WT", names(fold_changes_subset)) & 
#                       !grepl("KO_supp", names(fold_changes_subset))
#                     )
#                     

results_dir_timeline <- 'timeline_genes/results/'

extract_genes_and_prep_for_plotting <- function(genes_to_extract, scaled_data){
  #Extract the relevant genes
  subset_data <- subset(scaled_data, row.names(scaled_data) %in% genes_to_extract$genes)
  subset_data$gene <- rownames(subset_data)
  
  #Gather to long format for plotting
  long_df <- gather(subset_data, timepoint, expression, KO_t.0:KO_t.24)
  
  #Define order of timepoints
  timepoint_order <- c("KO_t.0", "KO_t.4", "KO_t.10", "KO_t.16", "KO_t.24")
  
  # Convert the timepoint variable to a factor with the desired order
  long_df$timepoint <- factor(long_df$timepoint, levels = timepoint_order)
  
  return(long_df)
}


plot_specific_genes_on_timeline <- function(genes_to_plot, y, path, legend, plot_title, group, color){
  print(color)
  
  genes_plot <- ggplot(data=genes_to_plot, aes(x=timepoint, y=!!sym(y), group=!!sym(group), color = !!sym(color))) +
    geom_line()+
    geom_point() +
    theme(legend.position = c(.5, .5), legend.background = element_rect(fill = alpha("white", 0.6)))
  
  if (legend == FALSE){
    genes_plot <- genes_plot +
      theme(legend.position="none")
  }
  
  if (plot_title != ''){
    genes_plot <- genes_plot +
      ggtitle(plot_title)
  }
  ggsave(paste0(results_dir_timeline, path), device = "png")
}



