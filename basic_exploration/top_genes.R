library(pheatmap)
library(tidyverse)
library(ggplot2)
library(scales)

plot_dir <- 'basic_exploration/plots/top_genes/'
number_of_genes <- 20

process_pipeline <- function(diff_expression, label){

  diff_expression <- na.omit(diff_expression[c("log2FoldChange","external_gene_name")])

  
  ordered_by_foldchange <-  diff_expression[order(diff_expression$log2FoldChange, decreasing = TRUE),]

  #
  # # Then select the first 20 genes:
  top_genes_max_pos_foldchange <-  ordered_by_foldchange[1:number_of_genes,]

  top_genes_max_neg_foldchange <- tail(ordered_by_foldchange, n = number_of_genes)

  ggplot(top_genes_max_pos_foldchange, aes(x = reorder(external_gene_name, log2FoldChange), y = log2FoldChange)) +
    geom_bar(stat = "identity") +
    theme_classic()+
    xlab("Gene name") +
    ylab("Fold change") +
    labs(title = "Top 20 Positive Fold Change genes",subtitle = label) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave(paste0(plot_dir, "Positive/Top20Genes_PositiveFC_",label,".png"))

  ggplot(top_genes_max_neg_foldchange, aes(x = reorder(external_gene_name, log2FoldChange), y = log2FoldChange)) +
    geom_bar(stat = "identity") +
    theme_classic()+
    xlab("Gene name") +
    ylab("Fold change") +
    labs(title = "Top 20 Negative Fold Change genes",subtitle = label) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave(paste0(plot_dir, "Negative/Top20Genes_NegativeFC_",label,".png"))

}


preprocessed <- TRUE
loaded_data <- get_data(preprocessed)
mapply(process_pipeline, loaded_data[[1]], loaded_data [[2]])