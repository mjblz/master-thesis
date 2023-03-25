library(here)
library(hash)

library(ggplot2)

here::i_am("basic_exploration/significant_change_gene_count.R")

# This script counts the number of significantly changed (measured by padj) genes for each DE comparison and plots it.
# Furthermore it saves the list of those significantly changed genes in a condition name - list of genes hash.
#It depends on the utilities.R for loading the data

p_value_threshhold <- 0.01

#Global so other script can use it
#This hash holds the name of each comparison as the key and the list of significantly changed genes as a value
significant_genes_per_DE <- hash()

significant_gene_cnt_per_DE <- data.frame('Label'=character(0), 'Gene_Count'=numeric(0));

process_pipeline <- function(diff_expression, DE_name){
  print(DE_name)
  
  diff_expression <- na.omit(diff_expression[c("external_gene_name","padj")])
  #count_significantly_changed_genes(data, DE_name)
  
  #Filter against threshold
  filtered_below_threshold <- diff_expression[diff_expression$padj < p_value_threshhold, ]
  
  #get the gene names as list and store for later comparison, example t0_WT_vs_t0_KO : 61
  .set( significant_genes_per_DE, DE_name, as.list(filtered_below_threshold$external_gene_name))
  
  #Add label - count information to dataframe
  significant_gene_cnt_per_DE <<- significant_gene_cnt_per_DE %>% 
    add_row(Label = DE_name, Gene_Count = as.numeric(nrow(filtered_below_threshold)))
  #significant_gene_cnt_per_DE[nrow(significant_gene_cnt_per_DE) + 1, ] <<- list(current_DE, as.numeric(nrow(filtered_below_threshold)))
  
}

preprocessed <- FALSE
loaded_data <- get_lab_DE_data(preprocessed)
#loaded_data <- get_own_DE_data(preprocessed)
mapply(process_pipeline, loaded_data[[1]], loaded_data [[2]])


#Plotting
ggplot(significant_gene_cnt_per_DE, aes(x=Label, y=Gene_Count, fill=Label)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = Gene_Count), vjust = 0) +
  scale_x_discrete(guide = guide_axis(angle = 90))+
  #labs(title="Count of significantly differentiated genes with padj < 0.01 for all DE", y = "Gene Count") +
  theme_classic()
ggsave("basic_exploration/plots/significant_changed_genes_count.png", width = 10, height = 6, dpi = 300, units = "in", device='png')
#ggsave("basic_exploration/plots/own/significant_changed_genes_count_original_labels.png", width = 10, height = 6, dpi = 300, units = "in", device='png')







