library(reshape2)
#COMPARE AMOUNT OF COMMON GENES BETWEEN different DE

# Depends on the output of SignificantChangeGeneCount.R, specifically the significant_genes_per_DE (hash: DE - significant_genes_per_DE) and significant_gene_cnt_per_DE (DE - count)
# This script calculcates the amount of significantly changed genes two DE comparisons have in common. For example t0_WT_vs_to_K0 and t0_WT_supp_vs_t0_KO 


#Get all the labels
all_DE_labels <- get_all_labels()

#Prepare result dataframe
#Each DE will be compared with each other DE and the amount of common genes will be counted
DE_vs_DE_common_genes <- data.frame(matrix(NA, nrow = length(all_DE_labels), ncol = length(all_DE_labels)))
colnames(DE_vs_DE_common_genes) <- all_DE_labels
row.names(DE_vs_DE_common_genes) <- all_DE_labels



for(i in 1:nrow(DE_vs_DE_common_genes)) {       
  for(j in 1:ncol(DE_vs_DE_common_genes)) {  
    
    #Extract the names of the current comparison to be done
    row_name <- rownames(DE_vs_DE_common_genes)[i]
    col_name <- colnames(DE_vs_DE_common_genes)[j]
    
    #Get gene lists from hash and do pairwise comparison of differences
    if (!row_name == col_name){
      #Get the list of the significantly changed genes for the DE comparison
      comp1 <- unlist(significant_genes_per_DE[[row_name]])
      comp2 <- unlist(significant_genes_per_DE[[col_name]])
      #Count the common genes between two and set the value in the df
      DE_vs_DE_common_genes[i,j] <- length(comp1[comp1 %in% comp2])
    }
    
    
  }
}

#Rewrite to percent
#Make a copy which will have percent calculations instead of counts
DE_vs_DE_common_genes_perc <-DE_vs_DE_common_genes

# Calculate which percent of the total genes they have in common
get_percentage <- function(common_gene_count,name_of_comparison){
  total_gene_count <- significant_gene_cnt_per_DE[significant_gene_cnt_per_DE$Label == name_of_comparison,]$Gene_Count
  (common_gene_count/total_gene_count)*100
  #Get the total count for that specific column/comparison

}

DE_vs_DE_common_genes_perc <- mapply(FUN=get_percentage, common_gene_count=DE_vs_DE_common_genes_perc, name_of_comparison=names(DE_vs_DE_common_genes_perc))
row.names(DE_vs_DE_common_genes_perc) <- all_DE_labels
DE_vs_DE_common_genes_perc[is.na(DE_vs_DE_common_genes_perc)] <- 100

#Plot Heatmap
melted_matrix <- melt(DE_vs_DE_common_genes_perc)
ggplot(data = melted_matrix, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  #labs(title="Heatmap of percentage of genes in common across DE comparisons", x ="", y = "") +
  scale_x_discrete(guide = guide_axis(angle = 90))
ggsave("basic_exploration/plots/heatmap_common_genes.png", width = 10, height = 5, dpi = 300, units = "in", device='png')



