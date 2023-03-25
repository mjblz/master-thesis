#Overrepresentation analysis


library <- 'enrichGO'

plot_enrichment <- function(enriched, enriched_directional, label, ontology){
  
  plots <- list()
  plots[["sig"]] <- dotplot(enriched) + labs(title = "All significant genes")
  plots[["dir"]] <- dotplot(enriched_directional) + labs(title = "DE + direction of change")
  grid <- plot_grid(plotlist = plots)
  ggsave(plot=grid, paste0(plot_base_dir, "overrepresentation/", library, "/dotplot/", label , "/", ontology, ".png"), scaling = 0.5) 
  
  
}


get_common_sets <- function(enriched, enriched_directional, label, ontology){
  
  sign_enriched <- unlist(subset(enriched@result, p.adjust < 0.01, select = "Description"))
  sign_enriched_dir <- unlist(subset(enriched_directional@compareClusterResult, p.adjust < 0.01, select = "Description"))
  
  common_set <- intersect(sign_enriched_dir,sign_enriched)
  common_set_df <- as.data.frame(common_set)
  colnames(common_set_df) <- c(label)
  write.csv(common_set_df, paste0("enrichment_analysis/gene_sets/enrichGO/", ontology, "/", label, ".csv"), row.names=FALSE)
  
  return(common_set) 
  
}

run_enrichGO <- function(diff_expression, label){
  
  combined_common_sets <- list()
  
  #Get the direction
  diff_expression$direction <- NA
  diff_expression$direction[diff_expression$padj<0.01 & diff_expression$log2FoldChange > 0] <- "up"
  diff_expression$direction[diff_expression$padj<0.01 & diff_expression$log2FoldChange < 0] <- "down"
  
  #Subset to significant genes
  sig_genes <- subset(diff_expression, !is.na(direction))
  
  #Select only the genes actually detected in the assay to use for the analysis universe
  uni <- unique(na.omit(diff_expression$external_gene_name[ !is.na(diff_expression$log2FoldChange) ]))
  
  for (ontology in ontologies){
    #For all significantly changed genes
    enriched <- enrichGO(
      gene = sig_genes$external_gene_name,
      OrgDb         = organism,
      keyType       = 'SYMBOL',
      universe = uni,
      ont           = ontology,
      readable      = TRUE)
    
    #Directional Analysis
    enriched_directional <- compareCluster(
      external_gene_name ~ direction, 
      data = sig_genes, 
      universe = uni, 
      OrgDb = organism, 
      keyType = "SYMBOL", 
      ont = ontology, 
      fun = "enrichGO"
    )
    
    #Get the common significantly enriched sets of those two algorithms
   
    common_set <- get_common_sets(enriched, enriched_directional,label, ontology)
    combined_common_sets <- append(combined_common_sets, common_set)
    
    
    #plot_enrichment(enriched, enriched_directional,label, ontology)
 
  } #End of for loop
 
  combined_common_sets_df <- unique(as.data.frame(unlist(combined_common_sets)))
  colnames(combined_common_sets_df) <- c(label)
  write.csv(combined_common_sets_df, paste0("enrichment_analysis/gene_sets/enrichGO/", label, "_all_ont.csv"), row.names=FALSE)
}




preprocessed <- TRUE
loaded_data <- get_data(preprocessed)
mapply(run_enrichGO, loaded_data[[1]], loaded_data [[2]])

