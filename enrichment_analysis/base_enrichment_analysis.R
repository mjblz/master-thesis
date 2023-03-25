library(clusterProfiler)
library(enrichplot)
library("readxl")
library(ggplot2)
require(DOSE)
library(ragg)
library(tryCatchLog)
library(fgsea)

library("cowplot")

# HUMAN Organism

organism = "org.Hs.eg.db"
ontologies <- c("MF", "BP","CC")
plot_base_dir <- "enrichment_analysis/plots/"

#Installs the Human annotations
#BiocManager::install(organism, character.only = TRUE)

library(organism, character.only = TRUE)

#GSE Analysis

#Sorted: Should the gene list be sorted, if yes it is sort in a decreasing order
#Identifier: Which column should be used for the names of resulting gene list, this can either be ensemble Id (ensembl_gene_id) or Gene Name (external_gene_name)
#Ranking column: Column to use for the actual value that will be compared (either pvalue or log2foldchange make sense)
# Reference: https://github.com/YuLab-SMU/DOSE/wiki/how-to-prepare-your-own-geneList
prepare_gene_list <- function(diff_expression_data, sorted, ranking_column, identifier){
  
  gene_list <- diff_expression_data[[ranking_column]]
  names(gene_list) <- diff_expression_data[[identifier]]
  
  # omit any NA values 
  gene_list<-na.omit(gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  if (sorted){
    return(sort(gene_list, decreasing = TRUE))
  }
  return(gene_list)
}