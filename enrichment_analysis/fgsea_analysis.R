library(fgsea)
library(enrichplot)
library(ReactomePA)
library(tibble)
library(stringr)
library(dplyr)


# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# https://stephenturner.github.io/deseq-to-fgsea/ DEseq to fgsea
# Which rank to use: https://www.biostars.org/p/9526168/#9526475

run_fgsea <- function(diff_expression, label, pathway_db){
  # ranks are from lowest to highest
  #Using the stat column for ranking
  sorted <- TRUE
  gene_list <- prepare_gene_list(diff_expression, sorted, 'stat', 'external_gene_name')
  pathway_db <- gmtPathways('enrichment_analysis/go_symbols.gmt')
  
  fgsea_res <- fgsea(pathways = pathway_db, 
                       stats = gene_list)
    
  #Extract the ids of the significantly enriched gene sets
  sign_enriched <- subset(fgsea_res, padj < 0.01, select = "pathway")
  colnames(sign_enriched) <- c(label)
  sign_enriched <- sign_enriched %>% 
    mutate_all(str_replace_all, "^GO.._", "") %>% 
    mutate_all(str_replace_all, "_", " ") %>% 
    mutate_all(funs(tolower))
  write.csv(sign_enriched, paste0("enrichment_analysis/gene_sets/fgsea/", label, "_all_ont.csv"), row.names=FALSE)
    # 
    
}

#ACTUAL RUN
preprocessed <- TRUE
loaded_data <- get_data(preprocessed)
#GO gene sets from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
#pathway_db <- gmtPathways('enrichment_analysis/go_symbols.gmt')
mapply(run_fgsea, loaded_data[[1]], loaded_data [[2]])



# 
# preprocessed <- TRUE
# loaded_data <- get_data(preprocessed)
# 
# DE_index <- 13
# 
# diff_expression <- loaded_data[[1]][DE_index][[1]]
# 
# # Just in case you want the label as well
# label <- loaded_data[[2]][DE_index][[1]]
# 

# gene_list <- prepare_gene_list(diff_expression, TRUE, 'stat', 'external_gene_name')







# topPathwaysUp <- signif_results[ES > 0][head(order(pval), n=10), pathway]
# topPathwaysDown <- signif_results[ES < 0][head(order(pval), n=10), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# plotGseaTable(pathway[topPathways], gene_list, fgsea_res, 
#               gseaParam=0.5)

#nes: The Normalized Enrichment Score (NES) measures the degree of enrichment of the input gene list for each term. 
# NES is a signed value, positive values indicate that the term is enriched in the input gene list, and negative values indicate that the term is depleted in the input gene list.
