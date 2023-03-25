#This script runs gene set enrichment on all the differential expression files, they need to contain gene names and log2 fold change values.
#This data is typically produced by differential expression analysis tool such as DESeq 2.

#https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/


#GLOBAL VARIABLES
log.path <- 'enrichment_analysis/gse_log.txt'

library <- 'gseGO'

#FUNCTIONS


#Plotting
plot_enrichment <- function(enrichment, DE_label){

  write(paste('Plotting enrichment for ', DE_label), log.path, append=TRUE)
  
  # Dotplot
  dotplot(enrichment, showCategory=10, split=".sign", font.size=10) + facet_grid(.~.sign)+
    labs(title=paste("Enriched Pathways for differential expression",DE_label))
  ggsave(paste0(plot_base_dir, "gse/", library, "/dotplot/", paste0('p', config_pvalueCutoff), "/", DE_label ,"_dotplot.png"), scaling = 0.5, device = agg_png) 
  
}

build_error_msg <- function(DE_name, pvaluecutoff, errormsg){
  paste0(" Error in ", DE_name, " for pvalue cutoff at ", pvaluecutoff, ": ", toString(errormsg))
}

run_gsea <- function(diff_expression, label){
  
  tryCatchLog({
    sorted <- TRUE
    gene_list <- prepare_gene_list(diff_expression, sorted, 'log2FoldChange', 'ensembl_gene_id')
    
    gseGO_res <- gseGO(geneList=gene_list,
                       ont ="ALL", #BP", "MF", and "CC", Biological Process, Molecular Function, Cellular Component
                       keyType = "ENSEMBL",
                       maxGSSize = 800,
                       pvalueCutoff = 0.01,
                       verbose = TRUE,
                       OrgDb = organism)
    
    #Plot
    #plot_enrichment(gseGO_res, DE_name)
    
    #Extract the ids of the significantly enriched gene sets
    sign_enriched <- subset(gseGO_res@result, p.adjust < 0.01, select = "Description")
    colnames(sign_enriched) <- c(label)
    write.csv(sign_enriched, paste0("enrichment_analysis/gene_sets/gseGO/", label, "_all_ont.csv"), row.names=FALSE)
    # 
    #Trigger garbage collection
    gc()
    
  },
  error = function(errormsg) {
    msg <- build_error_msg(label, config_pvalueCutoff, errormsg )
    write(msg, log.path, append=TRUE)
  }
  )
}

#ACTUAL RUN
preprocessed <- TRUE
loaded_data <- get_data(preprocessed)
mapply(run_gsea, loaded_data[[1]], loaded_data [[2]])


########## ----------- ################

#Running enrichment only for one specific DE comparison (specifically for Ro), uncomment to use

########## ----------- ################

#This loads all the DE matrices

# preprocessed <- TRUE
# loaded_data <- get_data(preprocessed)
# # 
# # #Select which one you like by looking at the order of the available ones:
# # loaded_data[[2]]
# # 
# # #And select the index of the you like
# # #13 for example is t24_KO_vs_t0_KO
# DE_index <- 13
# 
# diff_expression <- loaded_data[[1]][DE_index][[1]]
# 
# # Just in case you want the label as well
# label <- loaded_data[[2]][DE_index][[1]]
# 
# gene_list <- prepare_gene_list(diff_expression, TRUE, 'log2FoldChange', 'ensembl_gene_id')
# 
# #This does the actual enrichment
# gse <- gseGO(geneList=gene_list,
#              ont ="ALL", #BP", "MF", and "CC", Biological Process, Molecular Function, Cellular Component
#              keyType = "ENSEMBL",
#              maxGSSize = 800,
#              pvalueCutoff = 0.01,
#              verbose = TRUE,
#              OrgDb = organism) #Organism needs to be run from the start of the script
# 
# #This just changes the ENSEMBL IDs to readable names
# readable_enrichment <- setReadable(gse, 'org.Hs.eg.db', 'ENSEMBL')
# #The core_enrichment contains each of the gene sets that was enriched
# readable_enrichment@result$core_enrichment
# 

