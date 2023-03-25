#Create one big matrix which has the significant fold changes, padj of all the genes and all the DE

preprocessed <- TRUE
loaded_data <- get_data(preprocessed)
#
data_list <- loaded_data[[1]]
labels <- loaded_data[[2]]
#
preprocessed_df <- Map(
  function(data, label)
    data %>%
    dplyr::select(external_gene_name, log2FoldChange, padj) %>%  #Select only fold change
    rename(!!paste0('log2FoldChange.', label) := log2FoldChange) %>%
    rename(!!paste0('padj.', label) := padj),
  data_list,
  labels
) #add the information to which DE it belongs

# #Merge into one
merged <- preprocessed_df %>%
  reduce(left_join, by = "external_gene_name") %>%
  na.omit(.) %>%
  data.frame(., row.names = 1)

write.csv(merged, "clustering/data/fold_changes.csv", row.names = TRUE)

fold_changes <- read.csv(file = 'clustering/data/fold_changes.csv'
                         , row.names = 1)
