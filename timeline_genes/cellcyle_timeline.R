#For All Cell Cycle Genes
all_genes <- read.csv(file = paste0('CC_checkpoint_genes.csv'), sep =';')

# group the dataframe by the "Phase" column
genes_grouped <- all_genes %>% group_by(Phase)

#extract the grouping keys for plotting
group_keys <- group_keys(genes_grouped)

scaled_data <- get_scaled_expression_data()

#For each phase group, get the expression data and prepare it for plotting
list_of_result_dataframes <- genes_grouped %>% group_map(~ extract_genes_and_prep_for_plotting(.,scaled_data))

for (i in 1:length(list_of_result_dataframes)){
  phase <- group_keys[[1]][i]
  plot_specific_genes_on_timeline(list_of_result_dataframes[[i]], 'expression', paste0(phase, '.png'), FALSE, phase, 'gene', 'gene')
}


##Redo the same thing but with clustering 