all_genes <- read.csv(file = paste0('pathways.csv'), sep =';')
all_genes<-subset(all_genes, group!="Stress" & group!="Differentiation")

genes_grouped <- all_genes %>% group_by(group)

#extract the grouping keys for plotting
group_keys <- group_keys(genes_grouped)
scaled_data <- get_scaled_expression_data()

#For each phase group, get the expression data and prepare it for plotting
list_of_result_dataframes <- genes_grouped %>% 
  group_map(~ extract_genes_and_prep_for_plotting(.,scaled_data))

#adding the group identifier
for (i in 1:length(list_of_result_dataframes)){
  group <- group_keys[[1]][i]
  list_of_result_dataframes[[i]]$group <- group
}

combined_list_of_result_dataframes <- do.call(rbind, list_of_result_dataframes)

#Add the average for each group
# combined_list_of_result_dataframes <- combined_list_of_result_dataframes %>%
#   group_by(group, timepoint) %>%
#   mutate(avg_expression = mean(expression))

#plot_specific_genes_on_timeline(combined_list_of_result_dataframes, 'avg_expression', 'pathways.png', TRUE, 'Pathways between t0 and t24', 'group', 'group')

genes_summary <- aggregate(expression ~ group + timepoint, data = combined_list_of_result_dataframes, FUN = function(x) c(mean = mean(x), sd = sd(x)))



ggplot(genes_summary, aes(x = timepoint, y = expression[, "mean"], group = group, color = group)) +
  #geom_line() +
  geom_smooth(size=0.5, se = FALSE) +
  geom_ribbon(aes(ymin = expression[, "mean"] - expression[, "sd"], ymax = expression[, "mean"] + expression[, "sd"], fill = group), alpha = 0.2, colour = NA) +
  labs(x = "Timepoint", y = "Gene expression") +
  theme_classic()

ggsave(paste0(results_dir_timeline, 'pathways-sd-smooth.png'), device = "png")
