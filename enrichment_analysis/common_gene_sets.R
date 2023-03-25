
input_dir <- c('enrichment_analysis/gene_sets/enrichGO/', 'enrichment_analysis/gene_sets/fgsea/',  'enrichment_analysis/gene_sets/gseGO/')

all_labels <- get_all_labels()

#Get the gene sets for all different algorithms and find which they have in common
for (label in all_labels){
  
  gene_sets <- list()
  for (dir in input_dir){
    filepath <- paste0(dir,label,'_all_ont.csv')
    gene_set <- list(read.csv(file = filepath))
    if (dim(gene_set[[1]])[1] != 0) {
      gene_sets <- append(gene_sets, gene_set)
    }
  }
  
  common_gene_sets <- Reduce(intersect, gene_sets)
  write.csv(common_gene_sets, paste0("enrichment_analysis/gene_sets/common/", label, ".csv"), row.names=FALSE)
}