# Load required libraries:
library(tinytex)
# tinytex::install_tinytex()
library(DESeq2)
library(here)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
library(dplyr)
library(magrittr)
library(ggplot2)


## VARIABLES
plot_dir <- 'differential_expression/plots/'
output_files_dir <- 'differential_expression/DE/original_labels/'
#output_files_dir <- 'differential_expression/DE/additional_labels/'


# Load table of counts (expression abundances for all experimental samples):
count_matrix <- read.csv("differential_expression/raw_counts.csv", row.names = 1, sep=";" )

#head(count_matrix)
dim(count_matrix)

# Load the annotation file:
annotations <-  read.csv("differential_expression/annot.csv", row.names = 1, sep=";")

#Preprocess it
#annotations <- annotations[c("Condition","Group")]
annotations <- annotations[c("Condition","Group","Batch")]

# By default, R will choose a reference level for factors based on alphabetical order. 
# Then, if you never tell the DESeq2 functions which level you want to compare against (e.g. which level represents the control group), 
# the comparisons will be based on the alphabetical order of the levels.
#annotations$Condition <- as.factor(annotations$Condition)
#annotations$Group <- as.factor(annotations$Group)
annotations$Batch <- as.factor(annotations$Batch)
annotations$Combined <- factor(paste0(annotations$Group, "_", annotations$Condition))
summary(annotations)

all(row.names(annotations) == colnames(count_matrix))
# If you get TRUE, they are in the same order, otherwise you need to reorder the tables to match each other. 

#dds" is short for "DESeqDataSet

# # Create DESeq2 input object (DESeqDataSet):
dds <-  DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = annotations,
                               design = ~ Batch + Combined)

# # The DESeqDataSet is an R object to store the read counts and the intermediate results based on your analysis.
# # counts() function can be used the expression levels from deseqmat object:
#
head(counts(dds))
#
# # FILTER GENES WITH LOW EXPRESSION:
#
# # In this step, we try to remove the genes with low expression level across all samples.
# # This will help reducing the data matrix dimension and improving the analysis running time.
low_genes <-  rowSums(counts(dds)) < 5
#
# Here, counts() function provides access to table of count from DESeqDataSet object.
# This way, we find the genes that include less than 5 reads for all experimental samples.

# Check how many of the genes have less than 5 reads in total:
table(low_genes)

# #Filter them out
dds <-  dds[!low_genes,]

#
# #QUALITY CHECK
# # In order to do so, the data needs to be transformed (check the variance stablization using DESeq2) and normalized:
var_stablized <-  vst(dds)
#

#Plotting
categories <- c("Condition", "Group", "Combined")
lapply(X = categories, FUN = function(x) {
  DESeq2::plotPCA(var_stablized, intgroup = x, ntop=10000) + geom_text(aes(label=name),vjust=2,size=4)
  ggsave(paste0(plot_dir, "PCA_", x, ".png"), scaling = 0.7)
})

# #
#
# # DIFFERENTIAL EXPRESSION TESTING
#
# # In DESeq2, a single function (DESeq()) is responsible for performing normalization, variance stablization, model fitting and testing.
# # Of course the are also separate functions to do these tasks as well.
#

res_de <- DESeq(dds)

#Get all possible combinations
#In case I should ever need all
# all_combined <- as.character(unique(annotations$Combined))
# combinations <- all_combined %>%
#   expand.grid(.,.) %>%
#   distinct() %>%
#   filter(Var1 != Var2)


#Recreate the DE we got from the lab
original <- TRUE
all_labels <- get_all_labels(original)
signif_genes_count <- data.frame(comparison = character(0), sign_count_padj_0.01 = numeric(0))

for (label in all_labels) {
  split_res <- unlist(strsplit(label, "_vs_"))
  comb1 <- split_res[1]
  comb2 <- split_res[2]
  comp <- paste0(comb1, '_vs_', comb2)
  de_comparison <- results(res_de, contrast=c("Combined", comb1, comb2))
  de_comparison$external_gene_name <- rownames(de_comparison)
  #write.csv(de_comparison, file = paste0(output_files_dir, comp, '.csv'))
  
  #add the count to a dataframe
  signif_genes_count <- structure(
    rbind(
      signif_genes_count, c(comp,sign_count_padj_0.01=nrow(subset(de_comparison, padj < 0.01)))
    ),
    .Names = names(signif_genes_count)
  )
  
  # # Lets assume you are interested to check the expression levels for 20 first genes with highest positive log fold change.
  # # We can check this using heatmaps.
  # # First lets find the target genes by sortung the DESeq2 results matrix based on the log fold change:
  ordered_res_foldchange <-  de_comparison[order(de_comparison$log2FoldChange, decreasing = TRUE),]
  
  #
  # # Then select the first 20 genes:
  top_genes_max_pos_foldchange <-  row.names(ordered_res_foldchange)[1:20]
  
  top_genes_max_neg_foldchange <- rev(row.names(ordered_res_foldchange))[1:20]
  
  
  # # Now we can plot the heatmap for the selected genes using the normalized and transformed expression levels:
  pheatmap(assay(var_stablized)[top_genes_max_pos_foldchange,], cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=FALSE,
           filename= paste0(plot_dir, "top20genes/Top20Genes_PositiveFC_",comp,".png"),
           width=50,
           height=20,
           fontsize = 30
  )
  # #
  pheatmap(assay(var_stablized)[top_genes_max_neg_foldchange,], cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=FALSE,
           filename= paste0(plot_dir, "top20genes/Top20Genes_NegativeFC_",comp,".png"),
           width=50,
           height=20,
           fontsize=30
  )
  
}

write.csv(signif_genes_count, file = paste0(plot_dir, 'significant_genes_count_padj0.01.csv'))
###########################################################