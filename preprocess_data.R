library(stringr)
library(dplyr)
library(here)
here::i_am("preprocess_data.R")

#Depends on the utilities.R for loading the data

#GLOBAL VARIABLES
out_directory <- 'data/preprocessed/'
out_directory_own <- 'differential_expression/DE/preprocessed_original_labels/'
relevant_columns <- c("external_gene_name", "pvalue", "padj", "ensembl_gene_id", "log2FoldChange", "stat")

#FUNCTIONS

preprocess_data <- function(data){
  data %>%
    #Keep only the columns we need for the analysis
    select(contains(relevant_columns)) %>%
    #Remove the MT1 genes from the dataset, because they appear in the top differential genes and don't really have anything to do with our experiment 
    filter(!str_detect(external_gene_name, "^MT1"))
  
} 

process_pipeline <- function(data, DE_name){
  write.csv(preprocess_data(data), 
            file=paste0(out_directory_own, DE_name, ".csv"),
            row.names=FALSE
            )
}

#get the original data
preprocessed <- FALSE
#loaded_data <- get_lab_DE_data(preprocessed)
loaded_data <- get_own_DE_data(preprocessed)

mapply(process_pipeline, loaded_data[[1]], loaded_data [[2]])

