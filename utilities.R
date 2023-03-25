library(here)
library(stringr)
library("readxl")
here::i_am("utilities.R")

data_directory_original <- 'data/original'
data_directory_own_original <- 'differential_expression/DE/original_labels'
data_directory_preprocessed <- 'data/preprocessed'
data_directory_own_preprocessed <- 'differential_expression/DE/preprocessed_original_labels'



get_scaled_count_data_for_group_of_interest <- function(group_of_interest){
  counts <- read.csv(file = 'clustering/data/normalized_counts.csv', row.names = 1, sep=';', dec=',')
  relevant_columns <- get_relevant_columns(group_of_interest, counts)
  
  #Get all available timepoints for that group of interest
  timepoints <- unique(stringr::str_extract(names(relevant_columns), "(?<=KO_t.)[0-9]+"))
  
  #collapse replicates
  collapsed_replicates <- collapse_replicates(timepoints, group_of_interest, relevant_columns)
  
  #scale data
  scaled_data <- as.data.frame(scale_data(collapsed_replicates))
  return(scaled_data)
}

#Select only those columns that belong to that group
get_relevant_columns <- function(group_of_interest, raw_counts){
  selected_columns <- raw_counts
  selected_columns <- selected_columns %>%  
    dplyr::select(contains(group_of_interest, ignore.case = TRUE))
  return(selected_columns)
}

# there are technical duplicates for each time point. for our purpose
# # we just need one value per time point. It is sensible to compute the 
# # mean value for each gene at each time point.
collapse_replicates <- function(timepoints, group_of_interest, data){
  #Collapse all columns that belong to the same replicates into one 
  for (timepoint in timepoints){
    label <- paste0(group_of_interest,timepoint)
    
    #Get the column names of the replicates for one timepoint
    columns_to_group <- grep(label, names(data), value = TRUE)
    
    #calculate the mean over the rows
    data[label] <- rowMeans(data[columns_to_group])
    
    #Remove the columns that we are currently collapsing and leave everything else
    data <- data[,!(names(data) %in% columns_to_group)]
  }
  return(data)
}


##This function returns a dataframe which has a fixed color for each of the comparisons 
get_colors_for_comparisons <- function(){
  my_colors <- c("red", "blue", "green", "yellow", "orange", "purple")  
}

extract_DE_name_from_filepath <- function(filepath, preprocessed){
  if(preprocessed){
    #Naming pattern is different for the preprocessed files
    str_match(filepath, ".*preprocessed/\\s*(.*?)\\s*.csv")[,2]
  }else{
    str_match(filepath, ".*DE_Results_\\s*(.*?)\\s*.xlsx")[,2]
  }
  
}

extract_own_DE_name_from_filepath <- function(filepath, preprocessed){
  if(preprocessed){
    #Naming pattern is different for the preprocessed files
    str_match(filepath, ".*preprocessed_original_labels/\\s*(.*?)\\s*.csv")[,2]
  }else{
    str_match(filepath, ".*original_labels/\\s*(.*?)\\s*.csv")[,2]
  }
  
}

additional_labels <- list("t10_KO_vs_t4_KO",
                          "t16_KO_vs_t10_KO",
                          "t24_KO_vs_t16_KO"
                          )

lab_labels <- list("t0_WT_supp_vs_t0_KO",
                   "t0_WT_vs_t0_KO",
                   "t0_WT_vs_t0_WT_supp" ,
                   "t10_KO_supp_vs_t0_KO",
                   "t10_KO_vs_t0_KO",
                   "t16_KO_supp_vs_t0_KO",
                   "t16_KO_supp_vs_t16_KO",
                   "t16_KO_vs_t0_KO",
                   "t16_WT_vs_t0_WT",
                   "t16_WT_vs_t0_WT_supp",
                   "t16_WT_vs_t16_KO",
                   "t16_WT_vs_t16_KO_supp",
                   "t24_KO_vs_t0_KO",
                   "t4_KO_supp_vs_t0_KO", 
                   "t4_KO_vs_t0_KO")

#This function returns all the DE comparisons we have as a list
get_all_labels <- function(){
  return(c(lab_labels)) 
}

# Returns a list containing a list the loaded dataframes and a list of the DE names in the same order
#
# Data can be accessed like this:
# loaded_data <- get_data(preprocessed)
# dataframes <- loaded_data[[1]]
# DE_names <- loaded_data [[2]]

get_lab_DE_data <- function(preprocessed){
  
  input_directory <- NULL
  pattern <- NULL
  
  #Use preprocessed or original data
  if (preprocessed) {
    input_directory <- data_directory_preprocessed
    filetype <- '.csv'
    assign("read_function", read.csv)
  } else {
    input_directory <- data_directory_original
    filetype <- '.xlsx'
    assign("read_function", read_excel)
  }
  
  #Retrieve file paths from the right directory
  filepaths <- list.files(path=input_directory, pattern=filetype, all.files=TRUE,
                         full.names=TRUE)
  print(filepaths)
  
  DE_names <- lapply(filepaths, extract_DE_name_from_filepath, preprocessed=preprocessed)
  data <- lapply(filepaths, read_function)
  return(list(data, DE_names))
}

get_own_DE_data <- function(preprocessed){
  input_directory <- NULL
  #Use preprocessed or original data
  if (preprocessed) {
    input_directory <- data_directory_own_preprocessed
  } else {
    input_directory <- data_directory_own_original
  }
  
  #Retrieve file paths from the right directory
  filepaths <- list.files(path=input_directory, pattern='csv', all.files=TRUE,
                          full.names=TRUE)
  print(filepaths)
  
  DE_names <- lapply(filepaths, extract_own_DE_name_from_filepath, preprocessed=preprocessed)
  data <- lapply(filepaths, read.csv)
  return(list(data, DE_names))
}
