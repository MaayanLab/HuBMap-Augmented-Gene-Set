
library(tidyverse)
library(edgeR)
library(qpgraph)
library(DESeq2)

args <- commandArgs(TRUE)

# args[1] is directory of files to CTF normalize
directory <- args[1]

print("start")
#create output dir
output_dir <- paste0(directory, "_CTF_normalized")
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List files in dir to normalize
dir_files <- list.files(directory, full.names = TRUE)

# Loop through files to normalize
for (filename in dir_files){
    # Read in file, convert to matrix for norm factor function
    count_data <- read_delim(filename, 
                            delim = "\t", 
                            col_names = T, 
                            col_types = cols(.default = "d", gene = "c")) %>% 
    column_to_rownames("gene") %>% 
    as.matrix()
    # Get file basename
    file_basename <- basename(filename)
    print(paste0("Working on ", file_basename))
    # Library size = col sums of matrix
    lib_size <- base::colSums(count_data)
    # TMM normalization factors
    norm_factors <- calcNormFactors(object = count_data, lib.size = lib_size, method = "TMM")
    # Divide with norm factors here
    CTF_normalized <- sweep(count_data, 2, norm_factors, "/")
    # Convert Matrix to Data Frame
    CTF_normalized %>% 
      as.data.frame() -> CTF_normalized_df 
    # Log Transform using asinh
    asinh(CTF_normalized_df) -> CTF_normalized_log_df
    # Calculate correlation
    CTF_normalized_log_t_df <- t(CTF_normalized_log_df)
    rownames(CTF_normalized_log_t_df) <- NULL
    print(paste0("Calculating correlation on ", file_basename))
    CTF_normalized_log_t_cor_matrix <- cor(CTF_normalized_log_t_df)
    # Convert to dataframe
    CTF_normalized_log_t_cor_df <- as.data.frame(CTF_normalized_log_t_cor_matrix)
    # Write to data frame
    print(paste0("Writing ", file_basename))
    write_tsv(CTF_normalized_log_t_cor_df,
          paste0(file_basename),
          col_names = T)
}