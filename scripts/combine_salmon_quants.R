# run script from directory that contains 'metagene_salmon/' directory
# script requires gene_id to add as a suffix for TXNAME
# run as: Rscript: combine_salmon_quants.R [gene_id]
options(dplyr.summarise.inform = FALSE)

args <- commandArgs(trailingOnly = TRUE)
gene_id <- args[1]


suppressWarnings({
# get sample ids from within directory
sample_names <- c(list.files("metagene_salmon"))
#print(sample_names)

import_quants_function <- function(sample_id) {
  # get paths of file to import
  pathtofile <- paste0("metagene_salmon/", sample_id, "/quant.sf")
  df <- read.table(pathtofile, header=T)
  # keep only txname and read counts columns
  df <- df[, c("Name", "NumReads")]
  # replace all "." characters with a "_" in case IsoVis breaks
  sample_id <- gsub("\\.","_",sample_id)
  # rename cols to TXNAME and the sample id
  colnames(df) <- c("TXNAME", paste0(sample_id))
  rownames(df) <- df[,1]
  df[,1] <- NULL

  return(df)
}

# import all files and combine dfs
all_samples <- lapply(sample_names, import_quants_function)
combined_counts <- do.call("cbind", all_samples)

# remove version numbers
rownames(combined_counts) <- gsub("\\..*","", rownames(combined_counts))
# add gene id suffix
rownames(combined_counts) <- paste0(rownames(combined_counts), "_", gene_id)
# get order of columns
new_col_order <- c("TXNAME", paste0(colnames(combined_counts)))
# convert row names to column called TXNAME
combined_counts$TXNAME <- rownames(combined_counts)
# order df
combined_counts <- combined_counts[, new_col_order]
# remove rownames
#rownames(combined_counts) <- NULL



# export file
# change name
write.csv(combined_counts, paste0(gene_id, "_discoAnt_counts.csv"), quote = FALSE, row.names = FALSE)
#write.csv(combined_counts, paste0(out_path, "/", gene_id, "_discoAnt_counts.csv"), quote = FALSE, row.names = FALSE)

})