#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

options(dplyr.summarise.inform = FALSE)

option_list = list(
  make_option(c("-e", "--ens_var"), type="character", default=NULL,
              help="ENS_ID variable"),
  make_option(c("-m", "--read_min"), type="character", default=NULL,
              help="Read count threshold"),
  make_option(c("-o", "--output_path"), type="character", default=NULL,
              help="output path")
)

suppressWarnings({

  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)

  gene_id <- (opt$ens_var)
  read_count_min <- as.numeric(opt$read_min)
  outdir <- (opt$output_path)

  # get sample ids from within directory
  sample_names <- c(list.files(paste0(outdir,"/salmon_quants/")))

  import_quants_function <- function(sample_id) {
    # get paths of file to import
    pathtofile <- paste0(outdir,"/salmon_quants/", sample_id, "/quant.sf")
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

  combined_counts <- combined_counts %>% 
      mutate(total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
      dplyr::filter(total > read_count_min) %>%
      select(-total)

  # older Bambu versions name novel isoforms with "BambuTx' whereas new versions name them 'tx'
  combined_counts$TXNAME <- gsub('BambuTx', 'tx', combined_counts$TXNAME)

  # export file
  write.csv(combined_counts, paste0(outdir, "/", gene_id, "_discoAnt_counts.csv"), quote = FALSE, row.names = FALSE)

})