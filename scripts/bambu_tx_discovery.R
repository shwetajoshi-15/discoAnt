#!/usr/bin/env Rscript

suppressWarnings({
  
  suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(bambu)
  })
  options(dplyr.summarise.inform = FALSE)
  
  option_list = list(
    make_option(c("-b", "--bam"), type="character", default=NULL,
                help="merged bam file"),
    make_option(c("-f", "--fasta"), type="character", default=NULL,
                help="genome reference fasta"),
    make_option(c("-t", "--annotation"), type="character", default=NULL,
                help="reference annotations gtf"),
    make_option(c("-o", "--output_dir"), type="character", default=NULL,
                help="path to output directory")
    
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  if (is.null(opt$bam)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  }
  
  test.bam <- (opt$bam)
  fa.file <- (opt$fasta)
  gtf.file <- (opt$annotation)
  
  bambuAnnotations <- prepareAnnotations(gtf.file)
  
  # NOTES
  # NDR of 1 also means no false-discovery adjustment, do we really want to return low confidence isoforms?
  # min.readFractionByGene is extremely low. Currently outputs isoforms with 0.1% relative abundace, so 1 read in 1000.
  # have tried adjusting min.exonDistance, but dpesn't seem to do anything
  
  se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file, NDR=1, opt.discovery=list(min.readFractionByGene=0.001))
  
  writeBambuOutput(se, opt$output_dir)
  
  # only outputs one row?
  metadata_bambu <- data.frame(mcols(se))

  write.csv(metadata_bambu, file.path(opt$output_dir, "bambu_metadata.csv"))
  
})