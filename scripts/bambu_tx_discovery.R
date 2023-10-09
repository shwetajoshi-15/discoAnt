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
  
  se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file, opt.discovery = list(max.txNDR=1, min.readFractionByGene=0.001))
  #se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file, NDR=1, opt.discovery=list(min.readFractionByGene=0.001))
  
  
  writeBambuOutput(se, opt$output_dir)
  
  #UC <- assays(se)$uniqueCounts
  #CPM <- assays(se)$CPM
  RC <- assays(se)$counts
  
  #write.table(UC, file.path(opt$output_dir, "uniqueCounts.txt"), sep = "\t")
  #write.table(CPM, file.path(opt$output_dir, "CPM.txt"), sep = "\t")
  write.table(RC, file.path(opt$output_dir, "bambu_read_counts.txt"), sep = "\t")
  
})