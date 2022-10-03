#!/usr/bin/env Rscript
suppressMessages(library(optparse, warn.conflicts = F, quietly = T))
suppressMessages(library(dplyr, warn.conflicts = F, quietly = T))
suppressMessages(library(bambu, warn.conflicts = F, quietly = T))

# optimised for bambu v2.0.6

option_list = list(
 make_option(c("-b", "--bam"), type="character", default=NULL,
     help="merged bam file"),
 make_option(c("-f", "--fasta"), type="character", default=NULL,
     help="genome reference fasta"),
 make_option(c("-t", "--annotation"), type="character", default=NULL,
     help="reference annotations gtf"),
 make_option(c("-r", "--rfbg"), type="character", default=NULL,
     help="readFractionByGene generated based on total number of reads per barcode and read support for each transcript"),
 make_option(c("-o", "--output_dir"), type="character", default=NULL,
     help="path to output directory")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bam)){
  print_help(opt_parser)
  stop("Alignments (BAM) must be provided", call.=FALSE)
}

if (is.null(opt$rfbg)){
  print_help(opt_parser)
  stop("readFractionByGene value not found", call.=FALSE)
}

test.bam <- (opt$bam)
fa.file <- (opt$fasta)
gtf.file <- (opt$annotation)

bambuAnnotations <- prepareAnnotations(gtf.file)

se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file, opt.discovery = list(max.txNDR=1, min.readfractionByGene=opt$rfbg, fitReadClassModel = FALSE))

writeBambuOutput(se, opt$output_dir)

UC <- assays(se)$uniqueCounts
CPM <- assays(se)$CPM

write.table(UC, file.path(opt$output_dir, "uniqueCounts.txt"), sep = "\t")
write.table(CPM, file.path(opt$output_dir, "CPM.txt"), sep = "\t")

