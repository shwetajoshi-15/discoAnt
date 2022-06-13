# discoAnt
- version 04052022
- The repository is under-construction.

## Prepare FASTA files in a folder
All the sample files should be in FASTA format (.fa)

## Setting up discoAnt

1. git clone repository
2. cd discoAnt
3. Update discoAnt_params.txt with - Gene name, start and end coordinates and strands
4. Update discoAnt_params.txt with - path to FASTA folder and discoAnt folder
5. When running the pipeline for the first time - ```bash discoAnt_setup.sh```

#### Refer to the sample submission scripts for running the pipeline on a HPC cluster.

## What does discoAnt do?

1. Quality control for the fasta files (Under construction)  
1.a. Generates a file with no. of reads  \
1.b. average read length  \
1.c. read length in each barcode/sample  
2. Align the sample fasta files to reference genome
3. Filter out off-target alignments (based on the gene start and end site provided)
4. Merge alignments from all samples
5. Correct and collapse transcripts
6. Align reads to the metagene and generate counts/TPM  
7. Generate stacked bar plots representing transcript usage



