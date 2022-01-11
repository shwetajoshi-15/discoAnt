# discoAnt
- version 10012022
- The repository is under-construction.

## Prepare FASTA files in a folder
All the sample files should be in FASTA format (.fa)

## Setting up discoAnt

1. git clone repository
2. cd discoAnt
3. Update discoAnt_params.txt with - Gene name, start and end coordinates and strands
4. Update discoAnt_params.txt with - path to FASTA folder and discoAnt folder
5. When running the pipeline for the first time - ```bash discoAnt_setup.sh```
6. ```conda activate discoAnt.env```
7. ```cd /path/to/folder/discoAnt/programs/cDNA_Cupcake```
8. ```python setup.py build```
9. ```python setup.py install```

#### Refer to the sample submission scripts for running the pipeline on a HPC cluster.

## What does discoAnt do?

1. Quality control for the fasta files  
1.a. Generates a file with no. of reads  \
1.b. average read length  \
1.c. read length in each barcode/sample  
2. Align the sample fasta files to reference genome
3. Filter out off-target alignments (based on the gene start and end site provided)
4. Merge alignments from all samples
5. Correct (based on reference transcript annotations) and collapse transcripts
6. Classify the transcripts based on known annotations 
7. Align the samples fasta files to the metagene 
8. Generate counts/TPM based on the alignments to the metagene



