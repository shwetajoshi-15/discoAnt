# discoAnt
The repository is under-construction.

## Prepare FASTA files in a folder
All the sample files should be in FASTA format (.fa)

## Setting up discoAnt

1. git clone repository
2. cd discoAnt
3. Update discoAnt_params.txt with - Gene name, start and end coordinates and strands
4. Update discoAnt_params.txt with - path to FASTA folder and discoAnt folder
5. ```conda env create -f discoAnt.env```
6. ```conda activate discoAnt.env```
7. When running the pipeline for the first time - ```bash discoAnt_setup.sh```
8. cd /path/to/folder/discoAnt/programs/cDNA_Cupcake
9. python setup.py build
10. python setup.py install

#### Refer to the sample submission scripts for running the pipeline on a HPC cluster.

## What does discoAnt do?

1. Quality control for the fasta files 
1.a. Generates a file with no. of reads
1.b. average read length
1.c. read length in each barcode/sample 
2. Align the sample fasta files to reference genome
3. Filter out off-target alignments
4. Merge alignments from all samples
5. Correct the alignments based on reference genome
6. Generate transcript files (GTF) for the merged alignments
7. Classify the transcripts based on known annotations 
8. Generating a metagene
8.a. Create a metagene based on known and novel exons
8.b. Create a GTF file for the metagene
9. Align the samples fasta files to the metagene
10. Generate counts based on the alignments to the metagene
11. Merge the transcript annotations and counts in "results"
12. Generate a heatmap and PCA plot (under construction)


