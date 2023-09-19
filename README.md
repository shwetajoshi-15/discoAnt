# discoAnt
- version 07082022
- tested on conda v4.12.0

## Updates in v2.1
- optional primer site based filter
- report: includes statistics for the FASTA files and number of known and novel transcripts
- edits to the transcript file output for compatibility with IsoMix

## Prepare FASTA files in a folder
All the sample files (pass-reads only) should be in FASTA format (.fa)

## Setting up discoAnt

1. Downloading the workflow with ```git clone```
  
2. Downloading relevant reference file and creating a conda environment \
  ```cd discoAnt``` \
  ```bash discoAnt_setup.sh```
 
3. Running a test script to check installation \
  ```conda activate discoAnt``` \
  ```bash discoAnt_v2_SIRV_test.sh```
  
4. Update discoAnt_params.txt with gene info and paths to relevant folders before running the workflow with the data set of your choice \
  ```conda activate discoAnt.env``` \
  ```bash discoAnt_v2.1.sh```

#### Refer to the sample submission scripts for running the pipeline on a HPC cluster.

## Workflow

![This is an image](https://github.com/shwetajoshi-15/discoAnt/blob/main/discoAnt_workflow.png)



