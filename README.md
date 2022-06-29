# discoAnt
- version 28062022
- bambu v2.0.0
- SQANTI3 v4.2

## Prepare FASTA files in a folder
All the sample files should be in FASTA format (.fa)

## Setting up discoAnt

1. Downloading the workflow
  ```git clone repository``` \
  
2. Downloading relevant reference file and creating a conda environment
  ```cd discoAnt``` \
  ```bash discoAnt_setup.sh```
  
3. Running a test script to check installation \
  ```conda activate discoAnt.env``` \
  ```bash discoAnt_v2_SIRV.sh```
  
4. Update discoAnt_params.txt with gene info and paths to relevant folders before running the workflow with the data set of your choice
  ```conda activate discoAnt.env``` \
  ```bash discoAnt_v2.sh```

#### Refer to the sample submission scripts for running the pipeline on a HPC cluster.




