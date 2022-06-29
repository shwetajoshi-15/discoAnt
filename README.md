# discoAnt
- version 04052022
- bambu v2.2.0
- SQANTI3 v4.2

## Prepare FASTA files in a folder
All the sample files should be in FASTA format (.fa)

## Setting up discoAnt

1. git clone repository
2. cd discoAnt
3. When running the pipeline for the first time - ```bash discoAnt_setup.sh```
4. Once the setup/download is complete, run the test script - 
  ```conda activate discoAnt.env```
  ```bash discoAnt_v2_SIRV.sh```
5. Update discoAnt_params.txt with Gene info and paths to relevant folders before running the workflow with the data set of your choice


#### Refer to the sample submission scripts for running the pipeline on a HPC cluster.




