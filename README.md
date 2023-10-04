# discoAnt
- version 04102023_JG

## Setting up discoAnt

1. Downloading the workflow with ```git clone```
  
2. Downloading relevant reference file and creating a conda environment \
  ```cd discoAnt``` \
  ```bash scripts/discoAnt_setup.sh```
 
3. Running a test script to check installation \
  ```conda activate discoAnt``` \
  ```bash discoAnt_v2.1 sirv_test_data/SIRV_params.txt```
  
4. Update discoAnt_params.txt with gene info and paths to relevant folders before running the workflow with the data set of your choice \
  ```conda activate discoAnt.env``` \
  ```bash discoAnt_v2.1 example_parameters.txt```

#### Refer to the sample submission scripts for running the pipeline on a HPC cluster.

## Workflow

![This is an image](https://github.com/shwetajoshi-15/discoAnt/blob/main/discoAnt_workflow.png)



