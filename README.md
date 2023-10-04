# discoAnt
- v04102023_JG

## Setting up discoAnt

1. Download the workflow with \
   `git clone`
  
3. Downloading relevant reference file and creating a conda environment \
  `cd discoAnt` \
  `bash scripts/discoAnt_setup`
 
 Note: depending on your operating system you may need to edit the setup file to either 'conda activate' or 'source activate'

3. Running a test script to check installation \
  `conda activate discoAnt` \
  `bash discoAnt_main sirv_test_data/SIRV_params.txt`
  
4. Update example_parameters.txt with required sample information and paths to relevant folders before running the workflow with the data set of your choice \
  `conda activate discoAnt.env` \
  `bash discoAnt_main example_parameters.txt`

#### Refer to the sample submission scripts for running the pipeline on a HPC cluster

## Workflow

![This is an image](https://github.com/shwetajoshi-15/discoAnt/blob/main/discoAnt_workflow.png)



