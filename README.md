# discoAnt
- v04102023_JG

## Setting up discoAnt

1. Download the workflow with \
   `git clone`
  
2. Setup conda environment and dependencies \
  `cd discoAnt` \
  `bash scripts/discoAnt_setup`
 
 Note: depending on your operating system you may need to edit the setup file to either 'conda activate' or 'source activate'

3. Test discoAnt (from within the discoAnt directory) to check installation \
   `cd discoAnt` \
   `conda activate discoAnt` \
   `bash discoAnt_main sirv_test_data/sirv_params.txt`
  
5. Update example_parameters.txt with required sample information and paths to relevant folders \
  `conda activate discoAnt` \
  `bash discoAnt_main example_parameters.txt`

### Refer to the sample submission scripts folder for running the pipeline on an HPC cluster

## Workflow

![This is an image](https://github.com/shwetajoshi-15/discoAnt/blob/main/discoAnt_workflow.png)



