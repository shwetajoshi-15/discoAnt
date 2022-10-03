# discoAnt
- version 03102022
- tested on conda v4.12.0

## Prepare FASTA files in a folder
All the sample files (pass-reads only) should be in FASTA format (.fa)

## Setting up discoAnt

1. Download the workflow with ```git clone```
  
2. Download relevant reference file and create a conda environment \
  ```cd discoAnt``` \
  ```bash discoAnt_setup.sh```

3. Set up programs for SQANTI3 \
  ```conda activate discoAnt``` \
  ```chmod +x programs/SQANTI3-4.2/utilities/gtfToGenePred``` \
  ```cd programs/cDNA_Cupcake``` \
  ```python setup.py build``` \
  ```python setup.py install``` \
  ```conda deactivate```

3. Run a test script to check installation \
  ```conda activate discoAnt``` \
  ```bash discoAnt_v2_SIRV_test.sh```
  
4. Update discoAnt_params.txt with gene info and paths to relevant folders before running the workflow with the data set of your choice \
  ```conda activate discoAnt``` \
  ```bash discoAnt_v2.1.sh```

#### Refer to the sample submission scripts for running the pipeline on a HPC cluster.

## Workflow

![This is an image](https://github.com/shwetajoshi-15/discoAnt/blob/discoAnt_2022/workflow.png)



