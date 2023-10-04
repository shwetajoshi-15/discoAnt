#!/bin/bash

#source discoAnt_params_v2.1.txt
conda env create -f discoAnt.yml

## Installing programs
mkdir programs

cd programs
git clone https://github.com/ConesaLab/SQANTI3.git
chmod +x programs/SQANTI3/utilities/gtfToGenePred

echo "Installing cDNA_Cupcake..."
source activate discoAnt

git clone https://github.com/youyupei/cDNA_Cupcake.git
cd cDNA_Cupcake
python setup.py build
python setup.py install
conda deactivate 

