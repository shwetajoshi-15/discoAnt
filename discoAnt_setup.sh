#!/bin/bash

source discoAnt_params_v2.1.txt
conda env create -f discoAnt.yml

## Installing programs
ROOT_DIR=$PWD
mkdir -p $PROGRAMS

git clone https://github.com/ConesaLab/SQANTI3.git $PROGRAMS/SQANTI3
chmod +x $PROGRAMS/SQANTI3/utilities/gtfToGenePred

echo "Installing cDNA_Cupcake..."
conda activate discoAnt

cd $PROGRAMS
git clone https://github.com/youyupei/cDNA_Cupcake.git
cd cDNA_Cupcake
python setup.py build
python setup.py install
conda deactivate 

