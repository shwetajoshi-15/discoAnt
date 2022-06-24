#!/bin/bash
#SBATCH --partition=physical
#SBATCH --nodes=1
#SBATCH --account="abcd1234"
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=username@email.com
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#time in days-hours:mins:sec
#SBATCH --time=0-02:0:00

source activate discoAnt.env

bash discoAnt/discoAnt_v2.sh
