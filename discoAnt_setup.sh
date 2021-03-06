#!/bin/bash

source discoAnt_params.txt
conda env create -f discoAnt.yml

## Installing programs

mkdir -p $PROGRAMS

wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v4.2.tar.gz -P $PROGRAMS/
tar -xzf $PROGRAMS/v4.2.tar.gz
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P $PROGRAMS/SQANTI3-4.2/utilities/

cd $PROGRAMS
git clone https://github.com/Magdoll/cDNA_Cupcake.git


## GENCODE genome reference and annotation

mkdir -p $REF_HG38

        wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.p13.genome.fa.gz -P $REF_HG38/
        wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz -P $REF_HG38/

## Editing the reference and annotation files

        gunzip $REF_HG38/GRCh38.p13.genome.fa.gz
        gunzip $REF_HG38/gencode.v35.annotation.gtf.gz

        cut -d " " -f 1 $REF_HG38/GRCh38.p13.genome.fa > $REF_HG38/GRCh38.p13.genome_edit.fa

## Downlaoding reference files for SQANTI3 annotation

        wget http://reftss.clst.riken.jp/datafiles/3.3/human/refTSS_v3.3_human_coordinate.hg38.bed.gz -P $REF_HG38/
        wget https://raw.githubusercontent.com/Magdoll/images_public/master/SQANTI2_support_data/human.polyA.list.txt -P $REF_HG38/
        wget https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz --no-check-certificate -P $REF_HG38/

        gunzip $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed.gz 
        gunzip $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed.gz
        
        sed -i 's/^/chr/' $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed
        

        

