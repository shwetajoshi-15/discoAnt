#/bin/bash!

## This script performs the following steps
## 1. Quality control for the fasta files
## 1.a. Generates a file with no. of reads
## 1.b. average read length
## 1.c. read length in each barcode/sample
## 2. Align the sample fasta files to reference genome
## 3. Filter out off-target alignments
## 4. Merge alignments from all samples
## 5. Correct and collapse transcripts
## 6. Classify the transcripts based on known annotations
## 7. Align the samples fasta files to the metagene
## 8. Generate counts/TPM based on the alignments to the metagene

source discoAnt_params.txt
export PYTHONPATH="$PROGRAMS/cDNA_Cupcake/sequence/:$PYTHONPATH"

echo "Making folders"

mkdir -p $DISCOANT/$GENE
mkdir -p $DISCOANT/$GENE/fasta_stats
mkdir -p $DISCOANT/$GENE/minimap2
mkdir -p $DISCOANT/$GENE/minimap2_target
mkdir -p $DISCOANT/$GENE/flair
mkdir -p $DISCOANT/$GENE/flair/$GENE
mkdir -p $DISCOANT/$GENE/flair/collapse
mkdir -p $DISCOANT/$GENE/sqanti3
mkdir -p $DISCOANT/$GENE/flair_metagene_minimap2
mkdir -p $DISCOANT/$GENE/flair_metagene_counts
mkdir -p $DISCOANT/$GENE/plots

##########                                          ##########
##########  1. Quality control for the fasta files  ##########
##########                                          ##########



##########                                                      ##########
########## 2. Align the sample fasta files to reference genome  ##########
##########                                                      ##########


echo "minimap2 - Mapping fasta files to genome"

	if [[ ! -f $DISCOANT/$GENE/minimap2/"$GENE"_sorted_pri_align.COMPLETED ]]
	then

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	minimap2 -ax splice $REF_HG38/GRCh38.p13.genome_edit.fa $FASTA/${base}.fa > $DISCOANT/$GENE/minimap2/${base}.sam
	samtools view -S -h -b $DISCOANT/$GENE/minimap2/${base}.sam | samtools sort - > $DISCOANT/$GENE/minimap2/${base}_sorted.bam
	samtools view -h -F 2308 $DISCOANT/$GENE/minimap2/${base}_sorted.bam | samtools sort - > $DISCOANT/$GENE/minimap2/${base}_pri_sorted.bam
	done
	
	touch $DISCOANT/$GENE/minimap2/"$GENE"_sorted_pri_align.COMPLETED

	else
	echo "Alignments to the genome are present in $DISCOANT/$GENE/minimap2"
	fi

##########                                      ##########
########## 3. Filter out off-target alignments  ##########
##########                                      ##########

echo "Extracting alignments to the target region/gene"

	if [[ ! -f $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_align.COMPLETED  ]]
	then

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"
	
	samtools index $DISCOANT/$GENE/minimap2/${base}_pri_sorted.bam
	samtools view -h $DISCOANT/$GENE/minimap2/${base}_pri_sorted.bam "$CHR:$GENE_START-$GENE_END" > $DISCOANT/$GENE/minimap2_target/${base}_pri_tar_sorted.bam
	done
	
	touch $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_align.COMPLETED

	else
	echo "Filtered primary alignments are present in $DISCOANT/$GENE/minimap2_target"
	fi

##########                                      ##########
########## 4. Merge alignments from all samples ##########
##########                                      ##########

echo "Merging minimap2 primary alignments"

	if [[ ! -f $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_merged_align.COMPLETED ]]
	then

	samtools merge -f $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.bam $DISCOANT/$GENE/minimap2_target/*_pri_tar_sorted.bam
	samtools sort $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.bam > $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam
	samtools index $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam
	samtools view -h -o $DISCOANT/$GENE/minimap2_target/"$GENE"_merged.sam $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam

	touch $DISCOANT/$GENE/minimap2_target/"$GENE"_sorted_pri_tar_merged_align.COMPLETED

	else
	echo "Merged alignments are present in $DISCOANT/$GENE/minimap2_target"
	fi

##########                                       ##########
########## 5. Transcript correction and collapse ##########
##########                                       ##########

echo "FLAIR - correct and collapse"

	cat $FASTA/*.fa > $ALL/all.fa

	python $PROGRAMS/flair/bin/bam2Bed12.py -i $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bam > $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bed
	python $PROGRAMS/flair/flair.py correct -q $DISCOANT/$GENE/minimap2_target/"$GENE"_merged_sorted.bed -c $REF_HG38/chrom_sizes.txt -f $REF_HG38/gencode.v35.annotation.gtf -g $REF_HG38/GRCh38.p13.genome_edit.fa -o $DISCOANT/$GENE/flair/correct
	
	python $PROGRAMS/flair/flair.py collapse -r $ALL/all.fa -q $DISCOANT/$GENE/flair/correct_all_corrected.bed -f $REF_HG38/gencode.v35.annotation.gtf -g $REF_HG38/GRCh38.p13.genome_edit.fa -o $DISCOANT/$GENE/flair/collapse/default --temp_dir $DISCOANT/$GENE/flair
	python $PROGRAMS/flair/flair.py collapse -r $ALL/all.fa -q $DISCOANT/$GENE/flair/correct_all_corrected.bed -f $REF_HG38/gencode.v35.annotation.gtf -g $REF_HG38/GRCh38.p13.genome_edit.fa -o $DISCOANT/$GENE/flair/collapse/default_500 --temp_dir $DISCOANT/$GENE/flair -s 500
	python $PROGRAMS/flair/flair.py collapse -r $ALL/all.fa -q $DISCOANT/$GENE/flair/correct_all_corrected.bed -f $REF_HG38/gencode.v35.annotation.gtf -g $REF_HG38/GRCh38.p13.genome_edit.fa -o $DISCOANT/$GENE/flair/collapse/default_100 --temp_dir $DISCOANT/$GENE/flair -s 100


##########                                       ##########
########## 5. Gffcompare - annotate transcripts  ##########
##########                                       ##########

gffcompare 

##########                          ##########
########## 6. Transcript annotation ##########
##########                          ##########

echo "Annotated transcripts are present in $DISCOANT/$GENE/sqanti3"

	if [[ ! -f $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti3.COMPLETED ]]
	then
	python $PROGRAMS/SQANTI3-1.3/sqanti3_qc.py \
	--gtf $DISCOANT/$GENE/flair/collapse/default_500.isoforms.gtf \
	$REF_HG38/gencode.v35.annotation.gtf $REF_HG38/GRCh38.p13.genome_edit.fa \
	--cage_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
	--polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt \
	-d $DISCOANT/$GENE/sqanti3 -o "$GENE"_500
	
	python $PROGRAMS/SQANTI3-1.3/sqanti3_qc.py \
	--gtf $DISCOANT/$GENE/flair/collapse/default_100.isoforms.gtf \
	$REF_HG38/gencode.v35.annotation.gtf $REF_HG38/GRCh38.p13.genome_edit.fa \
	--cage_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
	--polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt \
	-d $DISCOANT/$GENE/sqanti3 -o "$GENE"_100

	touch $DISCOANT/$GENE/sqanti3/"$GENE"_sqanti3.COMPLETED
	
	else
	echo "Annotated transcripts are present in $DISCOANT/$GENE/sqanti3"
	fi

##########                              ##########
########## 7. Alignment to the metagene ##########
##########                              ##########

echo "Mapping the samples to the metagene (transcriptome) generated by flair"

	if [[ ! -f $DISCOANT/$GENE/flair_metagene_minimap2/"$GENE"_pri_sorted_align.COMPLETED ]]
	then

	for filename in $DISCOANT/$GENE/minimap2/*_pri_sorted.bam
	do
	base=$(basename $filename _pri_sorted.bam)
	echo "On sample : $base"

	minimap2 -ax map-ont $DISCOANT/$GENE/flair/collapse/default_500.isoforms.fa $FASTA/${base}.fa > $DISCOANT/$GENE/flair_metagene_minimap2/${base}.sam
	samtools view -S -h -b $DISCOANT/$GENE/flair_metagene_minimap2/${base}.sam | samtools sort - > $DISCOANT/$GENE/flair_metagene_minimap2/${base}_sorted.bam
	samtools view -h -F 2308 $DISCOANT/$GENE/flair_metagene_minimap2/${base}_sorted.bam | samtools sort - > $DISCOANT/$GENE/flair_metagene_minimap2/${base}_pri_sorted.bam
	samtools index $DISCOANT/$GENE/flair_metagene_minimap2/${base}_pri_sorted.bam

	done
	touch $DISCOANT/$GENE/flair_metagene_minimap2/"$GENE"_pri_sorted_align.COMPLETED

	else
	echo "Alignments to the metagene are present in $DISCOANT/$GENE/flair_metagene_minimap2"
	fi


##########                                                    ##########
########## 8. Quantification based on the metagene with Flair ##########
##########                                                    ##########

echo "Quantifying transcripts based on primary alignments to the metagene"

	python $PROGRAMS/flair/flair.py quantify -r $DISCOANT/reads_manifest_"$GENE".tsv -i $DISCOANT/$GENE/flair/collapse/default_500.isoforms.fa -o $DISCOANT/$GENE/flair_metagene_counts/default_500 --tpm

##########                                       ##########
########## 9. Plotting isoform usage with flair ##########
##########                                       ##########

python $PROGRAMS/flair/bin/plot_isoform_usage.py $DISCOANT/$GENE/flair/collapse/default_500.isoforms.bed $DISCOANT/$GENE/flair_metagene_counts/default $GENE_ID $DISCOANT/$GENE/plots/"$GENE"_isoform_plots
