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

source /data/gpfs/projects/punim0826/discoAnt/discoAnt_params.txt

export PYTHONPATH="$PROGRAMS/cDNA_Cupcake/sequence/:$PYTHONPATH"

echo "Making folders"

mkdir -p $DISCOANT/"$GENE"
mkdir -p $DISCOANT/"$GENE"/minimap2
mkdir -p $DISCOANT/"$GENE"/bambu
mkdir -p $DISCOANT/"$GENE"/bambu_metagene_salmon

##########                                                    ##########
########## 1. Aligning sample fasta files to reference genome ##########
##########                                                    ##########


echo "minimap2 - Mapping fasta files to genome"


        for filename in $FASTA/*.fa
        do
        base=$(basename $filename .fa)
        echo "On sample : $base"

        minimap2 -ax splice --splice-flank=no $REF_HG38/GRCh38.p13.genome_edit.fa $FASTA/${base}.fa > $DISCOANT/"$GENE"/minimap2/${base}.sam
        samtools view -S -h -b $DISCOANT/"$GENE"/minimap2/${base}.sam | samtools sort - > $DISCOANT/"$GENE"/minimap2/${base}_sorted.bam
        samtools view -h -F 2308 $DISCOANT/"$GENE"/minimap2/${base}_sorted.bam | samtools sort - > $DISCOANT/"$GENE"/minimap2/${base}_pri_sorted.bam
        done

samtools merge -f $DISCOANT/"$GENE"/minimap2/"$GENE"_pri_merged.bam $DISCOANT/"$GENE"/minimap2/*_pri_sorted.bam

##########                                                       ##########
########## 2.a. Correcting and collapsing transcripts with bambu ##########
##########                                                       ##########

Rscript $SCRIPTS/bambu_tx_discovery.R -b $DISCOANT/"$GENE"/minimap2_target/"$GENE"_pri_merged_sorted.bam \
-f $REF_HG38/GRCh38.p13.genome_edit.fa \
-t $REF_HG38/gencode.v35.annotation.gtf \
-o $DISCOANT/"$GENE"/bambu

## Extracting transcripts belonging to the Gene of Interest 

cat $DISCOANT/"$GENE"/bambu/counts_transcript.txt | grep "$GENE_ID" | awk '{ if ($3 >= 1) print}' > $DISCOANT/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1.txt
cat $DISCOANT/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1.txt | awk '{ print $1 }' > $DISCOANT/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1_transcript.txt
cat $DISCOANT/"$GENE"/bambu/extended_annotations.gtf | grep -wf $DISCOANT/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1_transcript.txt > $DISCOANT/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.gtf

##########                                                               ##########
########## 2.b. Creating a transcriptome based on the bambu transcripts  ##########
##########                                                               ##########

gffread -w $DISCOANT/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.fa -g $REF_HG38/GRCh38.p13.genome_edit.fa $DISCOANT/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.gtf
salmon index -t $DISCOANT/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.fa -i $DISCOANT/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1 -k 31

##########                                                              ##########
########## 2.c. Re-aligning and quantifying filtered bambu transcripts  ##########
##########                                                              ##########

for filename in $FASTA/*.fa
do
base=$(basename $filename .fa)
echo "On sample : $base"
        
salmon quant -i $DISCOANT/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1 -l A \
-r $FASTA/${base}.fa -o $DISCOANT/"$GENE"/bambu_metagene_salmon/${base}
        
done


##########                           ##########
########## 3. Annotating transcripts ##########
##########                           ##########

gffcompare -r $REF_HG38/gencode.v35.annotation.gtf \
-o $DISCOANT/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1 $DISCOANT/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.gtf

$PROGRAMS/SQANTI3-1.3/sqanti3_qc.py \
--gtf $DISCOANT/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.gtf \
$REF_HG38/gencode.v35.annotation.gtf $REF_HG38/GRCh38.p13.genome_edit.fa \
--cage_peak $REF_HG38/refTSS_v3.1_human_coordinate.hg38.bed \
--polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/polyA.list \
-o $DISCOANT/bambu
