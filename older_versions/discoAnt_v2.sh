
source discoAnt_params.txt

export PYTHONPATH="$PROGRAMS/cDNA_Cupcake/sequence/:$PYTHONPATH"

echo "Making folders"

mkdir -p $RESULTS
mkdir -p $RESULTS/"$GENE"
mkdir -p $RESULTS/"$GENE"/minimap2
mkdir -p $RESULTS/"$GENE"/bambu
mkdir -p $RESULTS/"$GENE"/bambu_metagene_salmon

##########                                                    ##########
########## 1. Aligning sample fasta files to reference genome ##########
##########                                                    ##########


echo "minimap2 - Mapping fasta files to genome"


        for filename in $FASTA/*.fa
        do
        base=$(basename $filename .fa)
        echo "On sample : $base"

        minimap2 -ax splice --splice-flank=yes $REF_HG38/GRCh38.primary_assembly.genome_edit.fa $FASTA/${base}.fa > $RESULTS/"$GENE"/minimap2/${base}.sam
        samtools view -S -h -b $RESULTS/"$GENE"/minimap2/${base}.sam | samtools sort - > $RESULTS/"$GENE"/minimap2/${base}_sorted.bam
        samtools view -h -F 2308 $RESULTS/"$GENE"/minimap2/${base}_sorted.bam | samtools sort - > $RESULTS/"$GENE"/minimap2/${base}_pri_sorted.bam
        done

samtools merge -f $RESULTS/"$GENE"/minimap2/"$GENE"_pri_merged.bam $RESULTS/"$GENE"/minimap2/*_pri_sorted.bam
samtools merge -f $RESULTS/"$GENE"/minimap2/"$GENE"_merged.bam $RESULTS/"$GENE"/minimap2/*_sorted.bam

samtools index $RESULTS/"$GENE"/minimap2/"$GENE"_pri_merged.bam
samtools index $RESULTS/"$GENE"/minimap2/"$GENE"_merged.bam

##########                                                       ##########
########## 2.a. Correcting and collapsing transcripts with bambu ##########
##########                                                       ##########

Rscript $SCRIPTS/bambu_tx_discovery.R -b $RESULTS/"$GENE"/minimap2/"$GENE"_pri_merged.bam \
-f $REF_HG38/GRCh38.primary_assembly.genome_edit.fa \
-t $REF_HG38/gencode.v41.annotation.gtf \
-o $RESULTS/"$GENE"/bambu

## Extracting transcripts belonging to the Gene of Interest 

cat $RESULTS/"$GENE"/bambu/counts_transcript.txt | grep "$GENE_ID" | awk '{ if ($3 >= 1) print}' > $RESULTS/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1.txt
cat $RESULTS/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1.txt | awk '{ print $1 }' > $RESULTS/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1_transcript.txt
cat $RESULTS/"$GENE"/bambu/extended_annotations.gtf | grep -wf $RESULTS/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1_transcript.txt > $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.gtf

## Editing filtered GTF and counts file for isomix compatibility

sed 's/tx./tx/g' $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.gtf > $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1_isomix.gtf
cat $RESULTS/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1.txt | sed 's/tx./tx/g' | awk '{ print $1"\t"$3}' > $RESULTS/"$GENE"/bambu/counts_transcript_"$GENE_ID"_count_1_isomix.txt


##########                                                               ##########
########## 2.b. Creating a transcriptome based on the bambu transcripts  ##########
##########                                                               ##########

gffread -w $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.fa -g $REF_HG38/GRCh38.primary_assembly.genome_edit.fa $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.gtf
salmon index -t $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.fa -i $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1 -k 31

##########                                                              ##########
########## 2.c. Re-aligning and quantifying filtered bambu transcripts  ##########
##########                                                              ##########

for filename in $FASTA/*.fa
do
base=$(basename $filename .fa)
echo "On sample : $base"
        
salmon quant -i $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1 -l A \
-r $FASTA/${base}.fa -o $RESULTS/"$GENE"/bambu_metagene_salmon/${base}
        
done


##########                           ##########
########## 3. Annotating transcripts ##########
##########                           ##########

gffcompare -r $REF_HG38/gencode.v41.annotation.gtf \
-o $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1 $RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.gtf

python $PROGRAMS/SQANTI3-4.2/sqanti3_qc.py \
$RESULTS/"$GENE"/bambu/extended_annotations_"$GENE_ID"_count_1.gtf \
$REF_HG38/gencode.v41.annotation.gtf $REF_HG38/GRCh38.primary_assembly.genome_edit.fa \
--cage_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
--polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt \
-d $RESULTS/"$GENE"/sqanti -o $GENE --report skip
