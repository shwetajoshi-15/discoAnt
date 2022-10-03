
source discoAnt_params.txt

export PYTHONPATH="$PROGRAMS/cDNA_Cupcake/sequence/:$PYTHONPATH"

echo "Making folders"

mkdir -p $RESULTS
mkdir -p $RESULTS/"$GENE"
mkdir -p $RESULTS/"$GENE"/bambu
mkdir -p $RESULTS/"$GENE"/filtered_transcripts
mkdir -p $RESULTS/"$GENE"/metagene_salmon


 echo "Reference annotation: hg38 and GENCODE v41"
                        
Rscript $SCRIPTS/bambu_tx_discovery_1.R -b $RESULTS/"$GENE"/minimap2/"$GENE"_pri_merged.bam \
-f $REF_HG38/GRCh38.primary_assembly.genome_edit.fa \
-t $REF_HG38/gencode.v41.annotation.gtf \
-r 0.01 \
-o $RESULTS/"$GENE"/bambu

awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s%%\n",$1,$2,100*$2/s}' $RESULTS/"$GENE"/bambu/uniqueCounts.txt $RESULTS/"$GENE"/bambu/uniqueCounts.txt > $RESULTS/"$GENE"/bambu/uniqueCounts_perc.txt
awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s%%\n",$1,$2,100*$2/s}' $RESULTS/"$GENE"/bambu/CPM.txt $RESULTS/"$GENE"/bambu/CPM.txt > $RESULTS/"$GENE"/bambu/CPM_perc.txt

echo "primer_site_based_filter = TRUE, filtering primer based on start and end coordinates provided and CPM threshold"

	## Create a 5bp window around the 3' sites of the primers
	((F1_end_5=F1_end-5))
	((R1_start_5=R1_start+5))

echo "Creating a primer window"

	echo $chr$'\t'$F1_end_5$'\t'$F1_end > $RESULTS/"$GENE"/bambu/F1_window_tmp.bed
	echo $chr$'\t'$R1_start$'\t'$R1_start_5 > $RESULTS/"$GENE"/bambu/R1_window_tmp.bed

echo "Filtering transcripts with bedtools based on primer window"

	bedtools subtract -A -a $RESULTS/"$GENE"/bambu/extended_annotations.gtf -b $RESULTS/"$GENE"/bambu/F1_window_tmp.bed | bedtools subtract -A -a - -b $RESULTS/"$GENE"/bambu/R1_window_tmp.bed > $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_no_primer_overlap.gtf
	bedtools intersect -wa -a $RESULTS/"$GENE"/bambu/extended_annotations.gtf -b $RESULTS/"$GENE"/bambu/F1_window_tmp.bed | bedtools intersect -wa -a - -b $RESULTS/"$GENE"/bambu/R1_window_tmp.bed > $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_with_primer_overlap.gtf

echo "Extracting transcripts (names) that overlap the primer window"

	cat $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_with_primer_overlap.gtf | awk '{ if ($3 == "transcript") print $12 }' | sed 's/"//g' | sed 's/;//g' > $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap.txt

echo "Extracting transcripts (names) that overlap the primer window and meet the CPM threshold"

	cat $RESULTS/"$GENE"/bambu/CPM.txt | grep -wf $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap.txt | awk '{ if ($2 >= 1) print }' > $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1.txt
	cat $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1.txt | awk '{ print $1 }' > $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt

echo "Filtering transcripts based on primer window overlap and CPM threshold"

	cat $RESULTS/"$GENE"/bambu/extended_annotations.gtf | grep -wf $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt > $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_primer_overlap_CPM_more_than_1.gtf
#	cat $RESULTS/"$GENE"/bambu/extended_annotations.gtf | grep -v -wf $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt > $RESULTS/"$GENE"/discarded_transcripts/extended_annotations_primer_overlap_CPM_more_than_1_discarded.gtf

	# Removing temporary files
	rm $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt

echo "Editing filtered GTF files for use with IsoMix"

	sed 's/tx./tx/g' $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_primer_overlap_CPM_more_than_1.gtf > $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered.gtf
	cat $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1.txt | sed 's/tx./tx/g' | awk '{ print $1"\t"$3}' > $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts.txt

echo "Creating a metatranscriptome based on the filtered transcripts"

	gffread -w $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered.fa -g $REF_HG38/GRCh38.primary_assembly.genome_edit.fa $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered.gtf
	salmon index -t $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered.fa -i $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered -k 31

echo "Quantifying transcripts (mapping based mode in salmon)"

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	salmon quant -i $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts -l A -r $FASTA/${base}.fa -o $RESULTS/"$GENE"/metagene_salmon/${base}

	done

echo "Annotating transcripts with gffcompare and SQANTI3"

gffcompare -r $REF_HG38/gencode.v41.annotation.gtf \
-o $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts_gffcomp $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts.gtf

python $PROGRAMS/SQANTI3-4.2/sqanti3_qc.py \
$RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts.gtf \
$REF_HG38/gencode.v41.annotation.gtf $REF_HG38/GRCh38.primary_assembly.genome_edit.fa \
--cage_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
--polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt \
-d $RESULTS/"$GENE"/sqanti -o $GENE --report skip

