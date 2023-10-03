
#/bin/bash!

echo \
"
██████  ██ ███████  ██████  ██████   █████  ███    ██ ████████ 
██   ██ ██ ██      ██      ██    ██ ██   ██ ████   ██    ██    
██   ██ ██ ███████ ██      ██    ██ ███████ ██ ██  ██    ██    
██   ██ ██      ██ ██      ██    ██ ██   ██ ██  ██ ██    ██    
██████  ██ ███████  ██████  ██████  ██   ██ ██   ████    ██ 
"

source discoAnt_params_v2.1.txt

export PYTHONPATH="$PROGRAMS/cDNA_Cupcake/sequence/:$PYTHONPATH"
export PYTHONPATH="$PROGRAMS/cDNA_Cupcake/:$PYTHONPATH"

echo "Making folders"

mkdir -p $RESULTS
mkdir -p $RESULTS/"$GENE"
mkdir -p $RESULTS/"$GENE"/minimap2
mkdir -p $RESULTS/"$GENE"/bambu
mkdir -p $RESULTS/"$GENE"/metagene_salmon
mkdir -p $RESULTS/"$GENE"/filtered_transcripts
mkdir -p $RESULTS/"$GENE"/discarded_transcripts/

##########                                                  ##########
########## 0. Counting the number of reads and downsampling ##########
##########                                                  ##########

# Downsampling

if [ "$downsampling" == TRUE ];
then

	echo "downsampling = TRUE, downsampling all barcodes"

	mkdir $FASTA/$reads_per_barcode_post_downsampling

	for filename in $FASTA/${base}*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	reformat.sh samplereadstarget=$reads_per_barcode_post_downsampling in=$FASTA/${base}.fa out=$FASTA/$reads_per_barcode_post_downsampling/${base}.fa

	done

else 

	echo "downsampling = FALSE, downsampling was not performed for this dataset"

fi


# Count number of reads in each barcode 

if [ "$downsampling" == TRUE ];
then

echo "downsampling = TRUE, creating a temporary file with number of pass reads in each barcode post-downsampling"

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	grep -c "^>" $FASTA/$reads_per_barcode_post_downsampling/${base}.fa >> $RESULTS/"$GENE"/no_of_reads_barcodewise_postdwnsmp_tmp.txt

	done

else

echo "downsampling = FALSE, creating a temporary file with number of pass reads in each barcode"

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	grep -c "^>" $FASTA/${base}.fa >> $RESULTS/"$GENE"/no_of_reads_barcodewise_tmp.txt

	done

fi

echo "Creating metrics for the report"

num_of_barcodes=$( ls -lh $FASTA | wc -l )
total_reads_pass=$(awk '{ sum += $1 } END { print sum }' $RESULTS/"$GENE"/no_of_reads_pass_barcodewise_tmp.txt)
total_reads_pass_post_dwnsmp=$(awk '{ sum += $1 } END { print sum }' $RESULTS/"$GENE"/no_of_reads_barcodewise_postdwnsmp_tmp.txt)
total_reads=$(awk '{ sum += $1 } END { print sum }' $RESULTS/"$GENE"/no_of_reads_barcodewise_tmp.txt)
read_supp_for_transcript=$(echo $total_reads \* $readFraction_for_bambu  | bc)

# Removing temporary files

rm $RESULTS/"$GENE"/no_of_reads_barcodewise_tmp.txt
rm $RESULTS/"$GENE"/no_of_reads_pass_barcodewise_tmp.txt
rm $RESULTS/"$GENE"/no_of_reads_barcodewise_postdwnsmp_tmp.txt

##########                                                    ##########
########## 1. Aligning sample fasta files to reference genome ##########
##########                                                    ##########


if [ "$downsampling" == TRUE ];
then

echo "Extracting chromossome from genome ref"

samtools faidx "${REF_GENOME_FN}" "${chr}" > REF_GENOME_FN_chr.fa

echo "Aligning pass reads to "${chr}""

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	if [ -z "$max_intron_length" ]; then
    max_intron_length="200k" # Set to the default value if max_intron_length is empty
	fi

	minimap2 -ax splice -G"${max_intron_length}"k --splice-flank=yes --eqx REF_GENOME_FN_chr.fa $FASTA/$reads_per_barcode_post_downsampling/${base}.fa > $RESULTS/"$GENE"/minimap2/${base}.sam
	samtools view -S -h -b $RESULTS/"$GENE"/minimap2/${base}.sam | samtools sort - > $RESULTS/"$GENE"/minimap2/${base}_sorted.bam
	samtools view -h -F 2308 $RESULTS/"$GENE"/minimap2/${base}_sorted.bam | samtools sort - > $RESULTS/"$GENE"/minimap2/${base}_pri_sorted.bam
	
	done

	samtools merge -f $RESULTS/"$GENE"/minimap2/"$GENE"_pri_merged.bam $RESULTS/"$GENE"/minimap2/*_pri_sorted.bam
	samtools merge -f $RESULTS/"$GENE"/minimap2/"$GENE"_merged.bam $RESULTS/"$GENE"/minimap2/*_sorted.bam

	samtools index $RESULTS/"$GENE"/minimap2/"$GENE"_pri_merged.bam
	samtools index $RESULTS/"$GENE"/minimap2/"$GENE"_merged.bam

fi

if [ "$downsampling" == FALSE ];
then

echo "Extracting chromossome from genome ref"

samtools faidx "${REF_GENOME_FN}" "${chr}" > REF_GENOME_FN_chr.fa

echo "Aligning pass reads to "${chr}""

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	if [ -z "$max_intron_length" ]; then
    max_intron_length="200k" # Set to the default value if max_intron_length is empty
	
	fi

	minimap2 -ax splice -G"${max_intron_length}"k --splice-flank=yes --eqx REF_GENOME_FN_chr.fa $FASTA/${base}.fa > $RESULTS/"$GENE"/minimap2/${base}.sam
	samtools view -S -h -b $RESULTS/"$GENE"/minimap2/${base}.sam | samtools sort - > $RESULTS/"$GENE"/minimap2/${base}_sorted.bam
	samtools view -h -F 2308 $RESULTS/"$GENE"/minimap2/${base}_sorted.bam | samtools sort - > $RESULTS/"$GENE"/minimap2/${base}_pri_sorted.bam

	done

	samtools merge -f $RESULTS/"$GENE"/minimap2/"$GENE"_pri_merged.bam $RESULTS/"$GENE"/minimap2/*_pri_sorted.bam
	samtools merge -f $RESULTS/"$GENE"/minimap2/"$GENE"_merged.bam $RESULTS/"$GENE"/minimap2/*_sorted.bam

	samtools index $RESULTS/"$GENE"/minimap2/"$GENE"_pri_merged.bam
	samtools index $RESULTS/"$GENE"/minimap2/"$GENE"_merged.bam
fi

rm REF_GENOME_FN_chr.fa

##########                                                       ##########
########## 2.a. Correcting and collapsing transcripts with bambu ##########
##########                                                       ##########

echo "Creating transcripts with bambu"

	Rscript $SCRIPTS/bambu_tx_discovery.R -b $RESULTS/"$GENE"/minimap2/"$GENE"_pri_merged.bam \
	-f $REF_GENOME_FN \
	-t $ANNA_GTF \
	-r $readFraction_for_bambu \
	-o $RESULTS/"$GENE"/bambu


awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s%%\n",$1,$2,100*$2/s}' $RESULTS/"$GENE"/bambu/uniqueCounts.txt $RESULTS/"$GENE"/bambu/uniqueCounts.txt > $RESULTS/"$GENE"/bambu/uniqueCounts_perc.txt
awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s%%\n",$1,$2,100*$2/s}' $RESULTS/"$GENE"/bambu/CPM.txt $RESULTS/"$GENE"/bambu/CPM.txt > $RESULTS/"$GENE"/bambu/CPM_perc.txt

##########                                                               ##########
########## 2.b. Filtering bambu transcripts and creating a transcriptome ##########
##########                                                               ##########

## Primer site + CPM based filter 

if [ "$primer_site_based_filter" == TRUE ];
then

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
	cat $RESULTS/"$GENE"/bambu/extended_annotations.gtf | grep -v -wf $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt > $RESULTS/"$GENE"/discarded_transcripts/extended_annotations_primer_overlap_CPM_more_than_1_discarded.gtf

	# Removing temporary files
	rm $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt

echo "Editing filtered GTF files for use with IsoMix"

	sed 's/tx./tx/g' $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_primer_overlap_CPM_more_than_1.gtf > $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered.gtf
	cat $RESULTS/"$GENE"/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1.txt | sed 's/tx./tx/g' | awk '{ print $1"\t"$3}' > $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts.txt

echo "Creating a metatranscriptome based on the filtered transcripts"

	gffread -w $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered.fa -g $REF_GENOME_FN $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered.gtf
	salmon index -t $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered.fa -i $RESULTS/"$GENE"/filtered_transcripts/extended_annotations_filtered -k 31

else

echo "primer_site_based_filter = FALSE, filtering primer based on CPM threshold only"
echo "Extracting transcripts (names) that meet the CPM threshold"

	cat $RESULTS/"$GENE"/bambu/CPM.txt |  awk '{ if ($2 >= 1) print }' > $RESULTS/"$GENE"/filtered_transcripts/transcripts_CPM_more_than_1.txt
	cat $RESULTS/"$GENE"/filtered_transcripts/transcripts_CPM_more_than_1.txt | awk '{ print $1 }' > $RESULTS/"$GENE"/filtered_transcripts/transcripts_CPM_more_than_1_list_tmp.txt

echo "Filtering transcripts based on CPM threshold"

	cat $RESULTS/"$GENE"/bambu/extended_annotations.gtf | grep -wf $RESULTS/"$GENE"/filtered_transcripts/transcripts_CPM_more_than_1_list_tmp.txt > $RESULTS/"$GENE"/bambu/extended_annotations_CPM_more_than_1.gtf
	cat $RESULTS/"$GENE"/bambu/extended_annotations.gtf | grep -v -wf $RESULTS/"$GENE"/filtered_transcripts/transcripts_CPM_more_than_1_list_tmp.txt > $RESULTS/"$GENE"/discarded_transcripts/extended_annotations_CPM_more_than_1_discarded.gtf
	
	# Removing temporary files 
	rm $RESULTS/"$GENE"/filtered_transcripts/transcripts_CPM_more_than_1_list_tmp.txt

echo "Editing filtered GTF files for use with IsoMix"

	sed 's/tx./tx/g' $RESULTS/"$GENE"/bambu/extended_annotations_CPM_more_than_1.gtf > $RESULTS/"$GENE"/bambu/extended_annotations_filtered.gtf
	cat $RESULTS/"$GENE"/filtered_transcripts/transcripts_CPM_more_than_1.txt | sed 's/tx./tx/g' | awk '{ print $1"\t"$3}' > $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts.txt

echo "Creating a metatranscriptome based on the filtered transcripts"

	gffread -w $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts.fa -g $REF_GENOME_FN $RESULTS/"$GENE"/bambu/extended_annotations_filtered.gtf
	salmon index -t $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts.fa -i $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts -k 31

echo "Creating metrics for the report"

	filtered_transcripts_known=$( cat $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts.txt | grep ENST | grep -v _ | wc -l )
	filtered_transcripts_novel=$( cat $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts.txt | grep tx | grep -v _ | wc -l )

fi


##########                                                              ##########
########## 2.c. Re-aligning and quantifying filtered bambu transcripts  ##########
##########                                                              ##########

if [ "$downsampling" == FALSE ];
then

echo "Quantifying transcripts (mapping based mode in salmon)"

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	salmon quant -i $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts -l A -r $FASTA/${base}.fa -o $RESULTS/"$GENE"/metagene_salmon/${base}

	done

else 

echo "Quantifying transcripts (mapping based mode in salmon)"

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample : $base"

	salmon quant -i $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts -l A -r $FASTA/$reads_per_barcode_post_downsampling/${base}.fa -o $RESULTS/"$GENE"/metagene_salmon/${base}
	done
fi


##########                           ##########
########## 3. Annotating transcripts ##########
##########                           ##########

echo "Annotating transcripts with gffcompare and SQANTI3"

gffcompare -r $ANNA_GTF \
-o $RESULTS/"$GENE"/filtered_transcripts/filtered_transcripts_gffcomp $RESULTS/"$GENE"/bambu/extended_annotations_filtered.gtf

python $PROGRAMS/SQANTI3/sqanti3_qc.py \
$RESULTS/"$GENE"/bambu/extended_annotations_filtered.gtf \
$ANNA_GTF $REF_GENOME_FN \
-d $RESULTS/"$GENE"/sqanti -o $GENE --report skip \
--CAGE_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
--polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt


##########           ##########
########## 3. Report ##########
##########           ##########


cat << EOT >> $RESULTS/"$GENE"/discoAnt_report.txt
`date`

Number of barcodes: $num_of_barcodes
Number of reads across barcodes: $total_reads_pass

Downsampling: $downsampling
Number of reads across barcodes post-downsampling: $total_reads_pass_post_dwnsmp

Read support for transcripts: 

Filters applied:
CPM threshold: ## recommended value is 1
Primer site based filter: $primer_site_based_filter

Known transcripts post filtering: $filtered_transcripts_known
Novel transcrupts post filtering: $filtered_transcripts_novel

EOT

##########                         ##########
########## 4. Genrate count matrix ##########
##########                         ##########

echo "
██████  ██ ███████  ██████  ██████   █████  ███    ██ ████████      ██████  ██████  ███    ███ ██████  ██      ███████ ████████ ███████ 
██   ██ ██ ██      ██      ██    ██ ██   ██ ████   ██    ██        ██      ██    ██ ████  ████ ██   ██ ██      ██         ██    ██      
██   ██ ██ ███████ ██      ██    ██ ███████ ██ ██  ██    ██        ██      ██    ██ ██ ████ ██ ██████  ██      █████      ██    █████   
██   ██ ██      ██ ██      ██    ██ ██   ██ ██  ██ ██    ██        ██      ██    ██ ██  ██  ██ ██      ██      ██         ██    ██      
██████  ██ ███████  ██████  ██████  ██   ██ ██   ████    ██         ██████  ██████  ██      ██ ██      ███████ ███████    ██    ███████ "                                                                                                                                        
