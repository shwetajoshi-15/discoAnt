
#/bin/bash!
PARAM_FILE=${1}
source $PARAM_FILE

echo \
"
██████  ██ ███████  ██████  ██████   █████  ███    ██ ████████ 
██   ██ ██ ██      ██      ██    ██ ██   ██ ████   ██    ██    
██   ██ ██ ███████ ██      ██    ██ ███████ ██ ██  ██    ██    
██   ██ ██      ██ ██      ██    ██ ██   ██ ██  ██ ██    ██    
██████  ██ ███████  ██████  ██████  ██   ██ ██   ████    ██ 
"

#source /data/scratch/users/yairp/discoAnt/params/RFTN2_RIC.params.txt

export PYTHONPATH="programs/cDNA_Cupcake/sequence/:$PYTHONPATH"
export PYTHONPATH="programs/cDNA_Cupcake/:$PYTHONPATH"

echo "Making folders"

mkdir -p $GENE
mkdir -p $GENE/minimap2
mkdir -p $GENE/bambu
mkdir -p $GENE/metagene_salmon
mkdir -p $GENE/filtered_transcripts
mkdir -p $GENE/discarded_transcripts/
mkdir -p $GENE/gffcomp_outs/
mkdir -p $GENE/temp_files/

#filter ref files 
grep ${ENSG_ID} ${ANNA_GTF} > $GENE/temp_files/filt_chr.gtf
ANNA_GTF_filt=$GENE/temp_files/filt_chr.gtf

echo "Extracting chromossome from genome references"

samtools faidx ${REF_GENOME_FN} ${chr} > $GENE/temp_files/filt_chr.fa

REF_GENOME_FN_filt=$GENE/temp_files/filt_chr.fa

### define max intron length
if [ -z "$max_intron_length" ]; then
    max_intron_length="200" # Set to the default value if max_intron_length is empty
fi

##########                                                  ##########
########## 0. Counting the number of reads and downsampling ##########
##########                                                  ##########

# Downsampling

if [ "$downsampling" == TRUE ];
then

	echo "Downsampling all barcodes"

	mkdir $FASTA/$reads_per_barcode_post_downsampling

	for filename in $FASTA/${base}*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample: $base"

	#to fix 
	#reformat.sh samplereadstarget=$reads_per_barcode_post_downsampling in=$FASTA/${base}.fa out=$FASTA/$reads_per_barcode_post_downsampling/${base}.fa

	done

else 

	echo ""

fi


# Count number of reads in each barcode 

if [ "$downsampling" == TRUE ];
then

echo "Creating a temporary file with number of pass reads in each barcode post-downsampling"

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample: $base"

	grep -c "^>" $FASTA/$reads_per_barcode_post_downsampling/${base}.fa >> $GENE/no_of_reads_barcodewise_postdwnsmp_tmp.txt

	done

else

echo "Creating a temporary file with number of pass reads in each barcode"

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample: $base"

	grep -c "^>" $FASTA/${base}.fa >> $GENE/temp_files/no_of_reads_barcodewise_tmp.txt

	done

fi

echo "Creating metrics for the report"

num_of_barcodes=$( ls -lh $FASTA | wc -l )
#total_reads_pass=$(awk '{ sum += $1 } END { print sum }' $GENE/no_of_reads_pass_barcodewise_tmp.txt)
#total_reads_pass_post_dwnsmp=$(awk '{ sum += $1 } END { print sum }' $GENE/no_of_reads_barcodewise_postdwnsmp_tmp.txt)
total_reads=$(awk '{ sum += $1 } END { print sum }' $GENE/temp_files/no_of_reads_barcodewise_tmp.txt)
read_supp_for_transcript=$(echo $total_reads \* $readFraction_for_bambu  | bc)


##########                                                    ##########
########## 1. Aligning sample fasta files to reference genome ##########
##########                                                    ##########


if [ "$downsampling" == TRUE ];
then

echo "Aligning pass reads to "${chr}""

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample: $base"

	minimap2 -ax splice -G"${max_intron_length}"k --splice-flank=yes --eqx $REF_GENOME_FN_filt $FASTA/$reads_per_barcode_post_downsampling/${base}.fa > $GENE/minimap2/${base}.sam
	samtools view -S -h -b $GENE/minimap2/${base}.sam | samtools sort - > $GENE/minimap2/${base}_sorted.bam
	samtools view -h -F 2308 $GENE/minimap2/${base}_sorted.bam | samtools sort - > $GENE/minimap2/${base}_pri_sorted.bam
	
	done

	samtools merge -f $GENE/minimap2/$GENE_pri_merged.bam $GENE/minimap2/*_pri_sorted.bam
	samtools merge -f $GENE/minimap2/$GENE_merged.bam $GENE/minimap2/*_sorted.bam

	samtools index $GENE/minimap2/$GENE_pri_merged.bam
	samtools index $GENE/minimap2/$GENE_merged.bam

fi

if [ "$downsampling" == FALSE ];
then

#echo "Extracting chromossome from genome ref"

#samtools faidx "${REF_GENOME_FN_filt}" "${chr}" > REF_GENOME_FN_filt_chr.fa

echo "Aligning pass reads to $chr"

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample: $base"

	minimap2 -ax splice -G"${max_intron_length}"k --splice-flank=yes --eqx $REF_GENOME_FN_filt $FASTA/${base}.fa > $GENE/minimap2/${base}.sam
	samtools view -S -h -b $GENE/minimap2/${base}.sam | samtools sort - > $GENE/minimap2/${base}_sorted.bam
	samtools view -h -F 2308 $GENE/minimap2/${base}_sorted.bam | samtools sort - > $GENE/minimap2/${base}_pri_sorted.bam

	done

	samtools merge -f $GENE/minimap2/$GENE_pri_merged.bam $GENE/minimap2/*_pri_sorted.bam
	samtools merge -f $GENE/minimap2/$GENE_merged.bam $GENE/minimap2/*_sorted.bam

	samtools index $GENE/minimap2/$GENE_pri_merged.bam
	samtools index $GENE/minimap2/$GENE_merged.bam
fi

rm $GENE/minimap2/*.sam

##########                                                       ##########
########## 2.a. Correcting and collapsing transcripts with bambu ##########
##########                                                       ##########

echo "Creating transcripts with bambu"

	Rscript $SCRIPTS/bambu_tx_discovery.R -b $GENE/minimap2/$GENE_pri_merged.bam \
	-f $REF_GENOME_FN_filt \
	-t $ANNA_GTF_filt \
	-r $readFraction_for_bambu \
	-o $GENE/bambu

awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s%%\n",$1,$2,100*$2/s}' $GENE/bambu/uniqueCounts.txt $GENE/bambu/uniqueCounts.txt > $GENE/bambu/uniqueCounts_perc.txt
awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s%%\n",$1,$2,100*$2/s}' $GENE/bambu/CPM.txt $GENE/bambu/CPM.txt > $GENE/bambu/CPM_perc.txt

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

	echo $chr$'\t'$F1_end_5$'\t'$F1_end > $GENE/bambu/F1_window_tmp.bed
	echo $chr$'\t'$R1_start$'\t'$R1_start_5 > $GENE/bambu/R1_window_tmp.bed

echo "Filtering transcripts with bedtools based on primer window"

	bedtools subtract -A -a $GENE/bambu/extended_annotations.gtf -b $GENE/bambu/F1_window_tmp.bed | bedtools subtract -A -a - -b $GENE/bambu/R1_window_tmp.bed > $GENE/filtered_transcripts/extended_annotations_no_primer_overlap.gtf
	bedtools intersect -wa -a $GENE/bambu/extended_annotations.gtf -b $GENE/bambu/F1_window_tmp.bed | bedtools intersect -wa -a - -b $GENE/bambu/R1_window_tmp.bed > $GENE/filtered_transcripts/extended_annotations_with_primer_overlap.gtf

echo "Extracting transcripts (names) that overlap the primer window"

	cat $GENE/filtered_transcripts/extended_annotations_with_primer_overlap.gtf | awk '{ if ($3 == "transcript") print $12 }' | sed 's/"//g' | sed 's/;//g' > $GENE/filtered_transcripts/transcripts_primer_overlap.txt

echo "Extracting transcripts (names) that overlap the primer window and meet the CPM threshold"

	cat $GENE/bambu/CPM.txt | grep -wf $GENE/filtered_transcripts/transcripts_primer_overlap.txt | awk '{ if ($2 >= 1) print }' > $GENE/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1.txt
	cat $GENE/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1.txt | awk '{ print $1 }' > $GENE/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt

echo "Filtering transcripts based on primer window overlap and CPM threshold"

	cat $GENE/bambu/extended_annotations.gtf | grep -wf $GENE/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt > $GENE/filtered_transcripts/extended_annotations_primer_overlap_CPM_more_than_1.gtf
	cat $GENE/bambu/extended_annotations.gtf | grep -v -wf $GENE/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt > $GENE/discarded_transcripts/extended_annotations_primer_overlap_CPM_more_than_1_discarded.gtf

	# Removing temporary files
	rm $GENE/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1_list_tmp.txt

echo "Editing filtered GTF files for use with IsoMix"

	sed 's/tx./tx/g' $GENE/filtered_transcripts/extended_annotations_primer_overlap_CPM_more_than_1.gtf > $GENE/${ENSG_ID}_discoAnt_isoforms.gtf
	cat $GENE/filtered_transcripts/transcripts_primer_overlap_CPM_more_than_1.txt | sed 's/tx./tx/g' | awk '{ print $1"\t"$3}' > $GENE/filtered_transcripts/filtered_transcripts.txt

echo "Creating a metatranscriptome based on the filtered transcripts"

	gffread -w $GENE/filtered_transcripts/extended_annotations_filtered.fa -g $REF_GENOME_FN_filt $GENE/${ENSG_ID}_discoAnt_isoforms.gtf
	salmon index -t $GENE/filtered_transcripts/extended_annotations_filtered.fa -i $GENE/filtered_transcripts -k 31

else

echo "primer_site_based_filter = FALSE, filtering primer based on CPM threshold only"
echo "Extracting transcripts (names) that meet the CPM threshold"

	cat $GENE/bambu/CPM.txt |  awk '{ if ($2 >= 1) print }' > $GENE/filtered_transcripts/transcripts_CPM_more_than_1.txt
	cat $GENE/filtered_transcripts/transcripts_CPM_more_than_1.txt | awk '{ print $1 }' > $GENE/filtered_transcripts/transcripts_CPM_more_than_1_list_tmp.txt

echo "Filtering transcripts based on CPM threshold"

	cat $GENE/bambu/extended_annotations.gtf | grep -wf $GENE/filtered_transcripts/transcripts_CPM_more_than_1_list_tmp.txt > $GENE/bambu/extended_annotations_CPM_more_than_1.gtf
	cat $GENE/bambu/extended_annotations.gtf | grep -v -wf $GENE/filtered_transcripts/transcripts_CPM_more_than_1_list_tmp.txt > $GENE/discarded_transcripts/extended_annotations_CPM_more_than_1_discarded.gtf
	
	# Removing temporary files 
	rm $GENE/filtered_transcripts/transcripts_CPM_more_than_1_list_tmp.txt

echo "Editing filtered GTF files for use with IsoMix"

	sed 's/tx./tx/g' $GENE/bambu/extended_annotations_CPM_more_than_1.gtf > $GENE/${ENSG_ID}_discoAnt_isoforms.gtf
	cat $GENE/filtered_transcripts/transcripts_CPM_more_than_1.txt | sed 's/tx./tx/g' | awk '{ print $1"\t"$3}' > $GENE/filtered_transcripts/filtered_transcripts.txt

echo "Creating a metatranscriptome based on the filtered transcripts"

	gffread -w $GENE/filtered_transcripts/filtered_transcripts.fa -g $REF_GENOME_FN_filt $GENE/${ENSG_ID}_discoAnt_isoforms.gtf
	salmon index -t $GENE/filtered_transcripts/filtered_transcripts.fa -i $GENE/filtered_transcripts -k 31

echo "Creating metrics for the report"
#not ENST
	filtered_transcripts_known=$( cat $GENE/filtered_transcripts/filtered_transcripts.txt | grep ENST | grep -v _ | wc -l )
	filtered_transcripts_novel=$( cat $GENE/filtered_transcripts/filtered_transcripts.txt | grep tx | grep -v _ | wc -l )

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
	echo "On sample: $base"

	salmon quant --quiet -i $GENE/filtered_transcripts -l A -r $FASTA/${base}.fa -o $GENE/metagene_salmon/${base}

	done

else 

echo "Quantifying transcripts (mapping based mode in salmon)"

	for filename in $FASTA/*.fa
	do
	base=$(basename $filename .fa)
	echo "On sample: $base"

	salmon quant --quiet -i $GENE/filtered_transcripts -l A -r $FASTA/$reads_per_barcode_post_downsampling/${base}.fa -o $GENE/metagene_salmon/${base}
	done
fi

##########                                ##########
########## 2.d. Generating count matrix   ##########
##########                                ##########


echo "Generating count matrix"
cd $GENE
Rscript $SCRIPTS/combine_salmon_quants.R $ENSG_ID
cd ..



##########                           ##########
########## 3. Annotating transcripts ##########
##########                           ##########

echo "Annotating transcripts with gffcompare and SQANTI3"

gffcompare -r $ANNA_GTF_filt \
	-o $GENE/gffcomp_outs/gffcomp $GENE/${ENSG_ID}_discoAnt_isoforms.gtf
mv $GENE/gffcomp.* $GENE/gffcomp_outs/ 


python programs/SQANTI3/sqanti3_qc.py \
	$GENE/${ENSG_ID}_discoAnt_isoforms.gtf \
	$ANNA_GTF_filt $REF_GENOME_FN_filt \
	-d $GENE/sqanti -o $GENE --report skip
	#--CAGE_peak $REF_HG38/refTSS_v3.3_human_coordinate.hg38.bed \
	#--polyA_peak $REF_HG38/atlas.clusters.2.0.GRCh38.96.bed --polyA_motif_list $REF_HG38/human.polyA.list.txt


##########           ##########
########## 3. Report ##########
##########           ##########


cat << EOT >> $GENE/${ENSG_ID}_discoAnt_report.txt
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

#removing temp files
rm -rf $GENE/temp_files

echo "
██████  ██ ███████  ██████  ██████   █████  ███    ██ ████████      ██████  ██████  ███    ███ ██████  ██      ███████ ████████ ███████ 
██   ██ ██ ██      ██      ██    ██ ██   ██ ████   ██    ██        ██      ██    ██ ████  ████ ██   ██ ██      ██         ██    ██      
██   ██ ██ ███████ ██      ██    ██ ███████ ██ ██  ██    ██        ██      ██    ██ ██ ████ ██ ██████  ██      █████      ██    █████   
██   ██ ██      ██ ██      ██    ██ ██   ██ ██  ██ ██    ██        ██      ██    ██ ██  ██  ██ ██      ██      ██         ██    ██      
██████  ██ ███████  ██████  ██████  ██   ██ ██   ████    ██         ██████  ██████  ██      ██ ██      ███████ ███████    ██    ███████ "                                                                                                                                        
