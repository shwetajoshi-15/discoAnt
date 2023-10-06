FASTA=/data/scratch/users/yairp/discoAnt/downsample_test/pass_subset
downsampling=TRUE
reads_per_barcode_post_downsampling=100

function downsampling_function() {
	if [ "$downsampling" == TRUE ]; then
		echo "Downsampling all barcodes"
		
		for subdir in "$FASTA"/*/
		do
			subdir=$(basename "$subdir") # Extract the subdirectory name
			mkdir -p $FASTA/reads_per_barcode_post_downsampling/
			
			## Concatenate all FASTQ files in the subdirectory into one
			cat $FASTA/$subdir/*.fastq > $FASTA/${subdir}_cat.fastq

			## Downsample the combined FASTQ file
			reformat.sh sample=$reads_per_barcode_post_downsampling \
			in=$FASTA/${subdir}_cat.fastq \
			out=$FASTA/reads_per_barcode_post_downsampling/$subdir.fastq

			# Clean up the combined file
			rm $FASTA/${subdir}_cat.fastq

			for filename in $FASTA/$subdir; do
				base=$(basename "$filename" .fa)
				echo "On sample $base"
				
				# If you still need to do something with individual files, you can do it here.
				# For example, you may want to move the original files to a different directory.
				
			done
		done
	fi
}

downsampling_function


