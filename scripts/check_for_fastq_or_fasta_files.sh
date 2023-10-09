FASTA=/data/scratch/users/yairp/discoAnt/downsample_test/pass_subset/dnzlfn


function check_for_fastq_or_fasta_files() {

    if [ -d "$FASTA" ]; then
        echo "Searching for FASTQ or FASTA files in $FASTA and its subdirectories..."
        
        # Use find to search for FASTQ or FASTA files
        if find "$FASTA" -type f \( -name "*.fastq"  -o -name "*.fastq.gz" -o -name "*.fasta" -o -name "*.fa" \) -print -quit | grep -q . ; then
            echo "FASTQ or FASTA files found. Continuing..."
            # Call the next function or perform the next task here
        else
            echo "No FASTQ or FASTA files found. Exiting."
            #$exit 1 
        fi
    else
        echo "Directory not found: $directory"
        #exit 1  # Quit the script with an error code (adjust as needed)
    fi
}


check_for_fastq_or_fasta_files