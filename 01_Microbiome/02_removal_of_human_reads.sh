################### 2. Removal of human reads ###################
###2-1. Run Bowtie2
#!/bin/bash

BOWTIE2_INDEX=/path/to/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as

INPUT_DIR=/path/to/TrimmomaticResults
OUTPUT_DIR=/path/to/Bowtie2Output

mkdir -p $OUTPUT_DIR

for file1 in $INPUT_DIR/*_R1_paired.fastq; do
    base=$(basename $file1 _R1_paired.fastq)    
    file2=$INPUT_DIR/${base}_R2_paired.fastq

    output_sam=$OUTPUT_DIR/${base}_mapped_and_unmapped.sam

    # Run Bowtie2
    bowtie2 --very-sensitive-local -p 24 --seed 99 \
        -x $BOWTIE2_INDEX \
        -1 $file1 \
        -2 $file2 \
        -S $output_sam

    echo "Bowtie2 alignment completed for $file1 and $file2"
done


###2-2. SAM to BAM with samtools
#!/bin/bash

SAMTOOLS=/path/to/samtools-1.20/samtools
INPUT_DIR=/path/to/Bowtie2Output
OUTPUT_DIR=/path/to/SAMtoBAM

mkdir -p $OUTPUT_DIR

threads=24

# Loop through all SAM files and convert to BAM
for samfile in $INPUT_DIR/*.sam; do
    base=$(basename $samfile .sam)  # Extract base name without extension
    bamfile=$OUTPUT_DIR/${base}.bam  # Set output BAM file name

    # Convert SAM to BAM
    $SAMTOOLS view -@ $threads -bS $samfile > $bamfile

    echo "Converted $samfile to $bamfile"
done


###2-3. Extract unmapped reads
#!/bin/bash

SAMTOOLS=/path/to/samtools-1.20/samtools
INPUT_DIR=/path/to/SAMtoBAM
OUTPUT_DIR=/path/to/UnmappedBAM

threads=24

# Loop through all BAM files and extract unmapped reads
for bamfile in $INPUT_DIR/*.bam; do
    # Extract base name without extension
    base=$(basename $bamfile .bam)
    # Set output BAM file name
    unmapped_bamfile=$OUTPUT_DIR/${base}_unmapped.bam

    # Extract unmapped reads
    $SAMTOOLS view -@ $threads -b -f 12 -F 256 $bamfile > $unmapped_bamfile

    echo "Extracted unmapped reads from $bamfile to $unmapped_bamfile"
done


###2-4. Sort the reads
#!/bin/bash

INPUT_DIR=/path/to/UnmappedBAM
OUTPUT_DIR=/path/to/SortedBAM

mkdir -p $OUTPUT_DIR

# Loop through all BAM files ending with _mapped_and_unmapped_unmapped.bam
for file in $INPUT_DIR/*_mapped_and_unmapped_unmapped.bam; do
  base=$(basename $file _mapped_and_unmapped_unmapped.bam)
  
  output_file=$OUTPUT_DIR/${base}_unmap.sorted.bam
  
  # Sort the BAM file by name
  /path/to/samtools-1.20/samtools sort -@ 48 -n $file -o $output_file
  
  echo "Sorted $file into $output_file"
done


###2-5. Convert sorted reads to fastq format
#!/bin/bash

INPUT_DIR=/path/to/SortedBAM
OUTPUT_DIR=/path/to/metagenome_fastq

mkdir -p $OUTPUT_DIR

# Loop through all BAM files ending with _mapped_and_unmapped_unmapped.bam
for file in $INPUT_DIR/*_unmap.sorted.bam; do

base=$(basename $file _unmap.sorted.bam)

output_file=$OUTPUT_DIR/${base}_metagenome.fastq

# Sorted BAM file to FASTQ format
/path/to/samtools-1.20/samtools bam2fq -@ 48 -n $file -o $output_file

echo "Sorted $file into $output_file"
done


###2-6. Compress fastq files
#!/usr/bin/env bash
cd /path/to/metagenome_fastq

pigz -p 4 M{009..049}_metagenome.fastq & 
pigz -p 4 M{050..099}_metagenome.fastq & 
pigz -p 4 M{100..150}_metagenome.fastq & 
pigz -p 4 M{151..180}_metagenome.fastq & 

wait
