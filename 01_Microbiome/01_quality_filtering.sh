################### 1. Quality filtering ###################
#!/bin/bash

TRIMMOMATIC_JAR=/path/to/Trimmomatic-0.39/trimmomatic-0.39.jar
ADAPTERS=/path/to/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa

INPUT_DIR=/path/to/RawData
OUTPUT_DIR=/path/to/TrimmomaticResults

mkdir -p $OUTPUT_DIR

for file1 in $INPUT_DIR/*_1.fastq.gz; do
    base=$(basename $file1 _1.fastq.gz)
    file2=$INPUT_DIR/${base}_2.fastq.gz

    paired1=$OUTPUT_DIR/${base}_R1_paired.fastq
    unpaired1=$OUTPUT_DIR/${base}_R1_unpaired.fastq
    paired2=$OUTPUT_DIR/${base}_R2_paired.fastq
    unpaired2=$OUTPUT_DIR/${base}_R2_unpaired.fastq

    # Trimmomatic
    java -jar $TRIMMOMATIC_JAR PE \
        $file1 $file2 \
        $paired1 $unpaired1 \
        $paired2 $unpaired2 \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:20 \
        MINLEN:105

    echo "Processed $file1 and $file2"
done