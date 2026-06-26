################### 3. Taxonomic and functional profiling ###################
###3-1. Taxonomic and functional profiling using HUMAnN3
#!/bin/bash

source /path/to/miniconda3/bin/activate biobakery3

INPUT_DIR="/path/to/QC_finished_cleaned_fastq_files"
OUTPUT_DIR="/path/to/humann3_results"
PROTEIN_DB="/path/to/uniref_db"
NUCLEOTIDE_DB="/path/to/chocophlan_db"
METAPHLAN_DB="/path/to/metaphlan_db"
THREADS=16
LOG_FILE="${OUTPUT_DIR}/humann.log"

mkdir -p "${OUTPUT_DIR}"

FILES=($(ls ${INPUT_DIR}/*.fastq.gz))
TOTAL_FILES=${#FILES[@]}

echo "HUMAnN processing started on $(date)" | tee "${LOG_FILE}"

for ((i=0; i<TOTAL_FILES; i++)); do
    INPUT_FILE="${FILES[$i]}"
    SAMPLE_NAME=$(basename "$INPUT_FILE" .fastq.gz)

    echo "Processing ${SAMPLE_NAME} (File $((i+1)) of ${TOTAL_FILES})..." | tee -a "${LOG_FILE}"

    /path/to/miniconda3/envs/biobakery3/bin/humann \
        --input "${INPUT_FILE}" \
        --output "${OUTPUT_DIR}/${SAMPLE_NAME}_humann" \
        --threads "${THREADS}" \
        --protein-database "${PROTEIN_DB}" \
        --nucleotide-database "${NUCLEOTIDE_DB}" \
        --metaphlan-options "--bowtie2db ${METAPHLAN_DB}" \
        >> "${LOG_FILE}" 2>&1

    echo "Finished ${SAMPLE_NAME} at $(date)" | tee -a "${LOG_FILE}"
done

echo "HUMAnN processing completed on $(date)" | tee -a "${LOG_FILE}"


###3-2. Merging feature tables of all samples using HUMAnN3
##3-2-1. Taxonomy feature table
merge_metaphlan_tables.py "/path/to/metaphlan_results"/*_metagenome_metaphlan_bugs_list.tsv \
    > "/path/to/metaphlan_results"/merged_taxa_table.txt


##3-2-2. Normalization of function feature table
#!/bin/bash

INPUT_DIR="/path/to/humann_output"
OUTPUT_DIR="/path/to/humann_results_pathabundance"

mkdir -p $OUTPUT_DIR

for file in $(find $INPUT_DIR -name "*_metagenome_pathabundance.tsv"); do
    sample=$(basename $file _metagenome_pathabundance.tsv)
    humann_renorm_table \
        --input $file \
        --output $OUTPUT_DIR/${sample}_pathabundance_relab.tsv \
        --units relab
done


##3-2-3. Merging of function feature table
humann_join_tables \
  --input /path/to/humann_results_pathabundance \
  --output /path/to/relab_merged/pathabundance_relab_merged.tsv \
  --file_name pathabundance_relab

#3-2-4. Generating stratified relative abundance bar plots for visualization of species-level contributions to metabolic pathways using HUMAnN3 [v.3.0]
humann_barplot \
  --input "pathabundance_relab_merged_with_metadata" \
  --last-metadata "last_metadata_name" \ # e.g., Group
  --focal-feature "focal_feature" \ # e.g., PWY-5973
  --focal-metadata "focal_metadata" \ # e.g., Group
  --sort sum metadata \
  --exclude-unclassified \
  --scaling logstack \
  --output "output"

##3-2-5. Stratification of function feature table
humann_split_stratified_table \
--input /path/to/relab_merged/pathabundance_relab_merged.tsv \
--output /path/to/relab_merged/pathabund_split
