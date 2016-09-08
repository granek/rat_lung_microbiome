source /opt/conda/bin/activate qiime1
# print_qiime_config.py -t

RAW_DATA_DIR="./raw_data/160614_McKenney_fastqs/"
READ1_FASTQ="$RAW_DATA_DIR/Undetermined_S0_L001_R1_001.fastq.gz"
READ2_FASTQ="$RAW_DATA_DIR/Undetermined_S0_L001_R2_001.fastq.gz"
INDEX_FASTQ="$RAW_DATA_DIR/Undetermined_S0_L001_I1_001.fastq.gz"

MAP_FILE="./notes_and_info/rat_lung_map.tsv"


WORKSPACE_DIR="./workspace"
SPLIT_FASTQ_DIR="$WORKSPACE_DIR/split_fastq"


# Generate a MAP file that QIIME will accept
QIIME_MAP_FILE="$WORKSPACE_DIR/rat_lung_qiime_map.tsv"
sed  '1s/^/#/' $MAP_FILE > $QIIME_MAP_FILE

validate_mapping_file.py --mapping_fp $QIIME_MAP_FILE --output_dir $WORKSPACE_DIR
# Split without filtering based on <http://benjjneb.github.io/dada2/faq.html>
for FASTQ in ${READ1_FASTQ} # ${READ2_FASTQ}
do
    split_libraries_fastq.py -r 999 -n 999 -q 0 -p 0.0001 \
			     --sequence_read_fps $FASTQ \
			     --output_dir $SPLIT_FASTQ_DIR \
			     --barcode_read_fps $INDEX_FASTQ \
			     --mapping_fps $QIIME_MAP_FILE \
			     --barcode_type golay_12 \
			     --rev_comp_mapping_barcodes \
			     --retain_unassigned_reads 
done

# split_sequence_file_on_sample_ids.py
#--------------------------------------------------------------------------------

# split_libraries_fastq.py OPTIONS
#---------------------------------
# --sample_ids
#     Comma-separated list of samples ids to be applied to all sequences, must be one per input file path (used when data is not multiplexed) [default: None]
# --store_demultiplexed_fastq
#     Write demultiplexed fastq files [default: False]
# --barcode_type
#     The type of barcode used. This can be an integer, e.g. for length 6 barcodes, or “golay_12” for golay error-correcting barcodes. Error correction will only be applied for “golay_12” barcodes. If data is not barcoded, pass “not-barcoded”. [default: golay_12]
# --max_barcode_errors
#     Maximum number of errors in barcode [default: 1.5]
# --phred_offset
#     The ascii offset to use when decoding phred scores (either 33 or 64). Warning: in most cases you don’t need to pass this value [default: determined automatically] 
