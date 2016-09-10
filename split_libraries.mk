# source /opt/conda/bin/activate qiime1
# print_qiime_config.py -t

.ONESHELL:
SHELL = /bin/bash
dir_guard=@mkdir -p $(@D)


RAW_DATA_DIR := raw_data/160614_McKenney_fastqs
WORKSPACE_DIR := workspace
SPLIT_FASTQ_BASE := $(WORKSPACE_DIR)/split_fastq
TAGGED_FASTQ_DIR := $(WORKSPACE_DIR)/tagged_fastq

#--------------------------------------------------


# READ1_FASTQ="$(RAW_DATA_DIR)/Undetermined_S0_L001_R1_001.fastq.gz"
# READ2_FASTQ="$RAW_DATA_DIR/Undetermined_S0_L001_R2_001.fastq.gz"
# INDEX_FASTQ="$RAW_DATA_DIR/Undetermined_S0_L001_I1_001.fastq.gz"

FASTQ_PREFIX := $(RAW_DATA_DIR)/Undetermined_S0_L001
FASTQ_SUFFIX := 001.fastq.gz
READ1_FASTQ := $(FASTQ_PREFIX)_R1_$(FASTQ_SUFFIX)
READ2_FASTQ := $(FASTQ_PREFIX)_R2_$(FASTQ_SUFFIX)
INDEX_FASTQ := $(FASTQ_PREFIX)_I1_$(FASTQ_SUFFIX)

#--------------------------------------------------

RAW_FASTQS := $(READ1_FASTQ) $(READ2_FASTQ)
TAGGED_FASTQS := $(addprefix $(TAGGED_FASTQ_DIR)/, $(notdir $(RAW_FASTQS:.fastq.gz=_tagged.fastq)))
SPLIT_FASTQS := $(addprefix $(WORKSPACE_DIR)/split_, $(addsuffix /COMPLETION_STAMP, $(notdir $(RAW_FASTQS:.fastq.gz=))))



all : $(TAGGED_FASTQS) $(SPLIT_FASTQS) run_rscript


test : 
	echo $(SPLIT_FASTQS)


#--------------------------------------------------
# Generate a MAP file that QIIME will accept
QIIME_MAP_FILE := $(WORKSPACE_DIR)/rat_lung_qiime_map.tsv
MAP_FILE := notes_and_info/rat_lung_map.tsv

$(QIIME_MAP_FILE) : $(MAP_FILE)
	sed  '1s/^/#/' $< > $@

#--------------------------------------------------
$(TAGGED_FASTQ_DIR)/%_tagged.fastq : $(RAW_DATA_DIR)/%.fastq.gz $(INDEX_FASTQ) $(QIIME_MAP_FILE)
	$(dir_guard)
	source /opt/conda/bin/activate qiime1
	echo $*
	split_libraries_fastq.py -r 999 -n 999 -q 0 -p 0.0001 \
		--sequence_read_fps $(word 1,$^) \
		--output_dir $(@D)/$* \
		--barcode_read_fps $(word 2,$^) \
		--mapping_fps $(word 3,$^) \
		--barcode_type golay_12 \
		--rev_comp_mapping_barcodes \
		--store_demultiplexed_fastq \
		--retain_unassigned_reads
	mv $(@D)/$*/seqs.fastq $@
	mv $(@D)/$*/split_library_log.txt $(@D)/$*_split_library_log.txt
	rm -r $(@D)/$*
#--------------------------------------------------
$(WORKSPACE_DIR)/split_%/COMPLETION_STAMP : $(TAGGED_FASTQ_DIR)/%_tagged.fastq
	source /opt/conda/bin/activate qiime1
	echo $*
	split_sequence_file_on_sample_ids.py -i $(word 1,$^) \
					 --file_type fastq \
					 --output_dir $(@D)
	touch $@
#--------------------------------------------------
TAXONOMY_DIR := $(WORKSPACE_DIR)/taxonomy_refs
SILVA_DB := $(TAXONOMY_DIR)/silva_nr_v123_train_set.fa.gz

run_rscript : $(SPLIT_FASTQS) $(SILVA_DB)
	# Rscript --no-restore process_fastq_to_counts.R --quality_plots 5 
	# Rscript --no-restore process_fastq_to_counts.R --filter_fastqs
	Rscript --no-restore process_fastq_to_counts.R


$(SILVA_DB) :
	$(dir_guard)
	wget -O $@_tmp "http://benjjneb.github.io/dada2/Training/silva_nr_v123_train_set.fa.gz"
	mv $@_tmp $@


# mkdir -p $SPLIT_FASTQ_BASE

# validate_mapping_file.py --mapping_fp $QIIME_MAP_FILE --output_dir $WORKSPACE_DIR
# Split without filtering based on <http://benjjneb.github.io/dada2/faq.html>
# for READ in R1 R2
# do
#     split_libraries_fastq.py -r 999 -n 999 -q 0 -p 0.0001 \
# 			     --sequence_read_fps ${FASTQ_PREFIX}_${READ}_001.fastq.gz \
# 			     --output_dir ${LABELED_FASTQ_BASE}_${READ} \
# 			     --barcode_read_fps ${FASTQ_PREFIX}_I1_001.fastq.gz \
# 			     --mapping_fps $QIIME_MAP_FILE \
# 			     --barcode_type golay_12 \
# 			     --rev_comp_mapping_barcodes \
# 			     --store_demultiplexed_fastq \
# 			     --retain_unassigned_reads
    
#     split_sequence_file_on_sample_ids.py -i ${LABEL_FASTQ_BASE}_${READ}/seqs.fastq \
# 					 --file_type fastq \
# 					 --output_dir ${SPLIT_FASTQ_BASE}_${READ}
# done

# split_sequence_file_on_sample_ids.py
# split_sequence_file_on_sample_ids.py -i ~/parker_rat_lung/workspace/split_fastq_R1/seqs.fastq --file_type fastq --output_dir ~/parker_rat_lung/workspace/split_by_id_R1
#--------------------------------------------------------------------------------

# split_libraries_fastq.py OPTIONS
#---------------------------------
# --sample_ids
#     Comma-separated list of samples ids to be applied to all sequences, must be one per input file path (used when data is not multiplexed) [default: None]
# --store_demultiplexed_fastq
#     Write demultiplexed fastq files [default: False]
# --max_barcode_errors
#     Maximum number of errors in barcode [default: 1.5]
# --phred_offset
#     The ascii offset to use when decoding phred scores (either 33 or 64). Warning: in most cases you don’t need to pass this value [default: determined automatically] 
