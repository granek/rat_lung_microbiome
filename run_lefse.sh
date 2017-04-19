#!/usr/bin/env bash
# Based on https://bitbucket.org/nsegata/lefse/raw/ee1653ca297801562a3118aa606d063085be20d8/example/run.sh


source /opt/conda/bin/activate qiime1

BOOL_INPUT="workspace/lefse/antibiotic_bool.tsv"
FACTOR_INPUT="workspace/lefse/antibiotic_factor.tsv"

# INPUT_FILE=${1:-"workspace/lefse/antibiotic_factor.tsv"}
BASENAME="${INPUT_FILE%.*}"

# Running the LEfSe commands with -h gives the list of available options

# format_input.py convert the input data matrix to the format for LEfSe.
#
# In this example we use the class information in the first line (-c 1)
# the subclass in the second line (-s 2) and the subject in the third (-u 3).
# If the subclass or the subject are not present in the data you need to set
# the value -1 for them.
# -o 1000000 scales the feature such that the sum (of the same taxonomic leve)
# is 1M: this is done only for obtaining more meaningful values for the LDA score

# lefse-format_input.py $INPUT_FILE  "${BASENAME}.in" -f c -u 1 -c 2 -s 3 -o 1000000  --output_table ${BASENAME}.tab

INPUT_FILE=$BOOL_INPUT
BASENAME="${INPUT_FILE%.*}"
lefse-format_input.py $INPUT_FILE  "${BASENAME}.in" -f c -u 1 -c 2 -s 3 -o 1000000  --output_table ${BASENAME}.tab

INPUT_FILE=$FACTOR_INPUT
BASENAME="${INPUT_FILE%.*}"
lefse-format_input.py $INPUT_FILE  "${BASENAME}.in" -f c -u 1 -c 2 -o 1000000  --output_table ${BASENAME}.tab

for INPUT_FILE in $BOOL_INPUT $FACTOR_INPUT; do
    BASENAME="${INPUT_FILE%.*}"
    # BASENAME=
    # BASENAME=$(basename "${INPUT_FILE}")
    # DIRNAME=$(dirname "${INPUT_FILE}")
    # run_lefse.py performs the actual statistica analysis
    #
    # Apply LEfSe on the formatted data producing the results (to be further processed
    # for visualization with the other modules). The option available
    # can be listed using the -h option 
    run_lefse.py "${BASENAME}.in" ${BASENAME}.res

    # plot_res.py visualizes the output
    #
    # Plot the list of biomarkers with their effect size
    # Severak graphical options are available for personalizing the output
    # ../plot_res.py hmp_aerobiosis_small.res hmp_aerobiosis_small.png
    lefse-plot_res.py --format pdf ${BASENAME}.res ${BASENAME}.pdf
    
    # plot_cladogram.py visualizes the output on a hierarchical tree
    #
    # Plot the representation of the biomarkers on the hierarchical tree
    # specified in the input data (using | in the name of the features)
    # In this case we will obtain the RDP taxonomy.
    # This is an early implementation of the module. I'm working on an improved version
    # that will be released independently from LEfSe
    lefse-plot_cladogram.py ${BASENAME}.res ${BASENAME}.cladogram.pdf --format pdf
    
    # Create a directory for storing the raw-data representation of the discovered biomarkers
    mkdir -p ${BASENAME}/biomarkers_raw_images
    
    # plot_features.py visualizes the raw-data features
    #
    # The module for exporting the raw-data representation of the features.
    # With the default options we will obtain the images for all the features that are
    # detected as biomarkers
    lefse-plot_features.py ${BASENAME}.in ${BASENAME}.res ${BASENAME}/biomarkers_raw_images/
done




