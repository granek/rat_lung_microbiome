#!/usr/bin/env bash
# Based on https://bitbucket.org/nsegata/lefse/raw/ee1653ca297801562a3118aa606d063085be20d8/example/run.sh

source /opt/conda/bin/activate qiime1

FORMAT="png"
INPUT_FILE=$1
BASENAME="${INPUT_FILE%.*}"

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
lefse-plot_res.py --format ${FORMAT} ${BASENAME}.res ${BASENAME}.${FORMAT}

# plot_cladogram.py visualizes the output on a hierarchical tree
#
# Plot the representation of the biomarkers on the hierarchical tree
# specified in the input data (using | in the name of the features)
# In this case we will obtain the RDP taxonomy.
# This is an early implementation of the module. I'm working on an improved version
# that will be released independently from LEfSe
lefse-plot_cladogram.py ${BASENAME}.res ${BASENAME}.cladogram.${FORMAT} --format ${FORMAT}

# Create a directory for storing the raw-data representation of the discovered biomarkers
mkdir -p ${BASENAME}_biomarkers_raw_images

# plot_features.py visualizes the raw-data features
#
# The module for exporting the raw-data representation of the features.
# With the default options we will obtain the images for all the features that are
# detected as biomarkers
lefse-plot_features.py ${BASENAME}.in ${BASENAME}.res ${BASENAME}_biomarkers_raw_images/

source deactivate


