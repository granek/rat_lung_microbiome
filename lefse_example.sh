

# https://bitbucket.org/biobakery/biobakery/wiki/lefse#rst-header-lefse-bitbucket
mkdir -p /home/rstudio/parker_rat_lung/workspace/lefse_workspace/
cd /home/rstudio/parker_rat_lung/workspace/lefse_workspace/
wget http://huttenhower.sph.harvard.edu/webfm_send/129
mkdir -p input tmp results output_figures
mv 129 input/hmp_aerobiosis_small.txt

source /opt/conda/bin/activate qiime1
lefse-format_input.py input/hmp_aerobiosis_small.txt tmp/hmp_aerobiosis_small.in -c 1 -s 2 -u 3 -o 1000000
run_lefse.py tmp/hmp_aerobiosis_small.in results/hmp_aerobiosis_small.res

lefse-plot_res.py results/hmp_aerobiosis_small.res output_figures/hmp_aerobiosis_small.png
lefse-plot_cladogram.py results/hmp_aerobiosis_small.res output_figures/hmp_aerobiosis_small.cladogram.png --format png

# Also see
# https://bitbucket.org/nsegata/lefse/src/54694b4b0d9e335ff1ecafff8db4f1e0cf7004da?at=default
# https://bitbucket.org/nsegata/metaphlan/wiki/MetaPhlAn_Pipelines_Tutorial (see Step 4: Taxonomic biomarker discovery with LEfSe)
