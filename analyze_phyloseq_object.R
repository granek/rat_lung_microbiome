#--------------------------------------------------
#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-f", "--filter_fastqs", action="store_true", default=FALSE,
                    help="Run filtering step (default is to assume it was already run)")
parser$add_argument("-q", "--quality_plots", type="integer", default=0,
                    help="Number of quality plots generate [default %(default)s]",
                    metavar="number")
parser$add_argument("--basedir", default=".",
                    help="Base directory for analysis [default %(default)s]",
                    metavar="DIR")
args <- parser$parse_args()
#--------------------------------------------------
basedir = args$basedir

map_file = file.path(basedir, "notes_and_info/rat_lung_map.tsv")
workdir = file.path(basedir, "workspace")

taxonomy_dir =file.path(workdir,"taxonomy_refs")
silva_ref = file.path(taxonomy_dir,"silva_nr_v123_train_set.fa.gz")

results_dir = file.path(workdir,"results")
psfile.prefix = file.path(results_dir, "rat_lung_ps")

##====================================================================
##====================================================================
#+ Setup: Load Libraries, include=FALSE
## print(.libPaths())
# library(dada2)
# library(ShortRead)
library(ggplot2)
library(phyloseq)
# library(dplyr)
# library(biom)
writeLines(capture.output(sessionInfo()), file.path(results_dir,"sessionInfo.txt"))
##====================================================================
loadPhyloseqFiles = function(otu_table_file,sample_data_file,tax_table_file){
  otab <- otu_table(read.csv(otu_table_file,row.names=1), taxa_are_rows=FALSE)
  taxtab <- tax_table(as.matrix(read.csv(tax_table_file,row.names=1)))

  sample.df = read.csv(sample_data_file,row.names=1)
  sample.df[is.na(sample.df)] <- c("")
  samdat = sample_data(sample.df)

  ps <- phyloseq(otab, samdat, taxtab)
  return(ps)
}
#==============================================================================
plot_richness(ps, x="animal", measures=c("Shannon", "Simpson"), color="antibiotic") + theme_bw()
plot_bar(ps, x="antibiotic", fill="Family") 

# output.files = outputPhyloseq(ps,seqid.map.df,psfile.prefix)
otu_table_file = output.files[[1]]
sample_data_file = output.files[[2]]
tax_table_file = output.files[[3]]

ps.loaded = loadPhyloseqFiles(otu_table_file,sample_data_file,tax_table_file)

print(proc.time() - ptm)
#==============================================================================
