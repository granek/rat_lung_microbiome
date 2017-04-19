#!/usr/bin/env Rscript

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Parse Commandline, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--RDS", default="workspace/antibiotic_ps.rds",
                    help="RDS containing phyloseq object",
                    metavar="DIR")
parser$add_argument("--outdir", default="workspace",
                    help="RDS containing phyloseq object",
                    metavar="DIR")
args <- parser$parse_args()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Libraries, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tibble"))
# suppressPackageStartupMessages(library("ggplot2"))
# suppressPackageStartupMessages(library("tidyr"))
# suppressPackageStartupMessages(library("stringr"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Phyloseq object from RDS, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
antibiotic_wcontrol_ps = readRDS(args$RDS)

#==============================================================================
#' # Generate BIOM file for LefSE
#+ Generate BIOM file for LefSE, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(biomformat)
otu.df = as.data.frame(t(as.matrix(otu_table(antibiotic_wcontrol_ps)))) %>%
  rownames_to_column("#OTU ID")

write.table(otu.df, file="antibiotic_wcontrol.otu.tsv", sep="\t", quote = FALSE,
            row.names = FALSE)
# source /opt/conda/bin/activate qiime1
# biom convert -i antibiotic_wcontrol.otu.tsv -o antibiotic_wcontrol.otu.biom --to-hdf5

# blah = read.delim("antibiotic_wcontrol.otu.tsv")


# otu.df = t(as.matrix(otu_table(antibiotic_wcontrol_ps)))
# sam = as.matrix(sample_data(antibiotic_wcontrol_ps))
sam.df = as.data.frame(sample_data(antibiotic_wcontrol_ps)) %>%
  rename("#SampleID" = SampleID)
write.table(sam.df, file="antibiotic_wcontrol.sam.tsv", sep="\t", quote = FALSE,
            row.names = FALSE)

# biom add-metadata -i antibiotic_wcontrol.otu.biom -o antibiotic_wcontrol.wsam.biom --sample-metadata-fp antibiotic_wcontrol.sam.tsv 
tax.df = as.data.frame(tax_table(antibiotic_wcontrol_ps)) %>%
  rownames_to_column("#OTUID") %>%
  mutate(taxonomy = paste(Kingdom,Phylum,Class,Order,Family,Genus, sep=";")) %>%
  select(1, taxonomy)
write.table(tax.df, file="antibiotic_wcontrol.tax.tsv", sep="\t", quote = FALSE,
            row.names = FALSE)

# biom add-metadata -i antibiotic_wcontrol.otu.biom -o antibiotic_wcontrol.wsam.biom --observation-metadata-fp antibiotic_wcontrol.tax.tsv --sample-metadata-fp antibiotic_wcontrol.sam.tsv
# lefse-format_input.py /home/rstudio/parker_rat_lung/antibiotic_wcontrol.biom antibiotic_wcontrol.in -biom_c antibiotic

# tax = as.matrix(tax_table(antibiotic_wcontrol_ps))
# # antibiotic_wcontrol.biom = make_biom(otu, sample_metadata=sam, observation_metadata = tax)
# antibiotic_wcontrol.biom = make_biom(otu)
# outfile = "antibiotic_wcontrol.biom"
# write_biom(antibiotic_wcontrol.biom, outfile)
# 
# # x = psmelt(antibiotic_wcontrol_ps)
# 
# y = read_biom(outfile)
# identical(antibiotic_wcontrol.biom, y)
# sessionInfo()
# z = read_biom("otu_table.biom")
# 

# biom convert -i /home/rstudio/parker_rat_lung/antibiotic_wcontrol.biom -o /home/rstudio/parker_rat_lung/antibiotic_wcontrol_hdf5.biom --to-hdf5


#'****************************************************************************
#' # Session Info
#--------------------------------------------------
#+ Session Info, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sessionInfo()
# writeLines(capture.output(sessionInfo()), file.path(results_dir,"analyze_phyloseq_object_sessionInfo.txt"))

