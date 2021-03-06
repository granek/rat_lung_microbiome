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
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("stringr"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Phyloseq object from RDS, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
antibiotic_wcontrol_ps = readRDS(args$RDS)
antibiotic_only_ps = subset_samples(antibiotic_wcontrol_ps,group %in% c("ASNLU", "APNLU", "AINLU"))

#==============================================================================
#' # Generate LefSE format directly
#+ Generate LefSE format directly, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Remove OTUs that do not show appear more than 5 times in more than half the samples
min_counts = 2
sample_proportion = 0.01
min_samples = ceiling(sample_proportion*nsamples(antibiotic_wcontrol_ps))
wh0 = genefilter_sample(antibiotic_wcontrol_ps, 
                        filterfun_sample(function(x) x >= min_counts), 
                        A=min_samples)
antibiotic_wcontrol.taxfilt.ps = prune_taxa(wh0, antibiotic_wcontrol_ps)
ntaxa(antibiotic_wcontrol_ps)
ntaxa(antibiotic_wcontrol.taxfilt.ps)

antibiotic_wcontrol.spread = psmelt(antibiotic_wcontrol.taxfilt.ps) %>% 
  mutate(SampleID = str_replace(SampleID, pattern="\\.", replacement="_")) %>%
  mutate(taxonomy = paste(Kingdom,Phylum,Class,Order,Family,Genus, sep="|")) %>%
  select(SampleID,OTU,Abundance,antibiotic_bool,antibiotic) %>%
  spread(OTU, Abundance) %>% 
  arrange(antibiotic)

lefse_outdir = file.path(args$outdir,"lefse")
dir.create(lefse_outdir, showWarnings = FALSE)

#' RepseqToTaxa:
#' Convert Repseq column names to Taxa column names in a spread data frame
#' The big problem here is that this needs to be done after all other 
#' manipulations to the dataframe, otherwise most functions will balk if there
#' are dataframe columns with identical names
#'
#' @param spread.df The dataframe generated from phyloseq object.
#' @param source.ps phyloseq object from which spread.df was derived.
RepseqToTaxa <- function(spread.df, source.ps) {
  tax.df = as.data.frame(tax_table(source.ps)) %>%
    rownames_to_column("repseq") %>%
    mutate(taxonomy = paste(Kingdom,Phylum,Class,Order,Family,Genus, sep="|")) %>%
    select(repseq, taxonomy)
  
  # need to preserve non-OTU column names (otherwise they get lost)
  colname_match = match(names(spread.df), tax.df$repseq)
  cols_to_keep = which(is.na(colname_match))
  colnames_to_keep = names(spread.df)[cols_to_keep]
  
  # replace repseqs with taxonomies
  names(spread.df) = tax.df$taxonomy[match(names(spread.df), tax.df$repseq)]
  # now reset the non-OTU column names
  names(spread.df)[cols_to_keep] = colnames_to_keep
  return(spread.df)
}

antibiotic_wcontrol.spread %>%
  select(-antibiotic_bool) %>%
  RepseqToTaxa(antibiotic_wcontrol.taxfilt.ps) %>%
  write.table(file=file.path(lefse_outdir, "antibiotic_factor.tsv"), 
            sep="\t", quote = FALSE,
            row.names = FALSE)

antibiotic_wcontrol.spread %>%
  RepseqToTaxa(antibiotic_wcontrol.taxfilt.ps) %>%
  write.table(file=file.path(lefse_outdir, "antibiotic_bool.tsv"), 
              sep="\t", quote = FALSE,
              row.names = FALSE)

non_oral = antibiotic_wcontrol.spread %>%
  filter(antibiotic_bool==TRUE) %>%
  select(-antibiotic_bool) %>%
  mutate(oral_bool = antibiotic=="oral") %>%
  select(SampleID,oral_bool,everything()) %>%
  RepseqToTaxa(antibiotic_wcontrol.taxfilt.ps) %>%
  write.table(file=file.path(lefse_outdir, "antibiotic_oral_bool.tsv"), 
              sep="\t", quote = FALSE,
              row.names = FALSE)


sessionInfo()
# writeLines(capture.output(sessionInfo()), file.path(results_dir,"analyze_phyloseq_object_sessionInfo.txt"))

