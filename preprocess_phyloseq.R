#!/usr/bin/env Rscript
#' ****************************************************************************
#+ Setup: Parse Commandline, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--basedir", default=".",
                    help="Base directory for analysis [default %(default)s]",
                    metavar="DIR")
args <- parser$parse_args()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Setup Paths, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
basedir = args$basedir

workdir = file.path(basedir, "workspace")
results_dir = file.path(workdir,"results")
figure_dir = file.path(workdir,"figures")
dir.create(figure_dir, showWarnings = TRUE)

phyloseq_rds.filename = file.path("results", "rat_lung_ps.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Libraries, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(dplyr))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Phyloseq object from RDS, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full_ps = readRDS(phyloseq_rds.filename)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Relevel so "none" treatment is first, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_data(full_ps)$antibiotic = relevel(sample_data(full_ps)$antibiotic, "none")
sample_data(full_ps)$left_aspiration = relevel(sample_data(full_ps)$left_aspiration, "none")
sample_data(full_ps)$right_aspiration = relevel(sample_data(full_ps)$right_aspiration, "none")
sample_data(full_ps)$sample_aspiration = relevel(sample_data(full_ps)$sample_aspiration, "none")

sample_data(full_ps)$antibiotic_bool <- factor(get_variable(full_ps, "antibiotic") != "none")
sample_data(full_ps)$aspiration_bool <- factor(get_variable(full_ps, "sample_aspiration") != "none")

saveRDS(full_ps, file.path(workdir, "full_relevelled_ps.rds"))

#------------------------------------------------------------------------------
#+ Extract subset of replicates with most counts in each pair, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bacteria_ps_both_reps = subset_taxa(full_ps,Kingdom=="Bacteria") # Drop non-bacteria taxa
total_count.df = data.frame(sample_data(bacteria_ps_both_reps), total_count=sample_sums(bacteria_ps_both_reps))

set.seed(1)
max_replicate = total_count.df %>% 
  select(SampleID, Description, group, total_count) %>%
  group_by(Description) %>%
  filter(rank(-total_count, ties.method="random")==1)

# check to be sure replicates were removed
max_replicate %>% select(Description) %>% duplicated() %>% any()

bacteria.ps = subset_samples(bacteria_ps_both_reps,SampleID %in% max_replicate$SampleID)

rm(full_ps, bacteria_ps_both_reps) # Clean up to be sure these aren't used

#------------------------------------------------------------------------------
#+ Extract "antibiotic" samples, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample_data(max_rep_bacteria_ps) %>% filter(group %in% ())
antibiotic_ps = subset_samples(bacteria.ps,group %in% c("NLU", "ASNLU", "APNLU", "AINLU"))
saveRDS(antibiotic_ps, file.path(workdir, "antibiotic_ps.rds"))


#'****************************************************************************
#' # Session Info
#--------------------------------------------------
#+ Session Info, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sessionInfo()
# writeLines(capture.output(sessionInfo()), file.path(results_dir,"analyze_phyloseq_object_sessionInfo.txt"))

