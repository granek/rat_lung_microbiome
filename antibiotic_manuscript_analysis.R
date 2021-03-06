#!/usr/bin/env Rscript


#'******************************************************************************
#' # Manuscript 1: Relationship of antibiotic modality on lung microbiome
#+ Goals, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 1. Samples
#'     + NLU1 through 9 (control)
#'     + ASNLU1 through 14 (subcutaneous antibiotic)
#'     + APNLU1 through 12 (oral antibiotic)
#'     + AINLU1 through 12 (IV antibiotic)
#' 2. Analyses
#'     + NMDS plots
#'         + Include no-antibiotic (NLU) samples
#'         + Excude no-antibiotic (NLU) samples
#'     + LEfSe
#'         + no-antibiotic vs yes-antibiotic (all administration methods combined)
#'         + compare among antibiotic administration method

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
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("stringr"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Phyloseq object from RDS, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
antibiotic_wcontrol_ps = readRDS(args$RDS)
antibiotic_only_ps = subset_samples(antibiotic_wcontrol_ps,group %in% c("ASNLU", "APNLU", "AINLU"))

#==============================================================================
#' # Remove rare taxa
#+ Remove rare taxa, include=FALSE
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

#==============================================================================
#' # Ordination plots
#+ Ordination plots: Pruning, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OrdintationPrep = function(in_ps, randseed=NULL, ...){
  set.seed(randseed)
  # Derived from https://joey711.github.io/phyloseq/plot_ordination-examples.html
  ## Remove OTUs that do not show appear more than 5 times in more than half the samples
  min_counts = 2
  sample_proportion = 0.01
  min_samples = ceiling(sample_proportion*nsamples(in_ps))
  wh0 = genefilter_sample(in_ps, 
                          filterfun_sample(function(x) x >= min_counts), 
                          A=min_samples)
  ps1 = prune_taxa(wh0, in_ps)
  ntaxa(in_ps)
  ntaxa(ps1)
  #' ## Pruning Taxa
  #' The data was pruned before ordination to remove rare taxa.
  #' To be included in the pruned dataset, a taxon must occur at least 
  #' `r min_counts` times in at least `r sample_proportion*100`% of samples 
  #' (`r min_samples` samples).
  #' Before pruning there are `r ntaxa(max_rep_bacteria_ps)` bacterial taxa 
  #' (all non-bacterial taxa have already been removed). After pruning there are 
  #' `r ntaxa(ps1)` bacterial taxa remaining.
  
  #--------------------------------------------------
  # Transform counts
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Check if any rows are all zero, because transform_sample_counts will generate NaNs from these
  ## which (apply(otu_table(ps1), 1, function(row) all(row ==0 )))
  
  ## Transform to even sampling depth.
  ps1 = transform_sample_counts(ps1, function(x) 1E6 * x/sum(x))
  
  ## Remove samples with NaN (samples with all zero rows generate NaN when transformed above because of division by zero)
  ps1 = prune_samples(complete.cases(otu_table(ps1)),ps1)
  
  # # Keep only the most abundant five phyla.
  # phylum.sum = tapply(taxa_sums(ps1), tax_table(ps1)[, "Phylum"], sum, na.rm=TRUE)
  # top5phyla = names(sort(phylum.sum, TRUE))[1:5]
  # ps1 = prune_taxa((tax_table(ps1)[, "Phylum"] %in% top5phyla), ps1)
  
  
  # ps1.ord <- ordinate(ps1, "NMDS", "bray", trace=2, sratmax=0.99999999999999, maxit = 300, sfgrmin = 1e-10)
  ps1.ord <- ordinate(ps1, "NMDS", "bray", trace=2, ...)
  ord.plot = plot_ordination(ps1, ps1.ord, type="samples")
  ps1.ord.data = ord.plot$data
  return(ps1.ord.data)
}
  
  



#' ## NMDS Plots
#' The goal of NMDS plots is to visualize relationships betwen datapoints.
#' In the NMDS plots below, each point represents a sample. The distance between
#' two points indicates how similar the two samples are in terms of bacterial 
#' community composition.
#---------------------------------------------------
#' ### NMDS Plot by Antibiotic
#' In this plot, points are colored by antibiotic treatment.
#' All datapoints are shown in each quadrant, but
#' in each quadrant only the points of interest are colored, other points are
#' shown in gray. The left quadrants highlight **Left Lung** samples, 
#' right quadrants **Right Lung** 
#' (right lungs never received aspiration treatments).  
#' The bottom quadrants highlight samples from lungs that were **aspirated**,
#' top quadrants highlight lungs that were **not aspirated**.
#' 
#' We can see that coordinate 1 (along the x-axis) does a reasonably good job of 
#' distinguishing aspirated lungs from unaspirated lungs, indicating that 
#' aspiration has a large effect on the bacterial community. However, 
#' differences in antibiotic treatment does not seem to have a large effect.
#+ NMDS plot: Antibiotic All Points, echo=FALSE

nsamples(antibiotic_wcontrol_ps)
nsamples(antibiotic_only_ps)


#' ## NMDS Plot For "Antibiotic" Samples with no antibiotic control
antibiotic_wcontrol_ord = OrdintationPrep(antibiotic_wcontrol_ps, randseed=1,
                                          sratmax=0.99999999999999999, maxit = 300, sfgrmin = 1e-12)

ggplot(antibiotic_wcontrol_ord, aes(NMDS1, NMDS2)) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(aes(color = antibiotic)) 


#' ## NMDS Plot For "Antibiotic" Samples without control
antibiotic_only_ord = OrdintationPrep(antibiotic_only_ps, randseed=1,
                                      sratmax=0.99999999999999999, maxit = 300, sfgrmin = 1e-12)
# antibiotic_only_ord = OrdintationPrep(antibiotic_only_ps)

ggplot(antibiotic_only_ord, aes(NMDS1, NMDS2)) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(aes(color = antibiotic)) 

# Add ellipses to ordination plots 
# https://github.com/joey711/phyloseq/issues/323
# stat_ellipse(type = "norm", linetype = 2, aes(color = antibiotic)) +
# stat_ellipse(type = "t", aes(color = antibiotic))

# Save plot to file
# ggsave(file=file.path(figure_dir,"antibiotic_nmds_bray.png"))
#==============================================================================
#' # Permanova
#+ Permanova, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stop("Does data need to be rarified before permanova?")
stop("Apply 'waste not' to preliminary analysis")
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#negative-numbers-in-my-transformed-data-table
# https://github.com/joey711/phyloseq/issues/689
# https://github.com/joey711/phyloseq/issues/184
# http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
# http://joey711.github.io/phyloseq-demo/Restroom-Biogeography.html
# https://rpubs.com/michberr/randomforestmicrobe
# 
# restroom0 = restroom
# restroom = prune_taxa(taxa_sums(restroom) > 0, restroom)
# Any empty samples? Apparently not, though the code for pruning “empty” samples is also shown. And procedes quickly since there is nothing in restroom to modify.
# 
# any(sample_sums(restroom) == 0)
# ## [1] FALSE
# restroom = prune_samples(sample_sums(restroom) > 0, restroom)
# What about the total reads per sample, and what does the distribution look like?
# 
# readsumsdf = data.frame(nreads = sort(taxa_sums(restroom), TRUE), sorted = 1:ntaxa(restroom), 
#                         type = "OTUs")
# readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(restroom), 
#                                                         TRUE), sorted = 1:nsamples(restroom), type = "Samples"))
# title = "Total number of reads"
# p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
# p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
# set.seed(28132)
# restroomR = rarefy_even_depth(restroom, sample.size = 500)
# And now an alternative to random subsampling, a simple proportional transformation, calling it restroomP.
# 
# restroomP = transform_sample_counts(restroom, function(x) 500 * x/sum(x))
# For sanity-check, let's replot the sample sums of each of these new data objects, to convince ourselves that all of the samples now sum to 500.

par(mfrow = c(1, 2))
title = "Sum of reads for each sample, restroomR"
plot(sort(sample_sums(restroomR), TRUE), type = "h", main = title, ylab = "reads", 
    ylim = c(0, 1000))
title = "Sum of reads for each sample, restroomP"
plot(sort(sample_sums(restroomP), TRUE), type = "h", main = title, ylab = "reads", 
    ylim = c(0, 1000))

# df = as(sample_data(restroomR), "data.frame")
# d = distance(restroomR, "bray")
# restroomadonis = adonis(d ~ SURFACE + GENDER + BUILDING + FLOOR, df)
# restroomadonis
# plot(restroomadonis$aov.tab)






#'******************************************************************************
#' # Further Analyses
#+ Todo List, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 1. Omnibus test to determine if there is difference among groups
#'     + ANOSIM
#'     + Mantel
#'     + PERMANOVA
#' 2. Taxa-specific test to determine what OTUs distinguish groups
#'     + LEfSE
#'     + Canonical correlation analysis
#'     + Dirichlet-multinomial regression

#--------------------------------------------------
#'****************************************************************************
#' # Session Info
#--------------------------------------------------
#+ Session Info, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sessionInfo()
# writeLines(capture.output(sessionInfo()), file.path(results_dir,"analyze_phyloseq_object_sessionInfo.txt"))

