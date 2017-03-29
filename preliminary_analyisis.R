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

phyloseq.rds = file.path("results", "rat_lung_ps.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Libraries, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(gridExtra))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Load Phyloseq object from RDS, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full_ps = readRDS(phyloseq.rds)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Setup: Relevel so "none" treatment is first, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_data(full_ps)$antibiotic = relevel(sample_data(full_ps)$antibiotic, "none")
sample_data(full_ps)$left_aspiration = relevel(sample_data(full_ps)$left_aspiration, "none")
sample_data(full_ps)$right_aspiration = relevel(sample_data(full_ps)$right_aspiration, "none")
sample_data(full_ps)$sample_aspiration = relevel(sample_data(full_ps)$sample_aspiration, "none")
levels(sample_data(full_ps)$antibiotic)
levels(sample_data(full_ps)$left_aspiration)
levels(sample_data(full_ps)$right_aspiration)
levels(sample_data(full_ps)$sample_aspiration)

# antibiotic_bool = 
sample_data(full_ps)$antibiotic_bool <- factor(get_variable(full_ps, "antibiotic") != "none")
# aspiration_bool = 
sample_data(full_ps)$aspiration_bool <- factor(get_variable(full_ps, "sample_aspiration") != "none")



#' ****************************************************************************
#' # How many taxa are from each kingdom?
#+ How many taxa are from each kingdom?, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tax_table(full_ps) %>% as.data.frame %>% 
  group_by(Kingdom) %>% 
  summarise(total_taxa = n()) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Number of Counts per Sample
#+ MinMaxFloatingBarplot, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MinMaxFloatingBarplot = function(ps,plot_file,plot_title=""){
  total_counts = as.data.frame(rowSums(otu_table(ps)))
  colnames(total_counts) = "totals"
  
  group_table = sample_data(ps) %>% select(group,Description) %>% unique
  
  min_max_counts = left_join(add_rownames(sample_data(ps)), add_rownames(total_counts)) %>% 
    select(rowname, Description, group, totals) %>%
    group_by(Description) %>% 
    summarise(min=min(totals),max=max(totals)) %>% 
    left_join(group_table)
  
  max_min_plot = ggplot(min_max_counts, 
                      aes(x=Description,ymin = `min`, ymax = `max`,color=group)) + 
    geom_linerange(stat = 'identity') + 
    xlab('Sample') + 
    ylab('Counts') + 
    theme(axis.ticks.x=element_blank(),           
          axis.text.x=element_blank(),           
          panel.background = element_blank()) + 
    ggtitle(plot_title)
  ggsave(file=plot_file, max_min_plot)
  
  print(paste("Lowest Maximum Value:", min(min_max_counts$max)))
  return(max_min_plot)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ MinMax Barplots by Kingdom, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_barplot = MinMaxFloatingBarplot(full_ps,
                                    file.path(figure_dir,"all_min_max_readcounts.pdf"),
                                    "All Counts")
eukaryota_barplot = MinMaxFloatingBarplot(subset_taxa(full_ps,Kingdom=="Eukaryota"),
                                          file.path(figure_dir,"eukaryota_min_max_readcounts.pdf"),
                                          "Eukaryota Counts")
bacteria_barplot = MinMaxFloatingBarplot(subset_taxa(full_ps,Kingdom=="Bacteria"),
                                         file.path(figure_dir,"bacteria_min_max_readcounts.pdf"),
                                         "Bacteria Counts")
archaea_barplot = MinMaxFloatingBarplot(subset_taxa(full_ps,Kingdom=="Archaea"),
                                        file.path(figure_dir,"archaea_min_max_readcounts.pdf"),
                                        "Archaea Counts")
na_barplot = MinMaxFloatingBarplot(subset_taxa(full_ps,is.na(Kingdom)),
                                   file.path(figure_dir,"na_min_max_readcounts.pdf"),
                                   "NA (Undetermined) Counts")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Eukaryote contamination note
#' A careful analysis of the phylogenetic makeup of data indicates that a large
#' portion of the reads are of Eukarytoic origin.  Most of these appear to be 
#' derived from the host (rat mitochondrian 16s rRNA).  These floating barplots
#' show the range of counts in each sample for: total counts, Eukaryote only, 
#' Bacteria only, Archea only, and Undetermined. The bar represents the read 
#' count range for each sample.  Since each sample was sequenced twice, the bar 
#' top gives the number of reads for the replicate with the most counts, and
#' bar bottom the counts for the replicate with the least counts.

#------------------------------------------------------------------------------
#+ MinMax Barplots by Kingdom Print, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(all_barplot)
print(eukaryota_barplot)
print(bacteria_barplot)
print(archaea_barplot)
print(na_barplot)
rm(all_barplot, eukaryota_barplot, bacteria_barplot, archaea_barplot, na_barplot)

#'*****************************************************************************
#+ Low count sample note
#' It is a common rule of thumb that for amplicon sequence analysis,
#' each sample must have a minimum of 5,000-10,000 reads.
#' However an analysis by 
#' [Kuczynski, et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2898070/)
#' indicates that this strict requirement is generally not necessary. They point out
#' that, as with any statistical analysis, the amount of data required to 
#' observe an effect is inversely proportional to the expected effect size.
#' This means that we should still be able to see large scale community changes,
#' despite the limited number of bacterial reads available for most samples.
#' 
#' The sample replicate from the duplicate pair with the most reads was used for
#' all further  analyses.

#------------------------------------------------------------------------------
#+ Extract subset of replicates with most counts in each pair, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bacteria_ps_both_reps = subset_taxa(full_ps,Kingdom=="Bacteria") # Drop non-bacteria taxa
# total_count.df = data.frame(sample_data(full_ps), total_count=sample_sums(full_ps))
total_count.df = data.frame(sample_data(bacteria_ps_both_reps), total_count=sample_sums(bacteria_ps_both_reps))
#   
#   
#   
#   as.data.frame(rowSums(otu_table(full_ps)))
# total_count.df = as.data.frame(rowSums(otu_table(full_ps)))
# mysums = sample_sums(full_ps)
# sum.df <- 
# colnames(sum.df)[length(colnames(sum.df))] <- "depth"
# 
# colnames(total_counts) = "totals"
set.seed(1)
max_replicate = total_count.df %>% 
  select(SampleID, Description, group, total_count) %>%
  group_by(Description) %>%
  filter(rank(-total_count, ties.method="random")==1)

# check to be sure replicates were removed
max_replicate %>% select(Description) %>% duplicated() %>% any()

bacteria.ps = subset_samples(bacteria_ps_both_reps,SampleID %in% max_replicate$SampleID)

rm(full_ps, bacteria_ps_both_reps) # Clean up to be sure these aren't used


#==============================================================================
#' # Alpha Diversity Plots
#' There does not seem to be a large treatment effect on measures of alpha 
#' diversity.  The unaspirated, no-antibiotic lungs seem to have a number of higher
#' diversity outliers, as might be expected for unperturbed communities.  However, 
#' this could be an artifact due to the fact that untreated lungs make up more 
#' samples than any other treatment group.
#+ Alpha Diversity Plots, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha_plot = plot_richness(bacteria.ps, x = "antibiotic", 
                           color = "sample_aspiration", 
                           measures = c("Chao1", "ACE", "Shannon", "InvSimpson"), nrow=2) + 
  geom_boxplot() + theme(panel.background = element_blank())
ggsave(file=file.path(figure_dir,"max_rep_alpha_diversity.pdf"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Alpha Diversity Plot Print, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(alpha_plot)

#==============================================================================
#==============================================================================
#' # Relative Abundance Plots
#+ Relative Abundance Plots, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bacteria_ps.rel  = transform_sample_counts(bacteria.ps, function(x) x / sum(x) )
bacteria_ps.rel.filt = filter_taxa(bacteria_ps.rel, function(x) var(x) > 1e-3, TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Relative Abundance Plots: Left Lung Samples, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p = ggplot(psmelt(subset_samples(bacteria_ps.rel.filt,lung=="left")), 
           aes_string(x = "animal", y = "Abundance", fill = "Genus"))
p = p + geom_bar(stat = "identity", position = "stack")
p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
p <- p + facet_grid(antibiotic~sample_aspiration)
p = p + ggtitle("Left Lungs")
ggsave(file=file.path(figure_dir,"left_lung_abundance.png"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Relative Abundance Plots Print: Left Lung Samples, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(p)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Relative Abundance Plots Setup: Right Lung Samples, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p = ggplot(psmelt(subset_samples(bacteria_ps.rel.filt,lung=="right")), 
           aes_string(x = "animal", y = "Abundance", fill = "Genus"))
p = p + geom_bar(stat = "identity", position = "stack")
p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
p <- p + facet_grid(antibiotic~left_aspiration)
p = p + ggtitle("Right Lungs (untreated)")
ggsave(file=file.path(figure_dir,"right_lung_abundance.png"))

# ### Observations from relative abundance plots
# I am surprised by the Rlung-whole_gastric-subq samples.
# I would have expected that they would look like the Llung-none-subq,
# but they look more like the Rlung-noantibiotic samples.
# 
# Is there cross-talk between the left and right lungs?
# 
# It is reasuring that the Llung-antibiotic samples look so similar despite 
# antibiotic mode-of-delivery 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#+ Relative Abundance Plots Print: Right Lung Samples, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(p)



# #---------------------------------------------------------------
# #' # BROKEN
# #' ## relative abundances between antibiotic treatments
# #+ Relative Abundance: Antibiotic Treatments, echo=FALSE
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot_bar(bacteria_ps.rel.filt, x="antibiotic", fill="Genus")


#==============================================================================
#' # Ordination plots
#+ Ordination plots: Pruning, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derived from https://joey711.github.io/phyloseq/plot_ordination-examples.html

min_counts = 2
sample_proportion = 0.01
min_samples = ceiling(sample_proportion*nsamples(bacteria.ps))
wh0 = genefilter_sample(bacteria.ps, 
                        filterfun_sample(function(x) x >= min_counts), 
                        A=min_samples)

bacteria.pruned.ps = prune_taxa(wh0, bacteria.ps)
ntaxa(bacteria.ps)
ntaxa(bacteria.pruned.ps)

#' ## Pruning Taxa
#' The data was pruned before ordination to remove rare taxa.
#' To be included in the pruned dataset, a taxon must occur at least 
#' `r min_counts` times in at least `r sample_proportion*100`% of samples 
#' (`r min_samples` samples).
#' Before pruning there are `r ntaxa(bacteria.ps)` bacterial taxa 
#' (all non-bacterial taxa have already been removed). After pruning there are 
#' `r ntaxa(bacteria.pruned.ps)` bacterial taxa remaining.

#--------------------------------------------------
#+ Ordination plots: Transform counts, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Check if any rows are all zero, because transform_sample_counts will generate NaNs from these
## which (apply(otu_table(bacteria.pruned.ps), 1, function(row) all(row ==0 )))

## Transform to even sampling depth.
bacteria.even.ps = transform_sample_counts(bacteria.pruned.ps, function(x) 1E6 * x/sum(x))

## Remove samples with NaN (samples with all zero rows generate NaN when transformed above because of division by zero)
bacteria.even.ps = prune_samples(complete.cases(otu_table(bacteria.even.ps)),bacteria.even.ps)

# # Keep only the most abundant five phyla.
# phylum.sum = tapply(taxa_sums(bacteria.even.ps), tax_table(bacteria.even.ps)[, "Phylum"], sum, na.rm=TRUE)
# top5phyla = names(sort(phylum.sum, TRUE))[1:5]
# bacteria.even.ps = prune_taxa((tax_table(bacteria.even.ps)[, "Phylum"] %in% top5phyla), bacteria.even.ps)

bacteria.even.nmds <- ordinate(bacteria.even.ps, "NMDS", "bray")
nmds.plot = plot_ordination(bacteria.even.ps, bacteria.even.nmds, type="samples")
bacteria.even.nmds.data = nmds.plot$data

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
ggplot(bacteria.even.nmds.data, aes(NMDS1, NMDS2)) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(data = transform(bacteria.even.nmds.data, aspiration_bool = NULL, lung = NULL), 
             color = "grey90") +
  geom_point(aes(color = antibiotic)) + 
  facet_grid(aspiration_bool~lung,labeller = "label_both")
#+ NMDS plot: Antibiotic All Points SAVE, include=FALSE
ggsave(file=file.path(figure_dir,"antibiotic_nmds_bray.png"))


#' ### NMDS Plot by Aspiration
#' In this plot, points are colored by the aspiration treatment to the **left lung**.
#' As with the above plot, all datapoints are shown in each quadrant, 
#' but only a subset is highlighted with color, the rest are gray. 
#' Again, the left and right quadrants highlight left and right lungs, respectively.
#' The top quadrants highlight rats that **recevied antibiotics**,
#' the bottom quadrants highlight rats did **not** receive antibiotics.
#' 
#' We can agains see that coordinate 1 (along the x-axis) does a reasonably good job of 
#' distinguishing aspirated lungs from unaspirated lungs. While the saline 
#' treated lungs are close to the other aspirated lungs, they form a tight cluster, 
#' which seems to indicate that they have distinict communities from the other 
#' aspirated lungs. Furthermore, the saline cluster is closest to the untreated
#' lungs, indicating that, of the aspiration treatments, the saline treatment
#' produces the smallest perturbation.
#' 
#' It also seems that the antibiotic treatment results in a larger difference 
#' between the unaspirated and whole gastric aspirated lungs, although this could
#' be an artifact of the larger number of antibiotic treated, unaspirated lungs.
#' 
#' Finally, it appears that left lung aspiration treatment does not affect right 
#' (unaspirated) lung community structure. All points are colored by the **left lung**
#' aspiration treatment, and there does not appear to be clustering in the right
#' lung samples based on the treatment to the left lung.  It is possible that there
#' are small effects that are swamped out by the much larger effects in the left lung.
#' A more careful analysis might be able to tease out finer effects.
#+ NMDS plot: Aspiration All Points, echo=FALSE
ggplot(bacteria.even.nmds.data, aes(NMDS1, NMDS2)) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(data = transform(bacteria.even.nmds.data, antibiotic_bool = NULL, lung = NULL), 
             color = "grey90") +
  geom_point(aes(color = left_aspiration)) + 
  scale_colour_brewer(palette="Set1",type="qual") +
  facet_grid(antibiotic_bool~lung,labeller = "label_both")
#+ NMDS plot: Aspiration All Points SAVE, include=FALSE
ggsave(file=file.path(figure_dir,"aspiration_nmds_bray.png"))

#----------------------------------------------------------------
#+ Permanova with bacteria.even.ps, include=FALSE

#' # Permanova

# Let’s see if Bray Curtis varied significantly by tooth class
# otus <- t(otu_table(psb_ab))
# psb_ab.bray <- vegdist(otus, method="bray")
# psb_ab.aov <- adonis(psb_ab.bray~sample_data(psb_ab)$Class)
# psb_ab.aov 


# otus <- t(otu_table(bacteria.even.ps))
otus <- otu_table(bacteria.even.ps)
bacteria.even.bray <- vegdist(otus, method="bray")

adonis(bacteria.even.bray ~ sample_data(bacteria.even.ps)$antibiotic_bool + sample_data(bacteria.even.ps)$aspiration_bool)

adonis(bacteria.even.bray ~ sample_data(bacteria.even.ps)$antibiotic + sample_data(bacteria.even.ps)$sample_aspiration)

adonis(bacteria.even.bray ~ sample_data(bacteria.even.ps)$sample_aspiration + sample_data(bacteria.even.ps)$antibiotic)

# bacteria.even.aov <- adonis(bacteria.even.bray ~ sample_data(bacteria.even.ps)$antibiotic)
# # ps1.aov <- adonis(bacteria.even.bray ~antibiotic_bool, data=sample_data(ps1))
# bacteria.even.aov 

#' # Check Sample Read Depth
#----------------------------------------------------------------
#+ plotCounts Definition, include=FALSE
plotCounts = function(physeq){
  sum.df <- data.frame(sample_data(physeq), sample_sums(physeq))
  colnames(sum.df)[length(colnames(sum.df))] <- "depth"
  
  grouped_boxplot = ggplot(sum.df, aes(antibiotic_bool, depth,group=antibiotic_bool)) + 
    geom_boxplot() +
    facet_grid(~sample_aspiration) +
    xlab("Antibiotic") +
    ylab("Bacterial Read Counts")
  
  sample_barplot = ggplot(sum.df, aes(animal, depth)) + 
    geom_bar(stat="identity") + 
    facet_wrap(~group)
  return(list(grouped_boxplot=grouped_boxplot, sample_barplot=sample_barplot))
}

#----------------------------------------------------------------
#+ plotCounts of bacteria.ps, echo=FALSE
#' # Check Raw Count Sample Read Depth
raw.plots = plotCounts(bacteria.ps)
#' ## Distribution by Treatment
print(raw.plots["grouped_boxplot"])
#' ## Per Sample Read Depth, Organized by Treatment Group
print(raw.plots["sample_barplot"])

#----------------------------------------------------------------
#+ plotCounts of bacteria.even.ps, echo=FALSE
#' # Check Raw Count Sample Read Depth
bacteria.even.plots = plotCounts(bacteria.even.ps)
#' ## Distribution by Treatment
print(bacteria.even.plots["grouped_boxplot"])
#' ## Per Sample Read Depth, Organized by Treatment Group
print(bacteria.even.plots["sample_barplot"])


#+ rlog_transformation, echo=FALSE
#--------------------------------------------------
#' # Try Regularized Log Transformation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Based on http://statweb.stanford.edu/~susan/papers/oralRR.html
#make a new object that we will transform
bacteria.rlog <- bacteria.pruned.ps

#add 1 to each sample count prior to transforming
# otu_table(bacteria.rlog) <- otu_table(bacteria.ps) + 1
# colnames(sample_data(bacteria.ps))

#convert the phyloseq object to a DESeq object
bacteria.rlog.dds = phyloseq_to_deseq2(bacteria.rlog, ~antibiotic+sample_aspiration)

# Use alternative method to calculate geometric means to avoid problem with zeroes
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(bacteria.rlog.dds), 1, gm_mean)
bacteria.rlog.dds = estimateSizeFactors(bacteria.rlog.dds, geoMeans = geoMeans)
bacteria.rlog.dds = DESeq(bacteria.rlog.dds, fitType="local")



#do a regularized log transformation
# bacteria.rlog.dds.rld  <- rlog(bacteria.rlog.dds , blind=FALSE, fitType="local")

#extract the otu counts from the R-log transformed object
# bacteria.rlog.dds.assay <- as.matrix(assay(bacteria.rlog.dds))
# bacteria.rlog.dds.counts <- as.matrix(counts(bacteria.rlog.dds,normalized=TRUE))
bacteria.rlog.dds.counts <- as.matrix(assay(bacteria.rlog.dds))


#exchange the otu tables
otu_table(bacteria.rlog) <- otu_table(bacteria.rlog.dds.counts, taxa_are_rows = TRUE)

#' ## Regularized Log Transformation
rlog.plots = plotCounts(bacteria.rlog)
#' ### Distribution by Treatment
print(rlog.plots["grouped_boxplot"])
#' ### Per Sample Read Depth, Organized by Treatment Group
print(rlog.plots["sample_barplot"])



#--------------------------------------------------
#+ Ordination plots: transform rlog to rlog-relative, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Check if any rows are all zero, because transform_sample_counts will generate NaNs from these
## which (apply(otu_table(ps1), 1, function(row) all(row ==0 )))

# ## Transform to even sampling depth.
# bacteria.rlog.even = transform_sample_counts(bacteria.rlog, function(x) 1E6 * x/sum(x))
#
# ## Remove samples with NaN (samples with all zero rows generate NaN when transformed above because of division by zero)
# bacteria.rlog.even = prune_samples(complete.cases(otu_table(bacteria.rlog.even)),bacteria.rlog.even)

#+  plotCounts: bacteria.rlog.rel, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bacteria.rlog.rel  = transform_sample_counts(bacteria.rlog, function(x) x / sum(x) )
# bacteria.rlog.rel.filt = filter_taxa(bacteria.rlog.rel, function(x) var(x) > 1e-3, TRUE)
# x = taxa_sums(bacteria.rlog.rel)


#' # Check Raw Count Sample Read Depth
#' ## Relative Abundance of Regularized Log Transformed

rel.plots = plotCounts(bacteria.rlog.rel)
#' ## Distribution by Treatment
print(rel.plots["grouped_boxplot"])
#' ## Per Sample Read Depth, Organized by Treatment Group
print(rel.plots["sample_barplot"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Do an ordination on Bray Curtis
# Use the regularized log transformed data frame
# Also use the relative abundance standardized data frame
#relative abundances
rel.pcoa <- ordinate(bacteria.rlog.rel, method="PCoA", distance="bray")
rel.pcoa.p1 <- plot_ordination(bacteria.rlog.rel, rel.pcoa, "samples", color="group") + ggtitle("Relative Abundance")

#log transformed
rlog.pcoa <- ordinate(bacteria.rlog, method="PCoA", distance="bray")
rlog.pcoa.p2 <- plot_ordination(bacteria.rlog, rlog.pcoa, "samples", color="group") + ggtitle("Rlog-transformed")

#plot them side by side
## Loading required package: gridExtra
grid.arrange(rel.pcoa.p1, rlog.pcoa.p2, ncol=2)


#-----------------------------------------------------------------

rel.nmds <- ordinate(bacteria.rlog.rel, method="NMDS", distance="bray")
rel.nmds.p1 <- plot_ordination(bacteria.rlog.rel, rel.nmds, "samples") + ggtitle("Relative Abundance")

#log transformed
rlog.nmds <- ordinate(bacteria.rlog, method="NMDS", distance="bray")
rlog.nmds.p2 <- plot_ordination(bacteria.rlog, rlog.nmds, "samples", color="group") + ggtitle("Rlog-transformed")

#plot them side by side
## Loading required package: gridExtra
grid.arrange(rel.nmds.p1, rlog.nmds.p2, ncol=2)

#'******************************************************************************
#' # Further Analyses
#+ Todo List, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 1. Summary relative abundance plots (i.e. abundance across all samples in a group)
#'     + treatment group
#'     + antibiotic treatment
#'     + lung treatment
#' 1. Identify taxa that distinguish groups (e.g. unaspirated vs gastric)
#' 1. Compare duplicates from each sample to determine how well min sample replicates max 
#' 1. Paired analysis of L and R lungs from same animal to look for communication
#' 1. Switch to using 
#- 1. Anything else?

#--------------------------------------------------
#'****************************************************************************
#' # Session Info
#--------------------------------------------------
#+ Session Info, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sessionInfo()
writeLines(capture.output(sessionInfo()), file.path(results_dir,"analyze_phyloseq_object_sessionInfo.txt"))

