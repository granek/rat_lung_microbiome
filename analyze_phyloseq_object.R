#!/usr/bin/env Rscript

#--------------------------------------------------
#+ Setup: Parse Commandline, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--basedir", default=".",
                    help="Base directory for analysis [default %(default)s]",
                    metavar="DIR")
args <- parser$parse_args()
#--------------------------------------------------
#+ Setup: Setup Paths, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
basedir = args$basedir

workdir = file.path(basedir, "workspace")
results_dir = file.path(workdir,"results")
figure_dir = file.path(workdir,"figures")
dir.create(figure_dir, showWarnings = TRUE)

phyloseq.rds = file.path("results", "rat_lung_ps.rds")

#--------------------------------------------------
#+ Setup: Load Libraries, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(phyloseq)
library(dplyr)
library(DESeq2)

#--------------------------------------------------
#+ Setup: Load Phyloseq object from RDS, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full_ps = readRDS(phyloseq.rds)

#==============================================================================
#==============================================================================
#' # How many taxa are from each kingdom?
#--------------------------------------------------
#+ How many taxa are from each kingdom?, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tax_table(full_ps) %>% as.data.frame %>% 
  group_by(Kingdom) %>% 
  summarise(total_taxa = n()) 

#' # Examine the Number of Counts per Sample
#--------------------------------------------------
#+ MinMaxFloatingBarplot, include=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#--------------------------------------------------
#+ MinMax Barplots by Kingdom, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MinMaxFloatingBarplot(full_ps,
                      file.path(figure_dir,"all_min_max_readcounts.pdf"),
                      "All Counts")
MinMaxFloatingBarplot(subset_taxa(full_ps,Kingdom=="Eukaryota"),
                      file.path(figure_dir,"eukaryota_min_max_readcounts.pdf"),
                      "Eukaryota Counts")
MinMaxFloatingBarplot(subset_taxa(full_ps,Kingdom=="Bacteria"),
                      file.path(figure_dir,"bacteria_min_max_readcounts.pdf"),
                      "Bacteria Counts")
MinMaxFloatingBarplot(subset_taxa(full_ps,Kingdom=="Archaea"),
                      file.path(figure_dir,"archaea_min_max_readcounts.pdf"),
                      "Archaea Counts")
MinMaxFloatingBarplot(subset_taxa(full_ps,is.na(Kingdom)),
                      file.path(figure_dir,"na_min_max_readcounts.pdf"),
                      "NA (Undetermined) Counts")

#==============================================================================
#==============================================================================
#==============================================================================
## ---------------------------------------------
## Extract subset of replicates with most counts in each pair
##---------------------------------------------
bacteria_ps = subset_taxa(full_ps,Kingdom=="Bacteria") # Drop non-bacteria taxa
total_counts = as.data.frame(rowSums(otu_table(full_ps)))
colnames(total_counts) = "totals"
max_replicate = left_join(add_rownames(sample_data(bacteria_ps)), add_rownames(total_counts)) %>% 
  select(rowname, Description, group, totals) %>%
  group_by(Description) %>% 
  top_n(n=1)

# Get rid of full_ps to be sure it isn't accidentally used
rm(full_ps)

# check to be sure replicates were removed
max_replicate %>% select(Description) %>% duplicated() %>% any()

max_rep_bacteria_ps = subset_samples(bacteria_ps,SampleID %in% max_replicate$rowname)

#==============================================================================
## ---------------------------------------------
## Alpha Diversity Plots
##---------------------------------------------
# plot_richness(max_rep_ps, x = "sample_aspiration", color = "antibiotic") + geom_boxplot()
plot_richness(max_rep_bacteria_ps, x = "antibiotic", color = "sample_aspiration", 
              measures = c("Chao1", "ACE", "Shannon", "InvSimpson"), nrow=2) + 
  geom_boxplot() +
  theme(panel.background = element_blank())
ggsave(file=file.path(figure_dir,"max_rep_alpha_diversity.pdf"))

# plot_richness(max_rep_ps, x = "antibiotic", color = "sample_aspiration", 
#               measures = c("Shannon")) + geom_boxplot() +
#   theme(panel.background = element_blank())

#==============================================================================

## ---------------------------------------------
## Abundance Plots
##---------------------------------------------
# tax_table(max_rep_ps) %>% as.data.frame %>% add_rownames() %>% filter(Kingdom == "Eukaryota") %>% head
# tax_table(max_rep_ps) %>% as.data.frame %>% add_rownames() %>% filter(Kingdom == "Archaea") %>% head
# tax_table(max_rep_ps) %>% as.data.frame %>% add_rownames() %>% filter(is.na(Kingdom)) %>% select(rowname) %>% head

max_rep_bacteria_ps.rel  = transform_sample_counts(max_rep_bacteria_ps, function(x) x / sum(x) )
max_rep_bacteria_ps.rel.filt = filter_taxa(max_rep_bacteria_ps.rel, function(x) var(x) > 1e-3, TRUE)


# plot_bar(max_rep_bacteria_ps.rel.filt, facet_grid=~antibiotic, fill="Genus")
# plot_bar(max_rep_bacteria_ps.rel.filt, facet_grid=antibiotic~sample_aspiration, fill="Genus")
# max_rep_bacteria_ps.rel.filt.left = subset_samples(bacteria_ps,lung=="left")
# mdf = psmelt(subset_samples(max_rep_bacteria_ps.rel.filt,lung=="left"))

## ---------------------------------------------
## Plot relative abundances in Left Lung Samples
## ---------------------------------------------
p = ggplot(psmelt(subset_samples(max_rep_bacteria_ps.rel.filt,lung=="left")), 
           aes_string(x = "animal", y = "Abundance", fill = "Genus"))
p = p + geom_bar(stat = "identity", position = "stack")
p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
p <- p + facet_grid(antibiotic~sample_aspiration)
p = p + ggtitle("Left Lungs")
print(p)
ggsave(file=file.path(figure_dir,"left_lung_abundance.png"))

## ---------------------------------------------
## Plot relative abundances in Right Lung Samples
## ---------------------------------------------
p = ggplot(psmelt(subset_samples(max_rep_bacteria_ps.rel.filt,lung=="right")), 
           aes_string(x = "animal", y = "Abundance", fill = "Genus"))
p = p + geom_bar(stat = "identity", position = "stack")
p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
p <- p + facet_grid(antibiotic~left_aspiration)
p = p + ggtitle("Right Lungs (untreated)")
print(p)
ggsave(file=file.path(figure_dir,"right_lung_abundance.png"))

#' ## Observations from relative abundance plots
#' I am surprised by the Rlung-whole_gastric-subq samples.
#' I would have expected that they would look like the Llung-none-subq,
#' but they look more like the Rlung-noantibiotic samples.
#' 
#' Is there cross-talk between the left and right lungs?
#' 
#' It is reasuring that the Llung-antibiotic samples look so similar despite 
#' antibiotic mode-of-delivery 

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
## Everything above here seems to be working
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

## ---------------------------------------------
## Plot relative abundances between antibiotic
## ---------------------------------------------
plot_bar(max_rep_bacteria_ps.rel.filt, x="antibiotic", fill="Genus")


# ggplot(psmelt(max_rep_bacteria_ps.rel.filt), aes(x="antibiotic", y=Abundance, fill="Genus")) + 
#   geom_bar(stat="identity") + 
#   facet_grid(antibiotic~sample_aspiration)




# plot_bar(max_rep_bacteria_ps.rel.filt,x="sample_aspiration", facet_grid="antibiotic", fill="Genus")
# plot_bar(max_rep_ps.rel.filt, "group", "Abundance", title=title)
# plot_richness(max_rep_ps, x = "antibiotic", color = "sample_aspiration", 
#               measures = c("Chao1", "ACE", "Shannon", "InvSimpson"), nrow=2) + 
#   geom_boxplot() +
#   theme(panel.background = element_blank())
# ggsave(file=file.path(figure_dir,"max_rep_alpha_diversity.pdf"))

#==============================================================================
# Ordination plots
# Derived from https://joey711.github.io/phyloseq/plot_ordination-examples.html

## Remove OTUs that do not show appear more than 5 times in more than half the samples
wh0 = genefilter_sample(max_rep_bacteria_ps, filterfun_sample(function(x) x > 5), A=0.5*nsamples(max_rep_bacteria_ps))
ps1 = prune_taxa(wh0, max_rep_bacteria_ps)

## Check if any rows are all zero, because transform_sample_counts will generate NaNs from these
## which (apply(otu_table(ps1), 1, function(row) all(row ==0 )))

## Transform to even sampling depth.
ps1 = transform_sample_counts(ps1, function(x) 1E6 * x/sum(x))

## Remove samples with NaN (samples with all zero rows generate NaN when transformed above because of division by zero)
ps1 = prune_samples(complete.cases(otu_table(ps1)),ps1)

# Keep only the most abundant five phyla.
phylum.sum = tapply(taxa_sums(ps1), tax_table(ps1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
ps1 = prune_taxa((tax_table(ps1)[, "Phylum"] %in% top5phyla), ps1)


#' NMDS Plot by Aspiration and Antibiotic
ps1.ord <- ordinate(ps1, "NMDS", "bray")
p2 = plot_ordination(ps1, ps1.ord, type="samples", color="sample_aspiration", shape="antibiotic")
p2 + geom_point(size=3)
# p2 + geom_polygon(aes(fill=sample_aspiration)) + geom_point(size=5) + ggtitle("samples")
# print(p2)
#==============================================================================


#==============================================================================
#==============================================================================
#' # Session Info
#--------------------------------------------------
#+ Session Info?, echo=FALSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sessionInfo()
writeLines(capture.output(sessionInfo()), file.path(results_dir,"analyze_phyloseq_object_sessionInfo.txt"))


