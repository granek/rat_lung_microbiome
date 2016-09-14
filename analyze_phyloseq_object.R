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
library(dplyr)
# library(biom)
writeLines(capture.output(sessionInfo()), file.path(results_dir,"analyze_phyloseq_object_sessionInfo.txt"))
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
list.files(results_dir)
setwd("./parker_rat_lung/")
# output.files = outputPhyloseq(ps,seqid.map.df,psfile.prefix)
otu_table_file = paste0(psfile.prefix,"_otu.csv")
sample_data_file = paste0(psfile.prefix,"_samdat.csv")
tax_table_file = paste0(psfile.prefix,"_tax.csv")

ps = loadPhyloseqFiles(otu_table_file,sample_data_file,tax_table_file)

plot_richness(ps, x="animal", measures=c("Shannon", "Simpson"), color="antibiotic") + theme_bw()
plot_bar(ps, x="antibiotic", fill="Family") 

#==============================================================================
# Ordination plots
# Derived from https://joey711.github.io/phyloseq/plot_ordination-examples.html

## Remove OTUs that do not show appear more than 5 times in more than half the samples
wh0 = genefilter_sample(ps, filterfun_sample(function(x) x > 5), A=0.5*nsamples(ps))
ps1 = prune_taxa(wh0, ps)

## Transform to even sampling depth.
ps1 = transform_sample_counts(ps1, function(x) 1E6 * x/sum(x))

# Keep only the most abundant five phyla.
phylum.sum = tapply(taxa_sums(ps1), tax_table(ps1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
ps1 = prune_taxa((tax_table(ps1)[, "Phylum"] %in% top5phyla), ps1)


# We will want to investigate a major prior among the samples, which is that some are human-associated microbiomes, and some are not. Define a human-associated versus non-human categorical variable:
  
# human = get_variable(GP1, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
# sample_data(GP1)$human <- factor(human)

## (2) Just samples

## Next, let’s plot only the samples, and shade the points by “SampleType” while also modifying the shape according to whether they are human-associated. There are a few additional ggplot2 layers added to make the plot even nicer…
ps1.ord <- ordinate(ps1, "NMDS", "bray")
p2 = plot_ordination(ps1, ps1.ord, type="samples", color="sample_aspiration", shape="antibiotic") 
p2 + geom_polygon(aes(fill=SampleType)) + geom_point(size=5) + ggtitle("samples")

#==============================================================================
dim(psmelt(ps))
## ---------------------------------------------
## Floating barplot for replicate Min/max 
##---------------------------------------------
total_counts = as.data.frame(rowSums(otu_table(ps)))
colnames(total_counts) = "totals"

group_table = sample_data(ps) %>% select(group,Description) %>% unique

min_max_counts = left_join(add_rownames(sample_data(ps)), add_rownames(total_counts)) %>% 
  select(rowname, Description, group, totals) %>%
  group_by(Description) %>% 
  summarise(min=min(totals),max=max(totals)) %>% 
  left_join(group_table)

ggplot(min_max_counts, aes(x=Description,ymin = `min`, ymax = `max`,color=group)) + 
  geom_linerange(stat = 'identity') +
  xlab('Sample') + 
  ylab('Counts') +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background = element_blank())

