suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("vegan"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("tibble"))

ShowFigure = function(plotfile, dummy_figure=FALSE){
  if (file.exists(plotfile)){
    fig_markdown = paste0("![](", plotfile, ")")
  } else if (dummy_figure) {
    # make a dummy figure if the real one doesn't exist so pandoc doesn't barf
    ggplot() + ggtitle(paste("placeholder for:", plotfile))
    ggsave(plotfile)
    fig_markdown = paste0("![](", plotfile, ")")
  } else {
    fig_markdown = paste("> Image doesn't exist:", plotfile)
  }
  cat(fig_markdown, fill=TRUE)
}

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

GenerateLefseOutput = function(ps,output_columns,outfile){
  base_columns = c("SampleID","OTU","Abundance")
  spread.df = psmelt(ps) %>% 
    mutate(SampleID = str_replace(SampleID, pattern="\\.", replacement="_")) %>%
    mutate(taxonomy = paste(Kingdom,Phylum,Class,Order,Family,Genus, sep="|")) %>%
    select(one_of(c(base_columns,output_columns))) %>%
    spread(OTU, Abundance)
  
  RepseqToTaxa(spread.df, ps) %>%
    write.table(file=outfile, 
                sep="\t", quote = FALSE,
                row.names = FALSE)
}


PruneTaxa <- function(ps, min_counts = 2, sample_proportion = 0.01) {
  min_samples = ceiling(sample_proportion*nsamples(ps))
  wh0 = genefilter_sample(ps, 
                          filterfun_sample(function(x) x >= min_counts), 
                          A=min_samples)
  
  ps.pruned = prune_taxa(wh0, ps)
  
  cat(paste(c("The data was pruned before ordination to remove rare taxa.",
              "To be included in the pruned dataset, a taxon must occur at least",
              min_counts,  "times in at least", sample_proportion*100, 
              "% of samples (", min_samples, "samples).\n", 
              "Before pruning there are", ntaxa(ps), "taxa, after pruning there are", 
              ntaxa(ps.pruned), "taxa remaining.")))
  return(ps)
}


TransformCounts <- function(ps) {
  ## Transform to even sampling depth.
  even.ps = transform_sample_counts(ps, function(x) 1E6 * x/sum(x))
  
  ## Remove samples with NaN (samples with all zero rows generate NaN when transformed above because of division by zero)
  even.ps = prune_samples(complete.cases(otu_table(even.ps)), even.ps)
  return(even.ps)
}


NMDSPlot <- function(ps,grouping="aspiration", nmds_seed=1) {
  set.seed(nmds_seed)
  ps.nmds <- ordinate(ps, "NMDS", "bray")
  ps.nmds.plot = plot_ordination(ps, ps.nmds, type="samples")
  nmds.ggplot = ggplot(ps.nmds.plot$data, aes(NMDS1, NMDS2)) +
    theme_classic() +
    geom_point(aes_string(color = grouping)) +
    annotate("text",x=-Inf,y=-Inf,hjust=0,vjust=0,
             label= paste("Stress:", ps.nmds$stress, 
                          "\nConverged:", ps.nmds$converged))
  
  cat(paste(c("- **Stress:**", ps.nmds$stress, 
              "\n**Converged:**",  ps.nmds$converged)))
  return(nmds.ggplot)
}

RunPermanova <- function(even.ps, samvar) {
  even.otus <- otu_table(even.ps)
  even.bray <- vegdist(even.otus, method="bray")
  even.adonis = adonis(even.bray ~ sample_data(even.ps)[[samvar]])
  
  even.beta <- betadisper(even.bray, sample_data(even.ps)[[samvar]])
  even.beta.permute = permutest(even.beta)
  
  return(list(adonis=even.adonis, beta=even.beta.permute))
}

