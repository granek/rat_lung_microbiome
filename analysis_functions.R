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
