lefse_sample = read.delim("workspace/lefse_workspace/input/hmp_aerobiosis_small.txt",
                          header = FALSE, stringsAsFactors=FALSE)
lefse_sample %>% filter(row_number()<=2) %>%
  column_to_rownames("V1") %>%
  t %>%
  as_data_frame %>%
  table
