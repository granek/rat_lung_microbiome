library(tidyr)
library(plyr); library(dplyr)

# Load Sample File
f <- file("160613_McKenney_MF_NJK_06012016_16SRV.TXT", open="r", encoding="UTF-16LE")
df <- read.table(f, sep='\t', dec=',', header=TRUE,comment="")

# Load Sample Key
sample_key <- read.csv("metadata_key.csv", header=TRUE)


# Split sample ID
sep_id = df %>% rename(SampleID=X.SampleID) %>%
  extract(col=SampleID,into=c("group","animal"),regex="([[:alpha:]]+)([[:digit:]]+)",remove = FALSE)

sep_id %>% mutate(antibiotic = mapvalues(group, 
                                         c("NLU",  "ASNLU", "AINLU", "APNLU"), 
                                         c("none", "subq", "iv", "oral"))) 
                    
                    
#                     SampleID)
  
  
plyr::mapvaluesx <- c("a", "b", "c")
mapvalues(x, c("a", "c"), c("A", "C"))

list(c("a",1),c("b",2))
unique(sep_id$group)
