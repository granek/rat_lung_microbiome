library(tidyr)
library(plyr); library(dplyr)

# Set up directory and files
info_dir = "notes_and_info"
raw_sample_map = file.path(info_dir, "160613_McKenney_MF_NJK_06012016_16SRV.TXT")
raw_sample_key = file.path(info_dir, "metadata_key.csv")
final_sample_map = file.path(info_dir, "rat_lung_map.tsv")

# need to jump through some hoops to open file with UTF-16LE encoding
f <- file(raw_sample_map, open="r", encoding="UTF-16LE")
df <- read.table(f, sep='\t', dec=',', header=TRUE,comment="")
close(f)

# Load Sample Key
sample_key <- read.csv(raw_sample_key, header=TRUE)

# Split sample ID
sep_id = df %>% rename(SampleID=X.SampleID) %>%
  extract(col=SampleID,into=c("group","animal"),regex="([[:alpha:]]+)([[:digit:]]+)",remove = FALSE)

sep_id$techrep = as.numeric(duplicated(sep_id$SampleID))+1

# sep_id %>% mutate(SampleID = paste(SampleID,techrep,sep="_")) %>% head

map_df = sep_id %>% mutate(antibiotic = mapvalues(group, 
                                         c("NLU",  "ASNLU", "AINLU", "APNLU"), 
                                         c("none", "subq", "iv", "oral"))) %>%
  mutate(SampleID = paste(SampleID,techrep,sep=".")) %>%
  select(SampleID, BarcodeSequence, LinkerPrimerSequence, 
         techrep, group, animal, antibiotic, 
         BarcodePlate, Well, 
         Description)




write.table(map_df,final_sample_map,quote=FALSE,row.names = FALSE,sep="\t")
                    
# unique(sep_id$group)

