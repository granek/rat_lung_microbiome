# Based on <http://benjjneb.github.io/dada2/R/tutorial.html>
##====================================================================
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

read1_dir = file.path(workdir,"split_Undetermined_S0_L001_R1_001")
read2_dir = file.path(workdir,"split_Undetermined_S0_L001_R2_001")
taxonomy_dir =file.path(workdir,"taxonomy_refs")
# greengenes_ref = file.path(taxonomy_dir,"ggtrain_97.fa")
silva_ref = file.path(taxonomy_dir,"silva_nr_v123_train_set.fa.gz")

results_dir = file.path(workdir,"results")
dir.create(results_dir, showWarnings = FALSE,recursive = TRUE)

filtered_fastq_dir = file.path(workdir, "filtered_fastqs")
filtered_fastqs_stamp = file.path(filtered_fastq_dir,"filtered_fastq_STAMP")
qual_plot_dir = file.path(workdir, "qual_plots")

# psfile.prefix is the base name to use for naming phyloseq output files
psfile.prefix = file.path(results_dir, "rat_lung_ps")

##====================================================================
##====================================================================
#+ Setup: Load Libraries, include=FALSE
## print(.libPaths())
library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(phangorn)
library(msa)
# library(dplyr)
# library(biom)
## sessionInfo()
writeLines(capture.output(sessionInfo()), file.path(results_dir,"process_fastq_to_counts_sessionInfo.txt"))
##====================================================================
#+ Filtering and Trimming, include=FALSE
# First we read in the file names for all the fastq files and do a little string
# manipulation to get lists of the forward and reverse fastq files in matched order:

findFastqs = function(fastqdir,fastq.suffix=".fq.gz",filter="",remove=""){

  fns <- list.files(fastqdir)
  fastqs <- fns[grepl(paste0(fastq.suffix,"$"), fns)]
  if (filter != ""){
        fastqs <- fastqs[grepl(filter, fastqs)]
    }
  if (remove != ""){
      print(paste("removing fastqs with",remove,"in name"))
        fastqs <- fastqs[!grepl(remove, fastqs)]
  }
  fastqs <- sort(fastqs) # Sort should keep them paired in order

  return(file.path(fastqdir,fastqs))
}

## ##+ Quality Score Visualization, include=FALSE
visualizeQuality = function(forward_fastqs, reverse_fastqs,plot_dir,fastq_ext=".fastq"){
  
  #+Examine quality profiles of forward and reverse reads
  # It is always important to look at your data. There are many ways to do this, 
  # but here we use a visualization from the ShortRead package.
  
  # Visualize the quality profile of the forward reads:
  # for(fnF in fnFs[1:1]) {
  for(fnF in forward_fastqs) {
    qqF <- qa(fnF)[["perCycle"]]$quality
    qqF.plot = ShortRead:::.plotCycleQuality(qqF, main="Forward")
    # trellis.device(device="pdf", 
    pdf(file=file.path(plot_dir,gsub(fastq_ext,"_R1_quality.pdf",basename(fnF))))
    print(qqF.plot)
    dev.off()
    # ggsave(file.path(fastqdir,paste0(fnF,"_quality.pdf")),plot=qqF.plot)
  }
  
  # Visualize the quality profile of the reverse reads:
  # for(fnR in fnRs[1:2]) {
  for(fnR in reverse_fastqs) {
    qqR <- qa(fnR)[["perCycle"]]$quality
    qqR.plot = ShortRead:::.plotCycleQuality(qqR, main="Reverse")
    pdf(file=file.path(plot_dir,gsub(fastq_ext,"_R2_quality.pdf",basename(fnR))))
    print(qqR.plot)
    dev.off()
  }
}


#+ Process And Analyze FASTQs, include=FALSE
filterFASTQs = function(fnFs,fnRs,filtered_fastq_dir,
                        fastq_ext=".fastq",
                        F_trimLeft=10,
                        R_trimLeft=10,
                        F_truncLen=240,
                        R_truncLen=200){
  #+ Perform filtering and trimming
  
  # The trimming parameters were decided by inspecting the quality profiles. 
  # The filtering parameters we’ll use are standard: maxN=0 (DADA2 requires no Ns), 
  # truncQ=2 (quality score 2 in Illumina means “stop using this read”) and maxEE=2. 
  # The maxEE parameter sets the maximum number of “expected errors” allowed in a read. 
  # Setting a threshold on expected errors is a better filter than simply 
  # averaging quality scores. We use the fastqPairedFilter function to jointly 
  # filter the forward and reverse reads.
  
  # Filter the forward and reverse reads:
  
  ## > gsub(".fq.gz","","r2.sample_97.fq.gz")
  # filtFs <- paste0(filtered_fastq_dir, "/", gsub(".fq.gz","",fnFs), "_filt.fastq.gz")
  # filtRs <- paste0(filtered_fastq_dir, "/", gsub(".fq.gz","",fnRs), "_filt.fastq.gz")
  ## filtFs <- paste0(tmpdir, sapply(strsplit(fnFs, "\\."), `[`, 1), "_filt.fastq.gz")
  ## filtRs <- paste0(tmpdir, sapply(strsplit(fnRs, "\\."), `[`, 1), "_filt.fastq.gz")
  for_dir = dirname(fnFs[1])
  rev_dir = dirname(fnRs[1])
  ## print(fnFs)
  ## print(fnRs)
  ## print(fastq_ext)
  for(raw_for in fnFs) {
      read_base = basename(raw_for)
      raw_rev = file.path(rev_dir,read_base)
      filt_for = file.path(filtered_fastq_dir,gsub(fastq_ext,"_R1_filt.fastq.gz",read_base))
      filt_rev = file.path(filtered_fastq_dir,gsub(fastq_ext,"_R2_filt.fastq.gz",read_base))
      print(c(raw_for, raw_rev, filt_for, filt_rev))
      print("-------------------------------------")
      fastqPairedFilter(c(raw_for, raw_rev), c(filt_for, filt_rev),
                        maxN=0, maxEE=2, truncQ=2, trimLeft=c(F_trimLeft, R_trimLeft), 
                        truncLen=c(F_truncLen,R_truncLen), compress=TRUE, verbose=TRUE)
      ## fastqPairedFilter(file.path(fastqdir, c(fnFs[i], fnRs[i])), c(filtFs[i], filtRs[i]),
      ##                   maxN=0, maxEE=2, truncQ=2, trimLeft=c(F_trimLeft, R_trimLeft), 
      ##                   truncLen=c(F_truncLen,R_truncLen), compress=TRUE, verbose=TRUE)
  }
  return(list(fnFs,fnRs))
}

processFastqs = function(filtFs,filtRs,for_suffix="_R1_filt.fastq.gz",rev_suffix="_R2_filt.fastq.gz"){
  #+Dereplication
  # In the dereplication step, all reads with identical sequences are combined
  # into “unique sequences” with a corresponding abundance, i.e. the number of
  # reads with that unique sequence. Dereplication is a part of most pipelines 
  # because it reduces computation time by eliminating repeated comparisons of
  # identical sequences.
  
  # Dereplication in the DADA2 pipeline has one crucial addition: DADA2 retains 
  # a summary of the quality information associated with each unique sequence. 
  # DADA2 constructs a “consensus” quality profile for each unique sequence by 
  # averaging the positional qualities from the dereplicated reads. These 
  # consensus quality profiles inform the error model of the subsequent denoising
  # step, significantly increasing DADA2’s accuracy.
  
  #+ Dereplicate the filtered fastq files:
  filtRs = gsub(for_suffix, rev_suffix, filtFs)
  derepFs <- lapply(filtFs, derepFastq, verbose=TRUE)
  derepRs <- lapply(filtRs, derepFastq, verbose=TRUE)
  
  # Name the derep-class objects by the sample names
  sam_names = gsub(for_suffix, "", basename(filtFs))
  print(sam_names)
  # stop("check sam_names")
  names(derepFs) <- sam_names
  names(derepRs) <- sam_names
  
  # Inspect the derep-class object returned by derepFastq:
  derepFs[[1]]
  
  #+ Sample Inference
  
  # We are now ready to apply DADA2’s core sample inference algorithm to the dereplicated sequences.
  
  # Perform joint sample inference and error rate estimation (takes a few minutes):
  dadaFs <- dada(derepFs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
  dadaRs <- dada(derepRs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
  
  ## # need to make lists of one, since we only have one sample
  ## dadaFs = list(dadaFs)
  ## dadaRs = list(dadaRs)
  
  # Inspecting the dada-class object returned by dada:
  dadaFs[[1]]
  
  # Visualize estimated error rates:
  plotErrors(dadaFs[[1]], "A", nominalQ=TRUE)
  
  # Identify chimeras
  # Identify chimeric sequences:
  bimFs <- sapply(dadaFs, isBimeraDenovo, verbose=TRUE,simplify = FALSE)
  bimRs <- sapply(dadaRs, isBimeraDenovo, verbose=TRUE,simplify = FALSE)
  print(unname(sapply(bimFs, mean)), digits=2)
  print(unname(sapply(bimRs, mean)), digits=2)
  
  # Merge paired reads
  # Merge the denoised forward and reverse reads:
  mergers <- mapply(mergePairs, dadaFs, derepFs, dadaRs, derepRs, SIMPLIFY=FALSE)
  head(mergers[[1]])
  
  # Remove chimeras:
  mergers.nochim <- mapply(function(mm, bF, bR) mm[!bF[mm$forward] & !bR[mm$reverse],], 
                           mergers, bimFs, bimRs, SIMPLIFY=FALSE)
  
  # stop("fix for no chimeras")
  
  head(getUniques(mergers.nochim[[1]]), n=2)
  head(mergers.nochim[[1]][,c("sequence", "abundance")], n=2)
  
  #+ Constructing the sequence table
  # Construct sequence table:
  seqtab <- makeSequenceTable(mergers.nochim)
  print("table")
  print(table(nchar(colnames(seqtab))))
  print("length")
  print(length(unique(substr(colnames(seqtab), 1, 230))))
  print("dim")
  print(dim(seqtab))
  return(seqtab)
}

makePhyloseq = function(seqtab,referenceFasta,map_file,random.seed=100){
  # Assign taxonomy to dadatype sequences
  set.seed(random.seed)
  print("Assigning taxonomy")
  seqtab.tax <- assignTaxonomy(seqtab, referenceFasta, minBoot=50)
  colnames(seqtab.tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")[seq(ncol(seqtab.tax))]
  unname(head(seqtab.tax))
  
  ## Load Map
  map.df = read.delim(map_file)
  samdat = sample_data(map.df)
  rownames(samdat) = samdat$SampleID
  
  # Make "otu_table"
  ##---------------------------------------
  # seqs <- colnames(seqtab)
  # otab <- otu_table(seqtab, taxa_are_rows=FALSE)
  ##---------------------------------------------------------------------
  ## block below derived from http://f1000research.com/articles/5-1492/v1
  seqs <- getSequences(seqtab)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  mult <- msa(seqs, method="ClustalW", type="dna", order="input")

  ## The phangorn package is then used to construct a phylogenetic tree. Here we first construct a neighbor-joining tree, and then fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree as a starting point.

  phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data=phang.align)
  ## negative edges length changed to 0!
  
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
  ##---------------------------------------------------------------------
  
  
  ps <- phyloseq(tax_table(seqtab.tax),
                 samdat,
                 otu_table(seqtab, taxa_are_rows=FALSE), 
                 phy_tree(fitGTR$tree))
  
  # 
  # # colnames(otab) <- paste0("Seq", seq(ncol(otab)))
  # # print(otab)
  # 
  # ##---------------------------------------
  # ## Make taxtable
  # ##---------------------------------------
  # # Assign taxonomy to dadatype sequences
  # set.seed(random.seed)
  # print("Assigning taxonomy")
  # seqtab.tax <- assignTaxonomy(colnames(seqtab), referenceFasta, minBoot=50)
  # 
  # seqtab.tax.list = strsplit(seqtab.tax,";") # split to list
  # 
  # # pad sequences that are missing fine scale annotation
  # max.len <- max(sapply(seqtab.tax.list, length))
  # seqtab.tax.list.padded <- lapply(seqtab.tax.list, function(x) {c(x, rep("", max.len - length(x)))})
  # 
  # # make 
  # taxtab <- tax_table(do.call(rbind,seqtab.tax.list.padded))
  # # problem in next line
  # rownames(taxtab) = colnames(otab)
  # colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")[seq(ncol(taxtab))]
  # 
  # ##---------------------------------------
  # ## Make taxtable
  # ##---------------------------------------
  # map.df = read.delim(map_file)
  # samdat = sample_data(map.df)
  # # sample_names(samdat) = paste0("sample_",map.df$SampleID)
  # ps <- phyloseq(otab, samdat, taxtab)
  # 
  # ##---------------------------------------
  # ## Make seqID map
  # ##---------------------------------------
  # seqid.map.df = data.frame(sequence = colnames(seqtab),
  #                           row.names = paste0("Seq", seq(ncol(otab))))
  # return(list(ps,seqid.map.df))
  return(ps)
}

# outputPhyloseq = function(ps,seqid.map.df,outfile.prefix){
outputPhyloseq = function(ps,outfile.prefix){
  otu_table_file = paste0(outfile.prefix,"_otu.csv")
  sample_data_file = paste0(outfile.prefix,"_samdat.csv")
  tax_table_file = paste0(outfile.prefix,"_tax.csv")
  
  write.csv(otu_table(ps), file=otu_table_file)
  write.csv(sample_data(ps), file=sample_data_file)
  write.csv(tax_table(ps), file=tax_table_file)
  ## seqid_map_file = paste0(outfile.prefix,"_seq.csv")
  ## write.csv(seqid.map.df, file=seqid_map_file)
  
  return(list(otu_table_file,sample_data_file,tax_table_file))
}

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
# Files to use
##-----------------------------------------------------------------------------------
# F_trimLeft=10;R_trimLeft=10;F_truncLen=240;R_truncLen=200
# F_trimLeft=10;R_trimLeft=10;F_truncLen=150;R_truncLen=150
# F_trimLeft=0;R_trimLeft=0;F_truncLen=251;R_truncLen=251
# F_trimLeft=0;R_trimLeft=0;F_truncLen=251;R_truncLen=251
F_trimLeft=10;R_trimLeft=10;F_truncLen=140;R_truncLen=140

##-----------------------------------------------------------------------------------
print(paste("Trimming Parameters:",
            "F_trimLeft:", F_trimLeft,
            "R_trimLeft:", R_trimLeft,
            "F_truncLen:", F_truncLen,
            "R_truncLen:", R_truncLen))
separator = paste0("\n", paste(replicate(70, "-"), collapse=""),"\n")



ptm <- proc.time()

fastq_end = ".fastq"

# fnFs <- fastqs[grepl("r1.sample_", fastqs)]
for_fastqs = findFastqs(read1_dir,fastq_end,remove="Unassigned")
rev_fastqs = findFastqs(read2_dir,fastq_end,remove="Unassigned")

# FilterFASTQs = FALSE
# SampleSubsetString = "sample_2"
SampleSubsetString = ""

if (args$quality_plots > 0) {
  dir.create(qual_plot_dir, showWarnings = FALSE,recursive = TRUE)
  visualizeQuality(for_fastqs[1:args$quality_plots], 
                   rev_fastqs[1:args$quality_plots],
                   qual_plot_dir,fastq_end)
}

## print(paste("FILTER?:", args$filter_fastqs))
## print(read1_dir)
## print(for_fastqs)

if (args$filter_fastqs || !file.exists(filtered_fastqs_stamp)) {
  print("RUNNING FILTERING!")
  dir.create(filtered_fastq_dir, showWarnings = FALSE,recursive = TRUE)
  filtered.fastqs = filterFASTQs(for_fastqs,rev_fastqs,filtered_fastq_dir,
                                 fastq_ext=fastq_end,
                                 F_trimLeft=F_trimLeft,R_trimLeft=R_trimLeft,
                                 F_truncLen=F_truncLen,R_truncLen=R_truncLen)
  cat(Sys.time(), file = filtered_fastqs_stamp)
} else {
  print("Using existing filtered reads!")
}
filtFs = findFastqs(filtered_fastq_dir,"_filt.fastq.gz",filter="_R1_")
filtRs = findFastqs(filtered_fastq_dir,"_filt.fastq.gz",filter="_R2_")

## filtFs = file.path(filtered_fastq_dir,filtered.fastqs[[1]])
## filtRs = file.path(filtered_fastq_dir,filtered.fastqs[[2]])

print(filtFs)
print(filtRs)

seqtab = processFastqs(filtFs,filtRs)
# Next two lines for stepping through makePhyloseq
referenceFasta = silva_ref
random.seed=100

# ps_and_seqid = makePhyloseq(seqtab,silva_ref,map_file)
# ps = ps_and_seqid[[1]]
# seqid.map.df = ps_and_seqid[[2]]
ps = makePhyloseq(seqtab,silva_ref,map_file)
## plot_richness(ps, x="animal", measures=c("Shannon", "Simpson"), color="antibiotic") + theme_bw()
## plot_bar(ps, x="antibiotic", fill="Family") 
saveRDS(ps, paste0(psfile.prefix,".rds"))


# output.files = outputPhyloseq(ps,seqid.map.df,psfile.prefix)
output.files = outputPhyloseq(ps,psfile.prefix)
otu_table_file = output.files[[1]]
sample_data_file = output.files[[2]]
tax_table_file = output.files[[3]]
write("Phyloseq object saved as following files:",file="")
write(paste("otu_table_file:\t", otu_table_file),file="")
write(paste("sample_data_file:\t", sample_data_file),file="")
write(paste("tax_table_file:\t", tax_table_file), file="")



## ps.loaded = loadPhyloseqFiles(otu_table_file,sample_data_file,tax_table_file)
ps.loaded = readRDS(paste0(psfile.prefix,".rds"))

## # outputPhyloseq(ps.loaded,paste0(psfile.prefix,"_loaded"))
if (identical(ps,ps.loaded)) {
    print("YAY! Saved and reloaded ps is identical to original")
} else {
    print("*** ERROR: Saved and reloaded ps is DIFFERENT from original****")
}
## # phyloseqAnalysis(ps)
## # loadSampleData(map_file)

print(proc.time() - ptm)
#==============================================================================
