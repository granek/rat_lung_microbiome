# Based on <http://benjjneb.github.io/dada2/R/tutorial.html>
## [josh@hardac-login ~]$ srun --partition interactive --nodelist hardac-node01-1 Rscript ~/collabs/SeedLab/hartwell/dada/dada_for_simreads.R /data/gems/seed/hartwell/sim_reads_gz > /data/gems/seed/hartwell/sim_reads_gz/dada_results_sim_reads_gz.txt

##====================================================================
#+ Setup: Directories, include=FALSE
# Directories

args <- commandArgs(trailingOnly = TRUE)

# if (interactive()){
#   basedir<<-"XXXXXX"
# }
if (length(args) == 1){
  basedir<<-args[1]
} else {
  basedir<<-"."
}

workdir = "workspace"
fastqdir = file.path(workdir,"split_fastqs")
taxonomy_dir =file.path(workdir,"taxonomy_refs")
greengenes_ref = file.path(taxonomy_dir,"ggtrain_97.fa")
silva_ref = file.path(taxonomy_dir,"silva.bac.train.fa")

results_dir = file.path(workdir,"results")
filtered_fastq_dir = file.path(workdir, "filtered_fastqs")

dir.create(filtered_fastq_dir, showWarnings = FALSE,recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE,recursive = TRUE)

map_file = file.path(basedir, "notes_and_info/rat_lung_map.tsv")
psfile.prefix = file.path(results_dir, "mouse_csection_ps")

##====================================================================
testFilesAndDirs = function(path){
  if (!file.exists(path)){
    stop(paste("Path doesn't exist:", path))
    }
}

testFilesAndDirs(basedir)
testFilesAndDirs(fastqdir)
testFilesAndDirs(greengenes_ref)
testFilesAndDirs(silva_ref)
testFilesAndDirs(results_dir)
testFilesAndDirs(filtered_fastq_dir)
testFilesAndDirs(map_file)

##====================================================================
#+ Setup: Load Libraries, include=FALSE
print(.libPaths())
library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
# library(dplyr)
# library(biom)

sessionInfo()
##====================================================================



#+ Filtering and Trimming, include=FALSE
# First we read in the file names for all the fastq files and do a little string
# manipulation to get lists of the forward and reverse fastq files in matched order:

findFastqs = function(fastqdir,fastq.suffix=".fq.gz",filter=""){

  #+Sample Data
  fns <- list.files(fastqdir)
  # print(paste("fastqdir:", fastqdir))
  # print(paste("All FASTQs:", fns))
  fastqs <- fns[grepl(paste0(fastq.suffix,"$"), fns)]
  if (filter != ""){
        fastqs <- fastqs[grepl(filter, fastqs)]
  }
  fastqs <- sort(fastqs) # Sort should keep them paired in order

  fnFs <- fastqs[grepl("r1.sample_", fastqs)]
  fnRs <- fastqs[grepl("r2.sample_", fastqs)]
  return(list(fnFs,fnRs))
}

#+ Quality Score Visualization, include=FALSE
visualizeQuality = function(forward_fastqs, reverse_fastqs){
  
  #+Examine quality profiles of forward and reverse reads
  # It is always important to look at your data. There are many ways to do this, 
  # but here we use a visualization from the ShortRead package.
  
  # Visualize the quality profile of the forward reads:
  # for(fnF in fnFs[1:1]) {
  for(fnF in forward_fastqs) {
    qqF <- qa(paste0(fastqdir,"/",fnF))[["perCycle"]]$quality
    qqF.plot = ShortRead:::.plotCycleQuality(qqF, main="Forward")
    # trellis.device(device="pdf", 
    pdf(file=file.path(fastqdir,"/",paste0(fnF,"_quality.pdf")))
    print(qqF.plot)
    dev.off()
    # ggsave(file.path(fastqdir,paste0(fnF,"_quality.pdf")),plot=qqF.plot)
  }
  
  # Visualize the quality profile of the reverse reads:
  # for(fnR in fnRs[1:2]) {
  for(fnR in reverse_fastqs) {
    qqR <- qa(paste0(fastqdir,"/",fnR))[["perCycle"]]$quality
    qqR.plot = ShortRead:::.plotCycleQuality(qqR, main="Reverse")
    pdf(file=file.path(fastqdir,"/",paste0(fnR,"_quality.pdf")))
    print(qqR.plot)
    dev.off()
  }
}


#+ Process And Analyze FASTQs, include=FALSE
filterFASTQs = function(fnFs,fnRs,fastqdir,filtered_fastq_dir,
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
  filtFs <- paste0(filtered_fastq_dir, "/", gsub(".fq.gz","",fnFs), "_filt.fastq.gz")
  filtRs <- paste0(filtered_fastq_dir, "/", gsub(".fq.gz","",fnRs), "_filt.fastq.gz")
  # filtFs <- paste0(tmpdir, sapply(strsplit(fnFs, "\\."), `[`, 1), "_filt.fastq.gz")
  # filtRs <- paste0(tmpdir, sapply(strsplit(fnRs, "\\."), `[`, 1), "_filt.fastq.gz")
  for(i in seq_along(fnFs)) {
    fastqPairedFilter(file.path(fastqdir, c(fnFs[i], fnRs[i])), c(filtFs[i], filtRs[i]), 
                      maxN=0, maxEE=2, truncQ=2, trimLeft=c(F_trimLeft, R_trimLeft), 
                      truncLen=c(F_truncLen,R_truncLen), compress=TRUE, verbose=TRUE)
  }
  return(list(fnFs,fnRs))
}

processFastqs = function(filtFs,filtRs){
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
  derepFs <- lapply(filtFs, derepFastq, verbose=TRUE)
  derepRs <- lapply(filtRs, derepFastq, verbose=TRUE)
  
  # Name the derep-class objects by the sample names
  sam_names <- gsub("r1.", "", basename(filtFs))
  sam_names = gsub("_filt.fastq.gz", "", sam_names)
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
  ##---------------------------------------
  # Make "otu_table"
  ##---------------------------------------
  # seqs <- colnames(seqtab)
  otab <- otu_table(seqtab, taxa_are_rows=FALSE)
  colnames(otab) <- paste0("Seq", seq(ncol(otab)))
  # print(otab)
  
  ##---------------------------------------
  ## Make taxtable
  ##---------------------------------------
  # Assign taxonomy to dadatype sequences
  set.seed(random.seed)
  print("Assigning taxonomy")
  seqtab.tax <- assignTaxonomy(colnames(seqtab), referenceFasta, minBoot=50)
  
  seqtab.tax.list = strsplit(seqtab.tax,";") # split to list
  
  # pad sequences that are missing fine scale annotation
  max.len <- max(sapply(seqtab.tax.list, length))
  seqtab.tax.list.padded <- lapply(seqtab.tax.list, function(x) {c(x, rep("", max.len - length(x)))})

  # make 
  taxtab <- tax_table(do.call(rbind,seqtab.tax.list.padded))
  rownames(taxtab) = colnames(otab)
  colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[seq(ncol(taxtab))]
  
  ##---------------------------------------
  ## Make taxtable
  ##---------------------------------------
  map.df = read.delim(map_file)
  samdat = sample_data(map.df)
  sample_names(samdat) = paste0("sample_",map.df$SampleID)
  ps <- phyloseq(otab, samdat, taxtab)
  
  ##---------------------------------------
  ## Make seqID map
  ##---------------------------------------
  seqid.map.df = data.frame(sequence = colnames(seqtab),
                            row.names = paste0("Seq", seq(ncol(otab))))
  return(list(ps,seqid.map.df))
}

outputPhyloseq = function(ps,seqid.map.df,outfile.prefix){
  otu_table_file = paste0(outfile.prefix,"_otu.csv")
  sample_data_file = paste0(outfile.prefix,"_samdat.csv")
  tax_table_file = paste0(outfile.prefix,"_tax.csv")
  seqid_map_file = paste0(outfile.prefix,"_seq.csv")
  
  write.csv(otu_table(ps), file=otu_table_file)
  write.csv(sample_data(ps), file=sample_data_file)
  write.csv(tax_table(ps), file=tax_table_file)
  write.csv(seqid.map.df, file=seqid_map_file)
  
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
F_trimLeft=0;R_trimLeft=0;F_truncLen=251;R_truncLen=251
##-----------------------------------------------------------------------------------
print(paste("Trimming Parameters:",
            "F_trimLeft:", F_trimLeft,
            "R_trimLeft:", R_trimLeft,
            "F_truncLen:", F_truncLen,
            "R_truncLen:", R_truncLen))
separator = paste0("\n", paste(replicate(70, "-"), collapse=""),"\n")



ptm <- proc.time()
fastqs = findFastqs(fastqdir)
print(fastqs)


if (length(fastqs[[1]])==0){
  print(paste("No FASTQs in: ", fastqdir))
  next
}
for_fastqs = fastqs[[1]]
rev_fastqs = fastqs[[2]]

MakeQualityPlots = FALSE
FilterFASTQs = FALSE
# SampleSubsetString = "sample_2"
SampleSubsetString = ""

if (MakeQualityPlots) {
  visualizeQuality(for_fastqs, rev_fastqs)
}

if (FilterFASTQs) {
  filtered.fastqs = filterFASTQs(for_fastqs,rev_fastqs,fastqdir,filtered_fastq_dir,
                                 F_trimLeft=F_trimLeft,R_trimLeft=R_trimLeft,
                                 F_truncLen=F_truncLen,R_truncLen=R_truncLen)
} else {
  filtered.fastqs = findFastqs(filtered_fastq_dir,"_filt.fastq.gz",filter=SampleSubsetString)
}
filtFs = file.path(filtered_fastq_dir,filtered.fastqs[[1]])
filtRs = file.path(filtered_fastq_dir,filtered.fastqs[[2]])

print(filtered.fastqs)

seqtab = processFastqs(filtFs,filtRs)
ps_and_seqid = makePhyloseq(seqtab,silva_ref,map_file)

ps = ps_and_seqid[[1]]
seqid.map.df = ps_and_seqid[[2]]

output.files = outputPhyloseq(ps,seqid.map.df,psfile.prefix)
# otu_table_file = output.files[[1]]
# sample_data_file = output.files[[2]]
# tax_table_file = output.files[[3]]

# ps.loaded = loadPhyloseqFiles(otu_table_file,sample_data_file,tax_table_file)
# outputPhyloseq(ps.loaded,paste0(psfile.prefix,"_loaded"))

# identical(ps,ps.loaded)
# phyloseqAnalysis(ps)
# loadSampleData(map_file)


print(proc.time() - ptm)
#==============================================================================
