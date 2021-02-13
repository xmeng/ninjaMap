library("dada2"); packageVersion("dada2")
library("phyloseq"); packageVersion("phyloseq")

args <- commandArgs(trailingOnly = TRUE)
# set output paths
output_path  <- args[1]  #"/data/ampliseq/201014_16Sv4_Test2/DADA2_summary" #THIS SHOULD BE THE DADA2_summary folder
#output_path2 <- "/data/ampliseq/201014_16Sv4_Test2/DADA2_outputs" #THIS SHOULD BE THE DADA2_output folder

# create output directories
mainDir <- output_path
#dir.create( mainDir )
subDir <- "DADA2_summary"
dir.create(file.path(mainDir, subDir))

#mainDir <- "."
#subDir <- "DADA2_outputs"
#dir.create(file.path(mainDir, subDir))

# list samples
# the directory containing the fastq files
path <- args[2] # "/data/myBaseSpace/201014_16Sv4_Test2/fastq"
list.files(path)

#setwd("/data/ampliseq/Allison/data")
#list.files() # make sure what we think is here is actually here
## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
#samples is a file containing a list of samples
#samples <- scan("samples", what="character")

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME-XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_L001_R1_001.fastq.gz"), `[`, 1)

#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,180),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
#write.table(out, "Sample_stat.tsv", sep="\t", quote=F, col.names=NA)

# generate error model
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names


# Sample Inference    ="pseudo"
dadaFs <- dada(derep_forward, err=errF, pool=FALSE, multithread=TRUE)
dadaRs <- dada(derep_reverse, err=errR, pool=FALSE, multithread=TRUE)

dadaFs[[1]]

# merge PE reads
mergers <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
head (table(nchar(getSequences(seqtab))))


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# Track the number of reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, file.path(output_path,"Sample_stats.tsv"), sep="\t", quote=F, col.names=NA)


#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/data/ampliseq/silvaDB/silva_nr99_v138_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)

# no species level any more
#taxa <- addSpecies(taxa, "/data/ampliseq/silvaDB/silva_species_assignment_v138.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
#rownames(taxa.print) <- NULL
head(taxa.print)

 # giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

  # making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(output_path, "ASVs.fa" )  )

  # count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(output_path, "ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)

  # tax table:
  # creating table of taxonomy and setting any that are unclassified as "NA"
rownames(taxa) <- gsub(pattern=">", replacement="", x=asv_headers)
write.table(taxa, file.path(output_path, "ASVs_taxonomy_dada2.tsv"), sep = "\t", quote=F, col.names=NA)


###############################################################################
############## Output summary at each level of taxa    ########################
if (FALSE) {
  taxa <- assignTaxonomy(seqtab.nochim, "/data/ampliseq/silvaDB/silva_nr99_v138_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
  tax_tab <- tax_table(taxa)
  pseq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), tax_tab)


  otut<-otu_table(tax_glom(pseq,taxrank = "Phylum"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Phylum"))
  write.table(taxt,  file.path(output_path2,"Phylum_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path2,"Phylum_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  tax_ph<- paste(a1,a2,sep=";")
  write.table(otut, file.path(output_path,"Phylum_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)

  #otut2 <- cbind(otut, purity=round(rowMaxs(otut)/rowSums(otut)*100, 1))
  #write.table(otut2, file.path(output_path,"Phylum_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)

  otut<-otu_table(tax_glom(pseq,taxrank = "Class"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Class"))
  write.table(taxt,  file.path(output_path2,"Class_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path2,"Class_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  a3<-as.character(a[,3])
  tax_ph<- paste(a1,a2,a3,sep=";")
  write.table(otut, file.path(output_path,"Class_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)


  otut<-otu_table(tax_glom(pseq,taxrank = "Order"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Order"))
  write.table(taxt,  file.path(output_path2,"Order_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path2,"Order_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  a3<-as.character(a[,3])
  a4<-as.character(a[,4])
  tax_ph<- paste(a1,a2,a3,a4,sep=";")
  write.table(otut, file.path(output_path,"Order_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)


  otut<-otu_table(tax_glom(pseq,taxrank = "Family"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Family"))
  write.table(taxt,  file.path(output_path2,"Family_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path2,"Family_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  a3<-as.character(a[,3])
  a4<-as.character(a[,4])
  a5<-as.character(a[,5])
  tax_ph<- paste(a1,a2,a3,a4,a5,sep=";")
  write.table(otut, file.path(output_path,"Family_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)


  otut<-otu_table(tax_glom(pseq,taxrank = "Genus"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Genus"))
  write.table(taxt,  file.path(output_path2,"Genus_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path2,"Genus_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  a3<-as.character(a[,3])
  a4<-as.character(a[,4])
  a5<-as.character(a[,5])
  a6<-as.character(a[,6])
  tax_ph<- paste(a1,a2,a3,a4,a5,a6,sep=";")
  write.table(otut, file.path(output_path,"Genus_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)


  if (FALSE) {
      library(DECIPHER); packageVersion("DECIPHER")

      dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
      load("/data/ampliseq/silvaDB/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
      ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE) # use all processors

      ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
      asv_tax <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
      }))
      colnames(asv_tax) <- ranks
      rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

      #colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

      write.table(asv_tax, "ASVs_taxonomy_decipher.tsv", sep = "\t", quote=F, col.names=NA)


      #Removing likely contaminants
      library("phyloseq"); packageVersion("phyloseq")
      library(decontam); packageVersion("decontam")


      ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
                     tax_table(taxa))
      ps
  }
}
