# -----------------------------------------------------------------------------#
# Singapore and Malaysia Fungi Processing Raw Reads
# Author: Geoffrey Zahn
# Requirements: ITSx v 1.1b1
# -----------------------------------------------------------------------------#

# Prepare demultiplexed reads in Bash... ####
# Run ITSxpress on all fastq files 
# for i in ./Fwd/*.fastq.gz; do itsxpress --fastq $i  --outfile $i.FungalITS1.fastq.gz --region ITS1 --taxa Fungi -s --threads 10 --log $i.ITSx.log; done
# rename the files
    # rename -v -e 's/^(.{7}).*/$1.fastq.gz/' *.FungalITS1.fastq.gz
# No ITS1 found in REV reads


# Load packages ####
library(DESeq2)
library(DECIPHER)
library(decontam)
library(phangorn)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(gridExtra)
library(ggpubr)
library(dada2); packageVersion("dada2")
library(stringr)
library(purrr)

####################################################################################
############### NEED TO RUN DADA2 ON EACH ILLUMINA SET SEPARATELY! #################
####################################################################################



# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "./ITS1_Fastqs/fwd/MUX8046" 
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "fastq.gz"))
tail(fns)
sample.names <- sapply(strsplit(basename(fns), "_"), `[`, 1)


# visualize a couple of fwd read quality profiles to help you decide reasonable filtration parameters
plotQualityProfile(fns[3:4])


# Filter and trim ####
filts <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
# filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fns, filts, # fnRs, filtRs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=16) # On Windows set multithread=FALSE

out[order(out[,2]),]
no.pass = which(as.data.frame(out)$reads.out == 0)
passing_files = out[-no.pass,]


# learn error rates ####
# Since some samples had zero reads pass QC, reassign filtFs and filtRs
filts <- sort(list.files(filtpath, full.names = TRUE))

errF <- learnErrors(filts, multithread=TRUE, MAX_CONSIST = 20, nbases = 1e10)

# sanity check
plotErrors(errF, nominalQ=TRUE)


# Dereplication ####
derep <- derepFastq(filts, verbose=TRUE)

# Since some samples may have been removed (no reads passed QC), reassign sample.names
sample.names <- sapply(strsplit(basename(filts), "-"), `[`, 1)
names(derep) <- sample.names

# Sample inference ####
dadaFs <- dada(derep, err=errF, multithread=TRUE)

# Make a sequence table ####
seqtab <- makeSequenceTable(dadaFs)


# Remove Chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" and "seqtab.nochim" to remove missing reads
out = out[as.data.frame(out)$reads.out > 0,]
seqtab.nochim <- seqtab.nochim[sample.names,]
dada = dadaFs[sample.names]

# Track Reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[,1], sapply(dada, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$filter.loss = (track[,1]-track[,2])/track[,1]
write.csv(track, file = "./Output/MUX8046_read_counts_at_each_step.csv", row.names = TRUE)


# Save intermediate seqtable object
saveRDS(seqtab.nochim, "./Output/MUX8046_seqtab.nochim.RDS")

# Reload Point
seqtab.nochim = readRDS("./Output/MUX8046_seqtab.nochim.RDS")

# import metadata ####
meta = read.csv("./sample_mapping_data.csv")
meta = meta[meta$RunID == "MUX8046",]

row.names(meta) <- meta$IlluminaName
row.names(seqtab.nochim)

# reorder metadata
meta = meta[order(row.names(meta)),]

# Find controlsamples (extraction negatives) and clean ####
meta$controls <- meta$Species == "Blank"

# remove missing samples excluded due to poor QC 

row.names(seqtab.nochim) <- map(strsplit(row.names(seqtab.nochim), split = "_"),1)
good.samples = (row.names(meta) %in%  row.names(seqtab.nochim))
meta <- (meta[good.samples,])
rm(good.samples)


# find contaminants
# seqtab.nochim = readRDS(file = "./Output/clean_dada2_seqtable.RDS")
contams = isContaminant(seqtab.nochim, neg = meta$controls, normalize = TRUE)
table(contams$contaminant)  # No control samples passed QC...only 1 sequence made it through


# remove them
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$controls == FALSE,]
meta = meta[meta$controls == FALSE,]

# Remove all seqs with fewer than 100 nucleotides
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# Assign Taxonomy ####

# Save intermediate seqtab and re-load before assigning taxonomy to reduce virtual memory usage
saveRDS(seqtab.nochim, file = "./Output/MUX8046_clean_dada2_seqtable.RDS")
# rm(list = ls())
seqtab.nochim <- readRDS(file = "./Output/MUX8046_clean_dada2_seqtable.RDS")

taxa <- assignTaxonomy(seqtab.nochim, "./Taxonomy/UNITE_anthozoa_enhalus_halophila.fasta", multithread=20)

# Save intermediate files
write.csv(as.data.frame(seqtab.nochim), file = "./Output/MUX8046_SeqTable_no-chimera_no-contams_w_Taxonomy.csv", row.names = TRUE, quote = FALSE)
saveRDS(seqtab.nochim, file = "./Output/MUX8046_clean_dada2_seqtable_w_taxonomy.RDS")
saveRDS(taxa, file = "./Output/MUX8046_Taxonomy_from_dada2.RDS")
seqs <- getSequences(seqtab.nochim)

# Hand off to Phyloseq ####

seqtab.nochim = readRDS(file = "./Output/MUX8046_clean_dada2_seqtable_w_taxonomy.RDS")

taxa = readRDS(file = "./Output/MUX8046_Taxonomy_from_dada2.RDS")

unique(taxa[,2])

ps.MUX8046 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                       sample_data(meta), 
                       tax_table(taxa))


# Save RDS object for Phyloseq
saveRDS(ps.MUX8046, file = "./Output/MUX8046_clean_phyloseq_object.RDS")
