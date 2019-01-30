# -------------------------------------------------------------------------------#
# Fungal Community Analyses - Depends on "Process_Raw_Reads.R" and "plot_bar2.R
# Author: Geoffrey Zahn
# -------------------------------------------------------------------------------#

# Load packages
library(phyloseq)
library(DESeq2)
library(DECIPHER)
library(ggplot2); packageVersion("ggplot2")
library(gridExtra)
library(ggpubr)
library(stringr)
library(purrr)
library(dplyr)
library(vegan)
library(MASS)
library(ade4)
library(lme4)

#Load data ####
enhalus = readRDS("./Output/MUX8046_clean_phyloseq_object.RDS")

names(sample_data(enhalus))[5] <- "Source"

head(sample_data(enhalus))

# Rename from actual seqs to generic ESVs
    # First, save lookup table
    ESVcount = length(unlist(dimnames(enhalus@tax_table@.Data)[1]))
    ESV_ID_to_Sequence = data.frame(ESV = paste0("ESV_",1:ESVcount), Sequence = dimnames(enhalus@tax_table@.Data)[[1]])
    
    # Now, rename the phyloseq object attributes
dimnames(enhalus@tax_table@.Data)[[1]] <- paste0("ESV_",1:ESVcount)
dimnames(enhalus@otu_table@.Data)[[2]] <- paste0("ESV_",1:ESVcount)


# Clean up table ####
# Remove non-fungi
enhalus <- subset_taxa(enhalus, Kingdom == "k__Fungi")

# remove PNA-test
enhalus <- subset_samples(enhalus, Source != "Coral")

# remove empty and low-count ESVs
enhalus <- prune_taxa(taxa_sums(enhalus) >= 100, enhalus) # remove taxa with fewer than 100 occurences 

# remove newly-emptied samples
enhalus <- prune_samples(sample_sums(enhalus) != 0, enhalus)


# save levels of sam_data
source2 = levels(enhalus@sam_data$Source)
location2 = levels(enhalus@sam_data$Location)


# merge by island AND Structure ####
variable1 = as.character(get_variable(enhalus, "Source"))
variable2 = as.character(get_variable(enhalus, "Location"))
sample_data(enhalus)$NewPastedVar <- mapply(paste, variable1, variable2, collapse = "_", sep = "_")
enhalus_m = merge_samples(enhalus, "NewPastedVar")


# repair values of source and location
enhalus_m@sam_data$Source = rep(source2,each=6)
location2 = location2[c(1,2,3,6,4,5)]
enhalus_m@sam_data$Location = rep(location2,4)

# clean up new zeros
enhalus_m <- prune_samples(sample_sums(enhalus_m) != 0, enhalus_m)

# convert to relabund
enhalus_m = transformSampleCounts(enhalus_m, function(OTU) OTU/sum(OTU))


# heatmap ####
heatmap(as.matrix(as.data.frame(otu_table(enhalus_m))), col = gray.colors(100))


# Barplots ####
source("./plot_bar2.R")

# phylum relative abundance
plot_bar2(enhalus_m, fill = "Phylum", x= "Location") + facet_grid(~Source)
ggsave("./Output/enhalus_phylum_barplot.png", dpi = 300, width = 10, height = 9)

# class relative abundance
plot_bar2(enhalus_m, fill = "Class", x= "Location") + facet_grid(~Source)
ggsave("./Output/enhalus_class_barplot.png", dpi = 300, width = 12, height = 9)


# convert full phyloseq objects to relative abundance ####
enhalus = transformSampleCounts(enhalus, function(OTU) OTU/sum(OTU))

# Mantel Test ####
spatial.dist = dist(cbind(enhalus@sam_data$Lon_E, enhalus@sam_data$Lat_N))
comm.dist = vegdist(as.matrix(enhalus@otu_table))
mantel.test.enhalus = mantel.rtest(spatial.dist, comm.dist, nrepet = 999)

ggplot(mapping = aes(x=jitter(spatial.dist, amount = 1), y=comm.dist)) +
  geom_point(alpha=.05) + stat_smooth(method = "lm") +
  labs(x="Spatial Distance",y="Community Distance") + theme_bw()
ggsave("./Output/Enhalus_Mantel_Plot.png", dpi=300)


sink("./Output/mantel_tests.txt")
print("Enhalus")
print(mantel.test.enhalus)
sink(NULL)

# Weighted Metric multi-dimensional scaling ####
wcmd = (wcmdscale(vegdist(otu_table(enhalus)), 2, add = TRUE,eig = TRUE))
enhalus_wcmd = as.data.frame(scores(wcmd))

# Goodness of Fit for WCMD
sink("./Output/WCMD_Goodness-of-fit.txt")
print(noquote(" Goodness-of-fit statistic"))
print(noquote(""))
wcmd$GOF[1]
sink(NULL)

# Plotting
ggplot(enhalus_wcmd, aes(x=Dim1,y=Dim2,color=enhalus@sam_data$Source,shape=enhalus@sam_data$Country)) + 
  geom_point(size=3) + theme_bw() +
  labs(color="Source",shape="Country") +
  theme(legend.text = element_text(size=14), legend.title = element_text(size=16,face = "bold"))
ggsave("./Output/WCMD_enhalus_source-country.png", dpi = 300, height = 10, width = 12)

# Country is color and Source is shape
ggplot(enhalus_wcmd, aes(x=Dim1,y=Dim2,color=enhalus@sam_data$Country,shape=enhalus@sam_data$Source)) + 
  geom_point(size=3) + theme_bw() +
  labs(color="Region",shape="Structure") +
  theme(legend.text = element_text(size=14), legend.title = element_text(size=16,face = "bold"))
ggsave("./Output/WCMD_enhalus_source-country_2.png", dpi = 300, height = 10, width = 12)

# Sample Location (Site) is color
ggplot(enhalus_wcmd, aes(x=Dim1,y=Dim2,color=enhalus@sam_data$Location,shape=enhalus@sam_data$Source)) + 
  geom_point(size=3) + theme_bw() +
  labs(color="Site",shape="Structure") +
  theme(legend.text = element_text(size=14), legend.title = element_text(size=16,face = "bold"))
ggsave("./Output/WCMD_enhalus_site-structure.png", dpi=300, height = 10, width = 12)

# permANOVA ####
sink("./Output/PermANOVA_Tables.txt")
print("Enhalus PermANOVA Results Table")
adonis(otu_table(enhalus) ~ enhalus@sam_data$Location * enhalus@sam_data$Source)
sink(NULL)

