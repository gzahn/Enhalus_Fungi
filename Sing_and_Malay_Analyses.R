# --------------------------------------------------------------#
# Fungal Community Analyses - Depends on "Process_Raw_Reads.R"
# Author: Geoffrey Zahn
# --------------------------------------------------------------#

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
ps1 = readRDS("./Output/MUX8041_clean_phyloseq_object.RDS")
ps2 = readRDS("./Output/MUX8046_clean_phyloseq_object.RDS")

ps = merge_phyloseq(ps1,ps2)

names(sample_data(ps))[5] <- "Source"

head(sample_data(ps))

# Rename from actual seqs to generic ESVs
    # First, save lookup table
    ESVcount = length(unlist(dimnames(ps@tax_table@.Data)[1]))
    ESV_ID_to_Sequence = data.frame(ESV = paste0("ESV_",1:ESVcount), Sequence = dimnames(ps@tax_table@.Data)[[1]])
    
    # Now, rename the phyloseq object attributes
dimnames(ps@tax_table@.Data)[[1]] <- paste0("ESV_",1:ESVcount)
dimnames(ps@otu_table@.Data)[[2]] <- paste0("ESV_",1:ESVcount)


# Clean up table ####
# Remove non-fungi
ps <- subset_taxa(ps, Kingdom == "k__Fungi")

# remove PNA-test
ps <- subset_samples(ps, Source != "Coral")

# remove empty and low-count ESVs
ps <- prune_taxa(taxa_sums(ps) >= 100, ps) # remove taxa with fewer than 100 occurences 

# remove newly-emptied samples
ps <- prune_samples(sample_sums(ps) != 0, ps)
# Sanity check
which(sample_sums(ps) == 0)
summary(taxa_sums(ps))
summary(sample_sums(ps))

# cleanup objects
rm(ps1);rm(ps2)

# Explore data ####
levels(ps@sam_data@.Data[[4]])
plot(taxa_sums(ps))


# Separate by species ####
halophila = subset_samples(ps, Species == "Halophila ovalis")
enhalus = subset_samples(ps, Species == "Enhalus acoroides")


# save levels of sam_data
source = levels(halophila@sam_data$Source)
location = levels(halophila@sam_data$Location)
source2 = levels(enhalus@sam_data$Source)
location2 = levels(enhalus@sam_data$Location)


# merge by island AND Structure ####
variable1 = as.character(get_variable(halophila, "Source"))
variable2 = as.character(get_variable(halophila, "Location"))
sample_data(halophila)$NewPastedVar <- mapply(paste, variable1, variable2, collapse = "_", sep = "_")
halophila_m = merge_samples(halophila, "NewPastedVar")

variable1 = as.character(get_variable(enhalus, "Source"))
variable2 = as.character(get_variable(enhalus, "Location"))
sample_data(enhalus)$NewPastedVar <- mapply(paste, variable1, variable2, collapse = "_", sep = "_")
enhalus_m = merge_samples(enhalus, "NewPastedVar")


# repair values of source and location
halophila_m@sam_data$Source = rep(source, each=7)
halophila_m@sam_data$Location = rep(location,4)
enhalus_m@sam_data$Source = rep(source2,each=6)
location2 = location2[c(1,2,3,6,4,5)]
enhalus_m@sam_data$Location = rep(location2,4)

# clean up new zeros
halophila_m <- prune_samples(sample_sums(halophila_m) != 0, halophila_m)
enhalus_m <- prune_samples(sample_sums(enhalus_m) != 0, enhalus_m)

# convert to relabund
halophila_m = transformSampleCounts(halophila_m, function(OTU) OTU/sum(OTU))
enhalus_m = transformSampleCounts(enhalus_m, function(OTU) OTU/sum(OTU))


# heatmap ####
heatmap(as.matrix(as.data.frame(otu_table(halophila_m))), col = gray.colors(100))
heatmap(as.matrix(as.data.frame(otu_table(enhalus_m))), col = gray.colors(100))


# Barplots ####
source("./plot_bar2.R")

# phylum relative abundance
plot_bar2(halophila_m, fill = "Phylum", x= "Location") + facet_grid(~Source)
ggsave("./Output/halophila_phylum_barplot.png", dpi = 300, width = 10, height = 9)

plot_bar2(enhalus_m, fill = "Phylum", x= "Location") + facet_grid(~Source)
ggsave("./Output/enhalus_phylum_barplot.png", dpi = 300, width = 10, height = 9)

# class relative abundance
plot_bar2(halophila_m, fill = "Class", x= "Location") + facet_grid(~Source)
ggsave("./Output/halophila_class_barplot.png", dpi = 300, width = 12, height = 9)

plot_bar2(enhalus_m, fill = "Phylum", x= "Location") + facet_grid(~Source)
ggsave("./Output/enhalus_class_barplot.png", dpi = 300, width = 12, height = 9)


# convert full phyloseq objects to relative abundance ####
halophila = transformSampleCounts(halophila, function(OTU) OTU/sum(OTU))
enhalus = transformSampleCounts(enhalus, function(OTU) OTU/sum(OTU))

# Mantel Test ####
spatial.dist = dist(cbind(halophila@sam_data$Lon_E, halophila@sam_data$Lat_N))
comm.dist = vegdist(as.matrix(halophila@otu_table))
mantel.test.halophila = mantel.rtest(spatial.dist, comm.dist, nrepet = 999)

ggplot(mapping = aes(x=jitter(spatial.dist,amount=1), y=comm.dist)) +
  geom_point(alpha=.05) + stat_smooth(method = "lm") + 
  labs(x="Spatial Distance",y="Community Distance") + theme_bw()
ggsave("./Output/Halophila_Mantel_Plot.png", dpi=300)

spatial.dist = dist(cbind(enhalus@sam_data$Lon_E, enhalus@sam_data$Lat_N))
comm.dist = vegdist(as.matrix(enhalus@otu_table))
mantel.test.enhalus = mantel.rtest(spatial.dist, comm.dist, nrepet = 999)

ggplot(mapping = aes(x=jitter(spatial.dist, amount = 1), y=comm.dist)) +
  geom_point(alpha=.05) + stat_smooth(method = "lm") +
  labs(x="Spatial Distance",y="Community Distance") + theme_bw()
ggsave("./Output/Enhalus_Mantel_Plot.png", dpi=300)


sink("./Output/mantel_tests.txt")
print("Halophila")
print(mantel.test.halophila)
print(noquote(" "))
print("Enhalus")
print(mantel.test.enhalus)
sink(NULL)



# Metric multi-dimensional scaling ####

halophila_wcmd = as.data.frame(wcmdscale(vegdist(otu_table(halophila)), 2, add=TRUE))
ggplot(halophila_wcmd, aes(x=V1,y=V2,color=halophila@sam_data$Source,shape=halophila@sam_data$Country)) + 
  geom_point(size = 3) + theme_bw() +
  labs(color="Source",shape="Country")
ggsave("./Output/WCMD_halophila_source-country.png", dpi = 300, height = 10, width = 12)


enhalus_wcmd = as.data.frame(wcmdscale(vegdist(otu_table(enhalus)), 2, add = TRUE))
ggplot(enhalus_wcmd, aes(x=V1,y=V2,color=enhalus@sam_data$Source,shape=enhalus@sam_data$Country)) + 
  geom_point(size=3) + theme_bw() +
  labs(color="Source",shape="Country")
ggsave("./Output/WCMD_enhalus_source-country.png", dpi = 300, height = 10, width = 12)

# Country is color and Source is shape
ggplot(halophila_wcmd, aes(x=V1,y=V2,color=halophila@sam_data$Country,shape=halophila@sam_data$Source)) + 
  geom_point(size = 3) + theme_bw() +
  labs(color="Region",shape="Structure")
ggsave("./Output/WCMD_halophila_source-country_2.png", dpi = 300, height = 10, width = 12)


ggplot(enhalus_wcmd, aes(x=V1,y=V2,color=enhalus@sam_data$Country,shape=enhalus@sam_data$Source)) + 
  geom_point(size=3) + theme_bw() +
  labs(color="Region",shape="Structure")
ggsave("./Output/WCMD_enhalus_source-country_2.png", dpi = 300, height = 10, width = 12)


# permANOVA ####

sink("./Output/PermANOVA_Tables.txt")
print("Halophila")
adonis(otu_table(halophila) ~ halophila@sam_data$Location * halophila@sam_data$Source)
print("Enhalus")
adonis(otu_table(enhalus) ~ enhalus@sam_data$Location * enhalus@sam_data$Source)
sink(NULL)

