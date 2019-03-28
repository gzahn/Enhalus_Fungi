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
library(VennDiagram)
library(ggbiplot)


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

# Edit taxon names
enhalus_m@tax_table@.Data[,"Phylum"] <- as.character(gsub("p__","",enhalus_m@tax_table@.Data[,"Phylum"]))
enhalus_m@tax_table@.Data[,"Class"] <- as.character(gsub("c__","",enhalus_m@tax_table@.Data[,"Class"]))


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

library(ade4)
citation("ade4")
ggplot(mapping = aes(x=jitter(spatial.dist, amount = 1), y=comm.dist)) +
  geom_point(alpha=.05) + stat_smooth(method = "lm") +
  labs(x="Spatial Distance",y="Community Distance") + theme_bw()
ggsave("./Output/Enhalus_Mantel_Plot.png", dpi=300)

wcm
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

pc = princomp(t(decostand(otu_table(enhalus), method = "total", MARGIN = 1)),cor = TRUE,scores = TRUE)
pc.names = (which(abs(rowSums(pc$scores)) > 100))


ggbiplot(pc,labels.size = .1,) + lims(x=c(-5,5),y=c(-5,5)) 

# permANOVA ####
sink("./Output/PermANOVA_Tables.txt")
print("Enhalus PermANOVA Results Table")
adonis(otu_table(enhalus) ~ enhalus@sam_data$Location * enhalus@sam_data$Source)
sink(NULL)


# Heatmap ####
source("./heatmap_left.R")
hm.pal = c("#ffffff",RColorBrewer::brewer.pal(8,"OrRd"))

otus = as(otu_table(enhalus), "matrix")
otus.df = as.data.frame(otus)
orders = str_remove(tax_table(enhalus)[,3], "c__")
orders[is.na(orders)] <- "NA"
orders[orders == "NA"] <- "Unassigned"
names(otus.df) <- orders
df = as.data.frame(t(otus.df))
df$Order = row.names(df)
order.df = df %>% group_by(Order) %>% summarise_all(sum) %>% as.data.frame()
row.names(order.df) <- order.df$Order
order.df <- order.df[,-1]
order.df = order.df[which(row.names(order.df) != "Unassigned"),]
# reorder Cols by structure
meta <- as(sample_data(enhalus),"data.frame")
meta = meta %>% arrange(Source)
order.df = order.df[,meta$IlluminaName]
order.df = decostand(order.df,"total", MARGIN = 2)
order.matrix = as.matrix(order.df)

structures = meta$Source
str.cols = as.character(plyr::mapvalues(structures, from = levels(structures), to=c("#000000","#492c24","#3791b2","#397c1c")))


png("./Output/Heatmap_order.png",width = 16,height = 16,units = 'in',res = 300)
heatmap_left(t(order.matrix), col = hm.pal, Rowv = NA, Colv=NA, labRow = NA, RowSideColors = str.cols,
             margins = c(50,10), cexCol = 5)
dev.off()


# Venn_Diagram ####

pspa = transform_sample_counts(enhalus, function(abund) 1*(abund>0))  
lf.pa = subset_samples(pspa, Source == "Leaf")
rz.pa = subset_samples(pspa, Source == "Rhizome")
rt.pa = subset_samples(pspa, Source == "Root")
sl.pa = subset_samples(pspa, Source == "Soil")

# subset
area.lf = length(which(taxa_sums(lf.pa) > 0))
area.rz = length(which(taxa_sums(rz.pa) > 0))
area.rt = length(which(taxa_sums(rt.pa) > 0))
area.sl = length(which(taxa_sums(sl.pa) > 0))

# find areas
lf.pa = subset_taxa(lf.pa,taxa_sums(lf.pa) > 0)
rz.pa = subset_taxa(rz.pa,taxa_sums(rz.pa) > 0)
rt.pa = subset_taxa(rt.pa,taxa_sums(rt.pa) > 0)
sl.pa = subset_taxa(sl.pa,taxa_sums(sl.pa) > 0)

# find shared sets 
nlf_rz = sum(taxa_names(lf.pa) %in% taxa_names(rz.pa)) #n12
nlf_rt = sum(taxa_names(lf.pa) %in% taxa_names(rt.pa)) #n13
nlf_sl = sum(taxa_names(lf.pa) %in% taxa_names(sl.pa)) #n14
nrz_rt = sum(taxa_names(rz.pa) %in% taxa_names(rt.pa)) #n23
nrz_sl = sum(taxa_names(rz.pa) %in% taxa_names(sl.pa)) #n24
nrt_sl = sum(taxa_names(rt.pa) %in% taxa_names(sl.pa)) #n34

nlf_rz_rt = sum(taxa_names(lf.pa) %in% taxa_names(rz.pa) & taxa_names(lf.pa) %in% taxa_names(rt.pa))#n123
nlf_rz_sl = sum(taxa_names(lf.pa) %in% taxa_names(rz.pa) & taxa_names(lf.pa) %in% taxa_names(sl.pa))#n124
nlf_rt_sl = sum(taxa_names(lf.pa) %in% taxa_names(rt.pa) & taxa_names(lf.pa) %in% taxa_names(sl.pa))#n134
nrz_rt_sl = sum(taxa_names(rz.pa) %in% taxa_names(rt.pa) & taxa_names(rz.pa) %in% taxa_names(sl.pa))#n234
nlf_rz_rt_sl = sum(taxa_names(lf.pa) %in% taxa_names(rz.pa) & taxa_names(lf.pa) %in% taxa_names(rt.pa) & taxa_names(lf.pa) %in% taxa_names(sl.pa))#n1234  


dev.off()
png("./Output/VennDiagram_Source.png", res = 300, width = 2400,height = 1600)
draw.quad.venn(area.lf,area.rz,area.rt,area.sl,nlf_rz,nlf_rt,nlf_sl,nrz_rt,nrz_sl,nrt_sl,nlf_rz_rt,nlf_rz_sl,nlf_rt_sl,nrz_rt_sl,nlf_rz_rt_sl,
                 fill=c("#397c1c","#3791b2","#492c24","#000000"),
                 category = c("Leaf","Rhizome","Root","Soil"),
                 alpha = .45, cex = rep(2,15), cat.cex = rep(2,4), cat.dist = .22)
dev.off()

