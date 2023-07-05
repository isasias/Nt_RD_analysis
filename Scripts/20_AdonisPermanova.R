#### Adonis for Beta div ####

### Libraries

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("vegan")
library("dplyr")
library("tidyverse")
library("grid")
library("gridExtra") # not used

### Data
load("Psoil_filt.Rdata")
load("Psoil_rel.RData")

### Extract data frames from phyloseq object 

# count table: ASVs
SoilASVs <- as.data.frame(Psoil_filt@otu_table)

# metadata or Environment
SoilMeta <- read.csv("~/Grad_School/Maestria/Processed_data/metadata.csv")
rownames(SoilMeta) <- SoilMeta$Sample_ID

### Basic Adonis

Soil_adonis <- adonis2(SoilASVs ~ Treatment * Plant_Type, 
                      data=SoilMeta, permutations = 999,
                      method = "bray")

# No significant differences 

### Stratified data

# By plant type

Soil_adonis <- adonis2(SoilASVs ~ Treatment, 
                       strata = SoilMeta$Plant_Type, 
                       data=SoilMeta, permutations = 999, 
                       method = "bray")

# no differences 

# by treatment

Soil_adonis <- adonis2(SoilASVs ~ Plant_Type, 
                       strata = SoilMeta$Treatment, 
                       data=SoilMeta, permutations = 999,
                       method = "bray")
# no differences

### Plotting NMDS

SoilNMDS <- metaMDS(SoilASVs)
Soilscores <- scores(SoilNMDS)
data.scores <- as.data.frame(Soilscores[["sites"]])

ggplot(data=data.scores) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=SoilMeta$Plant_Type),level = 0.50) +
  geom_point(aes(x=NMDS1,y=NMDS2,shape=SoilMeta$Plant_Type,colour=SoilMeta$Plant_Type),size=4) + 
  theme_minimal() # using sites do not know why

# it seems there were no significant differences no matter the iterations

