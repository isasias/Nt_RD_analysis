#### Beta Diversity ITS ####


### Libraries ###

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("vegan")

### Data ###

load("ITS_filt.RData")

#### Distance calculation: Bray ####

bray_bdiv <- phyloseq::distance(ITS_filtered, 
                                method = "bray", 
                                type = "sample") 
### Ordination method ###
bray_ord <- ordinate(ITS_filtered,"NMDS", distance = bray_bdiv)

### PCoA plot ###

## Order treatments and plant type
ITS_filtered@sam_data$Treatment <- factor(ITS_filtered@sam_data$Treatment, 
                                          levels = c("Pi","Low P", "Phi", "Pi/Phi mix"), 
                                          ordered = TRUE)

ITS_filtered@sam_data$Plant_Type <- factor(ITS_filtered@sam_data$Plant_Type,
                                           levels= c("ptxDNt-36","Wild Type", "Bulk Soil"),
                                           ordered = TRUE)

## Graph

p_bray <- plot_ordination(ITS_filtered, bray_ord, "samples", 
                          color = "Plant_Type",
                          shape = "Treatment")+
  scale_color_manual(values = c("steelblue","lightpink","burlywood3"), name = "Plant Type")+
  theme_bw()+
  geom_point(size=3)

shape_names <- c("Pi","Low P","Phi", "Pi/Phi mix")
p.shape <- 15:(15 + length(shape_names) - 1)
names(p.shape) <- shape_names
p.shape["Treatments"] <- 4

p_bray <- p_bray +
  scale_shape_manual(values = p.shape)

p_bray <- p_bray +
  geom_polygon(aes(fill=Plant_Type),alpha = 0.2)+
  scale_fill_manual(values = c("steelblue","lightpink","burlywood3"), 
                    name = "Plant Type")

p_bray

#### Distance calculation: Jaccard ####

## Presence absence only 
jacc_bdiv <- phyloseq::distance(ITS_filtered, 
                                method = "jaccard", binary = TRUE, 
                                type = "sample") 
### Ordination method ###
jacc_ord <- ordinate(ITS_filtered,"PCoA", distance = jacc_bdiv)

### PCoA plot ###

## Graph

p_jacc <- plot_ordination(ITS_filtered, jacc_ord, "samples", 
                          color = "Plant_Type",
                          shape = "Treatment")+
  scale_color_manual(values = c("steelblue","lightpink","burlywood3"), name = "Plant Type")+
  theme_bw()+
  geom_point(size=3)

shape_names <- c("Pi","Low P","Phi", "Pi/Phi mix")
p.shape <- 15:(15 + length(shape_names) - 1)
names(p.shape) <- shape_names
p.shape["Treatments"] <- 4

p_jacc <- p_jacc +
  scale_shape_manual(values = p.shape)

p_jacc <- p_jacc +
  geom_polygon(aes(fill=Plant_Type),alpha = 0.2)+
  scale_fill_manual(values = c("steelblue","lightpink","burlywood3"), 
                    name = "Plant Type")

p_jacc

#### PERMANOVA analysis ####

### Extract data frames from phyloseq object ###

# count table: ASVs
SoilASVs <- as.data.frame(ITS_filtered@otu_table)

# metadata or Environment
SoilMeta <- read.csv("~/Grad_School/Maestria/Processed_data/metadata.csv")


### Basic Adonis

adonis2(SoilASVs ~ Treatment * Plant_Type, 
        data=SoilMeta, permutations = 999,
        method = "bray")
# no significant differences

adonis2(SoilASVs ~ Treatment * Plant_Type, 
        data=SoilMeta, permutations = 999,
        method = "cao")
# plant type marginal difference p = 0.096

adonis2(SoilASVs ~ Treatment * Plant_Type, 
        data=SoilMeta, permutations = 999,
        method = "jaccard", binary = TRUE) # binary means presence/ absence
# no significant difference

### Stratified data ###

## By plant type

adonis2(SoilASVs ~ Treatment, 
                       strata = SoilMeta$Plant_Type, 
                       data=SoilMeta, permutations = 999, 
                       method = "bray")
# no differences 

## By treatment

adonis2(SoilASVs ~ Plant_Type, 
                       strata = SoilMeta$Treatment, 
                       data=SoilMeta, permutations = 999,
                       method = "bray")
# no differences

### Plotting NMDS with scores:

# not exactly sure the differences between plotting methods

SoilNMDS <- metaMDS(SoilASVs)
Soilscores <- scores(SoilNMDS)
data.scores <- as.data.frame(Soilscores[["sites"]])

ggplot(data=data.scores) + 
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=SoilMeta$Plant_Type),level = 0.50) +
  geom_point(aes(x=NMDS1,y=NMDS2,shape=SoilMeta$Treatment,colour=SoilMeta$Plant_Type),size=4) + 
  theme_minimal() 


