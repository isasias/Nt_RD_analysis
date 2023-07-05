#### Diversity Analysis ####

#### Libraries ####
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("ape") # for phylogenetic trees not used yet
library("plyr"); packageVersion("plyr") # might not be necessary
library("cluster"); packageVersion("cluster") # not necessary it was used to calculate gap statistics which I ended up not using 
library("patchwork") # two visualize graphs together meant to be used with ggplot
library("betapart") # for beta diversity not used
library("dplyr") # tidyverse
library("remotes") # download remote libraries
library("devtools") # to install ranacapa
library("ranacapa") # for rarefaction curve 
library("vegan") # for rarefy function

#### Upload Phyloseq object ####

load("phylo_soil.Rdata")

#### Preprocessing ####

### Filter taxa with low abundance ###

# ASVs that are found in at least 3 different samples
PS_filtered <-  filter_taxa(phylo_soil, function(OTU) sum(OTU) > 3, TRUE)

# Remove unidentified taxa
PS_filtered <- subset_taxa(PS_filtered, !is.na(Phylum))

### Transform to relative abundance ###

PS_relative  <-  transform_sample_counts(PS_filtered, function(x) x / sum(x) )
# This divides the count of one species by the total number of species in that sample

### Agglomeration function: join ASVs with sames genus into one read ###

# PS_glom <- tax_glom(PS_filtered,taxrank = "Genus", NArm = FALSE) # 691 taxa
# taxa agglomeration will not be conducted until the analysis with the 
# original data is conducted for comparison


#### Alpha-diversity: Indices calculations ####

### PS_filtered ###

data_richness <- estimateR(PS_filtered@otu_table)  # calculate richness and Chao1 using vegan package

## Shannon
data_shannon <- diversity(PS_filtered@otu_table, index = "shannon")
shannon_evenness <- diversity(PS_filtered@otu_table, index = "shannon") / log(specnumber(phylo_soil@otu_table))

## Simpson
data_simpson <- diversity(PS_filtered@otu_table, index = "simpson")
simpson_evenness <- diversity(PS_filtered@otu_table, index = "simpson") / log(specnumber(phylo_soil@otu_table))

## Create table with indices
Alphadiv_filt <- cbind(phylo_soil@sam_data, t(data_richness), data_shannon, shannon_evenness, data_simpson, simpson_evenness )

### PS_relative ###

# data_richness <- estimateR(PS_relative@otu_table)  
# richness cannot be calculated with relative values
# same values than absolute abundances for diversity and evenness

### PS_glom ###

data_richness <- estimateR(PS_glom@otu_table)  # calculate richness and Chao1 using vegan package

## Shannon
data_shannon <- diversity(PS_glom@otu_table, index = "shannon")
shannon_evenness <- diversity(PS_glom@otu_table, index = "shannon") / log(specnumber(phylo_soil@otu_table))

## Simpson
data_simpson <- diversity(PS_glom@otu_table, index = "simpson")
simpson_evenness <- diversity(PS_glom@otu_table, index = "simpson") / log(specnumber(phylo_soil@otu_table))

## Create table with indices
Alphadiv_glom <- cbind(PS_glom@sam_data, t(data_richness), data_shannon, shannon_evenness, data_simpson, simpson_evenness )
# this is the one that I will be using for downstream analysis

save(PS_glom, file = "Phylo_glom.RData")
write.csv(Alphadiv_glom, "~/Grad_School/Maestria/Processed_data/Alphadiversity.csv")



#### Check for normality ####

### Filtered ###

## Histograms

# Observed and Chao (normal)
par(mfrow=c(1,2))
hist(Alphadiv_filt$S.obs) 
hist(Alphadiv_filt$S.chao1)

# Shannon and Simpson Diversity
par(mfrow=c(1,2))
hist(Alphadiv_filt$data_shannon) 
hist(Alphadiv_filt$data_simpson)

# Evenness
par(mfrow=c(1,2))
hist(Alphadiv_filt$shannon_evenness) 
hist(Alphadiv_filt$simpson_evenness) # skewed to the right 

## Shapiro tests
shapiro.test(Alphadiv_filt$S.obs) # yes
shapiro.test(Alphadiv_filt$S.chao1) # yes
shapiro.test(Alphadiv_filt$data_shannon) # yes
shapiro.test(Alphadiv_filt$data_simpson) # yes
shapiro.test(Alphadiv_filt$shannon_evenness) # yes
shapiro.test(Alphadiv_filt$simpson_evenness) # yes

### Agglomerated ###

## Histograms

# Observed and Chao (normal)
par(mfrow=c(1,2))
hist(Alphadiv_glom$S.obs)
hist(Alphadiv_glom$S.chao1)  # skewed to the left?

# Shannon and Simpson Diversity
par(mfrow=c(1,2))
hist(Alphadiv_glom$data_shannon) 
hist(Alphadiv_glom$data_simpson)

# Evenness
par(mfrow=c(1,2))
hist(Alphadiv_glom$shannon_evenness) 
hist(Alphadiv_glom$simpson_evenness) 

## Shapiro tests
shapiro.test(Alphadiv_glom$S.obs) # yes
shapiro.test(Alphadiv_glom$S.chao1) # yes, barely
shapiro.test(Alphadiv_glom$data_shannon) # yes
shapiro.test(Alphadiv_glom$data_simpson) # yes
shapiro.test(Alphadiv_glom$shannon_evenness) # yes
shapiro.test(Alphadiv_glom$simpson_evenness) # yes



#### Plot Alpha Diversity ####

### Set theme ###
theme_set(theme_linedraw())

### Bar plots ###

# Richness
ggplot(Alphadiv_glom, aes(x=Treatment, y=S.obs, fill = Plant_Type)) +
  geom_col() +
  labs(title= 'Richness', x= ' ', y= '') +
  scale_fill_manual(values = c("lightpink","palegreen3","steelblue"))

# Chao
ggplot(Alphadiv_glom, aes(x=Treatment, y=S.chao1, fill = Plant_Type)) +
  geom_col() +
  labs(title= 'Chao Index', x= ' ', y= '') +
  scale_fill_manual(values = c("lightpink","palegreen3","steelblue"))

# Shannon
P1 <- ggplot(Alphadiv_glom, aes(x=Treatment, y=data_shannon, fill = Plant_Type)) +
  geom_col() +
  labs(title= 'Diversity', x= ' ', y= '', tag = "A") +
  scale_fill_manual(values = c("lightpink","palegreen3","steelblue"))
P2 <- ggplot(Alphadiv_glom, aes(x=Treatment, y=shannon_evenness, fill = Plant_Type)) +
  geom_col() +
  labs(title= 'Evenness', x= ' ', y= '', tag = "B") +
  scale_fill_manual(values = c("lightpink","palegreen3","steelblue"))
(P1 | P2)

# Simpson
P3 <- ggplot(Alphadiv_glom, aes(x=Treatment, y=data_simpson, fill = Plant_Type)) +
  geom_col() +
  labs(title= 'Diversity', x= ' ', y= '', tag = "A") +
  scale_fill_manual(values = c("lightpink","palegreen3","steelblue"))
P4 <- ggplot(Alphadiv_glom, aes(x=Treatment, y=simpson_evenness, fill = Plant_Type)) +
  geom_col() +
  labs(title= 'Evenness', x= ' ', y= '', tag = "B") +
  scale_fill_manual(values = c("lightpink","palegreen3","steelblue"))
(P3 | P4)

# I do not like how the column plot look let's use plot richness function

### Final Plot for Alpha Diversity ###

PR <- plot_richness(PS_glom, x="Treatment", 
                    measures=c("Observed","Shannon","Simpson"), 
                    color="Plant_Type")+ 
  scale_color_manual(values = c("lightpink","palegreen3","steelblue"))
PR$layers <- PR$layers[-1]
PR + geom_point(size=4, alpha=0.7)

# plot_richness(psd5, color = "species", x = "species", 
#              measures = c("Observed", "Chao1", "Shannon")) + 
#  geom_boxplot(aes(fill = species), alpha=.7) + 
#  scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) + 
#  scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f"))

#### Beta diversity ####

### Jaccard index: Presence/ absence ###

jacc_bdiv <- distance(PS_glom, "jaccard",binary = TRUE) # same as using betadiver
jacc_ord <- ordinate(PS_glom,"NMDS", distance = jacc_bdiv)

## NMDS plots using jaccard

# Plant better with this index 
p_jacc <- plot_ordination(PS_glom, jacc_ord, "samples", 
                          color = "Plant_Type", shape = "Treatment")
p_jacc + geom_polygon(aes(fill=Plant_Type), alpha = 0.3) + 
  geom_point(size=2)  
# try biplot or split to add taxa into the pcoa

# Treatment
p_jacct <- plot_ordination(PS_glom, jacc_ord, "samples", 
                           color = "Treatment", shape = "Plant_Type")
p_jacct + geom_polygon(aes(fill=Treatment), alpha = 0.3) + 
  geom_point(size=2)  

### Bray-Curtis: dissimilarity method ###

# This index takes into consideration abundance  

bray_bdiv <- distance(PS_glom, "bray") 
bray_ord <- ordinate(PS_glom,"MDS", distance = bray_bdiv)

## NMDS plots using jaccard

# Plant Type
p_bray <- plot_ordination(PS_glom, bray_ord, "samples", 
                          color = "Plant_Type", shape = "Treatment")
p_bray + geom_polygon(aes(fill=Plant_Type), alpha = 0.3) + 
  geom_point(size=2)  

# Treatment better with this index
p_bt <- plot_ordination(PS_glom, bray_ord, "samples", 
                           color = "Treatment", shape = "Plant_Type")
p_bt + geom_polygon(aes(fill=Treatment), alpha = 0.3) + 
  geom_point(size=2)

### Other beta diversity indices: no graphs ###

sor_bdiv <- distance(PS_glom, "sor", binary = TRUE)
w_bdiv <- distance(PS_glom, "w")
jsd_bdiv <- distance(PS_glom, "jsd")

## Create distance data frame for ANOVA analysis

trt_groups <- as.factor(PS_glom@sam_data$Treatment)
pl_groups <- as.factor(PS_glom@sam_data$Plant_Type)

# Treatment groups
bdj <-betadisper(jacc_bdiv,trt_groups)
bdb <- betadisper(bray_bdiv,trt_groups)
bdw <- betadisper(w_bdiv,trt_groups)
bds <- betadisper(sor_bdiv,trt_groups)
bdd <- betadisper(jsd_bdiv,trt_groups)

betadiv_trt <- cbind(PS_glom@sam_data, bdj$distances, bdb$distances,
                     bdw$distances, bds$distances, bdd$distances)

# Plant groups
bdj <-betadisper(jacc_bdiv,pl_groups)
bdb <- betadisper(bray_bdiv,pl_groups)
bdw <- betadisper(w_bdiv,pl_groups)
bds <- betadisper(sor_bdiv,pl_groups)
bdd <- betadisper(jsd_bdiv,pl_groups)

betadiv_pl <- cbind(PS_glom@sam_data, bdj$distances, bdb$distances,
                     bdw$distances, bds$distances, bdd$distances)

rm(bdb,bdd,bdj,bds,bdw,presabs)

write.csv(betadiv_pl, "~/Grad_School/Maestria/Processed_data/b_diver_pl.csv")
write.csv(betadiv_trt, "~/Grad_School/Maestria/Processed_data/b_diver_trt.csv")

# Check distance to centoid using boxplots

## ANOVAS: maybe create tables with distances to better depict this 

# Missing this section just save the distances 
# Check permanova too

#### Rarefaction ####

# Rarefaction is a technique to assess expected species richness.
# Rarefaction allows the calculation of species richness for a given 
# number of individual samples, based on the construction of rarefaction 
# curves.

## Extract otu table for analysis

tab <- otu_table(PS_glom)
class(tab) <- "matrix"

spAbund <- rowSums(tab)
raremin <- min(rowSums(tab))

# Rarefy function

sRare <- rarefy(tab, raremin) #gives an "expected"rarefied" number
# of species (not obs) if only 15 individuals were present
# heat map

# Rarefaction curve 

rarecurve(tab, col = "black", label = FALSE, tidy = TRUE) # long running code
grare <- ggrare(PS_glom, step = 500, label = NULL, color = "Plant_Type")

#### Heatmaps ####

plot_heatmap(PS_glom, sample.label= "Treatment", 
             sample.order = "Treatment", distance = "bray")

taxtable <- as.data.frame(PS_glom@tax_table)

taxtable %>%
  count(Phylum)

### Subset most abundant phyla ###

## Actinobacteria

soil_actino <- subset_taxa(PS_glom, Phylum=="Actinobacteria")

soil_actino <- prune_taxa(names(sort(taxa_sums(soil_actino),TRUE)[1:20]), soil_actino)
soil_actino2 <- prune_taxa(names(sort(taxa_sums(soil_actino),TRUE)[50:100]), soil_actino)
soil_actino <- merge_phyloseq(soil_actino,soil_actino2)
soil_actino <- tax_glom(soil_actino,taxrank = "Family", NArm = FALSE)

plot_heatmap(soil_actino, sample.label = "Sample_ID", method = "NMDS",
             distance = "bray", taxa.label = "Family",
             taxa.order = "Family", sample.order = "Sample_ID",
             low="#66CCFF", high="#000033", na.value="white")

## Proteobacteria

soil_proteo <- subset_taxa(PS_glom, Phylum=="Proteobacteria")
soil_proteo <- subset_taxa(soil_proteo,!Family=="Unknown_Family")

soil_proteo <- prune_taxa(names(sort(taxa_sums(soil_proteo),TRUE)[1:50]), soil_proteo)
# soil_proteo <- subset_taxa(soil_proteo,Family== "Beijerinckiaceae")

soil_proteo <- tax_glom(soil_proteo,taxrank = "Family", NArm = TRUE)

plot_heatmap(soil_proteo, sample.label = "Sample_ID", method = "NMDS",
             distance = "bray", taxa.label = "Family",
             taxa.order = "Family", sample.order = "Sample_ID",
             low="#66CCFF", high="#000033", na.value="white")


#### Gap Statistics ####

# The purpose of this section is to know in how many clusters do my data
# separate into, that may help showing differences between treatments

### Compute Gap Statistics

pam1 <-  function(x, k) 
{list(cluster = pam(x,k, cluster.only=TRUE))}
x <-  phyloseq:::scores.pcoa(jacc_ord, display="species", physeq = PS_glom) # ordination made with MDS

gs_bray <-  clusGap(x[, 1:2], FUN=pam1, K.max = 8, B = 50) # K.max number of clusters

## Wrapper function

gap_statistic_ordination = function(ord, FUNcluster, type="species", physeq= PS_glom, K.max=8, axes=c(1:2), B=500, verbose=interactive(), ...){
  require("cluster")
  #   If "pam1" was chosen, use this internally defined call to pam
  if(FUNcluster == "pam1"){
    FUNcluster = function(x,k) list(cluster = pam(x, k, cluster.only=TRUE))     
  }
  # Use the scores function to get the ordination coordinates
  x = phyloseq:::scores.pcoa(ord, display=type, physeq = physeq)
  #   If axes not explicitly defined (NULL), then use all of them
  if(is.null(axes)){axes = 1:ncol(x)}
  #   Finally, perform, and return, the gap statistic calculation using cluster::clusGap  
  clusGap(x[, axes], FUN=FUNcluster, K.max=K.max, B=B, verbose=verbose, ...)
}

### Plot results

plot_clusgap = function(clusgap, title="Gap Statistic calculation results"){
  require("ggplot2")
  gstab = data.frame(clusgap$Tab, k=1:nrow(clusgap$Tab))
  p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
  p = p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim))
  p = p + ggtitle(title)
  return(p)
}

gs = gap_statistic_ordination(bray_ord, "pam1", B=50, verbose=FALSE)
print(gs, method="Tibs2001SEmax")

plot_clusgap(gs) ## not sufficient gaps



#### New graphs ####

## Subset phyla and agglomerate by class

# Actinobacteria
Act_soil <- subset_taxa(Psoil_rel, Phylum == "Actinobacteria")

Act_soil<- tax_glom(Act_soil,taxrank = "Class", NArm = FALSE)

Act_plot <-plot_composition(Act_soil, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_brewer(palette = "Paired",
                    name = "Class")+
  guides(x=guide_axis(angle = 35))+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

# Separate Rubro and actinobacteria classes

# Proteobacteria
Prot_soil <- subset_taxa(Psoil_rel, Phylum == "Proteobacteria")

Prot_soil<- tax_glom(Prot_soil,taxrank = "Class", NArm = FALSE)

Prot_plot <-plot_composition(Prot_soil, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_brewer(palette = "Paired",
                    name = "Class")+
  guides(x=guide_axis(angle = 35))+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

# Separate Gammaproteobacteria, and Alphaproteobacteria

# Bacteroidia Class
Bact_soil <- subset_taxa(Psoil_rel, Class == "Bacteroidia")

Bact_soil<- tax_glom(Bact_soil,taxrank = "Order", NArm = FALSE)

Bact_plot <-plot_composition(Bact_soil, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_brewer(palette = "Paired",
                    name = "Order")+
  guides(x=guide_axis(angle = 35))+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

# From this graph it might be useful orders Cytophagales and Chitinophagales

# Chloroflexi
Chloro_soil <- subset_taxa(Psoil_rel, Phylum == "Chloroflexi")

Chloro_soil<- tax_glom(Chloro_soil,taxrank = "Class", NArm = FALSE)

Chloro_plot <-plot_composition(Chloro_soil, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_brewer(palette = "Paired",
                    name = "Class")+
  guides(x=guide_axis(angle = 35))+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

# Low abundance taxa seems to have more differences

# Acidobacteria

Acido_soil <- subset_taxa(Psoil_rel, Phylum == "Acidobacteria")

Acido_soil<- tax_glom(Acido_soil,taxrank = "Class", NArm = FALSE)

Acido_plot <-plot_composition(Acido_soil, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_brewer(palette = "Paired",
                    name = "Class")+
  guides(x=guide_axis(angle = 35))+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

# All may be useful
```

#### Heatmaps 

```{r}
# Choose first sample and taxa

which(Psoil_filt@otu_table == max(Psoil_filt@otu_table), arr.ind = TRUE) # Sample D3, otu ZOVUZ0834K

sorder <- c("A1","B1","C1", "A2","B2","C2", "A3","B3","C3",
            "A4","B4","C4", "A5","B5","C5", "A6","B6","C6",
            "A7","B7","C7", "D1","D2","D3","D4")

# Original
plot_heatmap(Psoil_filt, sample.label = "Sample_ID",
             method = "NMDS", distance = "bray",
             taxa.order = names(sort(taxa_sums(Psoil_filt))),
             low="blue", high="red",sample.order = sorder)

# Upper subset
Upsoil <- prune_taxa(names(sort(taxa_sums(Psoil_filt),TRUE)[450:1000]), Psoil_filt)

plot_heatmap(Upsoil, sample.label = "Sample_ID",
             method = "NMDS", distance = "bray",
             taxa.order = names(sort(taxa_sums(Upsoil))),
             low="blue", high="red",sample.order = sorder)

# Lower subset
Lowsoil <- prune_taxa(names(sort(taxa_sums(Psoil_filt),TRUE)[5000:10973]), Psoil_filt)

plot_heatmap(Lowsoil, sample.label = "Sample_ID",
             method = "NMDS", distance = "bray",
             taxa.order = names(sort(taxa_sums(Lowsoil))),
             low="blue", high="red",sample.order = sorder)

### Agglomerated heatmaps by taxonomic hierarchy

# Phyla
Phyla_soil <- tax_glom(Psoil_filt,taxrank = "Phylum", NArm = FALSE)

plot_heatmap(Phyla_soil, sample.label = "Sample_ID",
             method = "NMDS", distance = "bray", 
             taxa.label = "Phylum",
             taxa.order = names(sort(taxa_sums(Phyla_soil))),
             low="blue", high="red",sample.order = sorder)

# From this initial graph I selected 4 phyla that show some differences betweeen treatments and not only on individual samples: Nanoarchaeaeota, GAL15, Euryarchaeota, Thermotogae

# Class

Class_soil <- subset_taxa(Psoil_filt, Phylum == c("Nanoarchaeaeota", "GAL15", "Euryarchaeota", "Thermotogae"))

plot_heatmap(Class_soil, sample.label = "Sample_ID",
             method = "NMDS", distance = "bray", 
             taxa.label = "Phylum",
             taxa.order = names(sort(taxa_sums(Class_soil))),
             low="blue", high="red",sample.order = sorder)

## Check plot by taxa

taxa_bdiv <- phyloseq::distance(Psoil_filt, method = "bray", type = "taxa") # check how it looks with taxa type
taxa_ord <- ordinate(Psoil_filt,"NMDS", distance = bray_bdiv)

p_bray <- plot_ordination(Psoil_filt, bray_ord, "samples", 
                          color = "Plant_Type",
                          shape = "Treatment")+
  scale_color_manual(values = c("steelblue","palegreen3","lightpink"), name = "Plant Type")+
  theme_bw()+
  geom_point(size=3)

shape_names <- c("Control","Low P","Phi", "Pi/Phi mix")
p.shape <- 15:(15 + length(shape_names) - 1)
names(p.shape) <- shape_names
p.shape["Treatments"] <- 4

p_bray <- p_bray +
  scale_shape_manual(values = p.shape)


p_bray <- p_bray +
  geom_polygon(aes(fill=Plant_Type),alpha = 0.2)+
  scale_fill_manual(values = c("steelblue","palegreen3","lightpink"), name = "Plant Type")

