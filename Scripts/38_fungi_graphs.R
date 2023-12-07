#### Fungi Taxonomical Composition Analysis ####

### Libraries ###

library(phyloseq)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(microbiome)
library(ggsignif) # has DOI
library(scales)

### Data ###

load("ITS_filt.RData")
Func_Fungi <- read.csv("~/Grad_School/Maestria/Processed_data/Fungal traits/Functional_Fungi.csv")
Simp_fungi <- read.csv("~/Grad_School/Maestria/Processed_data/Fungal traits/Simplified_Fungi.csv")

#### Phyla abundance ####

Phyla_fun <- tax_glom(ITS_filtered,taxrank = "Phylum", NArm = FALSE)


## Extract data from phyloseq object to use pheatmap

OTU_matrix <- as.data.frame(Phyla_fun@otu_table)
Tax_matrix <- as.data.frame(Phyla_fun@tax_table)

## Rename columns from otu_matrix with phylum from tax_matrix

colnames(OTU_matrix) <- Tax_matrix$Phylum
OTU_matrix <- OTU_matrix[,-4] # remove unknown phyla

Phyla_matrix <- as.matrix(t(OTU_matrix))

# arrange rows and columns

Phyla_matrix <- Phyla_matrix[order(Phyla_matrix[,1], 
                                   decreasing = TRUE),]

sorder <- c("A1","B1","C1","A2","B2","C2", 
            "A3","B3","C3", "A4","B4","C4", 
            "A5","B5","C5", 
            "A6","B6","C6", "A7","B7","C7", 
            "D1","D2","D3","D4")

Phyla_matrix <- Phyla_matrix[ , sorder]

phyla_rows <- c("Ascomycota","Chytridiomycota","Mucoromycota",
                "Basidiomycota","Mortierellomycota","Blastocladiomycota",
                "Glomeromycota","Rozellomycota","Olpidiomycota",
                "Kickxellomycota","Aphelidiomycota","Zoopagomycota",
                "Calcarisporiellomycota","Neocallimastigomycota","Basidiobolomycota")

row.names(Phyla_matrix) <- phyla_rows

#### Heatmap using pheatmap ####

## add breaks
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- seq(min(Phyla_matrix), max(Phyla_matrix), length.out = 30)
mat_breaks <- quantile_breaks(Phyla_matrix, n = 31)

pheatmap(Phyla_matrix,cluster_rows=FALSE, cluster_cols=FALSE, scale = "none",
         color = colorRampPalette(c("snow2","lightskyblue1","plum3","magenta4","steelblue4","black"))(25),
         breaks            = mat_breaks,
         drop_levels       = F,
         fontsize          = 10,
         gaps_col = c(3,6,9,12,15,18,21),
         labels_col = c (rep("36LP",3),rep("WTLP",3),  
                         rep("36Pi",3), rep("WTPi",3),  
                         rep("36Phi",3), 
                         rep("36PM",3), rep("WTPM",3),
                         "BSLP","BSPi","BSPhi","BSPM"))


#### Functional Heatmap ####

row.names(Func_Fungi) <- Func_Fungi[,1]
Func_Fungi <- Func_Fungi[,-1]

# change always to matrix
Func_Fungi <- data.matrix(Func_Fungi)

Func_Fungi <- Func_Fungi[order(Func_Fungi[,1], 
                               decreasing = TRUE),]

sorder <- c("A1","B1","C1","A2","B2","C2", 
            "A3","B3","C3", "A4","B4","C4", 
            "A5","B5","C5", 
            "A6","B6","C6", "A7","B7","C7", 
            "D1","D2","D3","D4")

Func_Fungi <- Func_Fungi[ , sorder]

## Quantile breaks ##

# function
quantile_breaks <- function(xs, n = 10 ) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# breaks object 20 breaks

mat_breaks <- seq(min(Func_Fungi), max(Func_Fungi), length.out = 10)
mat_breaks <- quantile_breaks(Func_Fungi, n = 21)

# heatmap

pheatmap(Func_Fungi,cluster_rows=FALSE, cluster_cols=FALSE, scale = "none",
         color = colorRampPalette(c("snow2", "lightskyblue1", "plum3","magenta4","steelblue4", "black"))(18),
         breaks            = mat_breaks,
         drop_levels       = F,
         fontsize          = 5,
         gaps_col = c(3,6,9,12,15,18,21),
         labels_col = c (labels=c(rep("36LP",3),rep("WTLP",3),  
                                  rep("36Pi",3), rep("WTPi",3),  
                                  rep("36Phi",3), 
                                  rep("36PM",3), rep("WTPM",3),
                                  "BSLP","BSPi","BSPhi","BSPM")),
         labels_row = c(1:97))

# check which ones show differences

### Simplified functions ###

### Energy source 

row.names(Simp_fungi) <- Simp_fungi[,1]
Simp_fungi <- Simp_fungi[,-1]

# change always to matrix
Simp_fungi <- data.matrix(Simp_fungi)

Simp_fungi <- Simp_fungi[order(Simp_fungi[,1], 
                         decreasing = TRUE),]

sorder <- c("A1","B1","C1","A2","B2","C2", 
            "A3","B3","C3", "A4","B4","C4", 
            "A5","B5","C5", 
            "A6","B6","C6", "A7","B7","C7", 
            "D1","D2","D3","D4")

Simp_fungi <- Simp_fungi[ , sorder]

##Quantile breaks

# function
quantile_breaks <- function(xs, n = 10 ) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# breaks object 10 breaks
mat_breaks <- seq(min(Simp_fungi), max(Simp_fungi), length.out = 13)
mat_breaks <- quantile_breaks(Simp_fungi, n = 10)

# heatmap
pheatmap(Simp_fungi,cluster_rows=FALSE, cluster_cols=FALSE, scale = "none",
         color = colorRampPalette(c("snow2", "lightskyblue1", "plum3","magenta4","steelblue4","black"))(9),
         breaks            = mat_breaks,
         drop_levels       = F,
         fontsize          = 12,
         gaps_col = c(3,6,9,12,15,18,21),
         labels_col = c (labels=c(rep("36LP",3),rep("WTLP",3),  
                                  rep("36Pi",3), rep("WTPi",3),  
                                  rep("36Phi",3), 
                                  rep("36PM",3), rep("WTPM",3),
                                  "BSLP","BSPi","BSPhi","BSPM")))


#### Taxonomic Composition ####

### Preview graphs

## Ascomycota
Asco <- subset_taxa(ITS_filtered, Class == "c__Sordariomycetes")
Asco_fun <- tax_glom(Asco,taxrank = "Genus", NArm = FALSE)
Asco_fun <- prune_taxa(names(sort(taxa_sums(Asco_fun),TRUE)[1:10]), Asco_fun)

plot_composition(Asco_fun)

# Fusarium equiseti most abundant species
# Fusarium most abundant genus

### Mucoromycota

Muco <- subset_taxa(ITS_filtered, Phylum == "p__Mucoromycota")
Muco <- tax_glom(Muco,taxrank = "Genus", NArm = FALSE)
plot_composition(Muco)

# arrhizus and muscae weird most abundant species 
# Rhizopus and Circinella most abundat genus

### Chytridiomycota
Chy <- subset_taxa(ITS_filtered, Phylum == "p__Chytridiomycota")
Chy <- tax_glom(Chy,taxrank = "Genus", NArm = FALSE) # only up to genus
plot_composition(Chy)

# Rhizophlyctis, Rhizophlyctidaceae_gen_Incertae_sedis, 
# Spizellomycetaceae_gen_Incertae_sedis most abundant genera



### Phyla with significant differences ###

## Prepare table

Phyla_fun <- tax_glom(ITS_filtered,taxrank = "Phylum", NArm = FALSE)

# Extract matrix
OTU_matrix <- as.data.frame(Phyla_fun@otu_table)
Tax_matrix <- as.data.frame(Phyla_fun@tax_table)
# Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Phylum
# Extract metadata 
metadata <- as.data.frame(Phyla_fun@sam_data)
# Join phyla with metadata 
Fungi_phyla <- cbind(metadata,OTU_matrix)

## Mucoromycota

ggplot(Fungi_phyla, aes(x = Plant_Type, y = p__Mucoromycota, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 17500, xmin= 0.6, xmax = 1.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 10000, xmin= 1.6, xmax = 2.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 10000, xmin= 2.6, xmax = 3.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  scale_y_continuous(breaks = pretty_breaks(n = 5))
  
## Chytridiomycota

ggplot(Fungi_phyla, aes(x = Plant_Type, y = p__Mucoromycota, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 17500, xmin= 0.6, xmax = 1.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 10000, xmin= 1.6, xmax = 2.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 10000, xmin= 2.6, xmax = 3.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  scale_y_continuous(breaks = pretty_breaks(n = 5))

## Zoopagomycota

ggplot(Fungi_phyla, aes(x = Plant_Type, y = p__Zoopagomycota, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 65, xmin= 0.6, xmax = 1.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 15, xmin= 1.6, xmax = 2.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 65, xmin= 2.6, xmax = 3.4,
              annotation = "AB", tip_length = 0.01,
              col= 1)+
  scale_y_continuous(breaks = pretty_breaks(n = 6))


### Fusarium ###

## Table 
Fusarium <- subset_taxa(ITS_filtered, Genus == "g__Fusarium")
Fusarium <- tax_glom(Fusarium,taxrank = "Species", NArm = FALSE)
# Extract matrix
OTU_matrix <- as.data.frame(Fusarium@otu_table)
Tax_matrix <- as.data.frame(Fusarium@tax_table)
# Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Species
# Extract metadata 
metadata <- as.data.frame(Fusarium@sam_data)
# Join phyla with metadata 
Fusarium <- cbind(metadata,OTU_matrix)

## verticillioides

ggplot(Fusarium, aes(x = Treatment, y = s__verticillioides, col = Plant_Type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Plant Type", 
                     values = c("steelblue","lightpink","burlywood3"))+
  xlab("Treatment") + ylab("Abundance")+
  geom_signif(y_position = 800, xmin= 0.6, xmax = 3.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 1100, xmin= 3.6, xmax = 4.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 2.5, y= 1020,  label="p=0.053",
           color="black", size= 4)
  
## tricinctum

ggplot(Fusarium, aes(x = Plant_Type, y = s__tricinctum, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 1000, xmin= 1.6, xmax = 3.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 250, xmin= 0.6, xmax = 1.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 1, y= 760,  label="p<0.01",
           color="black", size= 4)


### Mucoromycota: species ###
# Not useful 

## Table
Muco <- subset_taxa(ITS_filtered, Phylum == "p__Mucoromycota")
Muco <- tax_glom(Muco,taxrank = "Species", NArm = FALSE)
## Extract matrix
OTU_matrix <- as.data.frame(Muco@otu_table)
Tax_matrix <- as.data.frame(Muco@tax_table)
## Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Species
## Extract metadata 
metadata <- as.data.frame(Muco@sam_data)
## Join phyla with metadata 
Muco <- cbind(metadata,OTU_matrix)


## Microsporus

ggplot(Muco, aes(x = Plant_Type, y = s__microsporus, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 130, xmin= 0.6, xmax = 1.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 60, xmin= 1.6, xmax = 2.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 60, xmin= 2.6, xmax = 3.4,
              annotation = "AB", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 2.5, y= 90,  label="p<0.05",
           color="black", size= 4)

## muscae

ggplot(Muco, aes(x = Plant_Type, y = s__muscae, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")

### Mucoromycota: genus ###

## Table
Muco <- subset_taxa(ITS_filtered, Phylum == "p__Mucoromycota")
Muco <- tax_glom(Muco,taxrank = "Genus", NArm = FALSE)
# Extract matrix
OTU_matrix <- as.data.frame(Muco@otu_table)
Tax_matrix <- as.data.frame(Muco@tax_table)
# Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Genus
# Extract metadata 
metadata <- as.data.frame(Muco@sam_data)
# Join phyla with metadata 
Muco <- cbind(metadata,OTU_matrix)


## Apophysomyces
ggplot(Muco, aes(x = Treatment, y = g__Apophysomyces, col = Plant_Type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Plant Type", 
                     values = c("steelblue","lightpink","burlywood3"))+
  xlab("Treatment") + ylab("Abundance")+
  geom_signif(y_position = 100, xmin= 0.6, xmax = 1.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 37.5, xmin= 1.6, xmax = 2.4,
              annotation = "AB", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 25, xmin= 2.6, xmax = 3.4,
              annotation = "C", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 62.5, xmin= 3.6, xmax = 4.4,
              annotation = "BC", tip_length = 0.01,
              col= 1)

## Absidia
ggplot(Muco, aes(x = Plant_Type, y = g__Absidia, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 300, xmin= 0.6, xmax = 1.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 50, xmin= 1.6, xmax = 2.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 150, xmin= 2.6, xmax = 3.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  scale_y_continuous(breaks = pretty_breaks(n = 6))

## Mucor: maybe 

ggplot(Muco, aes(x = Plant_Type, y = g__Mucor, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")

### Chytridiomycota: genus ###

## Table
Chy <- subset_taxa(ITS_filtered, Phylum == "p__Chytridiomycota")
Chy <- tax_glom(Chy,taxrank = "Genus", NArm = FALSE)
# Extract matrix
OTU_matrix <- as.data.frame(Chy@otu_table)
Tax_matrix <- as.data.frame(Chy@tax_table)
# Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Genus
# Extract metadata 
metadata <- as.data.frame(Chy@sam_data)
# Join phyla with metadata 
Chy <- cbind(metadata,OTU_matrix)

## Rhizophlyctis

ggplot(Chy, aes(x = Plant_Type, y = g__Rhizophlyctis, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 8600, xmin= 0.6, xmax = 1.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 2000, xmin= 1.6, xmax = 3.4,
              annotation = "B", tip_length = 0.01,
              col= 1)

## Spizellomyces

ggplot(Chy, aes(x = Plant_Type, y = g__Spizellomyces, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 100, xmin= 0.6, xmax = 1.4,
              annotation = "A", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 75, xmin= 1.6, xmax = 2.4,
              annotation = "B", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 180, xmin= 2.6, xmax = 3.4,
              annotation = "AB", tip_length = 0.01,
              col= 1)

#### Functional classification graphs ####

### Soil Saprotroph ###

## Rename variables
colnames(Simp_fungi)[6] <- "Soil_sapro"
colnames(Norm_SF)[6] <- "Soil_sapro"
Simp_fungi$Plant_Type[Simp_fungi$Plant_Type=="No plant"] <- "Bulk Soil"


ggplot(Simp_fungi, aes(x = Plant_Type, y = Soil_sapro, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                    values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 13000, xmin= 1.6, xmax = 3.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 23000, xmin= 0.5, xmax = 2.5,
              annotation = "p=0.047", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 25000, xmin= 0.5, xmax = 3.5,
              annotation = "p=0.096", tip_length = 0.01,
              col= 1)

### Dung Saprotroph ###

colnames(Simp_fungi)[9] <- "Dung_sapro"
colnames(Norm_SF)[9] <- "Dung_sapro"

Simp_fungi$Plant_Type <- factor(Simp_fungi$Plant_Type, levels = c("ptxDNt-36","Wild Type", "Bulk Soil"), ordered = TRUE)

ggplot(Simp_fungi, aes(x = Plant_Type, y = Dung_sapro, col = Treatment)) +
  geom_boxplot(fill = "snow2",lwd = 0.7) +
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Treatment", 
                     values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab("Abundance")+
  geom_signif(y_position = 10000, xmin= 1.6, xmax = 3.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 23000, xmin= 0.5, xmax = 2.5,
              tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 23000, xmin= 0.5, xmax = 3.5,
              annotation = "p<0.05", tip_length = 0.01,
              col= 1) # one line added in paint

### Root Associated ###

Root_sums$Plant_Type[Root_sums$Plant_Type=="No plant"] <- "Bulk Soil"

## Add row to improve format 

# Add Row to the DataFrame with out changing data types
Root_cat[nrow(Root_cat) + 1,] <- list("E1", "WTPhi", 
                                      "Wild Type","Phi",NA,0)

ggplot(Root_cat, aes(x = Treatment, y = Root_sums, fill = Plant_Type)) +
  geom_boxplot() +
  theme_bw(base_size = 15)+
  scale_fill_manual(name = "Plant Type", 
                     values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Abundance")+
  geom_signif(y_position = 7500, xmin= 2.6, xmax = 3.4,
              annotation = "p<0.01", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 4100, xmin= 0.6, xmax = 1.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 5000, xmin= 3.6, xmax = 4.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  scale_y_continuous(breaks = pretty_breaks(n = 6))

# Check percentage of root associated taxa?
