#### Thesis graphs ####

### Libraries

library(dplyr)
library(tidyverse) # Necessary
library(ggplot2)
library(scales) # customize scales for graphs
library(phyloseq)
library(RColorBrewer)
library(reshape) # melt function has DOI
library(pheatmap)
library(ggforce)
library(ggsignif) # has DOI
library(stringr) # put text in two lines str_wrap function
library(wesanderson)
library(ggbreak)
library(microbiome)

### Import data

Graphs_OG <- read.csv("~/Grad_School/Maestria/Processed_data/RSratio_table.csv") # before uploading I added a 0 in biomass for WT Phi
Data_log <- read.csv("~/Grad_School/Maestria/Processed_data/Biomass_logtransf.csv")

load("Psoil_filt.Rdata")
load("Psoil_rel.RData")
alpha_div <- read.csv("~/Grad_School/Maestria/Processed_data/Alphadiversity.csv")
metabolism <- read.csv("~/Grad_School/Maestria/Processed_data/Metabolic_groups.csv")
energys <- read.csv("~/Grad_School/Maestria/Processed_data/ENERGYSOURCE.filtered.csv")

Exud_Top <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Top_exudates.csv")
Exud_OG <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Citmal_OG.csv")
Exud_stats <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Citmal_stats.csv")
Exud_final <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Final_exudates.csv")

#### Differences in plant growth between WT and ptxDNt-36 ####

Graphs_OG$Treatment <- factor(Graphs_OG$Treatment, levels = c("Control","Low P", "Phi", "Pi/Phi mix"), ordered = TRUE)
Graphs_OG$Plant_Type <- factor(Graphs_OG$Plant_Type, levels = c("ptxDNt-36","Wild Type", "No Plant"), ordered = TRUE)

### Root biomass

ggplot(Graphs_OG, aes(x = Treatment, y = Roots*1000, fill = Plant_Type)) +
  geom_boxplot()+
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Weight (mg)") +
  scale_x_discrete(labels = c("Pi","Low P", "Phi", "Pi/Phi mix"))+
  scale_y_continuous(limits = c(0, 300))+ 
  annotate(geom = "text", x = 1, y= 110,  label="abc        bcd",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2, y= 110,  label="d           d",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 3, y= 110,  label="abc         e",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 4, y= 210,  label="a           d",
           color="black", size= 4.5)

### Shoot biomass

ggplot(Graphs_OG, aes(x = Treatment, y = Shoots*1000, fill = Plant_Type)) +
  geom_boxplot()+
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Weight (mg)") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi mix"))+
  scale_y_continuous(limits = c(0, 1700))+
  annotate(geom = "text", x = 1, y= 770,  label="ab          ab",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2, y= 770,  label="c            bc",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 3, y= 770,  label="bc           d",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 4, y= 1270,  label="a              c",
           color="black", size= 4.5)

### Total biomass

ggplot(Graphs_OG, aes(x = Treatment, y = (Roots + Shoots)*1000, fill = Plant_Type)) +
  geom_boxplot()+
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Weight (mg)") +
  scale_x_discrete(labels = c("Pi","Low P", "Phi", "Pi/Phi mix"))+
  scale_y_continuous(limits = c(0, 2000))+
  geom_signif(y_position = 500, # low P 
              xmin = 1.6, 
              xmax = 2.4,
              annotation = "NS", 
              tip_length = 0.03)+
  geom_signif(y_position = 1250, # Pi and Phi
              xmin = c(0.6,2.6), 
              xmax = c(1.4,3.4),
              annotation = c("NS","***"), 
              tip_length = 0.03)+
  geom_signif(y_position = 2000, # Pi/Phi mix
              xmin = 3.6, 
              xmax = 4.4,
              annotation = "***", 
              tip_length = 0.03)+
  annotate(geom = "text", x = 1, y= 1550,
           label="a",
           color="black", size= 5.5)+
  annotate(geom = "text", x = 3, y= 1550,
           label="a",
           color="black", size= 5.5)+
  annotate(geom = "text", x = 2, y= 1550,
           label="b",
           color="black", size= 5.5)
  

### Root structure ###

# Root Area

Rootsgraphs <- Graphs_OG %>%
  filter(!Graphs_OG$Plant_Type == "No Plant")

Rootsgraphs$Treatment <- factor(Rootsgraphs$Treatment, 
                                levels = c("Control","Low P", "Phi", "Pi/Phi mix"), 
                                ordered = TRUE)

ggplot(Rootsgraphs, aes(x = Treatment, y = root_area, fill = Plant_Type)) +
  geom_boxplot()+
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab(expression(Area ~ cm^2)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi mix"))+
  scale_y_continuous(limits = c(0, 3800))+
  annotate(geom = "text", x = 1, y= 3600,  label="a           abc",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2, y= 3600,  label="bc            bc",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 3, y= 3600,  label="ab              ",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 4, y= 3600,  label="a              bc",
           color="black", size= 4.5)

## Length

ggplot(Rootsgraphs, aes(x = Treatment, y = root_length, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_minimal(base_size = 20)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Length (cm)") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

ggplot(Rootsgraphs, aes(x = Treatment, y = root_length, fill = Plant_Type)) +
  geom_boxplot()+
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Length (cm)") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi mix"))+
  scale_y_continuous(limits = c(0, 700))+
  annotate(geom = "text", x = 1, y= 680,  label="ab           ab",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2, y= 680,  label="c            c",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 3, y= 680,  label="bc              ",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 4, y= 680,  label="ab              c",
           color="black", size= 4.5)

## Tip count 

tips_sdev <- Data_roots %>%
  group_by(Treatment, Plant_Type)%>%
  summarise(dv_tips= sd(root_tip_count),
            mean_tips = mean(root_tip_count))

tips <- as.data.frame(tips_sdev) 

tips <- tips %>%
  add_row(Treatment= "Phi", Plant_Type = "Wild Type", 
          dv_tips = 0, mean_tips = 0)

ggplot(tips, aes(x = Treatment, y = mean_tips, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Number of Tips") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))+
  scale_y_continuous(limits = c(0, 350))+
  geom_errorbar(aes(ymin=mean_tips-dv_tips, ymax=mean_tips+dv_tips), width=.2,
                position=position_dodge(.9)) 

## Lateral roots

Lat_roots <- Rootsgraphs %>%
  add_column(lat_roots= c(21,36,32.5,20.5,	40,	68,	59,
                          61,49.5,53,	26,	29,	83,	81,
                          22,	25,	35,	28.5,	40,	51,	47,	41,
                          43,	55,	24.5,	21,	89,	93, 29,	38,	
                          23,	23,	47.5,	38,	48.5,	56.5,	55,	
                          54,	14,	31,	82,	104,0,0,0,0,0,0))

latr_sdev <- Lat_roots %>%
  group_by(Treatment, Plant_Type)%>%
  summarise(dv_lat= sd(lat_roots),
            mean_lat = mean(lat_roots))

ggplot(latr_sdev, aes(x = Treatment, y = mean_lat, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Number of Lateral Roots") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))+
  scale_y_continuous(limits = c(0, 105))+
  geom_errorbar(aes(ymin=mean_lat-dv_lat, ymax=mean_lat+dv_lat), width=.2,
                position=position_dodge(.9)) 

### Pi biomass ###

ggplot(Data_roots, aes(x = Treatment, y = biomass_Pi, fill = Plant_Type)) +
  geom_boxplot() +
  theme_bw(base_size = 15)+
  scale_fill_manual(values = c("steelblue","lightpink"), name = "Plant Type")+
  xlab("Treatment") + ylab(str_wrap("P concentration mg Pi/g tissue", width = 15)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))+
  scale_y_continuous(limits = c(0, 5.2))


#### Comparison of rhizosphere microbial diversity and composition ####

### Alpha diversity 

Psoil_filt@sam_data$Plant_Type[Psoil_filt@sam_data$Plant_Type == "ptxD-36"] <- "ptxDNt-36"
Psoil_filt@sam_data$Plant_Type <- factor(Psoil_filt@sam_data$Plant_Type, levels = c("ptxDNt-36", "Wild Type", "No plant"))


Shannon <- plot_richness(Psoil_filt, x="Treatment", 
                         measures=c("Observed","Shannon","Simpson","Chao1"), 
                         color="Plant_Type")+
  scale_fill_manual(breaks=c("ptxDNt-36","Wild Type","No plant"),
                    labels=c("ptxDNt-36","Wild Type","Bulk soil"))+
  scale_color_manual(values = c("steelblue","lightpink","burlywood3"), 
                     name = "Plant Type",
                     labels= c("ptxDNt-36","Wild Type","Bulk soil")) +
  theme_bw()

Shannon$layers <- Shannon$layers[-1]
Shannon <- Shannon + geom_point(size=4, alpha=0.7)
Shannon

## Pielou

ggplot(alpha_div, aes(x=Treatment, y=pielou, color = Plant_Type)) +
  geom_point(size = 4, alpha=0.7) +
  scale_color_manual(values = c("steelblue","burlywood3", "lightpink"))+
  theme_bw()

### Beta diversity ###

Psoil_filt@sam_data$Treatment[Psoil_filt@sam_data$Treatment=="Control"] <- "Pi"
bray_bdiv <- phyloseq::distance(Psoil_filt, method = "bray", type = "sample") # check how it looks with taxa type
bray_ord <- ordinate(Psoil_filt,"NMDS", distance = bray_bdiv)

# ordination methods: 
# c("RDA", "NMDS", "MDS", "PCoA")

# good graphs with bray: RDA, NMDS, MDS= PCoA
# jsd: RDA NMDS more or less
# jaccard = bray

p_bray <- plot_ordination(Psoil_filt, bray_ord, "samples", 
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

### Phyla heatmap ###

Phyla_soil <- tax_glom(Psoil_filt,taxrank = "Phylum", NArm = FALSE)


## Extract data from phyloseq object to use pheatmap

OTU_matrix <- as.data.frame(Phyla_soil@otu_table)
Tax_matrix <- as.data.frame(Phyla_soil@tax_table)

## Rename columns from otu_matrix with phylum from tax_matrix

colnames(OTU_matrix) <- Tax_matrix$Phylum
colnames(OTU_matrix)[32] <- "SARS324 Clade"

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

## Create heatmap using pheatmap

# add breaks too

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- seq(min(Phyla_matrix), max(Phyla_matrix), length.out = 10)
mat_breaks <- quantile_breaks(Phyla_matrix, n = 21)

pheatmap(Phyla_matrix,cluster_rows=FALSE, cluster_cols=FALSE, scale = "none",
         color = colorRampPalette(c("snow2","lightskyblue1", "plum3","magenta4","steelblue4", "black"))(20),
         breaks            = mat_breaks,
         drop_levels       = F,
         fontsize          = 10,
         main              = "Microbial Abundance by phyla",
         gaps_col = c(3,6,9,12,15,18,21),
         labels_col = c (rep("36LP",3),rep("WTLP",3),  
                         rep("36Pi",3), rep("WTPi",3),  
                         rep("36Phi",3), 
                         rep("36PM",3), rep("WTPM",3),
                         "BSLP","BSPi","BSPhi","BSPM"))

#### Relative abundance top phyla ####

# Subset top phyla
Top_phyla <- tax_glom(Psoil_rel,taxrank = "Phylum", NArm = FALSE)

Top_phyla <- prune_taxa(names(sort(taxa_sums(Top_phyla),TRUE)[1:5]), Top_phyla)

sorder <- c("A1","B1","C1","A2","B2","C2",
            "A3","B3","C3","A4","B4","C4",
            "A5","B5","C5",
            "A6","B6","C6","A7","B7","C7",
            "D1","D2","D3","D4")
# Top phyla barplot

plot_composition(Top_phyla, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling2"),name = "Phylum", 
                    labels = c("Actinobacteria", "Proteobacteria", "Bacteroidota","Chloroflexi","Acidobacteria"))+
  guides(x=guide_axis(angle = 90))+
  scale_x_discrete( labels=c(rep("36LP",3),rep("WTLP",3),  
                             rep("36Pi",3), rep("WPi",3),  
                             rep("36Phi",3), 
                             rep("36PM",3), rep("WTPM",3),
                             "BSLP","BSPi","BSPhi","BSPM"))+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")+
  scale_y_continuous(breaks = pretty_breaks(n = 10),labels = scales::percent)


### Subset Top Genus ###

Rubro <- subset_taxa(Psoil_rel, Genus == "Rubrobacter") # 86
Pseuda <- subset_taxa(Psoil_rel, Genus == "Pseudarthrobacter") # 7
Micro <- subset_taxa(Psoil_rel, Genus == "Microvirga") # 54
Cell <- subset_taxa(Psoil_rel, Genus == "Cellulomonas") # 18
Sphin <- subset_taxa(Psoil_rel, Genus == "Sphingomonas") # 64
MN <- subset_taxa(Psoil_rel, Genus == "MND1") # 51
Agro <- subset_taxa(Psoil_rel, Genus == "Agromyces") # 6

## Separate in high and low
GenusH <- merge_phyloseq(Rubro,Pseuda,Micro)
GenusH <- tax_glom(GenusH,taxrank = "Genus", NArm = FALSE) 

GenusL <- merge_phyloseq(Sphin,MN,Cell)
GenusL <- tax_glom(GenusL,taxrank = "Genus", NArm = FALSE) 

# High plot
plot_composition(GenusH, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw(base_size = 18)+
  scale_fill_manual(values = wes_palette("IsleofDogs1"),name = "Genus", 
                    labels = c("Rubrobacter", "Pseudarthrobacter", 
                               "Microvirga"))+
  guides(x=guide_axis(angle = 90))+
  scale_x_discrete(labels=c(rep("36LP",3),rep("WTLP",3),  
                            rep("36Pi",3), rep("WPi",3),  
                            rep("36Phi",3), 
                            rep("36PM",3), rep("WTPM",3),
                            "BSLP","BSPi","BSPhi","BSPM"))+
  scale_y_continuous(breaks = pretty_breaks(n = 10),
                     labels = scales::percent)+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

# Low Genus
plot_composition(GenusL, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw(base_size = 18)+
  scale_fill_manual(values = wes_palette("Chevalier1"),name = "Genus", 
                    labels = c("Sphingomonas", "Nitrosomonas", 
                               "Cellulomonas"))+
  guides(x=guide_axis(angle = 90))+
  scale_x_discrete( labels=c(rep("36LP",3),rep("WTLP",3),  
                             rep("36C",3), rep("WTC",3),  
                             rep("36Phi",3), 
                             rep("36PM",3), rep("WTPM",3),
                             "NPLP","NPC","NPPhi","NPPM"))+
  scale_y_continuous(breaks = pretty_breaks(n = 10),
                     labels = scales::percent)+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

#### Functional Heat map ####

row.names(metabolism) <- metabolism[,1]
metabolism <- metabolism[,-1]

# change always to matrix
metabolism <- data.matrix(metabolism)

metabolism <- metabolism[order(metabolism[,1], 
                                   decreasing = TRUE),]

sorder <- c("A1","B1","C1","A2","B2","C2", 
            "A3","B3","C3", "A4","B4","C4", 
            "A5","B5","C5", 
            "A6","B6","C6", "A7","B7","C7", 
            "D1","D2","D3","D4")

metabolism <- metabolism[ , sorder]

row.names(metabolism) <- c("Ammonia Oxidizer", "Sulfate Reducer",
                           "Nitrite Reducer","Dehalogenation",
                           "Sulfide Oxidizer", "Nitrogen Fixation",
                           "Atrazine Metabolism", "Chitin degradation",
                           "Chlorophenol degrading", "Streptomycin Prducer",
                           "Xylan Degrader", "Arom. Hydrocarb. Degrader",
                           "Sulfur Metabolizer","Carbon Fixation",
                           "Sulfur Oxidizer", "Ligning degrader",
                           "Stores Polyhydroxybutyrate","Gramicidin Producer",
                           "Sulfur Reducer","Carbon Monoxide Oxidizer",
                           "Biomass Degrader", "Propionate Metabolism",
                           "Sugars Fermentor", "Dinitrogen Fixing",
                           "Napthalene Degrader","Celullose Degrader")

##Quantile breaks

# function
quantile_breaks <- function(xs, n = 10 ) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# breaks object 20 breaks
mat_breaks <- seq(min(metabolism), max(metabolism), length.out = 10)
mat_breaks <- quantile_breaks(metabolism, n = 21)

# heatmap
pheatmap(metabolism,cluster_rows=FALSE, cluster_cols=FALSE, scale = "none",
         color = colorRampPalette(c("snow2", "lightskyblue1", "plum3","magenta4","steelblue4", "black"))(20),
         breaks            = mat_breaks,
         drop_levels       = F,
         fontsize          = 12,
         gaps_col = c(3,6,9,12,15,18,21),
         labels_col = c (labels=c(rep("36LP",3),rep("WTLP",3),  
                                  rep("36Pi",3), rep("WTPi",3),  
                                  rep("36Phi",3), 
                                  rep("36PM",3), rep("WTPM",3),
                                  "BSLP","BSPi","BSPhi","BSPM")))

### Energy source 

row.names(energys) <- energys[,1]
energys <- energys[,-1]

# change always to matrix
energys <- data.matrix(t(energys))

energys <- energys[order(energys[,1], 
                               decreasing = TRUE),]

sorder <- c("A1","B1","C1","A2","B2","C2", 
            "A3","B3","C3", "A4","B4","C4", 
            "A5","B5","C5", 
            "A6","B6","C6", "A7","B7","C7", 
            "D1","D2","D3","D4")

energys <- energys[ , sorder]

##Quantile breaks

# function
quantile_breaks <- function(xs, n = 10 ) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# breaks object 20 breaks
mat_breaks <- seq(min(energys), max(energys), length.out = 10)
mat_breaks <- quantile_breaks(energys, n = 11)

# heatmap
pheatmap(energys,cluster_rows=FALSE, cluster_cols=FALSE, scale = "none",
         color = colorRampPalette(c("snow2", "lightskyblue1", "plum3","magenta4","steelblue4","black"))(10),
         breaks            = mat_breaks,
         drop_levels       = F,
         fontsize          = 12,
         gaps_col = c(3,6,9,12,15,18,21),
         labels_col = c (labels=c(rep("36LP",3),rep("WTLP",3),  
                                  rep("36Pi",3), rep("WTPi",3),  
                                  rep("36Phi",3), 
                                  rep("36PM",3), rep("WTPM",3),
                                  "BSLP","BSPi","BSPhi","BSPM")))

#### Root exudates ####

## Citrate and malate

# rename variables

Exud_OG$Plant_type[Exud_OG$Plant_type=="No Plant"] <- "Bulk Soil"
Exud_stats$Plant_type[Exud_stats$Plant_type=="No Plant"] <- "Bulk Soil"

# order factors

Exud_OG$Treatment <- factor(Exud_OG$Treatment, levels = c("Pi","Low P", "Phi", "Pi/Phi mix"), ordered = TRUE)
Exud_stats$Treatment <- factor(Exud_stats$Treatment, levels = c("Pi","Low P", "Phi", "P Mix"), ordered = TRUE)


# citrate
ggplot(Exud_OG, aes(x = Plant_type, y = Citrate_sum, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_bw(base_size = 15,)+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), 
                    name = "Treatment",
                    labels = c("Pi", "Low P", "Phi","Pi/Phi mix"))+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_sqrt()+
  annotate(geom = "text", x = 1, y= 3e+08,  label="a",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2, y= 2.5e+09,  label="b",
         color="black", size= 4.5)+
  annotate(geom = "text", x = 3, y= 4.2e+09,  label="b",
           color="black", size= 4.5)
  
  

# malate
ggplot(Exud_OG, aes(x = Treatment, y = Malate_sum, fill = Plant_type)) +
  geom_boxplot() +
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  scale_y_sqrt()+
  annotate(geom = "text", x = 1.1, y= 4e+09,  label="a",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2.1, y= 1.5e+10,  label="b",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 3.18, y= 8.3e+09,  label="ab",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 4.1, y= 4e+09,  label="a",
             color="black", size= 4.5)
    
## Main exudates

row.names(Exud_Top) <- Exud_Top$X
Exud_Top <- Exud_Top[,-1]

### Stacked plot 

## Calculating relative abundance 

Ex_relHi <- Exud_Top%>%
  select(Exud_7,Exud_19,Exud_26,Exud_27,Exud_34)

Ex_relC18 <- Exud_Top%>%
  select(Exud_5,Exud_11,Exud_14,Exud_39,Exud_47)

# Sums
Totals <- c()
for(i in 1:25){
  newv = sum(Ex_relC18[i,])
  Totals <- c(Totals, newv)
}
Totals <- c()
for(i in 1:25){
  newv = sum(Ex_relHi[i,])
  Totals <- c(Totals, newv)
}

# Relative abundance for loop
for(j in 1:5){
  for (i in 1:25) {
    Ex_relC18[i,j]= Ex_relC18[i,j]/Totals[i]
  }
}

for(j in 1:5){
  for (i in 1:25) {
    Ex_relHi[i,j]= Ex_relHi[i,j]/Totals[i]
  }
}
# Check
sum(Ex_relC18[23,])
sum(Ex_relHi[23,])

Ex_relC18 <- Ex_relC18/2
Ex_relHi <- Ex_relHi/2

# Add plant type and tratment

Ex_rel <- cbind(Ex_relC18,Ex_relHi)
Ex_rel <- Ex_rel[,c(4,8,3,7,5,2,9,10,6,1)]
Ex_rel <- Ex_rel[c(1,2,3,13,14,15,
                   4,5,6,16,17,18,
                   7,8,9,
                   10,11,12,19,20,21,
                   22,23,24,25),]

Exud_Top <- Exud_Top[c(1,2,3,13,14,15,
                   4,5,6,16,17,18,
                   7,8,9,
                   10,11,12,19,20,21,
                   22,23,24,25),]

Ex_rel <-  Ex_rel %>%
  add_column(Plant_type = Exud_Top$Plant_type)%>%
  add_column(Treatment = Exud_Top$Treatment)%>%
  add_column(Sample = row.names(Ex_rel))

Erel_melt <- melt(Ex_rel)

#Bar Plot
ggplot(Erel_melt, aes(x= Sample, y=value))+
  geom_bar(stat='identity', aes(fill= variable))+
  theme_bw()+
  scale_fill_brewer(palette = "Set3",name = "Exudates",
                    labels = c("3,5-Dimethoxy-4-hydroxy \n cinnamic acid",
                               "(2S)-2-Isopropyl-3- \n oxosuccinate",
                               "Methylsuccinate","2,3-Dimethylmaleate",
                               "Fumaric diamide", "Nicotinamide", "D-(+)-Glucose", 
                               "Diethyl(2R,3R)-2-methyl-3- \n hydroxysuccinate", 
                               "Fumaric acid", "1,3,6-Octatriene")
  )+ 
  guides(x=guide_axis(angle = 90))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.height= unit(9, 'mm'))+
  scale_x_discrete("Sample", labels= c ("BSLP", "BSPi", "BSPhi","BSPM",
                                        rep("36LP",3), rep("WTLP",3),
                                        rep("36Pi",3), rep("WTPi",3),
                                        rep("36Phi",3),
                                        rep("36PM",3),rep("WTPM",3)))+
  scale_y_continuous(breaks = pretty_breaks(n = 10),labels = scales::percent)+
  ylab("Relative concentration")

## PCA plot

Exud_PCA <- prcomp(Exud_Top[,1:10], center=TRUE, scale. = TRUE)
PCAscores <- Exud_PCA[["x"]]  

PCAscores <- as.data.frame(PCAscores)%>%
  add_column(Plant_type = Exud_Top$Plant_type)%>%
  add_column(Treatment = Exud_Top$Treatment)

# change treatment names
PCAscores$Treatment[PCAscores$Treatment == "P Mix"] <- "Pi/Phi mix"
PCAscores$Treatment[PCAscores$Treatment == "Control"] <- "Pi"
PCAscores$Plant_type[PCAscores$Plant_type == "No Plant"] <- "Bulk Soil"


# change color order with change in names
Exu_ellipse <- ggplot(data = PCAscores, 
                      aes(x= PC1, y= PC2, color = Treatment, shape = Plant_type))+
  geom_point(size= 3)+
  theme_bw()+
  scale_color_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), 
                     name = "Treatment")

Exu_ellipse <- Exu_ellipse +  
  geom_mark_ellipse(aes(fill = Treatment,
                        color = Treatment))+
  scale_fill_manual(values = c("olivedrab3","maroon","plum4","cadetblue"),
                    name = "Treatment")

Exu_ellipse <- Exu_ellipse +
  labs(x = "PCA 1 (49%)",
       y = "PCA 2 (32.2%)")
Exu_ellipse

# Loading plot: to select interesting exudates

# centered not scaled
Exud_PCA <- prcomp(Exud_Top[,1:10], center=TRUE) # scale. = TRUE

PCAloadings <- Exud_PCA$rotation

rownames(PCAloadings)[rownames(PCAloadings) == "Exud_7"] <- "Fumaric Acid"
rownames(PCAloadings)[rownames(PCAloadings) == "Exud_27"] <- "D-(+)-Glucose"
rownames(PCAloadings)[rownames(PCAloadings) == "Exud_19"] <- "2_3-Dimethylmaleate"
rownames(PCAloadings)[rownames(PCAloadings) == "Exud_26"] <- "(2S)-2-Isopropyl-3-oxosuccinate"

PCAloadings <- as.data.frame(PCAloadings)

ggplot(data = PCAloadings, 
       aes(x= PC1, y= PC2))+
  geom_point(shape=16,color="darkblue",size=3)+
  geom_text(aes(label=ifelse(PC1>0.05,
                             rownames(PCAloadings),"")),
            hjust=-.035, vjust=-.78)+
  xlim(-0.05,.9)+ ylim (-.85,.65)+
  labs(x = "PCA 1 (49%)", y = "PCA 2 (32.2%)")+
  theme_bw()

### boxplots exudates

Exud_Top$Treatment[Exud_Top$Treatment == "P Mix"] <- "Pi/Phi mix"
Exud_Top$Treatment[Exud_Top$Treatment == "Control"] <- "Pi"
Exud_Top$Plant_type[Exud_Top$Plant_type == "No Plant"] <- "Bulk Soil"

Exud_Top$Treatment <- factor(Exud_Top$Treatment, levels = c("Pi","Low P", "Phi", "Pi/Phi mix"), ordered = TRUE)

ggplot(Exud_final, aes(x = Treatment, y = Exud_25, 
                     col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)

# Exud 5: as expected
ggplot(Exud_Top, aes(x = Treatment, y = Exud_5, 
                     col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Plant type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("1,3,6-Octatriene")+
  theme(plot.title = element_text(hjust = 0))+
  geom_signif(y_position = 8.12, xmin= 1.9, xmax = 2.4,
    annotation = "NS", tip_length = 0,
    col= 1)+
  scale_y_log10()

# Exud 11: different LowP and C
ggplot(Exud_Top, aes(x = Treatment, y = Exud_11, 
                     col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Plant type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("Nicotinamide")+
  theme(plot.title = element_text(hjust = 0))

#  Exud_27
ggplot(Exud_Top, aes(x = Treatment, y = Exud_27, 
                     col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Plant type", 
                     values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("D-(+)-Glucose")+
  theme(plot.title = element_text(hjust = 0))

# Exud_ 47

ggplot(Exud_Top, aes(x = Treatment, y = Exud_47, 
                     col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name = "Plant type", 
                     values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("Fumaric diamide")+
  theme(plot.title = element_text(hjust = 0))

#### TOC & MBC ####

Graphs_OG$Treatment[Graphs_OG$Treatment == "Control"] <- "Pi"
Graphs_OG$Plant_Type[Graphs_OG$Plant_Type == "No Plant"] <- "Bulk Soil"
Graphs_OG$Treatment <- factor(Graphs_OG$Treatment, levels = c("Pi","Low P", "Phi", "Pi/Phi mix"), ordered = TRUE)

ggplot(Graphs_OG, aes(x = Plant_Type, y = TOC_m, fill = Treatment)) +
  geom_boxplot() +
  theme_bw(base_size = 15)+
  scale_fill_manual(name = "Treatment", 
                    values = c("olivedrab3","maroon","cadetblue","plum4"))+
  xlab("Plant Type") + ylab(expression(C ~ g ~ Kg^-1)) +
  scale_x_discrete(labels = c("Bulk Soil","ptxDNt-36", "Wild Type"))+
  geom_signif(y_position = 10.2, 
              xmin = c(0.51, 1.51, 2.6), 
              xmax = c(1.4, 2.4, 3.45),
              annotation = c("b", "a","a"), 
              tip_length = 0.05)

  
## MBC

ggplot(Graphs_OG, aes(x = Treatment, y = MBC_m, fill = Plant_Type)) +
  geom_boxplot() +
  theme_bw(base_size = 15)+
  scale_fill_manual(name = "Plant type", 
                     values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab(expression(C ~ g ~ Kg^-1)) +
  scale_x_discrete(labels = c("Pi","Low P", "Phi", "Pi/Phi Mix"))+
  annotate(geom = "text", x = 0.75, y= 4,  label="c",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 1.75, y= 4,  label="c",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2.8, y= 4,  label="c",
           color="black", size= 4.5)+  
  annotate(geom = "text", x = 3.75, y= 4,  label="c",
           color="black", size= 4.5)+
  
  annotate(geom = "text", x = 1, y= 7.6,  label="bc",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2, y= 7.6,  label="bc",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2.25, y= 7.6,  label="ab",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 4.25, y= 7.6,  label="ab",
           color="black", size= 4.5)+
  
  annotate(geom = "text", x = 1.25, y= 10.2,  label="abc",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 3.2, y= 8.3,  label="bc",
           color="black", size= 4.5)+  
  annotate(geom = "text", x = 4, y= 10.2,  label="ab",
           color="black", size= 4.5)
  


                                                

  


