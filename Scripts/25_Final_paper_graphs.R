#### Paper graphs ####

### Libraries

library(dplyr)
library(tidyverse) # Necessary
library(ggplot2)
library(scales) # customize scales for graphs
library(phyloseq); packageVersion(phyloseq)
library(RColorBrewer)
library(reshape) # melt function has DOI
library(pheatmap)
library(ggforce)
library(ggsignif) # has DOI
library(stringr) # put text in two lines str_wrap function

### Import data

Graphs_OG <- read.csv("~/Grad_School/Maestria/Processed_data/RSratio_table.csv") # before uploading I added a 0 in biomass for WT Phi
Data_log <- read.csv("~/Grad_School/Maestria/Processed_data/Biomass_logtransf.csv")
load("Psoil_filt.Rdata")
metabolism <- read.csv("~/Grad_School/Maestria/Processed_data/Metabolic_groups.csv")
Exud_Top <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Top_exudates.csv")
Exud_OG <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Citmal_OG.csv")
Exud_stats <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Citmal_stats.csv")


#### Complementary metrics ####

### Root biomass

Graphs_OG$Treatment <- factor(Graphs_OG$Treatment, levels = c("Control","Low P", "Phi", "Pi/Phi mix"), ordered = TRUE)
Graphs_OG$Plant_Type <- factor(Graphs_OG$Plant_Type, levels = c("ptxDNt-36","Wild Type", "No Plant"), ordered = TRUE)


ggplot(Graphs_OG, aes(x = Treatment, y = Roots*1000, fill = Plant_Type)) +
  geom_boxplot()+
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Weight (mg)") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "P Mix"))+
  scale_y_continuous(limits = c(0, 300))+ 
  annotate(geom = "text", x = 1, y= 110,  label="abc    bcd",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 2, y= 110,  label="d       d",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 3, y= 110,  label="abc     e",
           color="black", size= 4.5)+
  annotate(geom = "text", x = 4, y= 210,  label="a       d",
           color="black", size= 4.5)

### Pi Biomass SUPP

ggplot(Graphs_OG, aes(x = Treatment, y = biomass_Pi, fill = Plant_Type)) +
  geom_boxplot() +
  theme_bw(base_size = 15)+
  scale_fill_manual(values = c("steelblue","lightpink"), name = "Plant Type")+
  xlab("Treatment") + ylab(str_wrap("P concentration mg Pi/g tissue", width = 15)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

### MBC

ggplot(Graphs_OG, aes(x = Treatment, y = MBC_m, fill = Plant_Type)) +
  geom_boxplot() +
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","burlywood3","lightpink"))+
  xlab("Treatment") + ylab(expression(MBC ~ g ~ Kg^-1 ~ soil)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))


#### Microbial Diversity ####

### Alpha Diversity

Psoil_filt@sam_data$Plant_Type <- factor(Psoil_filt@sam_data$Plant_Type, levels = c("ptxD-36", "Wild Type", "No plant"))

Shannon <- plot_richness(Psoil_filt, x="Treatment", 
                    measures=c("Shannon"), 
                    color="Plant_Type")+
  scale_fill_manual(breaks=c("ptxD-36","Wild Type","No plant"),
                    labels=c("ptxD-36","Wild Type","No plant"))+
  scale_color_manual(values = c("steelblue","lightpink","burlywood3"), 
                     name = "Plant Type") +
  theme_bw()
  
Shannon$layers <- Shannon$layers[-1]
Shannon <- Shannon + geom_point(size=4, alpha=0.7) 

### Beta Diversity

bray_bdiv <- phyloseq::distance(Psoil_filt, method = "bray", type = "sample") # check how it looks with taxa type
bray_ord <- ordinate(Psoil_filt,"NMDS", distance = bray_bdiv)

p_bray <- plot_ordination(Psoil_filt, bray_ord, "samples", 
                          color = "Plant_Type",
                          shape = "Treatment")+
  scale_color_manual(values = c("steelblue","lightpink","burlywood3"), name = "Plant Type")+
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
  scale_fill_manual(values = c("steelblue","lightpink","burlywood3"), 
                    name = "Plant Type")

### Functional Heat map

row.names(metabolism) <- metabolism[,1]
metabolism <- metabolism[,-1]

# change always to matrix
metabolism <- data.matrix(metabolism)

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
         fontsize          = 10,
         main              = "Community metabolic profiles",
         gaps_col = c(4,7,10,13,16,19,22,25),
         labels_col = c ("NPLP", "NPC", "NPPhi","NPPM",
                         rep("36LP",3),rep("36C",3),
                         rep("36Phi",3),rep("36PM",3),
                         rep("WTLP",3),rep("WTC",3),rep("WTPM",3)))

### Phyla heatmap

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

sorder <- c("D1","D2","D3","D4",
            "A1","B1","C1","A2","B2","C2","A3","B3","C3",
            "A4","B4","C4","A5","B5","C5","A6","B6","C6",
            "A7","B7","C7")
Phyla_matrix <- Phyla_matrix[ , sorder]

## Create heatmap using pheatmap

# add breaks too
mat_breaks <- seq(min(Phyla_matrix), max(Phyla_matrix), length.out = 10)
mat_breaks <- quantile_breaks(Phyla_matrix, n = 21)

pheatmap(Phyla_matrix,cluster_rows=FALSE, cluster_cols=FALSE, scale = "none",
         color = colorRampPalette(c("snow2","lightskyblue1", "plum3","magenta4","steelblue4", "black"))(20),
         breaks            = mat_breaks,
         drop_levels       = F,
         fontsize          = 10,
         main              = "Microbial Abundance by phyla",
         gaps_col = c(4,7,10,13,16,19,22,25),
         labels_col = c ("NPLP", "NPC", "NPPhi","NPPM",
                         rep("36LP",3),rep("36C",3),
                         rep("36Phi",3),rep("36PM",3),
                         rep("WTLP",3),rep("WTC",3),rep("WTPM",3)))


#### Root Exudates ####

row.names(Exud_Top) <- Exud_Top$X
Exud_Top <- Exud_Top[,-1]

### Stacked plot 

## Calculating relative abundance 

Ex_rel <- Exud_Top[,1:10]

# Sums
Totals <- c()
for(i in 1:25){
  newv = sum(Ex_rel[i,])
  Totals <- c(Totals, newv)
}

# Relative abundance for loop
for(j in 1:10){
  for (i in 1:25) {
    Ex_rel[i,j]= Ex_rel[i,j]/Totals[i]
  }
}

# Check
sum(Ex_rel[23,])

# Add plant type and tratment

Ex_rel <-  Ex_rel %>%
  add_column(Plant_type = Exud_Top$Plant_type)%>%
  add_column(Treatment = Exud_Top$Treatment)%>%
  add_column(Sample = row.names(Ex_rel))

Erel_melt <- melt(Ex_rel)

### Bar Plot: 

ggplot(Erel_melt, aes(x= Sample, y=value*100))+
  geom_bar(stat='identity', aes(fill= variable))+
  theme_bw()+
  scale_fill_brewer(palette = "Set3",name = "Exudates",
                    labels = c("1,3,6-Octatriene", "Fumaric acid",
                               "Nicotinamide", "Methylsuccinate", "2_3-Dimethylmaleate",
                               "(2S)-2-Isopropyl-3- \n oxosuccinate", "D-(+)-Glucose", 
                               "Diethyl(2R_3R)-2-methyl-3- \n hydroxysuccinate", 
                               "3_5-Dimethoxy-4-hydroxy \n cinnamic acid",
                               "N,N-bis[(S)-1-methoxy \n carbonylethyl"))+
  guides(x=guide_axis(angle = 35))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.height= unit(9, 'mm'))+
  scale_x_discrete("Sample", labels= c ("NPLP", "NPC", "NPPhi","NPPM",
          rep("36LP",3),rep("36C",3),
          rep("36Phi",3),rep("36PM",3),
          rep("WTLP",3),rep("WTC",3),rep("WTPM",3)))+
  ylab("Compound % per Sample")

### Exudate Heatmap 

## Transpose
Exud_trans <- as.data.frame(t(Exud_Top[1:10]))

## Change to data matrix

Exud_trans <- data.matrix(Exud_trans)
exud_breaks <- seq(min(Exud_trans), max(Exud_trans), length.out = 10)

## Change column names and order

samplenames <- c("A1","B1","C1","A2","B2","C2","A3","B3","C3",
            "A4","B4","C4","A5","B5","C5","A6","B6","C6",
            "A7","B7","C7","D1","D2","D3","D4")

colnames(Exud_trans) <- samplenames

sorder <- c("D1","D2","D3","D4",
            "A1","B1","C1","A2","B2","C2","A3","B3","C3",
            "A4","B4","C4","A5","B5","C5","A6","B6","C6",
            "A7","B7","C7")

Exud_trans <- Exud_trans[ , sorder]

## Breaks 
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

exud_breaks <- quantile_breaks(Exud_trans, n = 11)

## Final graph

pheatmap(Exud_trans, cluster_cols=FALSE, scale = "none",
         color = colorRampPalette(c("snow2","magenta4", "navyblue"))(10),
         breaks= exud_breaks,
         drop_levels = F,
         fontsize = 10,
         gaps_col = c(4,7,10,13,16,19,22,25),
         border_color      = "black",
         labels_col = c ("NPLP", "NPC", "NPPhi","NPPM",
                         rep("36LP",3),rep("36C",3),
                         rep("36Phi",3),rep("36PM",3),
                         rep("WTLP",3),rep("WTC",3),rep("WTPM",3)),
         labels_row = c("1,3,6-Octatriene", "Fumaric acid",
                        "Nicotinamide", "Methylsuccinate", "2_3-Dimethylmaleate",
                        "(2S)-2-Isopropyl-3- \n oxosuccinate", "D-(+)-Glucose", 
                        "Diethyl(2R_3R)-2-methyl-3- \n hydroxysuccinate", 
                        "3_5-Dimethoxy-4-hydroxy \n cinnamic acid",
                        "N,N-bis[(S)-1-methoxy \n carbonylethyl"))

### Distance plots 

# PCA plot

Exud_PCA <- prcomp(Exud_Top[,1:10], center=TRUE, scale. = TRUE)
PCAscores <- Exud_PCA[["x"]]  

PCAscores <- as.data.frame(PCAscores)%>%
  add_column(Plant_type = Exud_Top$Plant_type)%>%
  add_column(Treatment = Exud_Top$Treatment)

Exu_ellipse <- ggplot(data = PCAscores, 
                      aes(x= PC1, y= PC2, color = Treatment, shape = Plant_type))+
  geom_point(size= 3)+
  theme_bw()+
  scale_color_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), 
                     name = "Treatment")

Exu_ellipse <- Exu_ellipse +  
  geom_mark_ellipse(aes(fill = Treatment,
                        color = Treatment))+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"),
                    name = "Treatment")

Exu_ellipse <- Exu_ellipse +
  labs(x = "PCA 1 (49%)",
       y = "PCA 2 (32.2%)")

# Loading plot

# centered not scaled
Exud_PCA <- prcomp(Exud_Top[,1:10], center=TRUE)

PCAloadings <- Exud_PCA$rotation

rownames(PCAloadings)[rownames(PCAloadings) == "Exud_7"] <- "Fumaric Acid"
rownames(PCAloadings)[rownames(PCAloadings) == "Exud_27"] <- "D-(+)-Glucose"
rownames(PCAloadings)[rownames(PCAloadings) == "Exud_19"] <- "2_3-Dimethylmaleate"
rownames(PCAloadings)[rownames(PCAloadings) == "Exud_26"] <- "(2S)-2-Isopropyl-3-oxosuccinate"

ggplot(data = PCAloadings, 
       aes(x= PC1, y= PC2))+
  geom_point(shape=16,color="darkblue",size=3)+
  geom_text(aes(label=ifelse(PC1>0.05,
                             rownames(PCAloadings),"")),
            hjust=-.035, vjust=-.78)+
  xlim(-0.05,.9)+ ylim (-.85,.65)+
  labs(x = "PCA 1 (49%)", y = "PCA 2 (32.2%)")+
  theme_bw()

#### Citrate and malate graphs ####

# Preprocess data
row.names(Exud_OG) <- Exud_OG$X
Exud_OG <- Exud_OG[,-1]
sapply(Exud_OG,mode)

# Order data
Exud_OG$Treatment <- factor(Exud_OG$Treatment, levels = c("Control","Low P", "Phi", "P Mix"), ordered = TRUE)
Exud_stats$Treatment <- factor(Exud_stats$Treatment, levels = c("Control","Low P", "Phi", "P Mix"), ordered = TRUE)

### Cumulative citrate graph ###

# Squared scale
ggplot(Exud_OG, aes(x = Plant_type, y = Citrate_sum, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_bw(base_size = 15,)+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), 
                    name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_sqrt()

### Cumulative malate graph ###

# Squared scale
ggplot(Exud_OG, aes(x = Treatment, y = Malate_sum, fill = Plant_type)) +
  geom_boxplot() +
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  geom_signif(
    comparisons = list(c("Low P", "P Mix")),
    map_signif_level = TRUE)+
  scale_y_sqrt()

### Citric Acid SUPP

ggplot(Exud_OG, aes(x = Plant_type, y = E1, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_bw(base_size = 15)+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), 
                    name = "Treatment")+
  xlab("Plant Type") + ylab("Normalized Area")+
  scale_y_sqrt()

### Titanium III citrate SUPP

ggplot(Exud_OG, aes(x = Plant_type, y = E2, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= .6, color = "grey100")+
  theme_bw(base_size = 15)+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"),
                    name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+ 
  scale_y_continuous(trans = "log10")

### Tetrahomocitrate SUPP

ggplot(Exud_OG, aes(x = Plant_type, y = E3, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_bw(base_size = 15)+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+ 
  scale_y_continuous(trans = "log10")

### (2R_3S)-2_3-Dimethylmalate by treatment SUPP

ggplot(Exud_OG, aes(x = Treatment, y = E4, fill = Plant_type)) +
  geom_boxplot() +
  theme_bw(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  scale_y_sqrt()

### (R)-Malate SUPP

ggplot(Exud_OG, aes(x = Plant_type, y = E7, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_bw(base_size = 15)+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+ 
  scale_y_sqrt()
