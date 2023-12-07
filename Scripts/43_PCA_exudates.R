#### ROOT EXUDATE PCA PLOTS ####

### Libraries ###

library(dplyr)
library(tidyverse) # Necessary
library(ggplot2)
library(scales) # customize scales for graphs
library(RColorBrewer)
library(reshape) # melt function has DOI
library(ggforce)
library(ggsignif) # has DOI
library(stringr) # put text in two lines str_wrap function
library(wesanderson)
library(ggbreak)
library(microbiome)

#### Data upload and pre-processing ####

## Metadata ##
metadata <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Exud_metadata.csv")

## Data for PCA: filling NAs ##
Exudates <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Sum_metabolites_noNAs.csv")
Exudates <- as.data.frame(t(Exudates)) # transpose
colnames(Exudates) <- Exudates[1,]
Exudates <- Exudates[-1,]
i <- c(1:106) # Specify columns you want to change
Exudates[,i] <- apply(Exudates[,i],2, # Specify own function within apply
                       function(x) as.numeric(as.character(x)))
Exudates <- cbind(metadata, Exudates)

## Top exudate names ##

Exud_sums <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Normalized_exudates.csv")
row.names(Exud_sums) <- Exud_sums[,1]
Exud_sums <- Exud_sums[,-1]
Top_names <- colnames(Exud_sums[c(38,29,68,46,48,51,71,78,77,16,31,
                                  22,25,27,28,41,49,83,84,87,30,12,
                                  15,6,10,26,33,34,40,55,58,62,72,
                                  91,93)]) # 35 total


#### PCA plot: all metabolites ####

# remove NA values using na.omit removed all the first row NOT USE
# Remove columns with NA values
Exud_PCA <- Exudates%>%
  select_if(~ !any(is.na(.)))

# reduced to 99 exudates
PCAe <- prcomp(Exud_PCA[,3:101], center=TRUE, 
                   scale. = TRUE)

summary(PCAe) # to check percentages for PCA
# proportion of variance gives PCA percentages
# PC1 50.29% PC2 20.24%

PCAscores <- PCAe[["x"]]  

PCAscores <- as.data.frame(PCAscores)%>%
  add_column(Plant_type = metadata$Plant_type)%>%
  add_column(Treatment = metadata$Treatment)

### Plot

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
  labs(x = "PCA 1 (50.29%)",
       y = "PCA 2 (20.24%)")
Exu_ellipse

## Using all exudates creates more overlap

#### PCA plot: Top cases ####

Norm_exud <- Exudates[,Top_names]
Unnorm_exud <- Exudates[,c("Salicylate_HILN","Benzoate_HILN",
                           "Glyoxalate_HILN","Pyruvate_C18N")]
Exudates_f <- cbind(metadata,Norm_exud,Unnorm_exud) # reduced to 39

### PCA values ###
PCAc <- prcomp(Exudates_f[,3:41], center=TRUE, 
               scale. = TRUE)

summary(PCAc) # to check percentages for PCA
# proportion of variance gives PCA percentages
# PC1 65.36% PC2 15.06%

PCAscores <- PCAc[["x"]]  

PCAscores <- as.data.frame(PCAscores)%>%
  add_column(Plant_type = metadata$Plant_type)%>%
  add_column(Treatment = metadata$Treatment)

### Plot

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
  labs(x = "PCA 1 (65.36%)",
       y = "PCA 2 (15.06%)")
Exu_ellipse

## This one has better separation between groups
## Highest variation WTPM

#### PCA of significants ####

Top_exudates <- c("Phthalate_C18P","Chorismate_C18P","Maleic_C18N",
                  "Acetate_HILN", "Nicotinate_C18N","Malate_C18N",
                  "Mugineate_C18N","Valarate_C18N","Tartarate_C18P",
                  "Exud_47_C18P","Ferulate_C18P","Salicylate_HILN",
                  "Benzoate_HILN") 
# 11 normal compounds, 2 nonparametric

# filter by top exudates
Exudates_f <- Exudates[,Top_exudates]
Exudates_f <- cbind(metadata,Exudates_f)

## PCA

PCAf <- prcomp(Exudates_f[,3:15], center=TRUE, 
               scale. = TRUE)

summary(PCAf) # to check percentages for PCA
# proportion of variance gives PCA percentages
# PC1 57.61% PC2 28.65%

PCAscores <- PCAf[["x"]]  

PCAscores <- as.data.frame(PCAscores)%>%
  add_column(Plant_type = metadata$Plant_type)%>%
  add_column(Treatment = metadata$Treatment)

### Plot

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
  labs(x = "PCA 1 (57.61%)",
       y = "PCA 2 (28.65%)")
Exu_ellipse



