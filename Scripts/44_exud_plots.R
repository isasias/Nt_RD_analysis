#### Root exudate boxplots ####

### Libraries ###

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
library(patchwork)

#### Data upload and pre-processing ####

## Metadata ##
metadata <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Exud_metadata.csv")

## Data ##
Exudates <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Sum_metabolites.csv")
Exudates <- as.data.frame(t(Exudates)) # transpose
colnames(Exudates) <- Exudates[1,]
Exudates <- Exudates[-1,]
i <- c(1:106) # Specify columns you want to change
Exudates[,i] <- apply(Exudates[,i],2, # Specify own function within apply
                      function(x) as.numeric(as.character(x)))
Exudates <- cbind(metadata, Exudates)

### Not significant exudate set

Exud_sums <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Normalized_exudates.csv")
row.names(Exud_sums) <- Exud_sums[,1]
Exud_sums <- Exud_sums[,-1]
NS_names <- colnames(Exud_sums[c(3,67,54,69,70)])
NS_names <- c(NS_names,"Nicotinate_C18P","Vanillate_HILN")

NS_exudates <- Exudates[,NS_names]
NS_exudates <- cbind(metadata,NS_exudates)

### Top exudates

Top_names <- colnames(Exud_sums[c(38,29,68,46,48,51,71,78,77,16,31,
                                  22,25,27,28,41,49,83,84,87,30,12,
                                  15,6,10,26,33,34,40,55,58,62,72,
                                  91,93)]) # 35 total

Top_names <- c(Top_names,"Salicylate_HILN","Benzoate_HILN",
               "Glyoxalate_HILN","Pyruvate_C18N")


# filter by top exudates
Top_Exudates <- Exudates[,Top_names]
Top_Exudates <- cbind(metadata,Top_Exudates)

#### For loops graphs ####

### Not significant 

NS_exudates$Treatment <- factor(NS_exudates$Treatment, levels = c("Low P","Pi","Phi","Pi/Phi mix"), ordered = TRUE)

# plot without a loop
ggplot(NS_exudates, aes(x = Treatment, y = NS_exudates[,3], 
                       col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  ggtitle(colnames(NS_exudates[3]))

## Boxplot by treatment
for(i in 3:ncol(NS_exudates)){
  NSp <- ggplot(NS_exudates, aes(x = Treatment, y = NS_exudates[,i], 
                          col = Plant_type)) +
    geom_boxplot(fill = "snow2",lwd = 0.7)+
    scale_y_sqrt()+
    ggtitle(colnames(NS_exudates[i]))
  print(NSp) # need to add print with base R is not necessary
}

## By plant type

ggplot(NS_exudates, aes(x = Plant_type, y = NS_exudates[,3], fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_bw(base_size = 15)

for(i in 3:ncol(NS_exudates)){
  NSp <- ggplot(NS_exudates, aes(x = Plant_type, y = NS_exudates[,i], fill = Treatment)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
    theme_bw(base_size = 15)+
    ggtitle(colnames(NS_exudates[i]))
  print(NSp) 
}

## Conclusion: 
# Sinapinic, Fumarate, Vanillate, Maleic, Butanoate all HILN show the 
# same pattern WT and ptxd in LP with high concentrations and the other
# samples were lower with no differences. 36Phi had a high variation
# Exud_10 showed no difference between treatments and Pi close to BS 
# with high variation
# Nicotinate high rizhosphere effect but both plants similar concentration
# in al treatments except PM with ptxd lower than WT


### Top exudates ###

Top_Exudates$Treatment <- factor(Top_Exudates$Treatment, levels = c("Low P","Pi","Phi","Pi/Phi mix"), ordered = TRUE)

# plot without a loop
ggplot(Top_Exudates, aes(x = Treatment, y = Top_Exudates[,5], 
                        col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  scale_y_sqrt()+
  theme_bw()+
  ggtitle(colnames(Top_Exudates[5]))

## By treatment
for(i in 3:ncol(Top_Exudates)){
  NSp <- ggplot(Top_Exudates, aes(x = Treatment, y = Top_Exudates[,i], 
                                 col = Plant_type)) +
    geom_boxplot(fill = "snow2",lwd = 0.7)+
    scale_y_sqrt()+
    theme_bw()+
    ggtitle(colnames(Top_Exudates[i]))
  print(NSp)
}


#### Comparisons ####

## Chorismate N vs P

Chor_P <- ggplot(Top_Exudates, aes(x = Treatment, y = Chorismate_C18P, 
                         col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  scale_y_sqrt()+
  theme_bw()+
  ggtitle("Chorismate_C18P")

Chor_N <- ggplot(Top_Exudates, aes(x = Treatment, y = Chorismate_C18N, 
                                   col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  scale_y_sqrt()+
  theme_bw()+
  ggtitle("Chorismate_C18N")

(Chor_P | Chor_N)

## Succinate N vs P

Succ_P <- ggplot(Top_Exudates, aes(x = Treatment, y = Succinate_C18P, 
                                   col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  scale_y_sqrt()+
  theme_bw()+
  ggtitle("Succinate_C18P")

Succ_N <- ggplot(Top_Exudates, aes(x = Treatment, y = Succinate_C18N, 
                                   col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  scale_y_sqrt()+
  theme_bw()+
  ggtitle("Succinate_C18N")

(Succ_P | Succ_N)

## Salicylate HILN vs P

Sali_P <- ggplot(Top_Exudates, aes(x = Treatment, y = Salicylate_C18P, 
                                   col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  scale_y_sqrt()+
  theme_bw()+
  ggtitle("Salicylate_C18P")

Sali_HN <- ggplot(Top_Exudates, aes(x = Treatment, y = Salicylate_HILN, 
                                   col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  scale_y_sqrt()+
  theme_bw()+
  ggtitle("Salicylate_HILN")

(Sali_P | Sali_HN)

## Glutarate P vs N

Glut_P <- ggplot(Top_Exudates, aes(x = Treatment, y = Glutarate_C18P, 
                                   col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  scale_y_sqrt()+
  theme_bw()+
  ggtitle("Glutarate_C18P")

Glut_N <- ggplot(Top_Exudates, aes(x = Treatment, y = Glutarate_C18N, 
                                   col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  scale_y_sqrt()+
  theme_bw()+
  ggtitle("Glutarate_C18N")

(Glut_P | Glut_N)

#### Trend 1 ####

# Well fertilized vs. not well fertilized plants
# Both plants on Low P and WTPM in higher concentrations
# 9 exudates (both of them chorismate)

### Phthalate

ggplot(Top_Exudates, aes(x = Treatment, y = Phthalate_C18P, 
                         col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("Phthalate")+
  scale_y_sqrt()+
  geom_signif(y_position = 77500, xmin= 0.83, xmax = 1.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 44900, xmin= 1.83, xmax = 3.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 77500, xmin= 3.83, xmax = 4.4,
              annotation = "p<0.001", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 1.12, y= 6.9e+09,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4.25, y= 5.7e+09,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 2.6, y= 3e+09,  label="b",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4, y= 2e+09,  label="b",
           color="darkred", size= 4.5, fontface =2)
  
### Chorismate

ggplot(Top_Exudates, aes(x = Treatment, y = Chorismate_C18P, 
                         col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("Chorismate")+
  scale_y_sqrt()+
  geom_signif(y_position = 120000, xmin= 0.83, xmax = 1.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 76000, xmin= 1.83, xmax = 3.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 100000, xmin= 3.83, xmax = 4.4,
              annotation = "p<0.001", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 1.12, y= 1.7e+10,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4.25, y= 9.4e+09,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 2.6, y= 7.37e+09,  label="b",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4, y= 2e+09,  label="c",
           color="darkred", size= 4.5, fontface =2)

#### Trend 2 ####
  
# Similar to 1 but here only WTPM has significantly higher concentrations
# In most Low P is slightly higher than the rest of the treatments
# but not significant, and in others there is no difference
# 16 exudates T2 plus 4 with the variation

### Malate T2

ggplot(Top_Exudates, aes(x = Treatment, y = Malate_C18N, 
                         col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("Malate")+
  scale_y_sqrt()+
  geom_signif(y_position = 40500, xmin= 0.83, xmax = 1.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 25000, xmin= 1.83, xmax = 3.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 60000, xmin= 3.83, xmax = 4.4,
              annotation = "p<0.05", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 1.1, y= 2.1e+09,  label="ab",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4.25, y= 3.4e+09,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 2.6, y= 1.1e+09,  label="b",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4, y= 5e+08,  label="b",
           color="darkred", size= 4.5, fontface =2)

### Shikimate T2 variation

ggplot(Top_Exudates, aes(x = Treatment, y = Shikimate_C18POS, 
                         col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("Shikimate")+
  scale_y_sqrt()+
  geom_signif(y_position = 29000, xmin= 0.83, xmax = 3.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 40000, xmin= 3.83, xmax = 4.4,
              annotation = "p<0.001", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 2.1, y= 1e+09,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4, y= 5e+08,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4.25, y= 1.53e+09,  label="b",
           color="darkred", size= 4.5, fontface =2)

#### Trend 3 ####

# Only 3 exudates the ones selected for thesis 
# Phi intermediate exudation

### Nicotinamide

ggplot(Top_Exudates, aes(x = Treatment, y = Exud_11_C18P, 
                         col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("Nicotinamide")+
  scale_y_sqrt()+
  geom_signif(y_position = 25500, xmin= 0.83, xmax = 1.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 24000, xmin= 1.83, xmax = 3.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 22000, xmin= 3.83, xmax = 4.4,
              annotation = "p<0.05", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 1.1, y= 7.7e+08,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4.25, y= 4.4e+08,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 2.13, y= 2.1e+08,  label="b",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 4, y= 1.5e+08,  label="b",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 3.18, y= 5.2e+08,  label="ab",
           color="darkred", size= 4.5, fontface =2)
  
#### Trend 4 ####

# Only with HILN measurements 5 exudates
# No rhizosphere effect
# Both LP plants with high concentrations and the rest treatments low

### Acetate

ggplot(Top_Exudates, aes(x = Treatment, y = Acetate_HILN, 
                         col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("Acetate")+
  scale_y_sqrt()+
  geom_signif(y_position = 11000, xmin= 0.83, xmax = 1.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 5500, xmin= 1.5, xmax = 4.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 1.1, y= 1.4e+08,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 2.95, y= 5e+07,  label="b",
           color="darkred", size= 4.5, fontface =2)

#### No trends ####

# 4 exudates with their own exudation pattern

### Ferulate: only rizosphere effect

ggplot(Top_Exudates, aes(x = Treatment, y = Ferulate_C18P, 
                         col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("Ferulate")+
  scale_y_sqrt()

### Glucose

ggplot(Top_Exudates, aes(x = Treatment, y = Exud_27_HILN, 
                         col = Plant_type)) +
  geom_boxplot(fill = "snow2",lwd = 0.7)+
  theme_bw(base_size = 15)+
  scale_color_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  ggtitle("D-(+)-Glucose")+
  geom_signif(y_position = 2.5e+09, xmin= 0.83, xmax = 1.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  geom_signif(y_position = 5.3e+09, xmin= 1.5, xmax = 4.4,
              annotation = "NS", tip_length = 0.01,
              col= 1)+
  annotate(geom = "text", x = 1.1, y= 3e+09,  label="a",
           color="darkred", size= 4.5, fontface =2)+
  annotate(geom = "text", x = 2.95, y= 5.8e+09,  label="b",
           color="darkred", size= 4.5, fontface =2)

  
