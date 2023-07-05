#### Main exudates potential graphs ####

#### Libraries ####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(pheatmap)
library(ggforce)
library(reshape)

### Data upload and preprocessing

Exud_Top <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Final_exudates.csv")
E_norm <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Normalized_exudates.csv")
E_nonp <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/UnNorm_exudates.csv")

## Select relevant exudates

row.names(Exud_Top) <- Exud_Top$X
Exud_Top <- Exud_Top[,-1]
Exud_Top <- Exud_Top[,c(5,7,11,14,19,26,27,34,39,47)]

#### Stacked bar plot ####

## Transform to relative abundance 

Ex_rel <- Exud_Top

# Sums
Totals <- c()
for(i in 1:25){
  newv = sum(Ex_rel[i,])
  Totals <- c(Totals, newv)
}

## Calculating relative abundance
for(j in 1:10){
  for (i in 1:25) {
    Ex_rel[i,j]= Ex_rel[i,j]/Totals[i]
  }
}

sum(Ex_rel[23,])

## Overall abundance

Ex_total <- Exud_Top

# Total sum  
Supertotal <- sum(Totals)

# Comparative abundance
for(j in 1:10){
  for (i in 1:25) {
    Ex_total[i,j]= Ex_total[i,j]/Supertotal
  }
}

## Add plant type and tratment

# Relative
Ex_rel <-  Ex_rel %>%
  add_column(Plant_type = E_norm$Plant_type)%>%
  add_column(Treatment = E_norm$Treatment)%>%
  add_column(Sample = row.names(Ex_rel))

# Total
Ex_total <-  Ex_total %>%
  add_column(Plant_type = E_norm$Plant_type)%>%
  add_column(Treatment = E_norm$Treatment)%>%
  add_column(Sample = row.names(Ex_rel))


Erel_melt <- melt(Ex_rel)
Etot_melt <- melt(Ex_total)

### Bar Plot: 

### Over 1
ggplot(Erel_melt, aes(x= Sample, y=value*100))+
  geom_bar(stat='identity', aes(fill= variable))+
  theme_bw()+
  scale_fill_brewer(palette = "Set3",name = "Exudates")+
  guides(x=guide_axis(angle = 35))+
  theme(axis.ticks.x = element_blank())

### Exudation amount

ggplot(Etot_melt, aes(x= Sample, y=value*100))+
  geom_bar(stat='identity', aes(fill= variable))+
  theme_bw()+
  scale_fill_brewer(palette = "Set3",name = "Exudates")+
  guides(x=guide_axis(angle = 35))+
  theme(axis.ticks.x = element_blank())


#### Heatmap ####

## Transpose
Exud_trans <- as.data.frame(t(Exud_Top))

## Change to data matrix

Exud_trans <- data.matrix(Exud_trans)

# Basic graph no column clustering no scale
pheatmap(Exud_trans, cluster_cols=FALSE, scale = "none")

### Breaks 

exud_breaks <- seq(min(Exud_trans), max(Exud_trans), length.out = 10)

# breaks function
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

exud_breaks <- quantile_breaks(Exud_trans, n = 11)

### Final graph
pheatmap(Exud_trans, cluster_cols=FALSE, scale = "none",
         color = colorRampPalette(c("snow2","magenta4", "navyblue"))(10),
         breaks= exud_breaks,
                  drop_levels = F,
                  fontsize = 10,
         gaps_col = c(3,6,9,12,15,18,21,25),
         border_color      = "black")

#### PCA plots ####

# PCA plot: centered and scale

Exud_Top <-  Exud_Top %>%
  add_column(Plant_type = E_norm$Plant_type)%>%
  add_column(Treatment = E_norm$Treatment)

Exud_PCA <- prcomp(Exud_Top[,1:10], center=TRUE, scale. = TRUE)

Exu_plot <- autoplot(Exud_PCA, data = Exud_Top, 
         colour = "Plant_type", shape = "Treatment", size= 3)+
  theme_bw()+
  scale_color_manual(values = c("burlywood3","steelblue","lightpink"), name = "Plant Type")

shape_names <- c("Control","Low P","Phi", "P Mix")
p.shape <- 15:(15 + length(shape_names) - 1)
names(p.shape) <- shape_names
p.shape["Treatments"] <- 4

Exu_plot <- Exu_plot +
  scale_shape_manual(values = p.shape)

ggplotly(Exu_plot)
# stat_ellipse(geom = "polygon",
#             aes(fill = "Plant_type"), alpha= .2)

PCAscores <- Exud_PCA[["x"]]  

PCAscores <- as.data.frame(PCAscores)%>%
  add_column(Plant_type = E_norm$Plant_type)%>%
  add_column(Treatment = E_norm$Treatment)

PCAscores$Plant_type <- as.factor(PCAscores$Plant_type)


Exu_ellipse <- ggplot(data = PCAscores, 
       aes(x= PC1, y= PC2, color = Treatment, shape = Plant_type))+
  geom_point(size= 3)+
  theme_bw()+
  scale_color_manual(values = c("cornflowerblue","gold2","hotpink","mediumorchid"), 
                     name = "Treatment")
  
Exu_ellipse <- Exu_ellipse +  
  geom_mark_ellipse(aes(fill = Treatment,
                        color = Treatment))+
  scale_fill_manual(values = c("cornflowerblue","gold2","hotpink","mediumorchid"),
                                           name = "Treatment")

Exu_ellipse <- Exu_ellipse +
  labs(x = "PCA 1 (49%)",
       y = "PCA 2 (32.2%)")

ggplotly(Exu_ellipse)

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

#### Save top exudates ####

write.csv(Exud_Top, "~/Grad_School/Maestria/Processed_data/Exudates/Top_exudates.csv")




