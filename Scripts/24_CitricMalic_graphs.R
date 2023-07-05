#### Citrate and malate graphs ####

#### Libraries ####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)

#### Data upload and pre-processing ####

Exud_OG <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Citmal_OG.csv")
Exud_stats <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Citmal_stats.csv")
row.names(Exud_OG) <- Exud_OG$X
Exud_OG <- Exud_OG[,-1]
sapply(Exud_OG,mode)

# order factors

Exud_OG$Treatment <- factor(Exud_OG$Treatment, levels = c("Control","Low P", "Phi", "P Mix"), ordered = TRUE)
Exud_stats$Treatment <- factor(Exud_stats$Treatment, levels = c("Control","Low P", "Phi", "P Mix"), ordered = TRUE)

#### E1 Citric acid ####

## Original data: squared scale BEST

ggplot(Exud_OG, aes(x = Plant_type, y = E1, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), name = "Treatment")+
  xlab("Plant Type") + ylab("Normalized Area")+ 
  scale_y_sqrt()

## Original data: log scale show better difference between plant types

ggplot(Exud_OG, aes(x = Plant_type, y = E1, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Normalized Area")+ 
  scale_y_continuous(trans = "log10")

## Squared transformed: remove scale transformation

ggplot(Exud_stats, aes(x = Plant_type, y = E1_sqrt, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 4000,alpha= .6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

#### E2 Titanium III citrate ####

# Original data

ggplot(Exud_OG, aes(x = Plant_type, y = E2, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = .3,alpha= .6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+ 
  scale_y_continuous(trans = "log10")

## Log transformed: remove scale transformation same graph

ggplot(Exud_stats, aes(x = Plant_type, y = E2_log, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = .7,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

#### E3 Tetrahomocitrate ####

# Original data

ggplot(Exud_OG, aes(x = Plant_type, y = E3, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = .3,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+ 
  scale_y_continuous(trans = "log10")

## Log transformed: remove scale transformation same graph

ggplot(Exud_stats, aes(x = Plant_type, y = E3_log, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = .7,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")


#### E4 (2R_3S)-2_3-Dimethylmalate ####

### Plant type

## no scale: THIS ONE IS BETTER

ggplot(Exud_OG, aes(x = Plant_type, y = E4, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

# squared scale
ggplot(Exud_stats, aes(x = Plant_type, y = E4, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_sqrt()

# log scale: WORST
ggplot(Exud_stats, aes(x = Plant_type, y = E4, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_continuous(trans = "log10")

### Treatment

## no scale 

# dots: NOPE
ggplot(Exud_OG, aes(x = Treatment, y = E4, fill = Plant_type)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("burlywood3","steelblue","lightpink"), name = "Plant type")+
  xlab("Plant Type") + ylab("Concentration")

# boxplot: BETTER

ggplot(Exud_OG, aes(x = Treatment, y = E4, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  scale_y_sqrt()

## Log scale: NOPE

ggplot(Exud_OG, aes(x = Treatment, y = E4, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  scale_y_continuous(trans = "log10")

## Squared scale: BEST

ggplot(Exud_OG, aes(x = Treatment, y = E4, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  scale_y_sqrt()

#### E5 (2R_3S)-3-Isopropylmalate ####

## No scale: shows no differences between plant types but low p higher 
ggplot(Exud_OG, aes(x = Plant_type, y = E5, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

# squared scale: Same as OG
ggplot(Exud_stats, aes(x = Plant_type, y = E5, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_sqrt()

# log scale: NOPE
ggplot(Exud_stats, aes(x = Plant_type, y = E5, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_continuous(trans = "log10")

#### E6 D-(+)-Malic acid ####

### Plant type

## no scale: Low P the lowest

ggplot(Exud_OG, aes(x = Plant_type, y = E6, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

# squared scale: better showing differences between plant types
ggplot(Exud_stats, aes(x = Plant_type, y = E6_sqrt, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

# log scale: BEST showing differences but it is the lowest in LowP
ggplot(Exud_OG, aes(x = Plant_type, y = E6, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_continuous(trans = "log10")

### Treatment

## no scale 

# dots: NOPE but no plant the lowest
ggplot(Exud_OG, aes(x = Treatment, y = E6, fill = Plant_type)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("burlywood3","steelblue","lightpink"), name = "Plant type")+
  xlab("Plant Type") + ylab("Concentration")

# boxplot: BETTER

ggplot(Exud_OG, aes(x = Treatment, y = E6, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")

## Log scale: NOPE

ggplot(Exud_OG, aes(x = Treatment, y = E6, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  scale_y_continuous(trans = "log10")

## Squared scale: BEST

ggplot(Exud_OG, aes(x = Treatment, y = E6, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  scale_y_sqrt()

## Squared data: Same as squared scale but better numbers 

ggplot(Exud_stats, aes(x = Treatment, y = E6_sqrt, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")

#### E7 (R)-Malate ####

## no scale: NO

ggplot(Exud_OG, aes(x = Plant_type, y = E7, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

# squared scale: BEST showing differences between plant types
ggplot(Exud_stats, aes(x = Plant_type, y = E7_sqrt, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

# log scale: NOPE lowers Low P 
ggplot(Exud_OG, aes(x = Plant_type, y = E6, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_continuous(trans = "log10")

#### Citrate compounds cumulative ####

## no scale: NO

ggplot(Exud_OG, aes(x = Plant_type, y = Citrate_sum, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

# squared scale: BEST
ggplot(Exud_OG, aes(x = Plant_type, y = Citrate_sum, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 15,)+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), 
                    name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_sqrt()

## squared data
ggplot(Exud_stats, aes(x = Plant_type, y = Citrate_scum, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("olivedrab3","maroon","cadetblue","plum4"), 
                    name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

# log scale: Better in showing differences between treatments but no LowP
ggplot(Exud_OG, aes(x = Plant_type, y = Citrate_sum, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_continuous(trans = "log10")

#### Malate compounds cumulative ####

### Plant Type

## no scale: P Mix highest

ggplot(Exud_OG, aes(x = Plant_type, y = Malate_sum, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")

# squared scale: Better showing patterns
ggplot(Exud_OG, aes(x = Plant_type, y = Malate_sum, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_sqrt()

# log scale: NOPE
ggplot(Exud_OG, aes(x = Plant_type, y = Malate_sum, fill = Treatment)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("royalblue4","gold2","hotpink","mediumorchid"), name = "Treatment")+
  xlab("Plant Type") + ylab("Concentration")+
  scale_y_continuous(trans = "log10")


### Treatment

## no scale 

# dots: No plant the lowest good graph
ggplot(Exud_OG, aes(x = Treatment, y = Malate_sum, fill = Plant_type)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5,alpha= 0.6, color = "grey100")+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(values = c("burlywood3","steelblue","lightpink"), name = "Plant type")+
  xlab("Plant Type") + ylab("Concentration")

# boxplot: BETTER but no difference between treatments

ggplot(Exud_OG, aes(x = Treatment, y = Malate_sum, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")

## Log scale: NOPE just no plant very low

ggplot(Exud_OG, aes(x = Treatment, y = Malate_sum, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  scale_y_continuous(trans = "log10")

## Squared scale: Same as OG

ggplot(Exud_OG, aes(x = Treatment, y = Malate_sum, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")+
  scale_y_sqrt()

## Squared data: Same as squared scale but better numbers 

ggplot(Exud_stats, aes(x = Treatment, y = Malate_scum, fill = Plant_type)) +
  geom_boxplot() +
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("burlywood3","steelblue","lightpink"))+
  xlab("Treatment") + ylab("Concentration")


