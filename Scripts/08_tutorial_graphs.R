#### Graphs for LGAGs 2022 ####

### Libraries

library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)

### Import data

Graphs_OG <- read.csv("~/Grad_School/Maestria/Processed_data/RSratio_table.csv") # before uploading I added a 0 in biomass for WT Phi
Data_log <- read.csv("~/Grad_School/Maestria/Processed_data/Biomass_logtransf.csv")

### Plant Biomass ###

## Roots

Graphs_OG$Treatment <- factor(Graphs_OG$Treatment, levels = c("Control","Low P", "Phi", "Pi/Phi mix"), ordered = TRUE)

ggplot(Graphs_OG, aes(x = Treatment, y = Roots*1000, fill = Plant_Type)) +
  geom_boxplot()+
  theme_light(base_size = 15)+
  scale_fill_manual(name= "Plant Type", values = c("steelblue","lightpink"))+
  xlab("Treatment") + ylab("Weight (mg)") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))+
  scale_y_continuous(limits = c(0, 300))+ 
  annotate(geom = "text", x = 1, y= 110,  label="abc       bcd",
           color="black", size= 5)+
  annotate(geom = "text", x = 2, y= 110,  label="d           d",
           color="black", size= 5)+
  annotate(geom = "text", x = 3, y= 110,  label="abc       e",
           color="black", size= 5)+
  annotate(geom = "text", x = 4, y= 210,  label="a           d",
           color="black", size= 5)

# Check graph with log transformed data: did not work

# Data_log$Treatment <- factor(Data_log$Treatment, levels = c("Control","Low P", "Phi", "Pi/Phi Mix"), ordered = TRUE)
# ggplot(Data_log, aes(x = Treatment, y = log_root, fill = Plant_Type)) +
#   geom_boxplot()+
#   theme_minimal(base_size = 30,)+
#   scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
#   xlab("Treatment") + ylab("Weight (mg)") +
#   scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

## Shoot

ggplot(Graphs_OG, aes(x = Treatment, y = Shoots*1000, fill = Plant_Type)) +
  geom_boxplot()+
  theme_minimal(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab("Weight (mg)") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

## Root/Shoot ratio

ggplot(Graphs_OG, aes(x = Treatment, y = Root_shoot_ratio, fill = Plant_Type)) +
  geom_boxplot()+
  theme_minimal(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab("Ratio") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

# lower values mean more root development compared to shoots... however values are pretty similar between treatments

### Root architecture ### 

Rootsgraphs <- Graphs_OG %>%
  filter(!Graphs_OG$Plant_Type == "No plant")

## Area

Rootsgraphs$Treatment <- factor(Rootsgraphs$Treatment, levels = c("Control","Low P", "Phi", "Pi/Phi mix"), ordered = TRUE)

ggplot(Rootsgraphs, aes(x = Treatment, y = root_area/100, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_minimal(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab(expression(Area ~ dm^2)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

## Length

ggplot(Rootsgraphs, aes(x = Treatment, y = root_length, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_minimal(base_size = 20)+
  scale_fill_manual(name= "Plant Type", values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab("Length (cm)") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

## Tip count 

ggplot(Rootsgraphs, aes(x = Treatment, y = root_tip_count, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_minimal(base_size = 20)+
  scale_fill_manual(name= "Plant Type", values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab("Tips") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

### Soil parameters ###

## TOC

ggplot(Graphs_OG, aes(x = Treatment, y = TOC, fill = Plant_Type)) +
  geom_boxplot() +
  theme_classic(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Darjeeling2"))+ # limits = c ("36ptxD","Wild Type", "No Plant") did not work
  xlab("Treatment") + ylab(expression(C ~ g ~ Kg^-1)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

## MBC

ggplot(Graphs_OG, aes(x = Treatment, y = MBC, fill = Plant_Type)) +
  geom_boxplot() +
  theme_classic(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Darjeeling2"))+
  xlab("Treatment") + ylab(expression(C ~ g ~ Kg^-1)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

## Inorganic Pi

# Roots
ggplot(Graphs_OG, aes(x = Treatment, y = root_Pi, fill = Plant_Type)) +
  geom_boxplot( show.legend = FALSE) +
  theme_minimal(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab(expression(mg~ "Pi" ~ g^-1)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

# Shoots
ggplot(Graphs_OG, aes(x = Treatment, y = shoot_Pi, fill = Plant_Type)) +
  geom_boxplot( show.legend = FALSE) +
  theme_minimal(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab(expression(mg~ "Pi" ~ g^-1)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

# Biomass

ggplot(Graphs_OG, aes(x = Treatment, y = biomass_Pi, fill = Plant_Type)) +
  geom_boxplot() +
  theme_minimal(base_size = 15)+
  scale_fill_manual(values = c("steelblue","lightpink"), name = "Plant Type")+
  xlab("Treatment") + ylab(str_wrap("P concentration mg Pi/g soil", width = 15)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))
