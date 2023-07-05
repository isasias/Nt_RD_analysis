#### Graphs for LGAGs 2022 ####

### Libraries

library(dplyr)
library(tidyverse)
library(ggplot2)
library(wesanderson)

### Import data

Data_graphs <- read.csv("~/Grad_School/Maestria/Processed_data/RSratio_table2.csv")
Data_tba <- read.csv("~/Grad_School/Maestria/Processed_data/Biomass_weights_transf.csv")

### Plant Biomass ###

## Root

Data_graphs$Treatment <- factor(Data_graphs$Treatment, levels = c("Control","Low P", "Phi", "Pi/Phi Mix"), ordered = TRUE)

ggplot(Data_graphs, aes(x = Treatment, y = Roots*1000, fill = Plant_Type)) +
  geom_boxplot()+
  theme_minimal(base_size = 30,)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab("Weight (mg)") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

## Shoot

ggplot(Data_graphs, aes(x = Treatment, y = Shoots*1000, fill = Plant_Type)) +
  geom_boxplot()+
  theme_minimal(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab("Weight (mg)") +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

### Root architecture ### 

Rootsgraphs <- Data_graphs %>%
  filter(!Data_graphs$Plant_Type == "No plant")

## Area

ggplot(Rootsgraphs, aes(x = Treatment, y = root_area, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_minimal(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab(expression(Area ~ cm^2)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

# Transformed

ggplot(Rootsgraphs, aes(x = Treatment, y = root_transf_area, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_minimal(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab(expression(Area ~ cm^2)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

## Length

ggplot(Rootsgraphs, aes(x = Treatment, y = root_length, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_minimal(base_size = 20)+
  scale_fill_manual(name= "Plant Type", values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab("Length (cm)") +
  scale_x_discrete(labels = c("Low P", "Phi","Control", "Pi/Phi Mix"))

# Transformed

ggplot(Rootsgraphs, aes(x = Treatment, y = root_transf_length, fill = Plant_Type)) +
  geom_bar(stat="identity", position=position_dodge())+ # use position dodge to separate color fill into two bars
  theme_minimal(base_size = 20)+
  scale_fill_manual(name= "Plant Type", values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab("Length (cm)") +
  scale_x_discrete(labels = c("Low P", "Phi","Control", "Pi/Phi Mix"))

### Soil parameters ###

## TOC

TOCgraph <- Data_graphs%>%
  filter(!Data_graphs$Sample_ID== "WT_graph")

ggplot(TOCgraph, aes(x = Treatment, y = TOC, fill = Plant_Type)) +
  geom_boxplot() +
  theme_classic(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab(expression(C ~ g ~ Kg^-1)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))

## MBC

ggplot(Data_graphs, aes(x = Treatment, y = MBC, fill = Plant_Type)) +
  geom_boxplot( show.legend = FALSE) +
  theme_classic(base_size = 30)+
  scale_fill_manual(name= NULL, values = wes_palette("Moonrise3"))+
  xlab("Treatment") + ylab(expression(C ~ g ~ Kg^-1)) +
  scale_x_discrete(labels = c("Control","Low P", "Phi", "Pi/Phi Mix"))
