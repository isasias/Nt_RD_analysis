# Data analysis Preview for abstract #

## Libraries
library(dplyr)
library(tidyverse)

## Import data 
Exp1_data <- read.csv("~/Grad_School/Maestria/Nt_RD_analysis/Exp1_data.csv")
 
#### Function original vs. transformed root metrics ####

root_transformation <- function(root_og, root_tr) {
  x <- root_tr - root_og
  percentage_inc <- (x * 100)/ root_og
  return(percentage_inc)
}

#### Adding percentage of change to dataframe ####

### Length
Exp1_data <- Exp1_data %>%
  add_column (Pr_inc_length = root_transformation(Exp1_data$root_length,Exp1_data$root_transf_length))

Exp1_data$Pr_inc_length <- format(round (Exp1_data$Pr_inc_length,2), nsmall = 2) # define decimal places

### Area 
Exp1_data <- Exp1_data %>%
   add_column(Pr_inc_area = root_transformation(Exp1_data$root_area,Exp1_data$root_transf_area))
 
Exp1_data$Pr_inc_area <- format(round (Exp1_data$Pr_inc_area,2), nsmall = 2)
 
### Tip count
Exp1_data <- Exp1_data %>%
  add_column(Pr_dec_tip = abs(root_transformation(Exp1_data$root_tip_count,Exp1_data$root_transf_tip_count)))

Exp1_data$Pr_dec_tip <- format(round (Exp1_data$Pr_dec_tip,2), nsmall = 2)

#### Statistical Analyses ####

## t.test for root metrics

t.test(Exp1_data$root_length, Exp1_data$root_transf_length, alternative = c("two.sided"), mu=0, paired = T)
t.test(Exp1_data$root_area, Exp1_data$root_transf_area, alternative = c("two.sided"), mu=0, paired = T)
t.test(Exp1_data$root_tip_count, Exp1_data$root_transf_tip_count, alternative = c("two.sided"), mu=0, paired = T)

## t.test for root and shoot weights

### Filtering data for WT
WT_data <- Exp1_data[Exp1_data$Plant_Type == "Wild Type",]
WT_data <- WT_data[WT_data$Treatment == "Low P" | WT_data$Treatment == "Pi  optimal nutrition",]

### roots
t.test(WT_data$Roots_.g.~WT_data$Treatment)

### shoots
t.test(WT_data$Shoots_.g.~WT_data$Treatment)

### root morphology
t.test(WT_data$root_length~WT_data$Treatment)
t.test(WT_data$root_area~WT_data$Treatment)
summary(aov(WT_data$root_area~WT_data$Treatment))
summary(aov(WT_data$root_length~WT_data$Treatment))

## 2-way ANOVAs: See if plant type and treatment have statistical difference

### For root and shoot weights
summary(aov(Exp1_data$Roots_.g.~ Exp1_data$Treatment + Exp1_data$Plant_Type))
summary(aov(Exp1_data$Shoots_.g.~ Exp1_data$Treatment + Exp1_data$Plant_Type))

### For MBC and DOC
summary(aov(Exp1_data$DOC_.ug.g.1.~ Exp1_data$Treatment + Exp1_data$Plant_Type))
summary(aov(Exp1_data$MBC._.ug.g.1.~ Exp1_data$Treatment + Exp1_data$Plant_Type))

### For root metrics
summary(aov(Exp1_data$root_transf_length~ Exp1_data$Treatment + Exp1_data$Plant_Type))
summary(aov(Exp1_data$root_transf_area~ Exp1_data$Treatment + Exp1_data$Plant_Type))
summary(aov(Exp1_data$root_transf_tip_count~ Exp1_data$Treatment + Exp1_data$Plant_Type))

## t.test for differences between plant types 

### weights
t.test(Exp1_data$Roots_.g.~Exp1_data$Plant_Type)
t.test(Exp1_data$Roots_.g.~Exp1_data$Plant_Type)

### MBC & DOC Anova
summary(aov(Exp1_data$DOC_.ug.g.1.~Exp1_data$Plant_Type))
summary(aov(Exp1_data$MBC._.ug.g.1.~Exp1_data$Plant_Type))

### root metrics
t.test(Exp1_data$root_length~Exp1_data$Plant_Type)
t.test(Exp1_data$root_area~Exp1_data$Plant_Type)
t.test(Exp1_data$root_tip_count~Exp1_data$Plant_Type)

## Effects of P-mix
Pmix_data <- Exp1_data[Exp1_data$Treatment == "Pi/Phi mix",]

### t.tests
t.test(Pmix_data$Roots_.g.~Pmix_data$Plant_Type)
t.test(Pmix_data$Shoots_.g.~Pmix_data$Plant_Type)
t.test(Pmix_data$root_length~Pmix_data$Plant_Type)
t.test(Pmix_data$root_area~Pmix_data$Plant_Type)
t.test(Pmix_data$root_tip_count~Pmix_data$Plant_Type)

summary(aov(Pmix_data$MBC._.ug.g.1.~Pmix_data$Plant_Type))
summary(aov(Pmix_data$DOC_.ug.g.1.~Pmix_data$Plant_Type)) # significant
TukeyHSD(aov(Pmix_data$DOC_.ug.g.1.~Pmix_data$Plant_Type)) # only with no plant


#### Preliminary analysis ####

# Roots
summary(aov(Data_graphs$Roots ~ Data_graphs$Treatment + Data_graphs$Plant_Type))
TukeyHSD(aov(Data_graphs$Roots ~ Data_graphs$Treatment + Data_graphs$Plant_Type))

# Shoots
summary(aov(Data_graphs$Shoots ~ Data_graphs$Treatment + Data_graphs$Plant_Type))
TukeyHSD(aov(Data_graphs$Shoots ~ Data_graphs$Treatment + Data_graphs$Plant_Type))

## Log transformed 

# Roots 
summary(aov(Data_tba$log_root ~ Data_tba$Treatment + Data_tba$Plant_Type))
TukeyHSD(aov(Data_tba$log_root ~ Data_tba$Treatment + Data_tba$Plant_Type))

## Paired tests

RootsTrt <- Data_graphs %>%
  filter(Data_graphs$Treatment == "Pi  optimal nutrition")
summary(aov(RootsLP$Roots~ RootsLP$Plant_Type))

RootsPlt <- Data_graphs %>%
  filter(Data_graphs$Plant_Type == "36ptxD")
summary(aov(RootsPlt$Roots~ RootsPlt$Treatment))
TukeyHSD(aov(RootsPlt$Roots~ RootsPlt$Treatment))

