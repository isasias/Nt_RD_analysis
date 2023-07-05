#### Two-way ANOVA ####

#### Libraries ####

library(dplyr)
library(tidyverse)
library(car)

#### Import data ####

Data_OG <- read.csv("~/Grad_School/Maestria/Processed_data/RSratio_table.csv")
Data_log <- read.csv("~/Grad_School/Maestria/Processed_data/Biomass_logtransf.csv")

#### Checking data ####


## Compute mean and SD by group ##

# code 1
group_by(Data_plants, Plant_Type, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(log_root, na.rm = TRUE),
    sd = sd(log_root, na.rm = TRUE))

# or code 2
model.tables(root_anovaS, type="means", se = TRUE)

#### Plant Biomass ####

## Create table removing no plant treatment for analysis

Data_plants <- Data_log %>%
  filter(!Plant_Type == "No plant") %>% # remove row
  dplyr::select(-"X") # remove column

table(Data_plants$Plant_Type, Data_plants$Treatment)

### Roots ###

## ANOVA ##

# Type 1 and 2 drop the same results so it seems balanced
root_anova <- aov(log_root ~ Plant_Type * Treatment, 
                   data = Data_plants) # synergistic effect

### How to save R objects for next sessions ###
save(root_anova,file="savetest.RData") # it worked

root_anova <- "savetest.RData"
load(root_anova,envir = parent.frame(), verbose = FALSE)

# Test homogeneity of variance

plot(root_anova, 1) # check outliers
leveneTest(log_root ~ Plant_Type * Treatment, 
           data = Data_plants) # higher than 0.05 so it meets ANOVA assumption

# Save table
Anova(root_anova) # Treatments and interactions are significant 
write.csv(Anova(root_anova), "~/Grad_School/Maestria/Processed_data/Root_Anova.csv")

## Tukey HSD ##

rt_Tukey <- TukeyHSD(root_anova)

# Save tables 
write.csv(rt_Tukey$Plant_Type, "~/Grad_School/Maestria/Processed_data/rt_Tukey_Plant.csv")
write.csv(rt_Tukey$Treatment, "~/Grad_School/Maestria/Processed_data/rt_Tukey_Trt.csv")
write.csv(rt_Tukey$`Plant_Type:Treatment`, "~/Grad_School/Maestria/Processed_data/rt_Tukey_Ptrt.csv")

### Shoots ###

## ANOVA ##

# Type 1 and 2 drop the same results so it seems balanced
shoot_anova <- aov(log_shoot ~ Plant_Type * Treatment, 
                  data = Data_plants) # synergistic effect

# Test homogeneity of variance

plot(shoot_anova, 1) # check outliers
leveneTest(log_shoot ~ Plant_Type * Treatment, 
           data = Data_plants) # higher than 0.05 so it meets ANOVA assumption

# Save table
Anova(shoot_anova) # Treatments and interactions are significant 
write.csv(Anova(shoot_anova), "~/Grad_School/Maestria/Processed_data/Shoot_Anova.csv")

## Tukey HSD ##

sh_Tukey <- TukeyHSD(shoot_anova)

# Save tables 
write.csv(sh_Tukey$Plant_Type, "~/Grad_School/Maestria/Processed_data/sh_Tukey_Plant.csv")
write.csv(sh_Tukey$Treatment, "~/Grad_School/Maestria/Processed_data/sh_Tukey_Trt.csv")
write.csv(sh_Tukey$`Plant_Type:Treatment`, "~/Grad_School/Maestria/Processed_data/sh_Tukey_Ptrt.csv")

#### TOC and MBC ####

### TOC ###

## ANOVA ##

# Type 1 and 2 drop the same results so it seems balanced
TOC_anova <- aov(log_TOC ~ Plant_Type * Treatment, 
                   data = Data_log) # synergistic effect

# Test homogeneity of variance

plot(TOC_anova, 1) # check outliers
leveneTest(log_TOC ~ Plant_Type * Treatment, 
           data = Data_log) # higher than 0.05 so it meets ANOVA assumption

Anova(TOC_anova) # Treatments and interactions are not significant 

### MBC ###

## ANOVA ##

# Type 1 and 2 drop the same results so it seems balanced
MBC_anova <- aov(log_MBC ~ Plant_Type * Treatment, 
                 data = Data_log) # synergistic effect

# Test homogeneity of variance

plot(MBC_anova, 1) # check outliers
leveneTest(log_MBC ~ Plant_Type * Treatment, 
           data = Data_log) # higher than 0.05 so it meets ANOVA assumption

Anova(MBC_anova) # Treatments and interactions are not significant 

#### Root Architecture ####

# Remove no plant

Data_archi <- Data_OG %>% 
  filter(!Plant_Type == "No plant") %>% # remove row
  filter(!Sample_ID== "CC")

Data_archi <- Data_archi[1:42,]

### Root Length ###

## ANOVA ##

# Type 1 and 2 drop the same results so it seems balanced
rlength_anova <- aov(root_length ~ Plant_Type * Treatment, 
                   data = Data_archi) # synergistic effect

# Test homogeneity of variance

plot(rlength_anova, 1) # check outliers
leveneTest(root_length ~ Plant_Type * Treatment, 
           data = Data_archi) # higher than 0.05 so it meets ANOVA assumption

# Save table
Anova(rlength_anova) # Treatments and interactions are significant 
write.csv(Anova(rlength_anova), "~/Grad_School/Maestria/Processed_data/Rlength_Anova.csv")

## Tukey HSD ##

rlength_Tukey <- TukeyHSD(rlength_anova)

# Save tables 
write.csv(rlength_Tukey$Plant_Type, "~/Grad_School/Maestria/Processed_data/rlength_Tukey_Plant.csv")
write.csv(rlength_Tukey$Treatment, "~/Grad_School/Maestria/Processed_data/rlength_Tukey_Trt.csv")
write.csv(rlength_Tukey$`Plant_Type:Treatment`, "~/Grad_School/Maestria/Processed_data/rlength_Tukey_Ptrt.csv")

### Root Area ###

## ANOVA ##

# Type 1 and 2 drop the same results so it seems balanced
rarea_anova <- aov(root_area ~ Plant_Type * Treatment, 
                     data = Data_archi) # synergistic effect

# Test homogeneity of variance

plot(rarea_anova, 1) # check outliers
leveneTest(root_area ~ Plant_Type * Treatment, 
           data = Data_archi) # higher than 0.05 so it meets ANOVA assumption

# Save table
Anova(rarea_anova) # Treatments and interactions are significant 
write.csv(Anova(rarea_anova), "~/Grad_School/Maestria/Processed_data/Rarea_Anova.csv")

## Tukey HSD ##

rarea_Tukey <- TukeyHSD(rarea_anova)

# Save tables 
write.csv(rarea_Tukey$Plant_Type, "~/Grad_School/Maestria/Processed_data/rarea_Tukey_Plant.csv")
write.csv(rarea_Tukey$Treatment, "~/Grad_School/Maestria/Processed_data/rarea_Tukey_Trt.csv")
write.csv(rarea_Tukey$`Plant_Type:Treatment`, "~/Grad_School/Maestria/Processed_data/rarea_Tukey_Ptrt.csv")

#### Root Pi ####

bpi_anova <- aov(biomass_Pi ~ Plant_Type * Treatment, 
                   data = Data_archi) # synergistic effect

# Test homogeneity of variance

plot(bpi_anova, 1) # check outliers
leveneTest(biomass_Pi ~ Plant_Type * Treatment, 
           data = Data_archi) # higher than 0.05 so it meets ANOVA assumption

# Save table
Anova(bpi_anova) # Treatments and interactions are significant 
write.csv(Anova(bpi_anova), "~/Grad_School/Maestria/Processed_data/Bpi_Anova.csv")

## Tukey HSD ##

bpi_Tukey <- TukeyHSD(bpi_anova)
write.csv(bpi_Tukey$Treatment, "~/Grad_School/Maestria/Processed_data/bpi_Tukey_Trt.csv")
write.csv(bpi_Tukey$`Plant_Type:Treatment`, "~/Grad_School/Maestria/Processed_data/bpi_Tukey_Ptrt.csv")



