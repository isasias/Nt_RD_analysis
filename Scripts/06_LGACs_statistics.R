#### Statistical analysis LGACs 2022 ####

### Libraries

library(dplyr)
library(tidyverse)
library(agricolae)

### Import data

Data_graphs <- read.csv("~/Grad_School/Maestria/Processed_data/RSratio_table2.csv")
Data_tba <- read.csv("~/Grad_School/Maestria/Processed_data/Biomass_weights_transf.csv")

### Plant biomass

## OG data

rtmodel <- aov(data = Data_graphs,Roots ~ Treatment + Plant_Type)
summary(rtmodel)
rtout <- HSD.test(rtmodel, trt = c("Plant_Type", "Treatment"), alpha = 0.05)

stmodel <- aov(data = Data_graphs,Shoots ~ Treatment + Plant_Type)
summary(stmodel)
stout <- HSD.test(stmodel, trt = c("Plant_Type", "Treatment"), alpha = 0.05)
stout

#log_transformed

rtmodel_log <- aov(data = Data_tba, log_root ~ Treatment + Plant_Type)
summary(rtmodel_log)
rtoutl <- HSD.test(rtmodel_log, trt = c("Plant_Type", "Treatment"), alpha = 0.05)
rtoutl

stmodel_log <- aov(data = Data_tba, log_shoot ~ Treatment + Plant_Type)
summary(stmodel_log)
stoutl <- HSD.test(stmodel_log, trt = c("Plant_Type", "Treatment"), alpha = 0.05)
stoutl

### Root architecture

areamodel <- aov(data = Data_graphs,root_transf_area ~ Treatment + Plant_Type)
summary(areamodel)
area_out <- HSD.test(areamodel, trt = c("Plant_Type", "Treatment"), alpha = 0.05)
area_out

### TOC & MBC

TOCmodel <- aov(data = TOCgraph,TOC ~ Treatment + Plant_Type)
summary(TOCmodel)
TOCout <- HSD.test(rtmodel, trt = c("Plant_Type", "Treatment"), alpha = 0.05)
TOCout

MBCmodel <- aov(data = TOCgraph, MBC ~ Treatment + Plant_Type)
summary(MBCmodel)
MBCout <- HSD.test(rtmodel, trt = c("Plant_Type", "Treatment"), alpha = 0.05)
MBCout
