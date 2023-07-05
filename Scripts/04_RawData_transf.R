#### Adding processed data columns to csv ####

# Start creating new tables to avoid transforming my original raw files

### Libraries

library(dplyr)
library(tidyverse)
library(ggplot2)

### Import data

Raw_file <- read.csv("~/Grad_School/Maestria/Nt_RD_analysis/Raw_data_mastersheet.csv")


### Adding root to shoot ratio 

RtSt_file <- Raw_file %>%
  add_column(Root_shoot_ratio = Raw_file$Roots/Raw_file$Shoots)

hist(RtSt_file$Root_shoot_ratio) # root to shoot ratio also needs transformation
hist(log(RtSt_file$Root_shoot_ratio))
shapiro.test(log(RtSt_file$Root_shoot_ratio))

## Saving it to a new file

write.csv(RtSt_file, "~/Grad_School/Maestria/Processed_data/RSratio_table.csv")

### Log transformed values for shoots and roots

Biomass_transf <- Raw_file %>%
  dplyr::select(Sample_ID:Date) %>%
  add_column(log_root = abs(log(Raw_file$Roots))) %>%
  add_column(log_shoot = abs(log(Raw_file$Shoots))) %>%
  add_column(log_rs_ratio = abs(log(RtSt_file$Root_shoot_ratio)))%>%
  add_column(log_TOC = abs(log(Raw_file$TOC)))%>%
  add_column(log_MBC =(log(Raw_file$MBC+1)))

## Check normality with Shapiro test in log transformed

shapiro.test(Biomass_transf$log_root)
shapiro.test(Biomass_transf$log_shoot)
shapiro.test(Biomass_transf$log_rs_ratio) # All of them are normal
shapiro.test(Biomass_transf$log_TOC)
shapiro.test(Biomass_transf$log_MBC)
shapiro.test(Raw_file$MBC)

## Save into a new csv

write.csv(Biomass_transf, "~/Grad_School/Maestria/Processed_data/Biomass_logtransf.csv")
