#### Normality ####

# Description: This script will use normality tests on all my sets of data 
# to transform them if necessary for analysis

### Libraries

library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)

### Import data

Raw_data <- read.csv("~/Grad_School/Maestria/Nt_RD_analysis/Raw_data_mastersheet.csv")

### Histograms: need to show normal distribution

# Roots and shoots 
par(mfrow=c(1,2)) # data skewed to the left in both
hist(Raw_data$Roots) 
hist(Raw_data$Shoots) 

# TOC and MBC are normally distributed

par(mfrow=c(1,2))
hist(Raw_data$TOC) 
hist(Raw_data$MBC)
hist(log(Raw_data$MBC+1))
qqnorm(Raw_data$MBC)

# Root morphology
par(mfcol=c(2,2))
hist(Raw_data$root_transf_diam_mean) # yes
hist(Raw_data$root_transf_length) # maybe
hist(Raw_data$root_transf_area) # bimodal?
hist(Raw_data$root_transf_tip_count) # yes

### QQnorm: 1 to 1 slope

# Roots and shoots

par(mfrow=c(1,2)) # data skewed to the left in both
qqnorm(Raw_data$Roots) 
qqnorm(Raw_data$Shoots) 

# Root morphology
par(mfcol=c(2,2))
qqnorm(Raw_data$root_transf_diam_mean) # yes
qqnorm(Raw_data$root_transf_length) # yes 
qqnorm(Raw_data$root_transf_area) # maybe? small jump
qqnorm(Raw_data$root_transf_tip_count) # yes

### Shapiro test: p-value is bigger than 0.05 data are normally distributed

# Roots and shoots

shapiro.test(Raw_data$Roots) # no
shapiro.test(Raw_data$Shoots) # no

# TOC and MBC are normally distributed

shapiro.test(Raw_data$TOC) # yes
shapiro.test(Raw_data$MBC) # no

# Root morphology
shapiro.test(Raw_data$root_transf_diam_mean) # yes
shapiro.test(Raw_data$root_transf_length) # yes
shapiro.test(Raw_data$root_transf_area) # yes
shapiro.test(Raw_data$root_transf_tip_count) # yes

## Transformation roots and shoots

shootstr <- abs(log(Raw_data$Shoots))
hist(shootstr)
shapiro.test(shootstr)

rootstr <- abs(log(Raw_data$Roots))
hist(rootstr)
shapiro.test(rootstr)

## Pi roots
hist(Raw_data$root_Pi)
rpi <- abs(log(Raw_data$root_Pi))
shapiro.test(rpi)
shapiro.test(Raw_data$root_Pi)
