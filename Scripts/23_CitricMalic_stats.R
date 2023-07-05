#### Root exudate analysis ####

### Libraries ###

library(dplyr)
library(tidyverse)
library(car)


#### Data pre-processing ####

## Upload data
Exud_CM <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/main_rootexud_weightcorr.csv")

## Filter for analysis
CM_norm <-  Exud_CM %>%
  select(Exud_code,X36C1_conc:WTPM3_conc)

## Transpose 
CM_norm <- as.data.frame(t(CM_norm))

## Pre-process

colnames(CM_norm) <- CM_norm[1,] 
CM_norm <- CM_norm%>%
  filter(!row_number() %in% c(1))

## Change all to numeric variables 
CM_norm[,1:7] <- apply(CM_norm[,1:7],2,            # Specify own function within apply
                            function(x) as.numeric(as.character(x)))

# Check mode
sapply(CM_norm,mode)

## Add plant type and treatment
CM_norm <- CM_norm %>%
  add_column(Citrate_sum = CM_norm$E1+CM_norm$E2+CM_norm$E3) %>%
  add_column(Malate_sum = CM_norm$E4+CM_norm$E5+CM_norm$E7)%>% # do not add Malic acid E6
  add_column(Plant_type = c(rep("ptxD36",12),
                            rep("No Plant",4),
                            rep("Wild Type",9))) %>%
  add_column(Treatment = c(rep("Control",3), rep("Low P",3),
                           rep("Phi",3),rep("P Mix",3),
                           "Control", "Low P","Phi","P Mix",
                           rep("Control",3), rep("Low P",3),
                           rep("P Mix",3)))
  

#### Data normalization #### 

# Histogram check (repeated with log transformed too)
for (i in 1:9) { 
  hist(CM_norm[,i],
       main = i)
} 

# Shapiro (repeated with log transformed too)
for(i in 1:ncol(CM_norm)){
  shapiro <- shapiro.test(CM_norm[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} ## None of them are normal

## Log transform

CM_log <- CM_norm

for(i in 1:9){
  CM_log[,i] <- abs(log(CM_log[,i]))
} 

# Shapiro for log transformed that worked
shapiro.test(CM_log$E2)
shapiro.test(CM_log$E3)
shapiro.test(CM_log$Malate_sum)

## Scale transformation did not worked

## Square root transformation:  Y’ = √Y+0.05

CM_sqrt <- CM_norm

for(i in 1:9){
  CM_sqrt[,i] <- sqrt(CM_sqrt[,i]) + 0.05
}

## Shapiro test for square transformation

shapiro.test(CM_sqrt$E1)
shapiro.test(CM_sqrt$E6)
shapiro.test(CM_sqrt$E7)
shapiro.test(CM_sqrt$Citrate_sum)

## Non- parametric data

shapiro.test(CM_norm$E4)


### Saving processed tables ###

CM_log <- CM_log[,c(2,3,9)] # separate data for log transformation

# Transformed data table
Exud_stats <- CM_norm[,4:5]

Exud_stats <- Exud_stats %>%
  add_column(E2_log = CM_log$E2) %>%
  add_column(E3_log = CM_log$E3) %>%
  add_column(E1_sqrt = CM_sqrt$E1)%>%
  add_column(E6_sqrt = CM_sqrt$E6)%>%
  add_column(E7_sqrt = CM_sqrt$E7)%>%
  add_column(Citrate_scum = CM_sqrt$Citrate_sum)%>%
  add_column(Malate_logsum = CM_log$Malate_sum)%>%
  add_column(Plant_type = CM_norm$Plant_type)%>%
  add_column(Treatment = CM_norm$Treatment)
  
# Save Processed table 
write.csv(CM_norm, "~/Grad_School/Maestria/Processed_data/Exudates/Citmal_OG.csv")
write.csv(Exud_stats, "~/Grad_School/Maestria/Processed_data/Exudates/Citmal_stats.csv")



#### Statistical Analysis ####

Exud_stats <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Citmal_stats.csv")
row.names(Exud_stats) <- Exud_stats$X
Exud_stats <- Exud_stats[,2:8]
sapply(Exud_stats,class)

### Normalized data ###

## Levene test: all passed 

for(i in 1:9){
  Lev_ex <- leveneTest(Exud_stats[,i] ~ Plant_type * Treatment, 
                       data = Exud_stats) 
  levene <- ifelse(Lev_ex[["Pr(>F)"]]>0.05, "YES","NO")
  print(c(i,levene))
}

## ANOVA ##

Anova(aov(Malate_logsum~ Plant_type * Treatment, 
          data = Exud_stats))

# Tukey HSD 
CM_Tukey <- TukeyHSD(aov(Malate_logsum~ Plant_type * Treatment, 
                           data = Exud_stats))
CM_Tukey[["Plant_type"]]


### Non-parametric data ###

sapply(Exud_stats, class)

# Kruskal Wallis test repeated for each exudate
kruskal.test(E5~ Plant_type, data = Exud_stats)
kruskal.test(E5~ Treatment, data = Exud_stats)


## Wilcox test 

pairwise.wilcox.test(Exud_stats$E4, Exud_stats$Treatment,
                     p.adjust.method = "BH")

  
