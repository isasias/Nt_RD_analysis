#### Root exudate ANOVA ####

### Libraries ###

library(dplyr)
library(tidyverse)
library(car)

# for non-parametric
library(ggpubr)
library(rstatix)

### Data upload ###

Exudates_norm <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Normalized_exudates.csv")
row.names(Exudates_norm) <- Exudates_norm[,1]
Exudates_norm <- Exudates_norm[,-1]

Exudates_unnorm <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/UnNorm_exudates.csv")
row.names(Exudates_unnorm) <- Exudates_unnorm[,1]
Exudates_unnorm <- Exudates_unnorm[,-1]

#### Anova check ####

# Levene Test: all passed 

for(i in 3:ncol(Exudates_norm)){
  Lev_ex <- leveneTest(Exudates_norm[,i] ~ Plant_type * Treatment, 
                 data = Exudates_norm) 
  levene <- ifelse(Lev_ex[["Pr(>F)"]]>0.05, "YES","NO")
  print(c(i,levene))
}

# Repeat this for each exudate 
# (cannot find out how to make a for loop)
Anova(aov(Exud_59~ Plant_type * Treatment, 
                  data = Exudates_norm))

### Tukey HSD

Exud_Tukey <- TukeyHSD(aov(Exud_59~ Plant_type * Treatment, 
          data = Exudates_norm))
Exud_Tukey[["Treatment"]]

# Save significant files 
write.csv(Exud_Tukey$`Plant_type:Treatment`, "~/Grad_School/Maestria/Processed_data/Exud51_Tukey_Ptrt.csv")


#### Non-parametric data ####

for(i in 1:23){
  Lev_ex <- leveneTest(Exudates_unnorm[,i] ~ Plant_type * Treatment, 
                       data = Exudates_unnorm) 
  levene <- ifelse(Lev_ex[["Pr(>F)"]]>0.05, "YES","NO")
  print(c(i,levene))
} # all passed

# check variables type
sapply(Exudates_unnorm, class)


# Kruskal Wallis test repeated for each exudate
kruskal.test(Exud_48~ Plant_type, data = Exudates_unnorm)
kruskal.test(Exud_48~ Treatment, data = Exudates_unnorm)


## Wilcox test 

pairwise.wilcox.test(Exudates_unnorm$Exud_30, Exudates_unnorm$Treatment,
                     p.adjust.method ="none")
pairwise.wilcox.test(Exudates_unnorm$Exud_8, Exudates_unnorm$Plant_type,
                     p.adjust.method ="none")
