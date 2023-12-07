#### Fungi Functional Statistics ####

### Libraries ###

library(dplyr)
library(tidyverse)
library(car)

#### Data Upload ####

Func_Fungi <- read.csv("~/Grad_School/Maestria/Processed_data/Fungal traits/Functional_Fungi.csv")
Simp_fungi <- read.csv("~/Grad_School/Maestria/Processed_data/Fungal traits/Simplified_Fungi.csv")
metadata <- read.csv("~/Grad_School/Maestria/Processed_data/metadata.csv")

#### Simplified Categories ####

### Data Manipulation ###

Simp_fungi <- as.data.frame(t(Simp_fungi))
colnames(Simp_fungi) <- Simp_fungi[1,]
Simp_fungi <- Simp_fungi[-1,]

Simp_fungi <- cbind(metadata,Simp_fungi)
Simp_fungi <- Simp_fungi[,-1]

# Remove low abundance 
Simp_fungi <- Simp_fungi[,-c(22:25)]

# Change to numeric type
for(i in 5:ncol(Simp_fungi)){
  Simp_fungi[,i] <- as.numeric(Simp_fungi[,i])
}

### Test for Normality ###

## Shapiro test
for(i in 5:ncol(Simp_fungi)){
  shapiro <- shapiro.test(Simp_fungi[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}
# Plant pathogen and Animal parasite are the only normal (5 & 11)

## Histograms 
for (i in 5:ncol(Simp_fungi)) { 
  hist(Simp_fungi[,i],
       main = i)
} # Wood saprotroph almost normal 10

### Log transform data ###

Norm_SF <- Simp_fungi
for(i in 5:ncol(Norm_SF)){
    Norm_SF[,i] <- abs(log10(Norm_SF[,i]+1))
} 

## re-check normality

# Shapiro test
for(i in 5:ncol(Norm_SF)){
  shapiro <- shapiro.test(Norm_SF[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # 6,7,8,10,11,13,14,15,16,18,19,21,

# Histograms 
for (i in 5:ncol(Norm_SF)) { 
  hist(Norm_SF[,i],
       main = i)
} # almost normal distribution for 20, 25

## 9,12,17,22,23 not normal

### Final table ###

Simp_FF <- Simp_fungi[,c(1:5)]
Norm_FF <- Norm_SF[,c(6,7,8,10,11,13,14,15,16,18,19,21,20,25)]

Simp_FF <- cbind(Simp_FF,Norm_FF)

Simp_FF <- Simp_FF[,-c(1:4)]

colnames(Simp_FF) <- c(paste0("Cat_",as.character(1:15)))
Simp_FF <- cbind(Simp_fungi[,c(1:4)],Simp_FF)


### ANOVA ###

## Levene test

leveneTest(Cat_15 ~ Plant_Type * Treatment, 
           data = Simp_FF) # all passed

## ANOVA
Anova(aov(Cat_15 ~ Plant_Type * Treatment, 
          data = Simp_FF))
# only Cat 2: Soil Saprotroph marginal diff

## Tukey HSD Soil saprotrophs
Fun_Tukey <- TukeyHSD(aov(Cat_2 ~ Plant_Type * Treatment, 
                          data = Simp_FF))

# ptxDNt-36-No plant   0.04758828
# Wild Type-No plant   0.09616748

### Kruskal Wallis ###

# Rename columns
colnames(Simp_fungi)[c(5:21)]<- c(paste0("Cat_",as.character(1:17)))

# 9,12,17 not normal
# Dung, Mycoparasite, Nectar/Tap

##  per Plant Type

# Dung saprotroph: significant p-value = 0.02787
kruskal.test(Cat_5 ~ Plant_Type, data = Simp_fungi)

# Wilcox 
pairwise.wilcox.test(Simp_fungi$Dung_sapro, Simp_fungi$Plant_Type,
                     p.adjust.method ="none")

# ptxD significantly dif from WT and BS
# BS 0.030
# WT 0.028 

# Mycoparasite: not significant
kruskal.test(Cat_8 ~ Plant_Type, data = Simp_fungi)

# Nectar/Tap: not significant 
kruskal.test(Cat_13 ~ Plant_Type, data = Simp_fungi)

## per treatment: none of them significant 

kruskal.test(Cat_5 ~ Treatment, data = Simp_fungi) # NS
kruskal.test(Cat_8 ~ Treatment, data = Simp_fungi) # NS
kruskal.test(Cat_13 ~ Treatment, data = Simp_fungi) # NS

#### All root associated categories ####

Func_Fungi <- as.data.frame(t(Func_Fungi))
colnames(Func_Fungi) <- Func_Fungi[1,]
Func_Fungi <- Func_Fungi[-1,]

## Change to numeric type

for(i in 1:ncol(Func_Fungi)){
  Func_Fungi[,i] <- as.numeric(Func_Fungi[,i])
}

## Sum all root associated categories 

Root_sums <- rowSums(Func_Fungi[,c(11,12,22,36,40,48,50,53,55,
                                   63,66,73,77,81,86,88,94)])
Root_sums <- cbind(metadata,Root_sums)

### ANOVA ###

# Shapiro
shapiro.test(Root_sums$Root_sums) # normal

# Levene test 

leveneTest(Root_sums ~ Plant_Type * Treatment, 
           data = Root_sums) # passed

# ANOVA 

Anova(aov(Root_sums ~ Plant_Type * Treatment, 
          data = Root_sums))
# Plant_Type 0.01341 *

# Tukey

Fun_Tukey <- TukeyHSD(aov(Root_sums ~ Plant_Type * Treatment, 
                          data = Root_sums))

# Wild Type:Pi-ptxDNt-36:Pi 0.01869647
# Pi/Phi mix-Phi   0.05673901
# Wild Type-No plant  0.05288917

## Remove Bulk Soil for analysis ##

Root_cat <- Root_sums[c(1:21),]

# Shapiro
shapiro.test(Root_cat$Root_sums) # normal

# Levene test 

leveneTest(Root_sums ~ Plant_Type * Treatment, 
           data = Root_cat) # passed

## ANOVA 

Anova(aov(Root_sums ~ Plant_Type * Treatment, 
          data = Root_cat))

# Plant_Type 0.009793 **
# Treatment 0.057376 . 
# Plant_Type:Treatment 0.014647 *

## Tukey

Fun_Tukey <- TukeyHSD(aov(Root_sums ~ Plant_Type * Treatment, 
                          data = Root_cat))

# Pi/Phi mix-Phi 0.04824884

#### Paper stats ####

# normal with log
shapiro.test(log(Simp_fungi$`Arbuscular Mycorrhizal`))
hist(log(Simp_fungi$`Arbuscular Mycorrhizal`))
qqnorm(log(Simp_fungi$`Arbuscular Mycorrhizal`))

# normal with sqrt
shapiro.test(sqrt(Simp_fungi$Ectomycorrhiza))
hist(sqrt(Simp_fungi$Ectomycorrhiza))

# normal with log+1
Endophytes <- Simp_fungi$`Root Endophyte Dark Septate Soil Saprotroph`+Simp_fungi$`Root Endophyte`
shapiro.test(log(Endophytes+1))
hist(log(Endophytes+1))

#### Abundance stats ####

sum(Simp_fungi[,5:21])
colSums(Simp_fungi[,5:21]) # 58% plant pathogens

# saprotrophs 40%
