#### Final Exudate Stats ####

### Libraries ###

library(dplyr)
library(tidyverse)
library(car)

#### Data upload and pre-processing ####

## Metadata ##
metadata <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Exud_metadata.csv")

## Exud sums ##

Exud_sums <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Normalized_exudates.csv")
row.names(Exud_sums) <- Exud_sums[,1]
Exud_sums <- Exud_sums[,-1]

## Organic Acids ##

Org_acids <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Normalized_acids.csv")
row.names(Org_acids) <- Org_acids[,1]
Org_acids <- Org_acids[,-1]

## Non-parametric data ##

Nonpar_exud <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Untransform_exudates.csv")
row.names(Nonpar_exud) <- Nonpar_exud[,1]
Nonpar_exud <- Nonpar_exud[,-1]
Nonpar_exud <- cbind(metadata,Nonpar_exud)


Nonpar_acids <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Untransform_acids.csv")
row.names(Nonpar_acids) <- Nonpar_acids[,1]
Nonpar_acids <- Nonpar_acids[,-1]
Nonpar_acids <- cbind(metadata,Nonpar_acids)


#### Levene test loop ####

for(i in 3:ncol(Exud_sums)){
  Lev_ex <- leveneTest(Exud_sums[,i] ~ Plant_type * Treatment, 
                       data = Exud_sums) 
  levene <- ifelse(Lev_ex[["Pr(>F)"]]>0.05, "YES","NO")
  print(c(i,levene))
}

for(i in 3:ncol(Org_acids)){
  Lev_ex <- leveneTest(Org_acids[,i] ~ Plant_type * Treatment, 
                       data = Org_acids) 
  levene <- ifelse(Lev_ex[["Pr(>F)"]]>0.05, "YES","NO")
  print(c(i,levene))
}

## All passed 

#### ANOVA using for loop ####

root_pvalues <- data.frame(Plant = rep(NA, 94),
                           Treatment = rep(NA,94),
                           Pl_Trt = rep(NA,94))  

for(i in 3:ncol(Exud_sums)){
  Exudate <- Exud_sums[,i] 
  AR <- Anova(aov(Exudate ~ Plant_type * Treatment, 
                       data = Exud_sums))
  j <- i-2
  root_pvalues$Plant[j] <- AR$`Pr(>F)`[1]
  root_pvalues$Treatment[j] <- AR$`Pr(>F)`[2]
  root_pvalues$Pl_Trt[j] <- AR$`Pr(>F)`[3]
}

row.names(root_pvalues) <- colnames(Exud_sums[3:96])

## Save anova table 

write.csv(root_pvalues, "~/Grad_School/Maestria/Processed_data/Exudates/Root_exud_ANOVAs.csv")

#### Tukey HSD ####

Exud_Tukey <- TukeyHSD(aov(Exud_sums[,3]~ Plant_type * Treatment, 
                           data = Exud_sums))
Exud_Tukey[["Treatment"]]

### For plant type

for(i in 3:ncol(Exud_sums)){
  ET <- TukeyHSD(aov(Exud_sums[,i]~ Plant_type * Treatment, 
                     data = Exud_sums))
  print(c(i,colnames(Exud_sums[i])))
  print(ET[["Plant_type"]])
}

# No SD (2 exudates): 3,67
# No Rhizosphere effect WT (3 Exudates): 54,69,70

# Rhizosphere effect (59):4,5,7:9,11:15,17:21,23,24,27,31:33,35,39,
# 42:47,50,,52,53,57,59:61,63,64,66,68,71,73:76,78,80:82,84:90,92,
# 95,96

# All SD (30 exudates): 6,10,16*,22,25,26,28:30,34,36,37*,38,40,41,
# 48,49,51,55,56,58,62,65,72,77,79*,83,91,93,94
# * means marginal effect

### For treatment

for(i in 3:ncol(Exud_sums)){
  ET <- TukeyHSD(aov(Exud_sums[,i]~ Plant_type * Treatment, 
                     data = Exud_sums))
  print(c(i,colnames(Exud_sums[i])))
  print(ET[["Treatment"]])
}

# No SD (55 significant; 8 marginal): 3*,4,5,7:11,13,14,17:19,20*,
# 23:25,27,30*,32,34,36*37,39,42:45,47*,50,52:54,56,57,59,60,62:67,
# 72,73,75,76,79,80,82,83,84*,85,86,87*,88,89,91:96

# Pi:LP (15 signif; 8 marginal): 6,16*,22,26,28,29,31,33,38,40,
# 41*,46,48*,49*,51*,55,58,61*,68,71,77,81*,90

# PM:LP (12 significant; 1 marginal): 12,15,16,21,29,31*,35,38,46,
# 69,70,74,81

# Phi:LP (3 signif; 4 marginal): 29,31*,38,46,61*,77*,81*

# PM:Phi 78*

# PM:Pi (5 signif; 2 marginal): 41*,48,49*,51,68,71,78


# Special 11 Exudates (more than 2 signif; *marginal):
# 16*,29,31*,38,46,48*,51*,68,71*,77*,78*

#### Plant:Treatment for loops ####

### For loops for PM comparisons

pl_trt <- as.data.frame(ET[["Plant_type:Treatment"]])

pl_trt$`p adj`[66]

for(i in 3:ncol(Exud_sums)){
  Exudate <- Exud_sums[,i] 
  ET <- TukeyHSD(aov(Exudate ~ Plant_type * Treatment, 
                  data = Exud_sums))
  pl_trt <- as.data.frame(ET[["Plant_type:Treatment"]])
  WTPM_Pi <- ifelse(pl_trt$`p adj`[63]<0.05, "YES","NO")
  WT36_PM <- ifelse(pl_trt$`p adj`[66]<0.05, "YES","NO")
  PM_Pi36 <- ifelse(pl_trt$`p adj`[66]<0.05, "YES","NO")
  print(c(i,WTPM_Pi,WT36_PM,PM_Pi36))
}
# ND (64):3,4,5,7:9,11:21,23,24,31,32,35:37,39,42:47,50,52:54,56,
# 57,59:61,63:67,69,70,73:76,79,81,82,85,86,88:90,92,94:96

# All (14): 22,25,27,28,38,41,48,49,51,71,78,83,84,87

# WT:36 PM & 36 PM:Pi(14):6,10,26,29,30,33,34,40,55,58,62,72,77,91,93

# WT PM:Pi (2): 68,80

### For loops for ptxd Pi Phi comparisons

pl_trt$`p adj`[14]

for(i in 3:ncol(Exud_sums)){
  Exudate <- Exud_sums[,i] 
  ET <- TukeyHSD(aov(Exudate ~ Plant_type * Treatment, 
                     data = Exud_sums))
  pl_trt <- as.data.frame(ET[["Plant_type:Treatment"]])
  PM_Phi36 <- ifelse(pl_trt$`p adj`[44]<0.05, "YES","NO")
  Pi_Phi36 <- ifelse(pl_trt$`p adj`[41]<0.05, "YES","NO")
  LP_Phi36 <- ifelse(pl_trt$`p adj`[41]<0.05, "YES","NO")
  print(c(i,PM_Phi36,Pi_Phi36, LP_Phi36))
} # no significant differences in any of them

### Low P comparisons

pl_trt$`p adj`[27]

for(i in 3:ncol(Exud_sums)){
  Exudate <- Exud_sums[,i] 
  ET <- TukeyHSD(aov(Exudate ~ Plant_type * Treatment, 
                     data = Exud_sums))
  pl_trt <- as.data.frame(ET[["Plant_type:Treatment"]])
  WT36LP <- ifelse(pl_trt$`p adj`[12]<0.05, "YES","NO")
  Pi_LP36 <- ifelse(pl_trt$`p adj`[17]<0.05, "YES","NO")
  Pi_LPWT <- ifelse(pl_trt$`p adj`[27]<0.05, "YES","NO")
  print(c(i,WT36LP,Pi_LP36,Pi_LPWT))
} 
# no significant differences between WT and ptxd in any exudate
# 36_Pi:LP: 29 Chorismate_C18P
# LP:Pi on both plants: 38 Phthalate_C18P

### Low P Phi PM comparisons

pl_trt$`p adj`[14]

for(i in 3:ncol(Exud_sums)){
  Exudate <- Exud_sums[,i] 
  ET <- TukeyHSD(aov(Exudate ~ Plant_type * Treatment, 
                     data = Exud_sums))
  pl_trt <- as.data.frame(ET[["Plant_type:Treatment"]])
  LP_Phi36 <- ifelse(pl_trt$`p adj`[14]<0.05, "YES","NO")
  LP_PM36 <- ifelse(pl_trt$`p adj`[20]<0.05, "YES","NO")
  print(c(i,LP_Phi36, LP_PM36))
} 
# 29 and 38 significantly different on both
# 36_LP:PM: 12,15,22,30


#### Plant:Treatment particular cases ####

## Top Exudates 

Top_names <- colnames(Exud_sums[c(38,29,68,46,48,51,71,78,77,16,31,
                                22,25,27,28,41,49,83,84,87,30,12,
                                15,6,10,26,33,34,40,55,58,62,72,
                                91,93)]) # 35 total

# First 11 main ones 38 to 31
# Second 13: 22 to 87
# Third weird cases 30,12,15
# Last place 12: 6 to 93


### 38: Phthalate_C18P ###

ET <- TukeyHSD(aov(Exud_sums[,38]~ Plant_type * Treatment, 
                   data = Exud_sums))

# Save Tukey HSD
Tuk_Phthalate <- ET[["Plant_type:Treatment"]]
write.csv(Tuk_Phthalate, "~/Grad_School/Maestria/Processed_data/Exudates/Tuk_Phthalate.csv")


### 29: Chorismate_C18P_SQ ###

ET <- TukeyHSD(aov(Exud_sums[,29]~ Plant_type * Treatment, 
                   data = Exud_sums))

# Save Tukey HSD
Tuk_Choris <- ET[["Plant_type:Treatment"]]
write.csv(Tuk_Choris, "~/Grad_School/Maestria/Processed_data/Exudates/Tuk_Chorismate.csv")



## Maleic_C18N

ET <- TukeyHSD(aov(Maleic_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Maleic_C18N <- pl_trt[,4]

## 46: Acetate_HILN

ET <- TukeyHSD(aov(Exud_sums[,46]~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Acetate_HILN <- pl_trt[,4]

## 48: Nicotinate_C18N

ET <- TukeyHSD(aov(Exud_sums[,48]~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Nicotinate_C18N <- pl_trt[,4]

## 49: Malate_C18N

ET <- TukeyHSD(aov(Malate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Malate_C18N <- pl_trt[,4]

## Mugineate_C18N

ET <- TukeyHSD(aov(Mugineate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Mugineate_C18N <- pl_trt[,4]

## Valarate_C18N_L

ET <- TukeyHSD(aov(Valarate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Valarate_C18N <- pl_trt[,4]

## Tartarate_C18P

ET <- TukeyHSD(aov(Tartarate_C18P~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Tartarate_C18P <- pl_trt[,4]

## 16: Ex47 N,N'-bis[(S)-1-methoxycarbonylethyl]fumaric diamide

ET <- TukeyHSD(aov(Exud_sums[,16]~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Exud_47 <- pl_trt[,4]

## Ferulate_C18P RHIZOSPHERE EFFECT only

ET <- TukeyHSD(aov(Ferulate_C18P~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Ferulate_C18P <- pl_trt[,4]

## 22: Ex5 1,3,6-Octatriene

ET <- TukeyHSD(aov(Exud_sums[,22]~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Exud_5 <- pl_trt[,4]

## Maleic_C18P_SQ

ET <- TukeyHSD(aov(Maleic_C18P~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Maleic_C18P <- pl_trt[,4]

## Glutarate_C18N_SQ

ET <- TukeyHSD(aov(Glutarate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Glutarate_C18N <- pl_trt[,4]

## Asparate_C18N_SQ

ET <- TukeyHSD(aov(Asparate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Asparate_C18N <- pl_trt[,4]

## Butyrate_C18N

ET <- TukeyHSD(aov(Butyrate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Butyrate_C18N <- pl_trt[,4]

## Glutamate_C18N

ET <- TukeyHSD(aov(Glutamate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Glutamate_C18N <- pl_trt[,4]

## Shikimate_C18POS

ET <- TukeyHSD(aov(Shikimate_C18POS~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Shikimate_C18P <- pl_trt[,4]

## Coumarate_C18N

ET <- TukeyHSD(aov(Coumarate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Coumarate_C18N <- pl_trt[,4]

## 12: Ex27 D-(+)-Glucose

ET <- TukeyHSD(aov(Exud_sums[,12]~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Exud_27 <- pl_trt[,4]

## 15: Ex39 3_5-Dimethoxy-4-hydroxycinnamicacid

ET <- TukeyHSD(aov(Exud_sums[,15]~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Exud_39 <- pl_trt[,4]

## Glutarate_C18P_SQ

ET <- TukeyHSD(aov(Glutarate_C18P~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Glutarate_C18P <- pl_trt[,4]

## Chorismate_C18N_SQ

ET <- TukeyHSD(aov(Chorismate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Chorismate_C18N <- pl_trt[,4]

## Gallate_C18N

ET <- TukeyHSD(aov(Gallate_C18N~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Gallate_C18N <- pl_trt[,4]

## 55: Ex11 nicotinamide

ET <- TukeyHSD(aov(Exud_sums[,55]~ Plant_type * Treatment, 
                   data = Exud_sums))
pl_trt <- ET[["Plant_type:Treatment"]]
Exud_11 <- pl_trt[,4]

Tukey_pltrt <- cbind(Acetate_HILN,Asparate_C18N,Butyrate_C18N,
                     Chorismate_C18N,Coumarate_C18N,Exud_11,
                     Exud_27,Exud_39,Exud_47,Exud_5,Ferulate_C18P,
                     Gallate_C18N,Glutamate_C18N,Glutarate_C18N,
                     Glutarate_C18P,Malate_C18N,Maleic_C18N,
                     Mugineate_C18N,Nicotinate_C18N,Shikimate_C18P,
                     Tartarate_C18P,Valarate_C18N)

row.names(Tukey_pltrt) <- row.names(pl_trt)
write.csv(Tukey_pltrt, "~/Grad_School/Maestria/Processed_data/Exudates/Tukey_topex.csv")






#### Organic Acids ####

### Levene test loop ###

for(i in 3:ncol(Org_acids)){
  Lev_ex <- leveneTest(Org_acids[,i] ~ Plant_type * Treatment, 
                       data = Org_acids) 
  levene <- ifelse(Lev_ex[["Pr(>F)"]]>0.05, "YES","NO")
  print(c(i,levene))
}

### Anova for loop ###

acid_pvalues <- data.frame(Plant = rep(NA, 25),
                           Treatment = rep(NA,25),
                           Pl_Trt = rep(NA,25))  

for(i in 3:ncol(Org_acids)){
  Exudate <- Org_acids[,i] 
  AR <- Anova(aov(Exudate ~ Plant_type * Treatment, 
                  data = Org_acids))
  j <- i-2
  acid_pvalues$Plant[j] <- AR$`Pr(>F)`[1]
  acid_pvalues$Treatment[j] <- AR$`Pr(>F)`[2]
  acid_pvalues$Pl_Trt[j] <- AR$`Pr(>F)`[3]
}

row.names(acid_pvalues) <- colnames(Org_acids[3:27])

acid_pvalues <- acid_pvalues[-c(3,4),] 

Org_acids <- Org_acids[,c(row.names(acid_pvalues))]
Org_acids <- cbind(metadata,Org_acids)

#### Tukey HSD ####

### For plant type

for(i in 3:ncol(Org_acids)){
  ET <- TukeyHSD(aov(Org_acids[,i]~ Plant_type * Treatment, 
                     data = Org_acids))
  print(c(i,colnames(Org_acids[i])))
  print(ET[["Plant_type"]])
}

# No Rhizosphere effect: 5*

# Rhizosphere effect (16): 6:9,11:13,15,16,18:21,23:25

# All SD (6): 3,4,10,14,17,22

### For treatment

for(i in 3:ncol(Org_acids)){
  ET <- TukeyHSD(aov(Org_acids[,i]~ Plant_type * Treatment, 
                     data = Org_acids))
  print(c(i,colnames(Org_acids[i])))
  print(ET[["Treatment"]])
}

#### Plant:Treatment for loops ####

### For loops for PM comparisons

for(i in 3:ncol(Org_acids)){
  Exudate <- Org_acids[,i] 
  ET <- TukeyHSD(aov(Exudate ~ Plant_type * Treatment, 
                     data = Org_acids))
  pl_trt <- as.data.frame(ET[["Plant_type:Treatment"]])
  WTPM_Pi <- ifelse(pl_trt$`p adj`[63]<0.05, "YES","NO")
  WT36_PM <- ifelse(pl_trt$`p adj`[66]<0.05, "YES","NO")
  PM_Pi36 <- ifelse(pl_trt$`p adj`[66]<0.05, "YES","NO")
  print(c(i,WTPM_Pi,WT36_PM,PM_Pi36))
}
# Special cases 14,25
# WTPM:Pi 8,14,25
# PMWT:36 3,10,14,17,22,25
# 36PM:Pi 3,10,14,17,22,25

### For loops for ptxd Pi Phi comparisons

pl_trt$`p adj`[14]

for(i in 3:ncol(Org_acids)){
  Exudate <- Org_acids[,i] 
  ET <- TukeyHSD(aov(Exudate ~ Plant_type * Treatment, 
                     data = Org_acids))
  pl_trt <- as.data.frame(ET[["Plant_type:Treatment"]])
  PM_Phi36 <- ifelse(pl_trt$`p adj`[44]<0.05, "YES","NO")
  Pi_Phi36 <- ifelse(pl_trt$`p adj`[41]<0.05, "YES","NO")
  LP_Phi36 <- ifelse(pl_trt$`p adj`[41]<0.05, "YES","NO")
  print(c(i,PM_Phi36,Pi_Phi36, LP_Phi36))
} # no significant differences in any of them

### Low P comparisons

pl_trt$`p adj`[27]

for(i in 3:ncol(Org_acids)){
  Exudate <- Org_acids[,i] 
  ET <- TukeyHSD(aov(Exudate ~ Plant_type * Treatment, 
                     data = Org_acids))
  pl_trt <- as.data.frame(ET[["Plant_type:Treatment"]])
  WT36LP <- ifelse(pl_trt$`p adj`[12]<0.05, "YES","NO")
  Pi_LP36 <- ifelse(pl_trt$`p adj`[17]<0.05, "YES","NO")
  Pi_LPWT <- ifelse(pl_trt$`p adj`[27]<0.05, "YES","NO")
  print(c(i,WT36LP,Pi_LP36,Pi_LPWT))
} 
# no significant differences in any of them NO org acids

### Low P Phi PM comparisons

pl_trt$`p adj`[14]

for(i in 3:ncol(Org_acids)){
  Exudate <- Org_acids[,i] 
  ET <- TukeyHSD(aov(Exudate ~ Plant_type * Treatment, 
                     data = Org_acids))
  pl_trt <- as.data.frame(ET[["Plant_type:Treatment"]])
  LP_Phi36 <- ifelse(pl_trt$`p adj`[14]<0.05, "YES","NO")
  LP_PM36 <- ifelse(pl_trt$`p adj`[20]<0.05, "YES","NO")
  print(c(i,LP_Phi36, LP_PM36))
} 
# 5 LP:PM

#### Plant:Treatment particular cases ####

## 4: Benzoic HILN

ET <- TukeyHSD(aov(Org_acids[,4]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 7 significant 5 marginal

## 14: Glutamic_C18P

ET <- TukeyHSD(aov(Org_acids[,14]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 10 significant 6 marginal

## 25: Butyric_C18N

ET <- TukeyHSD(aov(Org_acids[,25]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 8 significant

## 5: Ex39

ET <- TukeyHSD(aov(Org_acids[,5]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 10 significant 6 marginal

## 7: Phthalic_C18N

ET <- TukeyHSD(aov(Org_acids[,7]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 15 significant 5 marginal

## 8: Acetic_C18N

ET <- TukeyHSD(aov(Org_acids[,8]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 21 significant 

## 10: Gallic_C18N

ET <- TukeyHSD(aov(Org_acids[,10]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 8 significant 

## 17: Nicotinic_C18P

ET <- TukeyHSD(aov(Org_acids[,17]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 26 significant 2 marginal

## 21: Succinic_C18N

ET <- TukeyHSD(aov(Org_acids[,21]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 21 significant 3 marginal

## 22: Tartaric_C18P

ET <- TukeyHSD(aov(Org_acids[,22]~ Plant_type * Treatment, 
                   data = Org_acids))
print(ET[["Plant_type:Treatment"]])

# 27 significant 5 marginal

#### Non-parametric data ####

### Kruskal Wallis: Exudates

nonpar_pvalues <- data.frame(Plant = rep(NA, 12),
                           Treatment = rep(NA,12))

for(i in 3:ncol(Nonpar_exud)){
  KTp <- kruskal.test(Nonpar_exud[,i] ~ Plant_type, 
                      data = Nonpar_exud)
  KTt <- kruskal.test(Nonpar_exud[,i] ~ Treatment, 
                      data = Nonpar_exud)
  j <- i-2
  nonpar_pvalues$Plant[j] <- KTp[["p.value"]]
  nonpar_pvalues$Treatment[j] <- KTt[["p.value"]]
}

row.names(nonpar_pvalues) <- colnames(Nonpar_exud[3:14])

# Write csv
write.csv(nonpar_pvalues, "~/Grad_School/Maestria/Processed_data/Exudates/Root_exud_Kruskal.csv")


### Wilcox test 

## For plant type

for(i in 3:ncol(Nonpar_exud)){
  WT <- pairwise.wilcox.test(Nonpar_exud[,i], 
                             Nonpar_exud$Plant_type,
                             p.adjust.method ="none")
  print(c(i,colnames(Nonpar_exud[i])))
  print(WT[["p.value"]])
}
# NoRS:WT 3,5,6,7,8,9,11
# RS 13,14
# All different 10
# ND 4,12


## For treatment

for(i in 3:ncol(Nonpar_exud)){
  WT <- pairwise.wilcox.test(Nonpar_exud[,i], 
                             Nonpar_exud$Treatment,
                             p.adjust.method ="none")
  print(c(i,colnames(Nonpar_exud[i])))
  print(WT[["p.value"]])
}

# Pi:LP 3:5,6*,7:9,11*,12
# PM:LP 3:5,6*,7:9,11*,12
# No difference 10,13,14
# Phi:LP 4,7,12
# Pi:PM 4*,5*,7*,12*
# Special cases to graph 4 & 7

### Kruskal Wallis: Organic Acids

nonpar_pvalues <- data.frame(Plant = rep(NA, 15),
                             Treatment = rep(NA,15))

for(i in 3:ncol(Nonpar_acids)){
  KTp <- kruskal.test(Nonpar_acids[,i] ~ Plant_type, 
                      data = Nonpar_acids)
  KTt <- kruskal.test(Nonpar_acids[,i] ~ Treatment, 
                      data = Nonpar_acids)
  j <- i-2
  nonpar_pvalues$Plant[j] <- KTp[["p.value"]]
  nonpar_pvalues$Treatment[j] <- KTt[["p.value"]]
}

### Wilcox test 

## For plant type

for(i in 3:ncol(Nonpar_acids)){
  WT <- pairwise.wilcox.test(Nonpar_acids[,i], 
                             Nonpar_acids$Plant_type,
                             p.adjust.method ="none")
  print(c(i,colnames(Nonpar_acids[i])))
  print(WT[["p.value"]])
}
# No rhizosphere effect WT 3:5,8,9,11,13,14
# Rhizosphere effect 6,10,12,15
# No difference 7,16,17

## For treatment

for(i in 3:ncol(Nonpar_acids)){
  WT <- pairwise.wilcox.test(Nonpar_acids[,i], 
                             Nonpar_acids$Treatment,
                             p.adjust.method ="none")
  print(c(i,colnames(Nonpar_acids[i])))
  print(WT[["p.value"]])
}

# Pi:LP 3*,4,5,7:9,10*,11*,13,14*,16
# PM:LP 3:5,7:9,11,13,14,16,17
# No difference 6,12,15
# Phi:LP 7,9*,10*,13,16
# Pi:PM 5*,7*,8*,16*,17

