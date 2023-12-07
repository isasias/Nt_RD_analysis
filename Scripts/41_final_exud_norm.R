#### Exudate Sums stats ####

### Libraries ###

library(dplyr)
library(tidyverse)
library(car)

#### Data upload and pre-processing ####

## Metadata ##
metadata <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Exud_metadata.csv")

## Exud sums ##

Exud_sums <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Sum_metabolites.csv")
Exud_sums <- as.data.frame(t(Exud_sums)) # transpose
colnames(Exud_sums) <- Exud_sums[1,]
Exud_sums <- Exud_sums[-1,]

# Add metadata
Exud_sums <- cbind(metadata, Exud_sums)

### Change to numeric variables ###

## Exudates

i <- c(3:108) # Specify columns you want to change
Exud_sums[,i] <- apply(Exud_sums[,i],2, # Specify own function within apply
                       function(x) as.numeric(as.character(x)))
# Check mode
sapply(Exud_sums,mode)

## Organic Acids ##

Org_acids <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/Org_acids_sums.csv")
Org_acids <- as.data.frame(t(Org_acids)) # transpose
colnames(Org_acids) <- Org_acids[1,]
Org_acids <- Org_acids[-1,]

# Add metadata
Org_acids <- cbind(metadata, Org_acids)


## Organic acids

i <- c(3:42) # Specify columns you want to change
Org_acids[,i] <- apply(Org_acids[,i],2, # Specify own function within apply
                       function(x) as.numeric(as.character(x)))
# Check mode
sapply(Org_acids,mode)

#### Test for normality: Exud list ####

### Shapiro Tests

# Normal data
for(i in 3:ncol(Exud_sums)){
  shapiro <- shapiro.test(Exud_sums[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # 29 51 54 59 62 65 71 80 84 90 91 93 97 99 105 108

# Log transformed
for(i in 3:ncol(Exud_sums)){
  shapiro <- shapiro.test(log(Exud_sums[,i]))
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # 11 69 72 85

# Square root transformed
for(i in 3:ncol(Exud_sums)){
  shapiro <- shapiro.test(sqrt(Exud_sums[,i]))
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # 3 7 13 14 16 18 19 21 22 24 30 31 33 35 38 39 41 43 45 47:50
# 53 58 66:68 70 73 75:79 81 86 89 92 95 96 102 107

### Final normalized table 

Ex_log <- log(Exud_sums[,c(11,69,72,85)])

Exud_sq <- sqrt(Exud_sums[,c(3,7,13,14,16,18,19,21,22,24,30,31,33,
                             35,38,39,41,43,45,47:50,53,58,66:68, 
                             70,73,75:79,81,86,89,92,95,96,102,107)])

Exud_unt <- Exud_sums[,c(1,2,27,29,51,54,59,62,65,71,80,
                     84,90,91,93,97,99,105,108)] # add 27 based on shapiro

Exud_norm <- cbind(Exud_unt,Ex_log,Exud_sq)

### Not normal table

Untransf_exud <- Exud_sums[,c(4,5,6,8,9,10,12,15,17,20,23,25,26,
                              28,32,34,36,37,40,42,44,46,52,
                              55,56,57,60,61,63,64,74,82,83,
                              87,88,94,98,100,101,103,104,106)]

#### Check particulars for normality ####

## Histograms
for (i in 1:ncol(Untransf_exud)) { 
  hist(sqrt(Untransf_exud[,i]),
       main = i)
}
# sqrt good graphs 1,8*,9*,10,14,16,17,20,24,28,29,31,32*,33,34,35*,37*,40
# log 1,4*,5*,6,11*,14,15,16*,17,19,20,22,23,24,27,28,29,33*,38,39,40,41
# no good graphs 2,3,12,13,18,21,25,26,30,36,42
# * good enough graphs and decent shapiro

## qqplot
for (i in 1:ncol(Untransf_exud)) { 
  qqnorm(sqrt(Untransf_exud[,i]),
       main = i)
} 

## Print shapiro

for(i in 1:ncol(Untransf_exud)){
  shapiro <- shapiro.test(log(Untransf_exud[,i]))
  print(c(i,shapiro[["p.value"]]))
} # Normal 27 already added to untransformed values
# Log 3,6,17,19,22,23,28,29,38 (remove 14,20 because of graph)
# sqrt 1,7,10,14,20,24,31,33,34,40:42
# remove 19,38 because log graph was better 
# remove 2,11,13,18,21,25,30,36 because of bad graphs

Log_ex <- Untransf_exud[,c(3,4,5,6,11,16,17,19,22,23,28,29,38)]

for(i in 1:ncol(Log_ex)){
  Log_ex[,i] <- log(Log_ex[,i])
} 

sqrt_ex <- Untransf_exud[,c(1,7,8,9,10,14,20,24,31:35,37,40:42)]

for(i in 1:ncol(sqrt_ex)){
  sqrt_ex[,i] <- sqrt(sqrt_ex[,i])
} 

# Final table 
Exud_norm <- cbind(Exud_norm,Log_ex,sqrt_ex)

# Untransform exudate
Untransf_exud <- Untransf_exud[,c(2,12,13,15,18,21,25,
                                  26,27,30,36,39)]

#### Save csvs ####

write.csv(Exud_norm, "~/Grad_School/Maestria/Processed_data/Exudates/Normalized_exudates.csv")
write.csv(Untransf_exud, "~/Grad_School/Maestria/Processed_data/Exudates/Untransform_exudates.csv")

#### Normalization Organic ACids ####

### Shapiro Tests

# Normal data
for(i in 3:ncol(Org_acids)){
  shapiro <- shapiro.test(Org_acids[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # 16 17 39:41

# Log transformed
for(i in 3:ncol(Org_acids)){
  shapiro <- shapiro.test(log(Org_acids[,i]))
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # 7 18 25 

# Square root transformed
for(i in 3:ncol(Org_acids)){
  shapiro <- shapiro.test(sqrt(Org_acids[,i]))
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # 4 14 15 24 35:38 41

## Normal acids

Norm_acids <- Org_acids[,c(1,2,16,17,39:41)]

Log_acids <- Org_acids[,c(7,18,25)]
for(i in 1:3){
  Log_acids[,i] <- log(Log_acids[,i])
}

sq_acids <- Org_acids[c(4,14,15,24,35:38)]
for(i in 1:8){
  sq_acids[,i] <- sqrt(sq_acids[,i])
}

Norm_acids <- cbind(Norm_acids,Log_acids,sq_acids)

# Untransformed table
Unt_acids <- Org_acids[,-c(1,2,16,17,39:41,7,18,25,4,14,15,24,
                           35:38,41)]

#### Check particulars for normality ####

## Histograms
for (i in 1:ncol(Unt_acids)) { 
  hist(sqrt(Unt_acids[,i]),
       main = i)
} # 21 18 
# Log 24 22 20 19 17 16 15 14 10 9 5 1
# sqrt 20 18 12 1

## qqplot
for (i in 1:ncol(Unt_acids)) { 
  qqnorm(sqrt(Unt_acids[,i]),
         main = i)
} # 21 18
# Log 24 22 17 16 14 10 9 8 5
# sqrt 22 20 18 

## Print shapiro

for(i in 1:ncol(Unt_acids)){
  shapiro <- shapiro.test(sqrt(Unt_acids[,i]))
  print(c(i,shapiro[["p.value"]]))
} # 21 18 
# Log 2 9 10* 11 14 15 16* 17* 20* 22 24
# sqrt 1* 3 4 5 7 9 11* 12* 20 23 24 


Lg_acids <- Unt_acids[,c(10,16,17,20)]
for(i in 1:4){
  Lg_acids[,i] <- log(Lg_acids[,i])
}

sq_acids <- Unt_acids[,c(1,11,12)]
for(i in 1:3){
  sq_acids[,i] <- sqrt(sq_acids[,i])
}

Norm_acids <- cbind(Norm_acids,Unt_acids[,c(18,21)],Lg_acids,sq_acids)
Unt_acids <- Unt_acids[,-c(21,18,10,16,17,20,1,11,12)]

### Save tables

write.csv(Norm_acids, "~/Grad_School/Maestria/Processed_data/Exudates/Normalized_acids.csv")
write.csv(Unt_acids, "~/Grad_School/Maestria/Processed_data/Exudates/Untransform_acids.csv")
