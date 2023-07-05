#### Root exudate normalization ####

### Libraries ###

library(dplyr)
library(tidyverse)
library(car)

### Data upload ###

targ_ex <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/target_exudates_weightcorr.csv")
group_ex <- read.csv("~/Grad_School/Maestria/Processed_data/Exudates/target_groupexu_weightcorr.csv")

### Pre-process table ###

Exud_table <-  targ_ex %>%
  select(Name,X36C1_conc:NPPM_corr)

## Transpose 
Exud_table <- as.data.frame(t(Exud_table))

colnames(Exud_table) <- c(paste0("Exud_",as.character(1:59)))
# change names to exudate number for an easier analysis

Exud_list <- as.data.frame(t(Exud_list))
# Exudate list 
Exud_list <- Exud_table[1,]

Exud_table <- Exud_table%>%
  filter(!row_number() %in% c(1))

## Change all to numeric variables 
i <- c(1:59)                                  # Specify columns you want to change
Exud_table[,i] <- apply(Exud_table[,i],2,            # Specify own function within apply
                      function(x) as.numeric(as.character(x)))
# Check mode
sapply(Exud_table,mode)

## Add Treatment and Plant Type columns 

Pl_type <- c("ptxD36","ptxD36","ptxD36","ptxD36",
             "ptxD36","ptxD36","ptxD36","ptxD36",
             "ptxD36","ptxD36","ptxD36","ptxD36",
             "Wild Type","Wild Type","Wild Type",
             "Wild Type","Wild Type","Wild Type",
             "Wild Type","Wild Type","Wild Type",
             "No Plant","No Plant","No Plant","No Plant")

Trt <- c("Control","Control","Control",
         "Low P","Low P","Low P","Phi","Phi","Phi",
         "P Mix","P Mix","P Mix",
         "Control","Control","Control",
         "Low P","Low P","Low P",
         "P Mix","P Mix","P Mix",
         "Control", "Low P","Phi","P Mix")

Exud_table <- Exud_table %>%
  add_column(Plant_type = Pl_type) %>%
  add_column(Treatment = Trt)

### Test for normality ###

for(i in 1:ncol(Exud_table)){
  shapiro <- shapiro.test(Exud_table[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}

# Filter normal and no normal metabolites 

Normal_Exud <- Exud_table %>%
  select(c(Exud_4,Exud_11,Exud_14,Exud_19,Exud_20,Exud_22,
           Exud_26,Exud_27,Exud_32,Exud_34,Exud_36:Exud_39,
           Exud_47,Exud_50,Exud_56,Exud_59)) %>%
  add_column(Plant_type = Pl_type) %>%
  add_column(Treatment = Trt)

Nnormal_Exud <- Exud_table %>%
  select(-c(Exud_4,Exud_11,Exud_14,Exud_19,Exud_20,Exud_22,
           Exud_26,Exud_27,Exud_32,Exud_34,Exud_36:Exud_39,
           Exud_47,Exud_50,Exud_56,Exud_59)) 

Exud_12,Exud_10,
Exud_33,Exud_5,Exud_28,Exud_51,Exud_53,Exud_54,
Exud_55,Exud_57,Exud_58,Exud_52,Exud_18,Exud_21,
Exud_29,Exud_31,Exud_45,Exud_49
### Normalize metabolites

## Log transformed
for(i in 1:43){
  Nnormal_Exud[,i] <- abs(log10(Nnormal_Exud[,i]+1))
} # worked for more or less half of the data

for(i in 1:ncol(Nnormal_Exud)){
  shapiro <- shapiro.test(Nnormal_Exud[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}

## Scale data

scale_Exud <- as.data.frame(abs(scale(Nnormal_Exud[,1:29])))

# Log transform?

for(i in 1:29){
  scale_Exud[,i] <- abs(log(scale_Exud[,i]))
} # worked for more or less half of the data

## Shapiro check

for(i in 1:ncol(scale_Exud)){
  shapiro <- shapiro.test(scale_Exud[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
  }

# Add new normal 

Normal_Exud <- Normal_Exud %>%
  add_column(Exud_5log = Nnormal_Exud$Exud_5) %>%
  add_column(Exud_28log = Nnormal_Exud$Exud_28)%>%
  add_column(Exud_51log = Nnormal_Exud$Exud_51)%>%
  add_column(Exud_53log = Nnormal_Exud$Exud_53)%>%
  add_column(Exud_54log = Nnormal_Exud$Exud_54)%>%
  add_column(Exud_55log = Nnormal_Exud$Exud_55)%>%
  add_column(Exud_57log = Nnormal_Exud$Exud_57)%>%
  add_column(Exud_58log = Nnormal_Exud$Exud_58)


## Check normal data with histograms


# Log transformed data
for (i in 1:30) { 
  hist(Nnormal_Exud[,i],
       main = i)
}

# scaled data
for (i in 1:29) { 
  hist(scale_Exud[,i],
       main = i)
}

## Exudates bordering in normality

# Exud 18: Add
hist(scale_Exud$Exud_18) 
shapiro.test(scale_Exud$Exud_18) # 0.01
qqnorm(scale_Exud$Exud_18)

# Exud 21: Add
hist(scale_Exud$Exud_21) 
shapiro.test(scale_Exud$Exud_21) # 0.02 
qqnorm(scale_Exud$Exud_21)

# Exud 29: Add
hist(scale_Exud$Exud_29) 
shapiro.test(scale_Exud$Exud_29) # 0.01 
qqnorm(scale_Exud$Exud_29)

# Exud 31: Add
hist(scale_Exud$Exud_31) 
shapiro.test(scale_Exud$Exud_31) # 0.02 
qqnorm(scale_Exud$Exud_31)

# Exud 45: Add
hist(scale_Exud$Exud_45) 
shapiro.test(scale_Exud$Exud_45) # 0.03 
qqnorm(scale_Exud$Exud_45)

# Exud 49: Add
hist(scale_Exud$Exud_49) 
shapiro.test(scale_Exud$Exud_49) # 0.02 
qqnorm(scale_Exud$Exud_49)

## Add exudates bordering on normal

Normal_Exud <- Normal_Exud %>%
  add_column(Exud_12 = Nnormal_Exud$Exud_12) %>% # 0.03
  add_column(Exud_10log = Nnormal_Exudlog$Exud_10)%>% # 0.02
  add_column(Exud_33log = Nnormal_Exudlog$Exud_33) # 0.01

Normal_Exud <- Normal_Exud %>%
  add_column(Exud_52log101 = Nnormal_Exud$Exud_52) # 0.01

Normal_Exud <- Normal_Exud %>%
  add_column(Exud_6sc = scale_Exud$Exud_6) %>% # 0.03
  add_column(Exud_8sc = scale_Exud$Exud_8)

Normal_Exud <- Normal_Exud %>%
  add_column(Exud_18sclog = scale_Exud$Exud_18) %>% # 0.01
  add_column(Exud_21sclog = scale_Exud$Exud_21) %>% # 0.02
  add_column(Exud_29sclog = scale_Exud$Exud_29) %>% # 0.01
  add_column(Exud_31sclog = scale_Exud$Exud_31) %>% # 0.02
  add_column(Exud_45sclog = scale_Exud$Exud_45) %>% # 0.03
  add_column(Exud_49sclog = scale_Exud$Exud_49) # 0.02

## Save target tables

Exud_list <- as.data.frame(t(Exud_list))

write.csv(Exud_list, "~/Grad_School/Maestria/Processed_data/Exudates/Exudate_list.csv")
write.csv(Normal_Exud, "~/Grad_School/Maestria/Processed_data/Exudates/Normalized_exudates.csv")
write.csv(Nnormal_Exud, "~/Grad_School/Maestria/Processed_data/Exudates/UnNorm_exudates.csv")
write.csv(Exud_table, "~/Grad_School/Maestria/Processed_data/Exudates/Final_exudates.csv")

