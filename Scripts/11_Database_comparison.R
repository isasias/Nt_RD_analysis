#### Phylogenetic Analysis Data Base Comparison ####

#### Libraries ####

library(dplyr)
library(tidyverse)
library(car)

#### Import Taxonomy data ####

Silvatax <- read.csv("~/Grad_School/Maestria/Processed_data/TaxaSilva_trial.csv")
RPDtax <- read.csv("~/Grad_School/Maestria/Processed_data/TaxaRPD_trial.csv")

#### Count NAs ####

## Phyla ##

sum(is.na(Silvatax$Phylum))
sum(is.na(RPDtax$Phylum))

## Class ##

sum(is.na(Silvatax$Class))
sum(is.na(RPDtax$Class))

## Order ##

sum(is.na(Silvatax$Order))
sum(is.na(RPDtax$Order))

## Family ##

sum(is.na(Silvatax$Family))
sum(is.na(RPDtax$Family))

## Genus ##

sum(is.na(Silvatax$Genus))
sum(is.na(RPDtax$Genus))

## Species ##

sum(is.na(Silvatax$Species))
sum(is.na(RPDtax$Species))

#### Merge to compare IDs ####

taxa_merge <- merge(RPDtax, Silvatax,
                          by = "X", all = TRUE)

### Phyla ###

taxa_merge$Compare <- ifelse (is.na(taxa_merge$Phylum_R), "RDP", 
                              ifelse(is.na(taxa_merge$Phylum), "Silva",
                                     ifelse(taxa_merge$Phylum == taxa_merge$Phylum_R, "Same", "Different")))
                              
Phyla_comp <- taxa_merge %>%
  count(Compare)

### Class ###

taxa_merge$Compare <- ifelse (is.na(taxa_merge$Class_R), "RDP", 
                              ifelse(is.na(taxa_merge$Class), "Silva",
                                     ifelse(taxa_merge$Class == taxa_merge$Class_R, "Same", "Different")))

Class_comp <- taxa_merge %>%
  count(Compare)

### Order ###

taxa_merge$Compare <- ifelse (is.na(taxa_merge$Order_R), "RDP", 
                              ifelse(is.na(taxa_merge$Order), "Silva",
                                     ifelse(taxa_merge$Order == taxa_merge$Order_R, "Same", "Different")))

Order_comp <- taxa_merge %>%
  count(Compare)

### Family ###

taxa_merge$Compare <- ifelse (is.na(taxa_merge$Family_R), "RDP", 
                              ifelse(is.na(taxa_merge$Family), "Silva",
                                     ifelse(taxa_merge$Family_R == taxa_merge$Family, "Same", "Different")))

Family_comp <- taxa_merge %>%
  count(Compare)

### Genus ###

taxa_merge$Compare <- ifelse (is.na(taxa_merge$Genus_R), "RDP", 
                              ifelse(is.na(taxa_merge$Genus), "Silva",
                                     ifelse(taxa_merge$Genus_R == taxa_merge$Genus, "Same", "Different")))

Genus_comp <- taxa_merge %>%
  count(Compare)

### Join Compare counts ###

taxa_count <- merge(Phyla_comp, Class_comp,
                    by = "Compare", all = TRUE)

taxa_count <- merge(taxa_count, Order_comp,
                    by = "Compare", all = TRUE)

taxa_count <- merge(taxa_count, Family_comp,
                    by = "Compare", all = TRUE)

taxa_count <- merge(taxa_count, Genus_comp,
                    by = "Compare", all = TRUE)

colnames(taxa_count) <- c("Compare", "Phyla","Class","Order","Family", "Genus")
write.csv(taxa_count, "~/Grad_School/Maestria/Processed_data/Taxa_count.csv")

#### Filter same IDs ####

## Filter all same IDs by phyla ##

Taxa_same <-  taxa_merge %>%
  filter(taxa_merge$Compare == "Same")

## Filter Same IDs by Class from new object ##

Taxa_same$Compare <- ifelse (is.na(Taxa_same$Class_R), "RDP", 
                              ifelse(is.na(Taxa_same$Class), "Silva",
                                     ifelse(Taxa_same$Class == Taxa_same$Class_R, "Same", "Different")))

Taxa_same <-  Taxa_same %>%
  filter(Taxa_same$Compare == "Same")

## Filter Order ##

Taxa_same$Compare <- ifelse (is.na(Taxa_same$Order_R), "RDP", 
                             ifelse(is.na(Taxa_same$Order), "Silva",
                                    ifelse(Taxa_same$Order_R == Taxa_same$Order, "Same", "Different")))

Taxa_same <-  Taxa_same %>%
  filter(Taxa_same$Compare == "Same")

## Filter Family ##

Taxa_same$Compare <- ifelse (is.na(Taxa_same$Family_R), "RDP", 
                             ifelse(is.na(Taxa_same$Family), "Silva",
                                    ifelse(Taxa_same$Family_R == Taxa_same$Family, "Same", "Different")))

Taxa_same <-  Taxa_same %>%
  filter(Taxa_same$Compare == "Same")

## Filter Genus ##

Taxa_same$Compare <- ifelse (is.na(Taxa_same$Genus_R), "RDP", 
                             ifelse(is.na(Taxa_same$Genus), "Silva",
                                    ifelse(Taxa_same$Genus_R == Taxa_same$Genus, "Same", "Different")))

Taxa_same <-  Taxa_same %>%
  filter(Taxa_same$Compare == "Same")

write.csv(Taxa_same, "~/Grad_School/Maestria/Processed_data/Taxa_same.csv")
# later removed repeated columns manually

#### Filter Different values ####

## Filter different IDs by phyla ##

Taxa_diff <-  taxa_merge %>%
  filter(taxa_merge$Compare == "Different")

write.csv(Taxa_diff, "~/Grad_School/Maestria/Processed_data/Taxa_diff.csv")

### Remove Same IDs from OG data ###

remove_IDs <- Taxa_same$X

Taxa_notsame <-  taxa_merge %>%
  slice(-remove_IDs) 

Taxa_notsame$Compare <- ifelse (is.na(Taxa_notsame$Phylum_R), "RDP", 
                              ifelse(is.na(Taxa_notsame$Phylum), "Silva", "Other"))

write.csv(Taxa_notsame, "~/Grad_School/Maestria/Processed_data/Taxa_notsame.csv")

#### Count analysis from different subsets ####

RPDNAs <- Taxa_notsame%>%
  filter(Compare == "RDP")

RPDNAs %>%
  count(Phylum)

## Check NAs on different IDs ##

sum(is.na(Taxa_diff$Class))
sum(is.na(Taxa_diff$Class_R))

sum(is.na(Taxa_diff$Order))
sum(is.na(Taxa_diff$Order_R))

sum(is.na(Taxa_diff$Family))
sum(is.na(Taxa_diff$Family_R))

sum(is.na(Taxa_diff$Genus))
sum(is.na(Taxa_diff$Genus_R))

### Compare IDs between databases ###

# Genus

Taxa_diff$Compare <- ifelse (is.na(Taxa_diff$Genus_R), "RDP", 
                             ifelse(is.na(Taxa_diff$Genus), "Silva",
                                    ifelse(Taxa_diff$Genus_R == Taxa_diff$Genus, "Same", "Different")))

count(Taxa_diff %>%
  filter(is.na(Genus_R)) %>%
  filter(is.na(Genus)))
