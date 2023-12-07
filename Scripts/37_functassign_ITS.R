#### Fungal traits ####

### Libraries ###

devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("traitecoevo/fungaltraits")
library(fungaltraits) ## the database
library(microeco) ## supplementary material to table
library(phyloseq)
library(tidyverse)
library(dplyr)

### Access the data base

Fungal_traits <- fungal_traits()
write.csv(Fungal_traits, "~/Grad_School/Maestria/Processed_data/Fungal traits/Fungal_traits.csv")

#### Using microeco package ####

### Phyloseq to microtable ###

load("ITS_filt.RData") # phyloseq2meco does not work 

### Data bases
FT <- fungi_func_FungalTraits # from supplementary material 
FUNGuild <- fungi_func_FUNGuild # not very useful

#### Genus functional assignation ####

### Remove unidentified genus 

ITS_genus <- subset_taxa(ITS_filtered, !is.na(Genus))

### Extract tax table 

Fun_genus <- as.data.frame(ITS_genus@tax_table)
write.csv(Fun_genus, "~/Grad_School/Maestria/Processed_data/Fungal traits/Fungal_genus.csv")

### Upload corrected table

Fun_genus <- read.csv("~/Grad_School/Maestria/Processed_data/Fungal traits/Fungal_genus.csv")
row.names(Fun_genus) <- Fun_genus[,1]
Fun_genus <- Fun_genus[,-1]

### Merge two data frames

# Rename Genus
colnames(FT)[6] = "Genus"

## Merge

Fun_fun <-merge(x = Fun_genus, y = FT, by = "Genus")

### Clean data frame

Fun_fun <- Fun_fun[,-(3:7)]
Fun_fun <- Fun_fun[,-4]

# Remove empty columns
empty_columns <- sapply(Fun_fun, function(x) all(is.na(x) | x == ""))
Fun_fun <- Fun_fun[,-(24:25)]

# save dataframe
write.csv(Fun_fun, "~/Grad_School/Maestria/Processed_data/Fungal traits/Fungi_funct.csv")

#### Join with counts ####

## Separate counts 

Fun_counts <- as.data.frame(ITS_genus@otu_table) # extract from phyloseq object
Fun_counts <- Fun_counts[,Fun_fun[,2]] # Filtered identified lifestyles

## Add primary and secondary lifestyles as rows 

Fun_counts[nrow(Fun_counts) + 1,] <- Fun_fun[,9]
Fun_counts[nrow(Fun_counts) + 1,] <- Fun_fun[,10]

# save dataframe
write.csv(Fun_counts, "~/Grad_School/Maestria/Processed_data/Fungal traits/Fungi_lifestlyes.csv")
