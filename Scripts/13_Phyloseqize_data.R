#### phyloseq-ize Data ####

#### Libraries ####

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("ape")

#### Data upload and modification ####

# **** Do not do this section if already have csv files ****

load("~/Grad_School/Maestria/Nt_RD_analysis/seqtabnoc.RData")
load("~/Grad_School/Maestria/Nt_RD_analysis/taxaSilva.RData")

### ASV table ###

## Change row names to sample names

sample_names <- c("A1","A2","A3","A4","A5","A6","A7",
                  "B1","B2","B3","B4","B5","B6","B7",
                  "C1","C2","C3","C4","C5","C6","C7",
                  "D1","D2","D3","D4")

row.names(seqtab.nochim) <- sample_names
  
## Generate random codes for columns

rando <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}
seqnames <- rando(15643)

## Rename columns
colnames(seqtab.nochim) <- seqnames

## Save as csv
write.csv(seqtab.nochim, "~/Grad_School/Maestria/Processed_data/ASVtable.csv")

### Taxa table ###

## Match sequence names
row.names(taxa) <- seqnames

## Save as csv
write.csv(taxa, "~/Grad_School/Maestria/Processed_data/taxa_table.csv")

# *****

### Download data for phyloseq function ###

Sample_data <- read.csv("~/Grad_School/Maestria/Nt_RD_analysis/Sample_data.csv")

## Name samples
row.names(Sample_data) <- sample_names

#### Create Phyloseq object ####

phylo_soil <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                       sample_data(Sample_data),
                       tax_table(taxa))

save(phylo_soil,file="phylo_soil.RData")


