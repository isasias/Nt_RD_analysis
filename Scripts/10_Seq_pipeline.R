#### DADA 2 pipeline for 16S ####

## Libraries

library(dada2); packageVersion("dada2")
library(tidyverse)
library(dplyr)

# BiocManager::install("FastqCleaner") # removing sequences with N bases
# library(FastqCleaner)

## Set path

path <- "~/Grad_School/Maestria/Raw_data/16s/Reads/filtered"
list.files(path)

## Extract sample names

fnFs <- sort(list.files(path, pattern="_1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Plot qualities

plotQualityProfile(fnFs[1:2]) # compared to tutorial ir seems that it does not need trimming
# However if trimming is needed cut to 180

plotQualityProfile(fnRs[1:2]) # same not trimming but if needed cut to 200


#### Filtering and trimming #### 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(5,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # this step removes n sequences too 

### Add filtered reads 

filtFs <- sort(list.files(path, pattern="_F_filt.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(path, pattern="_R_filt.fastq.gz", full.names = TRUE))

## Notes on the code

# trunQ=2 do not reduce it
# truncLen= I started modifying length based on error rate and quality score distribution graphs from QC report
# Decided to keep it on 160,200 because it shows the highest reads out 
# Reads are 225 bases in length
# maxEE= seems to be the parameter discarding more reads so I kept it at (5,2)


#### Error rates ####

errF <- learnErrors(filtFs, multithread=TRUE) # modify nbases for trial only 
errR <- learnErrors(filtRs, multithread=TRUE) # Time: 28 min and 32 min

### Save the objects just in case it crashes 
save(errF,file="errF.RData")
save(errR,file="errR.RData")

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# You need to calculate error rates to use them for sample inference
# Since the code depends on a parametric error model it uses this variable to 
# calculate the sample inference

## Do I need to dereplicate the sequences?

#### Sample inference ####

# Original code: no pooling
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) # multithread just has to do with processing speed
dadaRs <- dada(filtRs, err=errR, multithread=TRUE) 

# 24 min + 19 min

# Pooling
dadaFs_pool <- dada(filtFs, err=errF, multithread=TRUE, pool = TRUE) 
dadaRs_pool <- dada(filtRs, err=errR, multithread=TRUE, pool = TRUE) # have not run it yet

# Pseudo-pooling
dadaFs_pspool <- dada(filtFs, err=errF, multithread=TRUE, pool = "pseudo") 
dadaRs_pspool <- dada(filtRs, err=errR, multithread=TRUE, pool = "pseudo") # have not run it yet

save(dadaRs_pspool, file = "dadaRs_pspool.RData")

#### Merging forward and reverse  reads ####

# OG data
mergers <- mergePairs(dadaFs_pspool, filtFs, dadaRs_pspool, filtRs, verbose=TRUE) # just 6 minutes
# mergers <- mergePairs(dadaRs_pspool, filtRs, dadaRs, filtRs, verbose=TRUE, justConcatenate = TRUE)

save(mergers, file = "mergers.RData")

# Missing tests with pooled and pseudo-pooled data

#### Sequence table ####

seqtab <- makeSequenceTable(mergers)

# Inspect table
dim(seqtab) # number of samples times number os ASVs
table(nchar(getSequences(seqtab))) # sequence lengths need to be less than 450

#### Remove Chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
# Around 29 min

save(seqtab.nochim, file = "seqtabnoc.RData")
dim(seqtab.nochim) # 58% of my samples removed
sum(seqtab.nochim)/sum(seqtab) # abundances of chimeric variant account for only 8%

# Check if by truncating the samples this can be reduced

#### Track number of reads filtered through pipeline ####

getN <- function(x) sum(getUniques(x)) # create a function
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim") #order based on cbind
rownames(track) <- sample.names # kept for all data frames 
head(track) # show top of the data frame 

## Save track table
write.csv(track, "~/Grad_School/Maestria/Processed_data/Track_table_trial.csv")

#### Assign Taxonomy ####

### Using Silva version 138 ###

taxa <- assignTaxonomy(seqtab.nochim, "~/Grad_School/Maestria/Raw_data/16s/Reads/tax/silva_nr_v132_train_set.fa.gz", 
                       multithread=TRUE) # directory specific where is the data base
# 23 minutes approx

save(taxa, file = "taxaSilva.RData")

## Adding species

taxa <- addSpecies(taxa, "~/Grad_School/Maestria/Raw_data/16s/Reads/tax/silva_species_assignment_v132.fa.gz")
# 13 min

# Saving taxa file without the sequence

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.csv(taxa.print, "~/Grad_School/Maestria/Processed_data/TaxaSilva_trial.csv")

### Using RDP 18: to compare between data bases ###

taxaRPD <- assignTaxonomy(seqtab.nochim, "~/Grad_School/Maestria/Raw_data/16s/Reads/tax/rdp_train_set_18.fa.gz", multithread=TRUE) 
## 11 min

## Adding species

taxaRPD <- addSpecies(taxaRPD, "~/Grad_School/Maestria/Raw_data/16s/Reads/tax/rdp_species_assignment_18.fa.gz")
# 6 min
# Saving taxa file without the sequence

taxaR.print <- taxaRPD # Removing sequence rownames for display only
rownames(taxaR.print) <- NULL
head(taxaR.print)

write.csv(taxaR.print, "~/Grad_School/Maestria/Processed_data/TaxaRPD_trial.csv")
