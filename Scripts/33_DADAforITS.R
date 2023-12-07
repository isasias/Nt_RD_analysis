x#### ITS DADA workflow ####

# This workflow will begin with raw reads but there will be another 
# script using the effective files 

  
### Libraries ###

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

#### Setting working area ####

### Defining path ###

path <- "/Users/isasi/Documents/Grad_School/Maestria/Raw_data/ITS/Reads"
list.files(path)

## Forwards and reverse reads ##

fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))

#### Identify primers ####

# This section might not be necessary since primers have already been 
# removed from my samples 

FWD <- "GGAAGTAAAAGTCGTAACAAGG"
REV <- "GCTGCGTTCTTCATCGATGC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

### Filter ambiguous bases (N) ###

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

### Count primer hits ###

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


#### Remove Primers using cutadapt ####

# first install cutadapt

cutadapt <- "/Users/isasi/Documents/Grad_School/Maestria/Raw_data/ITS/Reads/cutadapt.exe"
system2(cutadapt, args = "--version") # Run shell commands from R
# this command drops cutadapt version so if it works then I can use it

### Create directories for filtered reads ###

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Reverse complement: not all of it was necessary check which one was not

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

### Run Cutadapt ###

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
} # to plot quality scores maybe add -m 1

### Re-check presence of primers ###

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

### Inspect final files ###

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))

# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

## Inspect read quality profiles ##

plotQualityProfile(cutFs[25])
plotQualityProfile(cutRs[25])

#### Filter and Trim ####

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, 
                     maxEE = c(5, 2), # modify maxEE to 5,2 like with 16s
                     truncQ = 2,
                     minLen = 50, rm.phix = TRUE, 
                     compress = TRUE, multithread = FALSE) # all values set to default 
head(out)

write.csv(out, "~/Grad_School/Maestria/Processed_data/ITS_reads.csv")


#### Error rates ####

errF <- learnErrors(filtFs, multithread = TRUE)
save(errF,file="errF_ITS.RData")

errR <- learnErrors(filtRs, multithread=TRUE)
save(errR,file="errR_ITS.RData")

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#### Sample inference ####

dadaFs_pspool <- dada(filtFs, err=errF, multithread=TRUE, pool = "pseudo") 
save(dadaFs_pspool, file = "dadaFs_ITS.RData")

dadaRs_pspool <- dada(filtRs, err=errR, multithread=TRUE, pool = "pseudo") 
save(dadaRs_pspool, file = "dadaRs_ITS.RData")

### Merge samples ###

mergers <- mergePairs(dadaFs_pspool, filtFs, 
                      dadaRs_pspool, filtRs, 
                      verbose=TRUE)
save(mergers, file = "mergers_ITS.RData")

#### Sequence table ####

seqtab <- makeSequenceTable(mergers)

dim(seqtab)
table(nchar(getSequences(seqtab)))

### Remove chimeras ###

tab_nochim_ITS <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
save(tab_nochim_ITS, file = "tab_nochim_ITS.RData")

table(nchar(getSequences(tab_nochim_ITS)))

#### Track reads through the pipeline ####

getN <- function(x) sum(getUniques(x))
out <- read.csv("~/Grad_School/Maestria/Processed_data/ITS_reads.csv")
rownames(out) <- out[,1] 
out <- out[,-1]

track <- cbind(out, sapply(dadaFs_pspool, getN), 
               sapply(dadaRs_pspool, getN), 
               sapply(mergers, getN),
               rowSums(tab_nochim_ITS))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
head(track)

write.csv(track, "~/Grad_School/Maestria/Processed_data/ITS_track.csv")

#### Assign taxonomy ####

unite.ref <- "~/Grad_School/Maestria/Nt_RD_analysis/UNITE_ITS/sh_general_release_dynamic_25.07.2023.fasta"
load("tab_nochim_ITS.RData")

sq40 <- getSequences(tab_nochim_ITS)[31:40]

taxa4 <- assignTaxonomy(sq40, unite.ref, 
                       multithread = TRUE, tryRC = TRUE)

taxaj <- rbind(taxa,taxa2)


# I did it by parts so it was able to run all the sequences

## Binding all parts together 

load("taxaj2.RData")
load("ITS_taxaj.RData")

ITS_taxa <- rbind(taxaj, taxaj2)

save(ITS_taxa, file = "ITS_taxa_final.RData")

taxa.print <- ITS_taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#### phyloseq-ize Data ####

### Libraries ###

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("ape")

### Upload data ###

### Data upload and modification ###

load("ITS_taxa_final.RData") # Taxa table
load("tab_nochim_ITS.RData") # count table 


## ASV table ##

## Change row names to sample names

sample_names <- c("A1","A2","A3","A4","A5","A6","A7",
                  "B1","B2","B3","B4","B5","B6","B7",
                  "C1","C2","C3","C4","C5","C6","C7",
                  "D1","D2","D3","D4")

row.names(tab_nochim_ITS) <- sample_names

## Generate random codes for columns

rando <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}
seqnames <- rando(7388)

## Rename columns
colnames(tab_nochim_ITS) <- seqnames

## Save as csv
write.csv(tab_nochim_ITS, "~/Grad_School/Maestria/Processed_data/ASV_ITS.csv")

### Taxa table ###

## Match sequence names
row.names(ITS_taxa) <- seqnames

## Save as csv
write.csv(ITS_taxa, "~/Grad_School/Maestria/Processed_data/ITS_taxa.csv")

### Download metadata ###

Sample_data <- read.csv("~/Grad_School/Maestria/Nt_RD_analysis/Sample_data.csv")

## Name samples
row.names(Sample_data) <- sample_names

#### Phyloseq object ####

phylo_ITS <- phyloseq(otu_table(tab_nochim_ITS, taxa_are_rows = FALSE),
                       sample_data(Sample_data),
                       tax_table(ITS_taxa))

save(phylo_ITS,file="phylo_ITS.RData")

### Filter taxa with low abundance ###

load("phylo_ITS.Rdata")

# Remove unidentified taxa
ITS_filtered <- subset_taxa(phylo_ITS, !is.na(Phylum))

# ASVs that are found in at least 3 different samples
# ITS_filtered <-  filter_taxa(ITS_filtered, function(OTU) sum(OTU) > 3, TRUE)
# removed because it reduces sample size too much

### Transform to relative abundance ###

ITS_rel  <-  transform_sample_counts(ITS_filtered, function(x) x / sum(x) )
# This divides the count of one species by the total number of species in that sample

save(ITS_filtered,file="ITS_filt.RData")
save(ITS_rel,file="ITS_rel.RData")
