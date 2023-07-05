library(dada2); packageVersion("dada2")

## Check path

path <- "miseqsopdata/MiSeq_SOP"
list.files(path)

## Reverse and forward 

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Plot qualities 

plotQualityProfile(fnFst[1:2])
plotQualityProfile(fnRst[1:2])

## Filter and trim

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
head(out)

## Error rates 

errF <- learnErrors(filtFs, multithread=FALSE) # in original example code it is true
errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)

## Sample inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE) 
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]] # Inspecting the returned dada-class object
# might be useful to pool the samples just in case there are some with few reads 

dadaRs[[1]]

## Merged paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

## Sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

## Remove chimeras 

# using left and right segments from abundant sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab) # divide sum no chimeras by sum original data to see % of chimeras on merged sequence reads
# 100% no chimeric sequences 

## Track number of reads filtered through pipeline

# the following commands are mostly related with data manipulation not that much with the package

getN <- function(x) sum(getUniques(x)) # create a function
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim") #order based on cbind
rownames(track) <- sample.names # kept for all data frames 
head(track) # show top of the data frame 

## Assign taxonomy

# Download the reference data base in the directory with the fastq files. Could download other references?

taxa <- assignTaxonomy(seqtab.nochim, "miseqsopdata/MiSeq_SOP/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE) # directory specific where is the data base
# for my data add it to the same folder

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Alternative done with Decipher package 

## Evaluate accuracy using a mock community

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

# This previous section might not be useful for my raw data 

#### Phyloseq ####

## libraries

BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

## Construct a data frame

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

## phyloseq object

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

## Rename taxa to a short string 

dna <- Biostrings::DNAStringSet(taxa_names(ps)) # create a new object with the taxa names
names(dna) <- taxa_names(ps) # names matching in both objects
ps <- merge_phyloseq(ps, dna)# join both objects one with complete names other with abbreviations
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) # paste0 concatenates vectors in order# code to illustrate the paste0() function 
ps

## Alpha diversity: within samples

plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
# No obvious systematic difference in alpha-diversity between early and late samples.

## Ordinate
  
# Transform data to proportions as appropriate for Bray-Curtis distances

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu)) # transform to proportions over the sum of otus
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray") # ordination using NMDS method
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS") #plot

## Plots

# bar plot

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20] # choose the top 20 using sort
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))# proportions
ps.top20 <- prune_taxa(top20, ps.top20) # remove unwanted taxa 
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
# facet wrap to attach several variables in one bar 
# free x so x axis in each graph can have different scales

