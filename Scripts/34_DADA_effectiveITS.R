#### ITS DADA with effective reads####

### Libraries ###

library(dada2)
packageVersion("dada2")

### Setting working area ###

### Defining path ###

path <- "/Users/isasi/Documents/Grad_School/Maestria/Raw_data/ITS/Effective_Reads/"
list.files(path)


### Effective reads preprocessing ###

ITSefs <- sort(list.files(path, pattern = "effective.fastq.gz", full.names = TRUE))

ITS.filtN <- file.path(path, "filtN", basename(ITSefs))

## Filter and trim ##
filterAndTrim(ITSefs, ITS.filtN, maxN = 0, minLen = 50, multithread = FALSE)

### Error rates ###
err_ITS <- learnErrors(ITS.filtN, multithread = TRUE)
save(err_ITS,file="err_effITS.RData")
plotErrors(err_ITS)

### Sample inference ###
dadaITS_pspool <- dada(ITS.filtN, err=err_ITS, multithread=TRUE, pool = "pseudo") 
save(dadaITS_pspool, file = "dadaFs_ITS.RData")

# No merge use directly dada?

#### Get sequence table ####
seqtab <- makeSequenceTable(dadaITS_pspool)
dim(seqtab)
table(nchar(getSequences(seqtab)))

### Remove chimeras ###

eff_nochim_ITS <- removeBimeraDenovo(seqtab, method="consensus", 
                                     multithread=TRUE, verbose=TRUE)
save(eff_nochim_ITS, file = "eff_nochim_ITS.RData")
table(nchar(getSequences(eff_nochim_ITS)))




