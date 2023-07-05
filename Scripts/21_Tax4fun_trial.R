#### Tax4Fun2 tutorial ####

library(Tax4Fun)
library(phyloseq)


## transpose count table

functional_counts <- read.csv("~/Grad_School/Maestria/Processed_data/functional_counts.csv", 
                              header=FALSE, row.names=1)
colnames(functional_counts) <- functional_counts[1,]
SoilsASVs<- functional_counts[-1,]

SoilsASVs <- as.data.frame(t(SoilsASVs))
SoilsASVs <- transform(SoilsASVs, A1 = as.numeric(A1), 
                       A2 = as.numeric(A2),A3 = as.numeric(A3),
                       A4 = as.numeric(A4),A5 = as.numeric(A5),
                       A6 = as.numeric(A6),A7 = as.numeric(A7),
                       B1 = as.numeric(B1),B2 = as.numeric(B2),
                       B3 = as.numeric(B3),B4 = as.numeric(B4),
                       B5 = as.numeric(B5),B6 = as.numeric(B6),
                       B7 = as.numeric(B7),C1 = as.numeric(C1),
                       C2 = as.numeric(C2),C3 = as.numeric(C3),
                       C4 = as.numeric(C4),C5 = as.numeric(C5),
                       C6 = as.numeric(C6),C7 = as.numeric(C7),
                       D2 = as.numeric(D2),D1 = as.numeric(D1),
                       D3 = as.numeric(D3),D4 = as.numeric(D4))

samplenames <- c("A1","A2","A3","A4","A5","A6","A7",
                 "B1","B2","B3","B4","B5","B6","B7",
                 "C1","C2","C3","C4","C5","C6","C7",
                 "D1","D2","D3","D4")

ASVmatrix <- data.matrix(SoilsASVs)


Soillist <- list(sampleNames = samplenames,
                 otuTable = ASVmatrix)
save(Soillist,file="soil_tax4fun.RData")

# Prediction
Funcs <- Tax4Fun(Soillist,
        folderReferenceData = "~/Grad_School/Maestria/Nt_RD_analysis/SILVA123/")

# Sample data
data(GN16SData)
Funcsample <- Tax4Fun(GN16SData,
                 folderReferenceData = "~/Grad_School/Maestria/Nt_RD_analysis/SILVA123/")


#### Functional profiles ####

FuncProfile <- Funcs[["Tax4FunProfile"]]
Prof_samp <- Funcsample[["Tax4FunProfile"]]

write.csv(FuncProfile, "~/Grad_School/Maestria/Processed_data/Functional_profiles.csv")

#### Heat maps ###

# Transpose matrix for heat map
TfunP <- t(FuncProfile)

# original heat map
heatmap(TfunP)

# no dendogram or reordering
heatmap(TfunP, Colv = NA, Rowv = NA)
