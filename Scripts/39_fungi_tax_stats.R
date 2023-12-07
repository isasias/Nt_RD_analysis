#### Fungi Statistics ####

### Libraries ###

library(dplyr)
library(tidyverse)
library(car)
library(phyloseq)


#### Data Upload ####

load("ITS_filt.RData")

#### Taxonomic Analysis ####

### Data Preparation ###

## Join phyla
Phyla_fun <- tax_glom(ITS_filtered,taxrank = "Phylum", NArm = FALSE)

## Extract matrix
OTU_matrix <- as.data.frame(Phyla_fun@otu_table)
Tax_matrix <- as.data.frame(Phyla_fun@tax_table)

## Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Phylum

## Extract metadata 
metadata <- as.data.frame(Phyla_fun@sam_data)

## Join phyla with metadata 
Fungi_phyla <- cbind(metadata,OTU_matrix)


### Test for normality ###

## Shapiro test

for(i in 5:ncol(Fungi_phyla)){
  shapiro <- shapiro.test(Fungi_phyla[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}
# only Ascomycota and Olpidiomycota normal

## Histograms

for (i in 5:ncol(Fungi_phyla)) { 
  hist(Fungi_phyla[,i],
       main = i)
} # incertae normal too

### Normalize phyla ###

## Log transform
Norm_Fungi <- Fungi_phyla
for(i in 6:19){
  Norm_Fungi[,i] <- abs(log10(Norm_Fungi[,i]+1))
}

# Recheck with shapiro
for(i in 5:ncol(Norm_Fungi)){
  shapiro <- shapiro.test(Norm_Fungi[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}

# histograms 
shapiro.test(Norm_Fungi$p__Mucoromycota)

## all normalized with log transform 

### Final table ###

Norm_Fungi <- Norm_Fungi[,-c(8,13,15,20,19,18)]

Fungi_stats <- cbind(Norm_Fungi,Fungi_phyla)
Fungi_stats <- Fungi_stats[,-c(15:21)]
Fungi_stats <- Fungi_stats[,-c(16:19)]
Fungi_stats <- Fungi_stats[,-c(17:23)]

write.csv(Fungi_stats, "~/Grad_School/Maestria/Processed_data/Fungal traits/Fungi_stats.csv")

#### ANOVA ####

## Levene Test
for(i in 5:ncol(Fungi_stats)){
  Lev_ex <- leveneTest(Fungi_stats[,i] ~ Plant_Type * Treatment, 
                       data = Fungi_stats) 
  levene <- ifelse(Lev_ex[["Pr(>F)"]]>0.05, "YES","NO")
  print(c(i,levene))
}

## Check separately for each phyla

Fungi_stats <- read.csv("~/Grad_School/Maestria/Processed_data/Fungal traits/Fungi_stats.csv")

Anova(aov(p__Mucoromycota~ Plant_Type * Treatment, 
          data = Fungi_stats))

### Tukey HSD ###

## Mucoromycota, Chytridiomycota, and Zoopagomycota

TukeyHSD(aov(p__Zoopagomycota~ Plant_Type * Treatment, 
                           data = Fungi_stats))

# Mucoromycota ptxd BS 
# Chytridiomycota both plant BS
# Zoopagomycota ptxd BS


#### Fusarium ####

Fusarium <- subset_taxa(ITS_filtered, Genus == "g__Fusarium")
Fusarium <- tax_glom(Fusarium,taxrank = "Species", NArm = FALSE)

## Extract matrix
OTU_matrix <- as.data.frame(Fusarium@otu_table)
Tax_matrix <- as.data.frame(Fusarium@tax_table)

## Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Species

## Extract metadata 
metadata <- as.data.frame(Fusarium@sam_data)

## Join phyla with metadata 
Fusarium <- cbind(metadata,OTU_matrix)

## Total fusarium counts

Sum_fus <- rowSums(Fusarium[,5:14])
Sum_fus <- cbind(metadata,Sum_fus)

# Shapiro test
shapiro.test(Sum_fus$Sum_fus) # normal

# Levene test
leveneTest(Sum_fus ~ Plant_Type * Treatment, 
           data = Sum_fus) # passed

# ANOVA
Anova(aov(Sum_fus ~ Plant_Type * Treatment, 
          data = Sum_fus)) # not significant

## Test for normality: separately ##

# Shapiro test

for(i in 5:ncol(Fusarium)){
  shapiro <- shapiro.test(Fusarium[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}

# Histograms 

for (i in 5:ncol(Fusarium)) { 
  hist(Fusarium[,i],
       main = i)
}

## Levene test
leveneTest(s__equiseti ~ Plant_Type * Treatment, 
           data = Fusarium)

## ANOVA
colnames(Fusarium)[6] <- "U"

Anova(aov(s__verticillioides~ Plant_Type * Treatment, 
           data = Fusarium))

# verticillioides marginal

TukeyHSD(aov(s__verticillioides~ Plant_Type * Treatment, 
             data = Fusarium))

## log transform data

Norm_fusa <- Fusarium

for(i in 5:14){
  Norm_fusa[,i] <- abs(log10(Norm_fusa[,i]+1))
}

# Recheck shapiro

for(i in 5:ncol(Norm_fusa)){
  shapiro <- shapiro.test(Norm_fusa[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # normalized 11,12,13


# Levene test
leveneTest(s__chlamydosporum ~ Plant_Type * Treatment, 
           data = Norm_fusa)

## ANOVA
colnames(Fusarium)
Anova(aov(s__tricinctum ~ Plant_Type * Treatment, 
          data = Norm_fusa))

# tricinctum only one significant 
# Plant Type 0.002714 **
# Plant_Type:Treatment 0.012077 *

TukeyHSD(aov(s__tricinctum~ Plant_Type * Treatment, 
             data = Norm_fusa))

#### Mucoromycota ####

### Species ###
Muco <- subset_taxa(ITS_filtered, Phylum == "p__Mucoromycota") # 3% of total species
Muco <- tax_glom(Muco,taxrank = "Species", NArm = FALSE)

Mucorales <- subset_taxa(Muco, Order == "o__Mucorales")

sum(Mucorales@otu_table)/sum(Muco@otu_table) # 99% pathogenic

Apop <- subset_taxa(ITS_filtered, Genus == "g__Apophysomyces") 
# all Apophysomyces jiangsuensis

## Extract matrix
OTU_matrix <- as.data.frame(Muco@otu_table)
Tax_matrix <- as.data.frame(Muco@tax_table)

## Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Species

## Extract metadata 
metadata <- as.data.frame(Muco@sam_data)

## Join phyla with metadata 
Muco <- cbind(metadata,OTU_matrix)

## Remove low abundance
Muco <- Muco[,-c(9:23)]

### Test for Normality ###

# Shapiro test
for(i in 5:ncol(Muco)){
  shapiro <- shapiro.test(Muco[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}

# Histograms 
for (i in 5:ncol(Muco)) { 
  hist(Muco[,i],
       main = i)
}

# None are normal 

### Log transform data ###

Norm_Muco <- Muco
for(i in 5:8){
  Norm_Muco[,i] <- abs(log10(Norm_Muco[,i]+1))
}

## re-check normality

# Shapiro test
for(i in 5:ncol(Norm_Muco)){
  shapiro <- shapiro.test(Norm_Muco[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # all but 7

# Histograms 
for (i in 5:ncol(Norm_Muco)) { 
  hist(Norm_Muco[,i],
       main = i)
} # almost normal distribution for 7

### ANOVA ###

## Levene test
colnames(Norm_Muco)
leveneTest(s__microsporus ~ Plant_Type * Treatment, 
           data = Norm_Muco) # all passed

## ANOVA
Anova(aov(s__microsporus ~ Plant_Type * Treatment, 
          data = Norm_Muco))
colnames(Norm_Muco)

# muscae significant, and microsporus marginal (0.05274)

### Tukey HSD ###

TukeyHSD(aov(s__microsporus ~ Plant_Type * Treatment, 
                          data = Norm_Muco))

# muscae plants with BS 
# microsporus ptxd with BS p= 0.0429144

## GRAPHS NOT USEFUL FOR SPECIES

### Genus ###
Muco <- subset_taxa(ITS_filtered, Phylum == "p__Mucoromycota")
Muco <- tax_glom(Muco,taxrank = "Genus", NArm = FALSE)

## Extract matrix
OTU_matrix <- as.data.frame(Muco@otu_table)
Tax_matrix <- as.data.frame(Muco@tax_table)

## Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Genus

## Extract metadata 
metadata <- as.data.frame(Muco@sam_data)

## Join phyla with metadata 
Muco <- cbind(metadata,OTU_matrix)

## Remove low abundance
Muco <- Muco[,-7]
Muco <- Muco[,-c(10,11)]

### Test for Normality ###

# Shapiro test
for(i in 5:ncol(Muco)){
  shapiro <- shapiro.test(Muco[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} 

# Histograms 
for (i in 5:ncol(Muco)) { 
  hist(Muco[,i],
       main = i)
}
# None are normal 

### Log transform data ###

Norm_Muco <- Muco

for(i in 5:9){
  Norm_Muco[,i] <- abs(log10(Norm_Muco[,i]+1))
}

## re-check normality

# Shapiro test
for(i in 5:ncol(Norm_Muco)){
  shapiro <- shapiro.test(Norm_Muco[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # all transform to normal

### ANOVA ###

## Levene test
colnames(Norm_Muco)
leveneTest(g__Absidia ~ Plant_Type * Treatment, 
           data = Norm_Muco) # all passed

## ANOVA
Anova(aov(g__Apophysomyces ~ Plant_Type * Treatment, 
          data = Norm_Muco))
colnames(Norm_Muco)

# Circinella, Apophysomyces (treatment p= 0.01024), Absidia significant 
# Mucor marginal

### Tukey HSD ###

TukeyHSD(aov(g__Apophysomyces ~ Plant_Type * Treatment, 
                          data = Norm_Muco))

# Circinella plants BS: NO GRAPH

## Apophysomyces Treatment Pi-Low P Pi-Phi  significant and Pi/Phi mix-Pi marginal
# Pi-Low P 0.0192837 Pi-Phi 0.0227198 Pi/Phi mix-Pi 0.0917214

# Absidia ptxd with WT and BS and  Bulk Soil:Pi-ptxDNt-36:Low P 0.03940615
# Mucor WT ptxD marginal p= 0.05399718


#### Chytridiomycota ####

Chy <- subset_taxa(ITS_filtered, Phylum == "p__Chytridiomycota")
Chy <- tax_glom(Chy,taxrank = "Genus", NArm = FALSE)

## Extract matrix
OTU_matrix <- as.data.frame(Chy@otu_table)
Tax_matrix <- as.data.frame(Chy@tax_table)

## Rename columns from otu_matrix with phylum from tax_matrix
colnames(OTU_matrix) <- Tax_matrix$Genus

## Extract metadata 
metadata <- as.data.frame(Chy@sam_data)

## Join phyla with metadata 
Chy <- cbind(metadata,OTU_matrix)

## Remove low abundance
Chy <- Chy[,-c(17:30)]
Chy <- Chy[,-11]

### Test for Normality ###

# Shapiro test
for(i in 5:ncol(Chy)){
  shapiro <- shapiro.test(Chy[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}

# Histograms 
for (i in 5:ncol(Chy)) { 
  hist(Chy[,i],
       main = i)
}

## None of them are normal

### Log transform data ###

Norm_Chy <- Chy
for(i in 5:15){
  Norm_Chy[,i] <- abs(log10(Norm_Chy[,i]+1))
}

## re-check normality

# Shapiro test
for(i in 5:ncol(Norm_Chy)){
  shapiro <- shapiro.test(Norm_Chy[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
} # 5, 12, 13, and 15 not normal

# Histograms 
for (i in 5:ncol(Norm_Chy)) { 
  hist(Norm_Chy[,i],
       main = i)
} # almost normal distribution for 13

### ANOVA ###

## Levene test
colnames(Norm_Chy)[10] <- "U"
leveneTest(g__Rhizophydiales_gen_Incertae_sedis ~ Plant_Type * Treatment, 
           data = Norm_Chy) # all passed

## ANOVA
Anova(aov(U ~ Plant_Type * Treatment, 
          data = Norm_Chy))
colnames(Norm_Chy)

# Rhizophlyctis NA.1 significant and Spizellomycetaceae marginal
# NA.1 not normal

### Tukey HSD ###

Fun_Tukey <- TukeyHSD(aov(NA.1 ~ Plant_Type * Treatment, 
                          data = Norm_Chy))

# Rhizophlyctis Plant BS
# Na.1 Plants BS; Pi-Low P marginal; Bulk Soil:Pi-Wild Type:Low P 0.03823750

### Kruskal Wallis for not normal data ###

kruskal.test(NA.1 ~ Plant_Type, 
             data = Norm_Chy)
colnames(Norm_Chy)

# Powellomyces marginal; Spizellomyces significant

## Wilcox
pairwise.wilcox.test(Norm_Chy$g__Spizellomyces, Norm_Chy$Plant_Type,
                     p.adjust.method ="none")

# Spizellomyces ptxD BS

#### Phyloseq stats analysis ####

Top_phyla <- tax_glom(ITS_filtered,taxrank = "Phylum", NArm = FALSE)
Top_phyla <- prune_taxa(names(sort(taxa_sums(Top_phyla),TRUE)[1:3]), Top_phyla)

sum(ITS_filtered@otu_table) # 2343951
sum(Top_phyla@otu_table) # 2165448

# Accounts for 92%

