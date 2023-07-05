#### Two-way ANOVA ####

#### Libraries ####

library(dplyr)
library(tidyverse)
library(car)
library(vegan)

#### Import data ####

Data_OG <- read.csv("~/Grad_School/Maestria/Processed_data/RSratio_table.csv")
Data_log <- read.csv("~/Grad_School/Maestria/Processed_data/Biomass_logtransf.csv")
alpha_div <- read.csv("~/Grad_School/Maestria/Processed_data/Alphadiversity.csv")
s_meta <- read.csv("~/Grad_School/Maestria/Processed_data/metadata.csv")
load("Psoil_filt.Rdata")

#### Plant Biomass ####

## Create table removing no plant treatment for analysis

Data_plants <- Data_log %>%
  filter(!Plant_Type == "No plant") %>% # remove row
  dplyr::select(-"X") # remove column

table(Data_plants$Plant_Type, Data_plants$Treatment)

### Roots ###

## ANOVA ##

root_anova <- aov(log_root ~ Plant_Type * Treatment, 
                  data = Data_plants) # synergistic effect
Anova(root_anova)

# Test homogeneity of variance

plot(root_anova, 1) # check outliers
leveneTest(log_root ~ Plant_Type * Treatment, 
           data = Data_plants) # higher than 0.05 so it meets ANOVA assumption

## Tukey HSD ##

rt_Tukey <- TukeyHSD(root_anova)

### Shoots ###

## ANOVA ##

shoot_anova <- aov(log_shoot ~ Plant_Type * Treatment, 
                   data = Data_plants) # synergistic effect
Anova(shoot_anova)

# Test homogeneity of variance

plot(shoot_anova, 1) # check outliers
leveneTest(log_shoot ~ Plant_Type * Treatment, 
           data = Data_plants) # higher than 0.05 so it meets ANOVA assumption

## Tukey HSD ##

sh_Tukey <- TukeyHSD(shoot_anova)

### Total biomass ###

biom_anova <- aov(log_root+log_shoot ~ Plant_Type * Treatment, 
                  data = Data_plants) # synergistic effect
Anova(biom_anova)
biom_Tukey <- TukeyHSD(biom_anova)

### Mean weight ###

# means by groups
Means_table <- Data_OG %>%
  group_by(Treatment, Plant_Type)%>%
  summarise(mean_roots= mean(Roots),
            mean_shoots= mean(Shoots))

Means_table <- Data_OG %>%
  group_by(Treatment, Plant_Type)%>%
  summarise(mean_length= mean(root_length),
            mean_area= mean(root_area))

#### Root structure ####

Data_roots <- Data_OG %>%
  filter(!Plant_Type == "No Plant") %>%
  filter(!is.na(root_area))

### Area ###

## ANOVA 
rtarea_anova <- aov(root_area ~ Plant_Type * Treatment, 
                   data = Data_roots) # synergistic effect

Anova(rtarea_anova)

# Test homogeneity of variance

plot(rtarea_anova, 1) # check outliers
leveneTest(root_area ~ Plant_Type * Treatment, 
           data = Data_roots) # higher than 0.05 so it meets ANOVA assumption

## Tukey HSD ##

rta_Tukey <- TukeyHSD(rtarea_anova)

### Length ###

## ANOVA 
rtl_anova <- aov(root_length ~ Plant_Type * Treatment, 
                    data = Data_roots) # synergistic effect

Anova(rtl_anova)

# Test homogeneity of variance

plot(rtl_anova, 1) # check outliers
leveneTest(root_length ~ Plant_Type * Treatment, 
           data = Data_roots) # higher than 0.05 so it meets ANOVA assumption

## Tukey HSD ##

rtl_Tukey <- TukeyHSD(rtl_anova)

### Tip count ###

## ANOVA 
rtti_anova <- aov(root_tip_count ~ Plant_Type * Treatment, 
                 data = Data_roots) # synergistic effect

Anova(rtti_anova)

# Test homogeneity of variance

plot(rtti_anova, 1) # check outliers
leveneTest(root_tip_count ~ Plant_Type * Treatment, 
           data = Data_roots) # higher than 0.05 so it meets ANOVA assumption

## Tukey HSD ##

rtti_Tukey <- TukeyHSD(rtti_anova)

## Error bars

tips_sdev <- Data_OG %>%
  group_by(Treatment, Plant_Type)%>%
  summarise(dv_tips= sd(root_tip_count))

### Lateral roots

Lat_roots[Lat_roots == 0] <- NA

rlat_anova <- aov(lat_roots ~ Plant_Type * Treatment, 
                  data = Lat_roots) # synergistic effect

Anova(rlat_anova)

# Test homogeneity of variance

plot(rlat_anova, 1) # check outliers
leveneTest(lat_roots ~ Plant_Type * Treatment, 
           data = Lat_roots) # higher than 0.05 so it meets ANOVA assumption

## Tukey HSD ##

rlat_Tukey <- TukeyHSD(rlat_anova)

### Biomas Pi ###

bpi_anova <- aov(shoot_Pi ~ Plant_Type * Treatment, 
                 data = Data_roots) # synergistic effect

# Test homogeneity of variance

plot(bpi_anova, 1) # check outliers
leveneTest(biomass_Pi ~ Plant_Type * Treatment, 
           data = Data_roots) # higher than 0.05 so it meets ANOVA assumption

Anova(bpi_anova) # Treatments and interactions are significant 
bpi_Tukey <- TukeyHSD(bpi_anova)

#### Alpha diversity ####

obs_anova <- aov(Observed ~ Plant_Type * Treatment, 
                 data = alpha_div)

plot(obs_anova, 1) # check outliers
leveneTest(Observed ~ Plant_Type * Treatment, 
           data = alpha_div) # higher than 0.05 so it meets ANOVA assumption

Anova(obs_anova) # not significant difference

### Chao 1

Chao1_anova <- aov(Chao1 ~ Plant_Type * Treatment, 
                   data = alpha_div)

plot(Chao1_anova, 1) # check outliers
leveneTest(Chao1 ~ Plant_Type * Treatment, 
           data = alpha_div) # higher than 0.05 so it meets ANOVA assumption

Anova(Chao1_anova) # not significant difference

### Shannon

Shannon_anova <- aov(Shannon ~ Plant_Type * Treatment, 
                     data = alpha_div)

plot(Shannon_anova, 1) # check outliers
leveneTest(Shannon ~ Plant_Type * Treatment, 
           data = alphadiv) # higher than 0.05 so it meets ANOVA assumption

Anova(Shannon_anova) # not significant difference

### Simpson

Simp_anova <- aov(Simpson ~ Plant_Type * Treatment, 
                  data = alpha_div)

plot(Simp_anova, 1) # check outliers
leveneTest(Simpson ~ Plant_Type * Treatment, 
           data = alphadiv) # higher than 0.05 so it meets ANOVA assumption

Anova(Simp_anova) # not significant difference

### Evenness Pielou

Pielou_anova <- aov(pielou ~ Plant_Type * Treatment, 
                    data = alphadiv)

plot(Pielou_anova, 1) # check outliers
leveneTest(pielou ~ Plant_Type * Treatment, 
           data = alphadiv) # higher than 0.05 so it meets ANOVA assumption

Anova(Pielou_anova) # not significant difference

#### PERMANOVA for Beta Diversity ####

### Extract data frames from phyloseq object 

# count table: ASVs
SoilASVs <- as.data.frame(Psoil_filt@otu_table)

# metadata or Environment
SoilMeta <- read.csv("~/Grad_School/Maestria/Processed_data/metadata.csv")
rownames(SoilMeta) <- SoilMeta$Sample_ID

adonis2(SoilASVs ~ Treatment * Plant_Type, 
        data=SoilMeta, permutations = 999,
        method = "bray")

adonis2(SoilASVs ~ Treatment * Plant_Type, 
        data=SoilMeta, permutations = 999,
        method = "cao")

# dataTables[["norm"]][["genus"]]

#### Phyloseq analyses ####

load("Psoil_filt.Rdata")

### Subset phyla

Bact <- subset_taxa(Psoil_filt, Phylum == "Bacteroidota")
Gemma <- subset_taxa(Psoil_filt, Phylum == "Gemmatimonadota")
Sulf <- subset_taxa(Psoil_filt, Phylum == "Desulfobacterota")
Arm <- subset_taxa(Psoil_filt, Phylum == "Armatimonadota")

Signf_phyla <- merge_phyloseq(Bact,Gemma,Sulf,Arm)
Signf_phyla <- tax_glom(Signf_phyla,taxrank = "Phylum", NArm = FALSE)

# Prepare data frame

Physig <- as.data.frame(Signf_phyla@otu_table)
colnames(Physig) <- c("Bacteroidota","Gemmatimonadota",
                      "Desulfobacterota","Armatimonadota")
row.names(s_meta) <- s_meta[,1]
Physig <- cbind(s_meta, Physig)

# Check for normality
hist(Physig$Bacteroidota)
hist(Physig$Gemmatimonadota) # yes
hist(Physig$Desulfobacterota)
hist(Physig$Armatimonadota)

shapiro.test(Physig$Bacteroidota)
shapiro.test(Physig$Gemmatimonadota)
shapiro.test(Physig$Desulfobacterota)
shapiro.test(Physig$Armatimonadota)

## ANOVAs

sulf_anova <- aov(Desulfobacterota ~ Plant_Type,  
                 data = Physig) # one way there is difference
Anova(sulf_anova)

Physig <- Physig %>%
  add_column(sc_bact = scale(Physig$Bacteroidota))
hist(Physig$log_bact)
  
bact_anova <- aov(sc_bact ~ Plant_Type,  
                  data = Physig) # no diff
Anova(bact_anova)

ggplot(Physig, aes(x = Treatment, y = Armatimonadota)) +
  geom_boxplot()

gem_anova <- aov(log(Gemmatimonadota) ~ Plant_Type,  
                  data = Physig) # marginal difference
Anova(gem_anova)

arm_anova <- aov(Armatimonadota ~ Treatment,  
                 data = Physig) # no
Anova(arm_anova)

#### Metabolism stats ####

row.names(metabolism) <- metabolism[,1]
metabolism <- metabolism[,-1]
met_stats <- as.matrix(t(metabolism))
met_stats <- cbind(s_meta,met_stats)

shapiro.test(met_stats$Sugars.fermentor)

biom_anova <- aov(Sugars.fermentor ~ Plant_Type* Treatment,  
                 data = met_stats) # no
Anova(biom_anova)

## Check loading plots for phylum from metagenassist

#### Citrate and malate ####

## citrate
leveneTest(Citrate_sum ~ Plant_type * Treatment, 
            data = Exud_OG)
shapiro.test(sqrt(Exud_OG$Citrate_sum))

cit_anova <- aov(sqrt(Citrate_sum) ~ Plant_type * Treatment,  
                  data = Exud_plants) # no
Anova(cit_anova)

## Malate
leveneTest(Malate_sum ~ Plant_type * Treatment, 
           data = Exud_OG)
shapiro.test(log(Exud_OG$Malate_sum))

mal_anova <- aov(Malate_sum ~ Plant_type * Treatment,  
                 data = Exud_OG) 
Anova(mal_anova)
mal_tuk <- TukeyHSD(mal_anova)

## Glucose

shapiro.test(Exud_Top$Exud_27)

glu_anova <- aov(Exud_27 ~ Plant_type * Treatment,  
                 data = Exud_Top) 
Anova(glu_anova)
glu_tuk <- TukeyHSD(glu_anova)


#### TOC & MBC ####

shapiro.test(Graphs_OG$TOC_m)
hist(Graphs_OG$TOC_m)

TOC_anova <- aov(TOC_m ~ Plant_Type * Treatment,  
                 data = Data_OG) 
Anova(TOC_anova)
TOC_tuk <- TukeyHSD(TOC_anova)

# MBC
shapiro.test(Graphs_OG$MBC_m)
hist(Graphs_OG$MBC_m)

MBC_anova <- aov(MBC_m ~ Plant_Type * Treatment,  
                 data = Data_OG) 
Anova(MBC_anova)
MBC_tuk <- TukeyHSD(MBC_anova)

