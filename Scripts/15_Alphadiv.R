#### Libraries ####
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("patchwork") # graph placement
library("dplyr")
library("remotes") # download remote libraries
library("devtools")
library("breakaway") # Chao1
library("microbiome") # data analysis and visualisation
library(tidyverse)

## Data
load("Psoil_filt.Rdata")
load("Psoil_rel.RData")

# Separate metadata for alpha diversity table

samplemeta <- as.data.frame(Psoil_filt@sam_data) 

### Estimating species richness

alphadiv <- estimate_richness(Psoil_filt, split = TRUE, measures = NULL)
chao_bunge(Psoil_filt, cutoff = 55)

alphadiv <- cbind(samplemeta,alphadiv) # join alpha indices to metadata

## Add Chao-Bunge to table 

# Richness
Chaob <- c(5962,5859,5984,5937,6019,
           6212,6078,5880,5618,5973,
           6085,5758,6214,6082,5979,
           5736,5766,6098,5859,5761,
           5352,5734,5961,5575,5868)

# Standard Error
Chaose <- c(53.59,44.1,52.23,53.06,55.52,
            56.34,54.68,57.12,55.98,48.61,
            59.12,60.9,51.6,48.08,58.74,
            56.14,64.13,53,53.13,62.83,
            64.32,56.86,59.87,59.87,57.07)

alphadiv <- alphadiv %>%
  add_column(Chaob) %>%
  add_column(Chaose)


### Adding Species evenness

Evenness <- evenness(Psoil_filt, index = 'all',zeroes = TRUE, detection = 0) 
# evennes from microbiome package

alphadiv <- cbind(alphadiv,Evenness)


write.csv(alphadiv, "~/Grad_School/Maestria/Processed_data/Alphadiversity.csv")


#### Alpha diversity plots

# This section will be done with the basic phyloseq function but also I create a graph comparing Chao1 and Chao-Bunge indices, and an evenness plot.

# Indices plot

Psoil_filt@sam_data$Plant_Type <- factor(Psoil_filt@sam_data$Plant_Type, levels = c("ptxD-36", "Wild Type", "No plant"))

PR <- plot_richness(Psoil_filt, x="Treatment", 
                    measures=c("Chao1"), 
                    color="Plant_Type")+
  scale_fill_manual(breaks=c("ptxD-36","Wild Type","No plant"),labels=c("ptxD-36","Wild Type","No plant"))+
  scale_color_manual(values = c("steelblue","lightpink","palegreen3"), name = "Plant Type") +
  theme_bw()

PR$layers <- PR$layers[-1]
PR <- PR + geom_point(size=4, alpha=0.7) 


# Chao1 and Chao-Bunge comparison

Chao1g <- ggplot(alphadiv, aes(x=Treatment, y=Chao1, color = Plant_Type)) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c("lightpink","palegreen3","steelblue"))+
  geom_errorbar(aes(ymin=Chao1-se.chao1, ymax=Chao1+se.chao1), width=.05)+
  labs(title= 'Chao1 Index', x= ' ', y= '') +
  theme_bw()

Chaobg <- ggplot(alphadiv, aes(x=Treatment, y=Chaob, color = Plant_Type)) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c("lightpink","palegreen3","steelblue"))+
  geom_errorbar(aes(ymin=Chaob-Chaose, ymax=Chaob+Chaose), width=.05)+
  labs(title= 'Chao-Bunge Index', x= ' ', y= '') +
  theme_bw()

(Chao1g | Chaobg)

# Evenness comparison

Shannoneven <- ggplot(alphadiv, aes(x=Treatment, y=pielou, color = Plant_Type)) +
  geom_point(size=4, alpha=0.7) +
  scale_color_manual(values = c("lightpink","palegreen3","steelblue"))+
  labs(title= 'Shannon evenness', x= ' ', y= '') +
  theme_bw()
Shannoneven <- Shannoneven + theme(legend.position = "none")


Simpsoneven <- ggplot(alphadiv, aes(x=Treatment, y=simpson, color = Plant_Type)) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c("lightpink","palegreen3","steelblue"))+
  labs(title= 'Simpson evenness', x= ' ', y= '') +
  theme_bw()
Simpsoneven <- Simpsoneven + theme(legend.position = "none")

(PR | Shannoneven)


#### Statistical analysis for Alpha Diversity

# I will check if there are any significant differences between treatments and plant type using ANOVAS.

### Observed 

obs_anova <- aov(Observed ~ Plant_Type * Treatment, 
                 data = alphadiv)

plot(obs_anova, 1) # check outliers
leveneTest(Observed ~ Plant_Type * Treatment, 
           data = alphadiv) # higher than 0.05 so it meets ANOVA assumption

Anova(obs_anova) # not significant difference

### Chao 1

Chao1_anova <- aov(Chao1 ~ Plant_Type * Treatment, 
                   data = alphadiv)

plot(Chao1_anova, 1) # check outliers
leveneTest(Chao1 ~ Plant_Type * Treatment, 
           data = alphadiv) # higher than 0.05 so it meets ANOVA assumption

Anova(Chao1_anova) # not significant difference

### Shannon

Shannon_anova <- aov(Shannon ~ Plant_Type * Treatment, 
                     data = alphadiv)

plot(Shannon_anova, 1) # check outliers
leveneTest(Shannon ~ Plant_Type * Treatment, 
           data = alphadiv) # higher than 0.05 so it meets ANOVA assumption

Anova(Shannon_anova) # not significant difference

### Simpson

Simp_anova <- aov(Simpson ~ Plant_Type * Treatment, 
                  data = alphadiv)

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

### Evenness Simpson

Evsim_anova <- aov(simpson ~ Plant_Type * Treatment, 
                   data = alphadiv)

plot(Evsim_anova, 1) # check outliers
leveneTest(simpson ~ Plant_Type * Treatment, 
           data = alphadiv) # higher than 0.05 so it meets ANOVA assumption

Anova(Evsim_anova) # not significant difference

