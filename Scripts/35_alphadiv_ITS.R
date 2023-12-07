#### Alpha Diversity ITS ####

### Libraries ###
library("phyloseq"); packageVersion("phyloseq")
library("breakaway") # chao-bunge
library("microbiome")
library("ggplot2"); packageVersion("ggplot2")
library("RColorBrewer")
library("tidyverse")
library(dplyr)
library(car)


### Load data ###

load("ITS_filt.RData")
#### Calculate alpha diversity using phyloseq ####

alphadiv <- estimate_richness(ITS_filtered, 
                              split = TRUE, 
                              measures = NULL)

# Separate metadata for alpha diversity table
samplemeta <- as.data.frame(ITS_filtered@sam_data) 
alphadiv <- cbind(samplemeta,alphadiv) # join alpha indices to metadata

## Chao-bunge index ##
chao1_bunge(ITS_filtered) # breakaway package

### Adding Species evenness ###

Evenness <- evenness(ITS_filtered, index = 'all',
                     zeroes = TRUE, detection = 0) # microbiome package
alphadiv <- cbind(alphadiv,Evenness)

# Save file
write.csv(alphadiv, "~/Grad_School/Maestria/Processed_data/Alphadiv_ITS.csv")

#### Alpha diversity plot ####

ITS_filtered@sam_data$Plant_Type <- factor(ITS_filtered@sam_data$Plant_Type, 
                                           levels = c("ptxDNt-36", "Wild Type", "Bulk Soil"))
ITS_filtered@sam_data$Treatment <- factor(ITS_filtered@sam_data$Treatment, 
                                          levels = c("Low P", "Pi", "Phi", "Pi/Phi mix"))

# plot richness with microbiome
PR <- plot_richness(ITS_filtered, x="Treatment", 
                           measures=c("Observed","Shannon","Simpson","Chao1"), 
                           color="Plant_Type")+
    scale_fill_manual(breaks=c("ptxDNt-36","Wild Type","Bulk Soil"),
                      labels=c("ptxDNt-36","Wild Type","Bulk soil"))+
    scale_color_manual(values = c("steelblue","lightpink","burlywood3"), 
                       name = "Plant Type",
                       labels= c("ptxDNt-36","Wild Type","Bulk soil")) +
    theme_bw()
  
PR$layers <- PR$layers[-1]
PR <- PR + geom_point(size=4, alpha=0.7)
PR
  
## Pielou

ggplot(alphadiv, aes(x=Treatment, y=pielou, color = Plant_Type)) +
  geom_point(size = 4, alpha=0.7) +
  scale_color_manual(values = c("burlywood3","steelblue","lightpink"))+
  theme_bw()


#### Statistical analysis ####

alphadiv <- read.csv("~/Grad_School/Maestria/Processed_data/Alphadiv_ITS.csv")
row.names(alphadiv) <- alphadiv[,1]
alphadiv <- alphadiv[,-1]

## Check data normalization ##

# histograms
for (i in 5:18) { 
  hist(alphadiv[,i],
       main = i)
}
# shapiro test
for(i in 5:18){
  shapiro <- shapiro.test(alphadiv[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}

shapiro.test(alphadiv$Shannon)

### Observed ###

## ANOVA

ad_anova <- aov(Shannon ~ Plant_Type * Treatment, 
                   data = alphadiv) # synergistic effect
Anova(ad_anova) # signifcant in plant:treatment

# Test homogeneity of variance
plot(ad_anova, 1) # check outliers
leveneTest(Observed ~ Plant_Type * Treatment, 
           data = alphadiv) # higher than 0.05 so it meets ANOVA assumption

## Tukey HSD 

ad_Tukey <- TukeyHSD(ad_anova) # bulk soil Pi/Phi is statistically different from several samples

### Chao1 ###

## Anova
ad_anova <- aov(Chao1 ~ Plant_Type * Treatment, 
                data = alphadiv) # synergistic effect
Anova(ad_anova) # significant in plant:treatment

# Test homogeneity of variance
plot(ad_anova, 1) # check outliers
leveneTest(Observed ~ Plant_Type * Treatment, 
           data = alphadiv) # higher than 0.05 so it meets ANOVA assumption

## Tukey HSD 

ad_Tukey <- TukeyHSD(ad_anova) # bulk soil Pi/Phi is statistically different from several samples

### Shannon ###

## Anova
ad_anova <- aov(Shannon ~ Plant_Type * Treatment, 
                data = alphadiv) # synergistic effect
Anova(ad_anova) # not significant

### Simpson ###

## Anova
ad_anova <- aov(Simpson ~ Plant_Type * Treatment, 
                data = alphadiv) # synergistic effect
Anova(ad_anova) # not significant 

### Pielou ###

## Anova
ad_anova <- aov(pielou ~ Plant_Type * Treatment, 
                data = alphadiv) # synergistic effect
Anova(ad_anova) # not significant

