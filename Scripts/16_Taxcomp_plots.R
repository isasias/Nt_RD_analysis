#### Taxonomic composition plots ####

### Libraries ###

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("vegan")
library("patchwork")
library("betapart")
library("dplyr")
library("remotes") # download remote libraries
library("devtools")
library("ranacapa")
library("breakaway")
library("microbiome") # data analysis and visualisation
library("wesanderson")
library("ggforce")

## Data
load("Psoil_filt.Rdata")
load("Psoil_rel.RData")

## Relative abundance 
# Subset top phyla
Top_phyla <- tax_glom(Psoil_rel,taxrank = "Phylum", NArm = FALSE)

Top_phyla <- prune_taxa(names(sort(taxa_sums(Top_phyla),TRUE)[1:5]), Top_phyla)

sample_labs <- c("WT LowP", "36ptXD LowP 1","WT Control",
                 "36ptXD Control", "36ptXD Phi", "WT Pmix",
                 "36ptXD Pmix", "No plant")

#"WT LowP 1","WT LowP 2","WT LowP 3",
#                "36ptXD LowP 1","36ptXD LowP 2","36ptXD LowP 3",
#               "WT Control 1"

sorder <- c("A1","B1","C1","A2","B2","C2","A3","B3","C3",
            "A4","B4","C4","A5","B5","C5","A6","B6","C6",
            "A7","B7","C7","D1","D2","D3","D4")
# Top phyla barplot

Phyla_plot <-plot_composition(Top_phyla, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling2"),name = "Phylum", labels = c("Actinobacteria", "Proteobacteria", "Bacteroidetes","Chloroflexi","Acidobacteria"))+
  guides(x=guide_axis(angle = 35))+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        axis.ticks.x = element_blank())


#### Phyla heatmap

```{r}
Phyla_soil <- tax_glom(Psoil_filt,taxrank = "Phylum", NArm = FALSE)

Heatmap <- plot_heatmap(Phyla_soil, sample.label = "Sample_ID",
                        method = "NMDS", distance = "bray", 
                        taxa.label = "Phylum",
                        taxa.order = names(sort(taxa_sums(Phyla_soil))),
                        low="rosybrown1", high="darkslateblue",
                        na.value ="oldlace",sample.order = sorder)+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

```
