#### Beta diversity trials ####

### Libraries ###
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("vegan")
library("patchwork")
library("betapart")
library("dplyr")
library("remotes") # download remote libraries
library("microbiome") # data analysis and visualisation
library("ggforce") # add polygons to pcas

## Data
load("Psoil_filt.Rdata")
load("Psoil_rel.RData")

bray_bdiv <- phyloseq::distance(Psoil_filt, method = "bray", type = "sample") # check how it looks with taxa type
bray_ord <- ordinate(Psoil_filt,"PCoA", distance = bray_bdiv)

p_bray <- plot_ordination(Psoil_filt, bray_ord, "samples", 
                          color = "Plant_Type",
                          shape = "Treatment")+
  scale_color_manual(values = c("steelblue","lightpink","palegreen3"), name = "Plant Type")+
  theme_bw()+
  geom_point(size=3)

shape_names <- c("Control","Low P","Phi", "Pi/Phi mix")
p.shape <- 15:(15 + length(shape_names) - 1)
names(p.shape) <- shape_names
p.shape["Treatments"] <- 4

p_bray <- p_bray +
  scale_shape_manual(values = p.shape)


p_bray <- p_bray +
  geom_polygon(aes(fill=Plant_Type),alpha = 0.2)+
  scale_fill_manual(values = c("steelblue","lightpink","palegreen3"), name = "Plant Type")
