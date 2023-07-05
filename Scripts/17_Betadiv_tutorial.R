#### Betapart and Vegan tutorial ####

### Libraries ###
library("ggplot2"); packageVersion("ggplot2")
library("vegan")
library("betapart")

#### Upload data ####

comm <- read.csv("communities.csv",row.names=1,sep=";")
data(package = "vegan") ## names of data sets in the package
data(dune) # Vegetation and Environment in Dutch Dune Meadows
str(dune) #a data frame of observations of 30 species at 20 sites

# sites are rows and species are columns

#### Beta diversity and Betapart ####

groups <- factor(c(rep(1,3), rep(2,3)), 
                 labels = c("undisturbed","disturbed"))

### Presence absecence table ###

# Some indices do not take into consideration abundance only a 
# presence or absence matrix

presabs<-ifelse(comm>0,1,0)

### Jaccard and Sorensen ###

distj<-beta.pair(presabs, index.family="jaccard")
dists<-beta.pair(presabs, index.family="sorensen")

### Plot and ANOVA ### 

bd<-betadisper(distj[[3]],groups) # object only useful when plotted 
bds<-betadisper(dists[[3]],groups)
plot(bds)

anova(bds)

### Including abundance ###

dist<-bray.part(comm) # how to plot it?

#### Distance (dissimilarity) measures on Vegan ####

par(mfrow = c(1, 2))
bray <-  vegdist(dune, "bray") 
gower <- vegdist(dune, "gower")
hist(bray, xlim = range(0.0,1.0))
hist(gower, xlim = range(0.0,1.0))

