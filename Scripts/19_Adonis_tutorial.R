#### Adonis function ####

### Libraries

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("vegan")
library("dplyr")
library("tidyverse")
library("grid")
library("gridExtra") # helps with plot display

### Base Tutorial ###

# sample data
data(dune)
data(dune.env)

#adonis(formula, data, permutations = 999, method = "bray",
#       strata = NULL, contr.unordered = "contr.sum",
#       contr.ordered = "contr.poly", parallel = getOption("mc.cores"), ...)

# Basic formula
sample_ad <- adonis(dune ~ Management*A1, data=dune.env, permutations=99)

## Using strata ##

dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
# From this line of code:
# expan.grid is to create a data frame with supplied vector factors
# gl chooses random numbers to fill in the vector
# factor chooses either or from the variables you provide

Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
# with function modifies a column from a data frame given an expression
# rnorm is for normal distribution

total <- Agropyron + Schizachyrium

dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium) 

mod <- metaMDS(Y) # to perform nonmetric multidimensional scale
plot(mod)

### Hulls show treatment NO3
with(dat, ordihull(mod, group=NO3, show="0"))
with(dat, ordihull(mod, group=NO3, show="10", col=3))

### Spider shows fields
with(dat, ordispider(mod, group=field, lty=3, col="red"))

# all the with code is adding hulls (or polygons) to the plot mod 
# basically the plot of the metaMDS
# ordihull adds the polygons, in this case for treatment
# and ordispider adds the fields in spider or line form 

### Correct hypothesis test (with strata)

adstrat <- adonis(Y ~ NO3, data=dat, strata=dat$field, perm=999) # not significant
# in this case you stratified (or separated by field) to do the analysis
# for my data I can try to stratified by treatment maybe?

### Incorrect (no strata)
adstratwrong <- adonis(Y ~ NO3, data=dat, perm=999)

#### Tutorial with explanation ####

# Creating fake data:  three sets of sites (30 sites, 10 species) 
# for each of three treatments.

set.seed(123456789)
num<-30
disp.a<-5
sites.a<-data.frame(sp.a=rnbinom(num,mu = 40, size = disp.a),
                    sp.b=rnbinom(num,mu = 60, size = disp.a),
                    sp.c=rnbinom(num,mu = 50, size = disp.a),
                    sp.d=rnbinom(num,mu = 70, size = disp.a),
                    sp.e=rnbinom(num,mu = 10, size = disp.a),
                    sp.f=rnbinom(num,mu = 180, size = disp.a),
                    sp.g=rnbinom(num,mu = 100, size = disp.a),
                    sp.h=rnbinom(num,mu = 80, size = disp.a),
                    sp.i=rnbinom(num,mu = 40, size = disp.a),
                    sp.j=rnbinom(num,mu = 50, size = disp.a))

disp.b<-50
sites.b<-data.frame(sp.a=rnbinom(num,mu = 40, size = disp.a),
                    sp.b=rnbinom(num,mu = 60, size = disp.b),
                    sp.c=rnbinom(num,mu = 50, size = disp.a),
                    sp.d=rnbinom(num,mu = 70, size = disp.b),
                    sp.e=rnbinom(num,mu = 10, size = disp.a),
                    sp.f=rnbinom(num,mu = 180, size = disp.a),
                    sp.g=rnbinom(num,mu = 100, size = disp.b),
                    sp.h=rnbinom(num,mu = 80, size = disp.a),
                    sp.i=rnbinom(num,mu = 40, size = disp.b),
                    sp.j=rnbinom(num,mu = 50, size = disp.a))

disp.c<-200
sites.c<-data.frame(sp.a=rnbinom(num,mu = 40, size = disp.a),
                    sp.b=rnbinom(num,mu = 60, size = disp.b),
                    sp.c=rnbinom(num,mu = 50, size = disp.c),
                    sp.d=rnbinom(num,mu = 70, size = disp.b),
                    sp.e=rnbinom(num,mu = 10, size = disp.c),
                    sp.f=rnbinom(num,mu = 180, size = disp.a),
                    sp.g=rnbinom(num,mu = 100, size = disp.b),
                    sp.h=rnbinom(num,mu = 80, size = disp.c),
                    sp.i=rnbinom(num,mu = 40, size = disp.b),
                    sp.j=rnbinom(num,mu = 50, size = disp.c))

all.sites<-rbind(sites.a,sites.b,sites.c)


trt<-rep(c("C","H","L"),each=nrow(sites.a))

# Running an NMDS

all.mds <- metaMDS(all.sites)

  data.scores <- as.data.frame(scores(all.mds))
  data.scores$site <- rownames(data.scores)  
  data.scores$grp<-trt

# how to create an NMDS plot using ggplot
ggplot(data=data.scores) + 
  stat_ellipse(aes(x=sites.NMDS1,y=sites.NMDS2,colour=trt),level = 0.50) +
  geom_point(aes(x=sites.NMDS1,y=sites.NMDS2,shape=trt,colour=trt),size=4) + 
  theme_minimal() # using sites do not know why

# running adonis

adon.results<-adonis2(all.sites ~ trt, method="bray",perm=999)
print(adon.results)

### Why do we see significant differences? ###

# uses centroid variance 

## Bray-Curtis distances between samples
dis <- vegdist(all.sites)

## Calculate multivariate dispersions
mod <- betadisper(dis, trt)
mod

## Visualizing the multivariate homogeneity of group dispersions

# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(grps=rownames(mod$centroids),data.frame(mod$centroids))
vectors<-data.frame(group=mod$group,data.frame(mod$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

# create the convex hulls of the outermost points
grp1.hull<-seg.data[seg.data$group=="C",1:3][chull(seg.data[seg.data$group=="C",2:3]),]
grp2.hull<-seg.data[seg.data$group=="H",1:3][chull(seg.data[seg.data$group=="H",2:3]),]
grp3.hull<-seg.data[seg.data$group=="L",1:3][chull(seg.data[seg.data$group=="L",2:3]),]
all.hull<-rbind(grp1.hull,grp2.hull,grp3.hull)

# display points: First points (black symbols) and the centroids (red symbols)
# Adding vector segments: geom_segment in ggplot
# Hulls: using geom_polygon

panel.a<-ggplot() +
  geom_polygon(data=all.hull[all.hull=="C",],aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[1:30,],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[1,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=16) + 
  geom_point(data=seg.data[1:30,], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=16) +
  labs(title="Control",x="",y="") +
  coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  theme_minimal() + 
  theme(legend.position="none")

panel.b<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="H",],aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[31:60,],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[2,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=17) + 
  geom_point(data=seg.data[31:60,], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=17) +
  labs(title="High",x="",y="") +
  coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  theme_minimal() + 
  theme(legend.position="none")

panel.c<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="L",],aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[61:90,],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) +
  geom_point(data=centroids[3,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=15) + 
  geom_point(data=seg.data[61:90,], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15) + 
  labs(title="Low",x="",y="") +
  coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  theme_minimal() + 
  theme(legend.position="none")

panel.d<-ggplot() + 
  geom_polygon(data=all.hull,aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data,aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[,1:3], aes(x=PCoA1,y=PCoA2,shape=grps),size=4,colour="red") + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2,shape=group),size=2) +
  labs(title="All",x="",y="") +
  coord_cartesian(xlim = c(-0.2,0.2), ylim = c(-0.25,0.2)) +
  theme_minimal() + 
  theme(legend.position="none")

grid.arrange(panel.a,panel.b,panel.c,panel.d,nrow=1)
