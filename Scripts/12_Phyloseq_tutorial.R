#### Phyloseq tutorial ####

### Libraries ###

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("ape")
library("plyr"); packageVersion("plyr")
library("cluster"); packageVersion("cluster")

## Set theme
theme_set(theme_bw())


#### Creating a phyloseq object ####

# Check how to put my files in this format:
# OTU table: abundances per sample of each sequence 
# tax_table: the taxa ID 
# Phy-tree: it seems to be the only one missing check software

### pretend OTU table ###

otumat <-  matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
otumat

### Add sample and OTU names ###

rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
otumat

### Pretend taxonomy table ###

taxmat <-  matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat

### Join it into a phyloseq object ###

## This can be done in one step with the phyloseq formula
OTU <-  otu_table(otumat, taxa_are_rows = TRUE)
TAX <-  tax_table(taxmat)

physeq <-  phyloseq(OTU, TAX)

## Check abundances 

plot_bar(physeq, fill = "Family")

### Add sample information ###

sampledata <-  sample_data(data.frame(
  Location = sample(LETTERS[1:4], size=nsamples(physeq), replace=TRUE),
  Depth = sample(50:1000, size=nsamples(physeq), replace=TRUE),
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE
))
sampledata

### Create a random tree ###

random_tree <-  rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

### Merge these new two objects to phyloseq object ### 

physeq1 <-  merge_phyloseq(physeq, sampledata, random_tree)
physeq1

### Or rebuild new object ###

physeq2 <-  phyloseq(OTU, TAX, sampledata, random_tree)

### Are they identical? ###
  
identical(physeq1, physeq2) # use this object to compare taxa tables 

#### Plots ####

### Tree-plots ###

plot_tree(physeq1, color="Location", 
          label.tips="taxa_names", 
          ladderize="left", plot.margin=0.3)

plot_tree(physeq1, color="Depth", 
          shape="Location", label.tips="taxa_names", 
          ladderize="right", plot.margin=0.3)

### Heat maps ###

plot_heatmap(physeq1)

plot_heatmap(physeq1, taxa.label="Phylum")

#### Loading sample data ####

data(GlobalPatterns)
data(esophagus)
data(enterotype)
data(soilrep)

example(enterotype, ask=FALSE)

#### Reviewing data ####

### Check data ###

ntaxa(GlobalPatterns)
nsamples(GlobalPatterns)
sample_names(GlobalPatterns)[1:5]
rank_names(GlobalPatterns)
sample_variables(GlobalPatterns)
otu_table(GlobalPatterns)[1:5, 1:5]
tax_table(GlobalPatterns)[1:5, 1:4]
phy_tree(GlobalPatterns)
taxa_names(GlobalPatterns)[1:10]

### Short phylogenetic trees ###

myTaxa <-  names(sort(taxa_sums(GlobalPatterns), decreasing = TRUE)[1:10])
ex1 <-  prune_taxa(myTaxa, GlobalPatterns)
plot(phy_tree(ex1), show.node.label = TRUE)


plot_tree(ex1, color = "SampleType", label.tips = "Phylum", 
          ladderize = "left", justify = "left" , size = "Abundance")

#### Pre-processing ####

## Add abundance and filter taxa

GPr  <-  transform_sample_counts(GlobalPatterns, function(x) x / sum(x) )
GPfr <-  filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

## Difference between prune and subset

GP.chl <-  subset_taxa(GlobalPatterns, Phylum=="Chlamydiae") # uses data from other tables in the phyloseq object
GP.chl <-  prune_samples(sample_sums(GP.chl)>=20, GP.chl) # only works with info from the abundances or ASV table

## Merging 

GP.chl.merged <-  merge_taxa(GP.chl, taxa_names(GP.chl)[1:5]) # merges the first 5 OTUs in the Chlamydiae-only dataset

## Agglomeration functions

gpsfbg <-  tax_glom(GP.chl, "Family")
plot_tree(gpsfbg, color="SampleType", shape="Class", size="abundance")

## Transform abundance values

transform_sample_counts(GP.chl, function(OTU) OTU/sum(OTU)) # fractional abundance

## Filter taxa with small abundance

# Remove taxa not seen more than 3 times in at least 20% of the samples. 
# This protects against an OTU with small mean & trivially large C.V.

GP <-  filter_taxa(GlobalPatterns, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

## Add human vs. non human categorical variable

sample_data(GP)$human <-  factor( get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue"))
# adding this variable to the sample data data frame in the phyloseq object

## Standardize abundances to the median sequencing depth

total <-  median(sample_sums(GP))
standf <-  function(x, t=total) round(t * (x / sum(x))) # function
gps <- transform_sample_counts(GP, standf)

## Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation

gpsf <-  filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)

## Subset to Bacteroidetes

gpsfb <-  subset_taxa(gpsf, Phylum=="Bacteroidetes") 
# use subset instead of prune because Phylum is from taxa table

#### Graphics with subset data ####

# Can use ggplot instead

title <- "plot_bar; Bacteroidetes-only"

plot_bar(gpsfb, "SampleType", "Abundance", title=title)
plot_bar(gpsfb, "SampleType", "Abundance", "Family", title=title)
plot_bar(gpsfb, "Family", "Abundance", "Family", 
         title=title, facet_grid="SampleType~.")

#### The distance function in phyloseq ####

# Using enterotype data: does not have a phylogenetic tree

### Preprocessing data ###

# Remove unassigned values 
enterotype <- subset_taxa(enterotype, Genus != "-1")

## Distance methods ##

dist_methods <- unlist(distanceMethodList)
print(dist_methods) # Jensen-Shannon Divergence, JSD and betadiversity

# These require tree
dist_methods[(1:3)]
# This is the user-defined method:
dist_methods["designdist"]

# Remove them from the vector
dist_methods <- dist_methods[-(1:3)]
dist_methods <- dist_methods[-which(dist_methods=="ANY")]

## for loop to test all methods ##

plist <- vector("list", length(dist_methods))

names(plist) <-  dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(enterotype, method=i)
  # Calculate ordination
  iMDS  <- ordinate(enterotype, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(enterotype, iMDS, color="SeqTech", shape="Enterotype")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

### Graphs ###

## Color by sequencing technology

df <-  ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p <-  ggplot(df, aes(Axis.1, Axis.2, color=SeqTech, shape=Enterotype))
p <-  p + geom_point(size=3, alpha=0.5)
p <-  p + facet_wrap(~distance, scales="free")
p <-  p + ggtitle("MDS on various distance metrics for Enterotype dataset")
p

## Color by Enterotype

df <-  ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p <-  ggplot(df, aes(Axis.1, Axis.2, color=Enterotype, shape=SeqTech))
p <-  p + geom_point(size=3, alpha=0.5)
p <-  p + facet_wrap(~distance, scales="free")
p <- p + ggtitle("MDS on various distance metrics for Enterotype dataset")
p

# Helpful to compare results depending on method

print(plist[["jsd"]])
print(plist[["jaccard"]])
print(plist[["bray"]])
print(plist[["gower"]])

      
# Jensen Shannon (JSD), jaccard, and Bray show similar results and seem to be the mos popular ones
# investigate gower, and MDS

#### Gap statistics ####

# The purpose of this section is to know in how many clusters do my data
# separate into, that may help showing differences between treatments

## Ordination

exord <-  ordinate(enterotype, method="MDS", distance="jsd")

### Gap ###

pam1 <-  function(x, k) 
  {list(cluster = pam(x,k, cluster.only=TRUE))}
## Pam is  a function short for partitioning around methods
# cluster.only= only the clustering will be computed and returned

x <-  phyloseq:::scores.pcoa(exord, display="sites") # scores for sites or species 
# no documentation for phyloseq:::scores.pcoa info based on scores function

# gskmn = clusGap(x[, 1:2], FUN=kmeans, nstart=20, K.max = 6, B = 500)
gskmn <-  clusGap(x[, 1:2], FUN=pam1, K.max = 6, B = 50) 
# x numeric matrix or dataframe
# B = Monte Carlo (“bootstrap”) samples
# clusGap goes together with maxSE 


## Wrapper functions

gap_statistic_ordination = function(ord, FUNcluster, type="sites", K.max=6, axes=c(1:2), B=500, verbose=interactive(), ...){
  require("cluster")
  #   If "pam1" was chosen, use this internally defined call to pam
  if(FUNcluster == "pam1"){
    FUNcluster = function(x,k) list(cluster = pam(x, k, cluster.only=TRUE))     
  }
  # Use the scores function to get the ordination coordinates
  x <-  phyloseq:::scores.pcoa(ord, display=type) # type from the input of the function
  if(is.null(axes)){axes = 1:ncol(x)}
  #   Finally, perform, and return, the gap statistic calculation using cluster::clusGap  
  clusGap(x[, axes], FUN=FUNcluster, K.max=K.max, B=B, verbose=verbose, ...)
}
## this is a function merging pam, scores, and clusGap or it can be done separately?


### Plot Gaps ###

## Function to plot methods 
plot_clusgap = function(clusgap, title="Gap Statistic calculation results"){
  require("ggplot2")
  gstab <-  data.frame(clusgap$Tab, k=1:nrow(clusgap$Tab)) # Tab is an element from he list obtained from clusGap
  p <-  ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
  p <-  p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim))
  p <-  p + ggtitle(title)
  return(p)
}

gs = gap_statistic_ordination(exord, "pam1", B=50, verbose=FALSE)
print(gs, method="Tibs2001SEmax")

# "Tibs2001SEmax", Tibshirani et al's recommendation of determining
# the number of clusters from the gap statistics

plot_clusgap(gs) # function from phyloseq

# Base function
plot(gs, main = "Gap statistic for the 'Enterotypes' data")
mtext("Looks like 4 clusters is best, with 3 and 5 close runners up.")

# These plots show how many clusters have the larger gap between them
# therefore, you shoose the number of cluesters more separated between them

