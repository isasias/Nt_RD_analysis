#### PCAs tutorial ####

### Libraries ###

library(ChemometricsWithR)
library(plyr)
library(car)
library(maptools)
library(rgeos)

library(plotly)
library(ggfortify)

#### Data ####

## Tutorial 1
data(wines, package = "kohonen")
vino <- wines
vint <- vintages

## Tutorial 2 
df <- iris[1:4]

#### Transform and scale data ####

# Log transformation and autoscaling to set all variables on same scale

## log transform, pseudo scale (not entirely necessary_
vino_log <- log(vino)

## Auto scale
vino_log_scale <- scale(vino_log)

#### PCA function ####

# Tutorial 1
vino_PCA <- prcomp(vino_log_scale, center=FALSE)
# Tutorial 2 
pca_res <- prcomp(df, scale. = TRUE)


#### PCA visualization ggplot2 tutorial 2 ####

p <- autoplot(pca_res) ## just to plot PCA... I will use this 
ggplotly(p) # ggplotly makes it interactive

p <- autoplot(pca_res, data = iris, # including the data frame to add species as color  
              colour = 'Species', label = TRUE, label.size = 3)
ggplotly(p)

# removing shape so only the numbers are visible and keeping labels 
p <- autoplot(pca_res, data = iris, colour = 'Species', shape = FALSE, label.size = 3)

ggplotly(p)

# Displaying eigen vectors or loadings 
p <- autoplot(pca_res, data = iris, 
              colour = 'Species', loadings = TRUE) # check if I can do a loadings plot only 
p <- autoplot(pca_res, data = iris, colour = 'Species',
              loadings = TRUE, loadings.colour = 'blue',
              loadings.label = TRUE, loadings.label.size = 3)

ggplotly(p)

## check how to do loading plots with ggplot2

### Loadings with ggplot2 ###

PCAloadings <- pca_res$rotation

# Scatter plot: it worked!!!

ggplot(data = PCAloadings, 
       aes(x= PC1, y= PC2, 
           label=rownames(PCAloadings)))+
  geom_point(shape=15,color="blue",size=6)+
  geom_text()

#### Base graphics tutorial 1 ####

PCAcolors <- c("#66c2a5","#fc8d62","#8da0cb")[as.integer(vint)]


PCAscores <- vino_PCA$x # scores are x in the list
PCAloadings <- vino_PCA$rotation # loadings are rotation

## Plot 1: PCA scores
plot(PCAscores[,1:2],  # x and y data
     pch=21,           # point shape
     col=PCAcolors,    # point border color
     bg=PCAcolors,     # point color
     cex=1.5,          # point size
     main="Scores"     # title of plot
)
legend("topright",                                # position of legend
       legend=levels(vint),                       # legend display
       pch=21,                                    # point shape
       pt.bg=c("#66c2a5","#fc8d62","#8da0cb"),    # point colors
       pt.cex=1.5,                                # point size
       col = c("#66c2a5","#fc8d62","#8da0cb")    # point border color
) # problem with legend too big 

## Plot 2 :Loading plot
plot(PCAloadings[,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1,               # point size
     main="Loadings"      # title of plot
)
text(PCAloadings[,1:2],             # sets position of labels
     labels=rownames(PCAloadings)   # print labels
) 

### Creating an ellipseplot() function ###

## Customize yaxis range to makes sure axis ticks cover data
## Axes ticks do not always cover data range in R plots - reviewer did not like!

plotat <- function(RANGE) {
  if(length(RANGE) != 2) stop("RANGE argument must have a length of 2")
  if(RANGE[1] > RANGE[2]) stop("First element in RANGE must be smaller than second element")
  prettyres <- pretty(sprintf("%.2f",RANGE[1]):sprintf("%.2f",RANGE[2]), 7)
  while((min(prettyres) < RANGE[1]) == FALSE) {
    prdiff <- prettyres[2] - prettyres[1]
    prettyres[length(prettyres) + 1] <- prettyres[1] - prdiff
    prettyres <- sort(prettyres)
  } 
  while((max(prettyres) > RANGE[2]) == FALSE) {
    prdiff <- prettyres[2] - prettyres[1]
    prettyres[length(prettyres) + 1] <- prettyres[length(prettyres)] + prdiff
    prettyres <- sort(prettyres)    
  }   
  plotticks <- as.numeric(sprintf("%.2f",prettyres))
  plotticks
}

# not a problem with ggplot?

## This is the code for the ellipse function: basically is a 
## very very long function to use in base R plots 

ellipseplot <- function(x, y, factr, 
                        elev=0.95, # Ellipse probability level
                        legpos=c("topright","topleft","bottomleft","bottomleft"), # Legend position
                        pcol=NULL, # manual addition of colors, must meet length of factors
                        cexsize=1, # point size
                        ppch=21, # Point type, must meet length of factors
                        legcexsize=2, # legend font size
                        legptsize=2, # legend point size
                        pbgcol=TRUE,
                        axissize=1, 
                        linewidth=1, 
                        font=1) {
  require(plyr)
  require(car)
  ## Set factor levels
  if(is.factor(factr)) {
    f <- factr
  } else {
    f <- factor(factr, levels=unique(as.character(factr)))
  }
  intfactr <- as.integer(f) # Set integer vector that matches factor levels
  # Checking to make sure length of ppch equals number of factor levels
  if((length(ppch) > 1 & length(unique(intfactr)) != length(ppch))) stop("Can only increase point shape if equal to factor levels")
  ## Get data for ellipses
  edf <- data.frame(LV1 = x, LV2=y, factr = f) # create data frame with data and factor
  ellipses <- dlply(edf, .(factr), function(x) {
    LV1 <- x[,1]
    LV2 <- x[,2]
    dataEllipse(LV1, LV2, levels=elev, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
  })
  ## Get range of x and y data
  xrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,1])), min(x), max(x))))
  yrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,2])), min(y), max(y))))
  
  
  ## Set colors for plots
  if(is.null(pcol) != TRUE) { # If colors are supplied by user
    ptcol <- pcol
    pgcol <- paste(pcol, "7e", sep="") # adds opaqueness
  } else { # Default
    pgcol <- c("#e41a1c7e","#377eb87e","#4daf4a7e","#984ea37e","#807f7d7e") # Defaults at 5 colors
    ptcol <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#807f7d") # For opaqueness
  }
  # Plotting graphic
  plot(x,y, type="n", xlab="", ylab="", main="", xlim=range(xrange), ylim=range(yrange), axes=FALSE)
  axis(1, at=xrange, labels=xrange, cex.axis=axissize,lwd=linewidth, font=font)
  axis(2, las=2, cex.axis=axissize,lwd=linewidth, font=font)
  box(lwd=linewidth, font=font)
  abline(h=0, v=0, col="gray", lty=2) # Adds lines at 0
  legpch <- c() # vector to collect legend pch data
  legcol <- c() # vector to collect legend col data
  ## Adds points, ellipse, and determines color specifications for legend 
  if(pbgcol==TRUE)  {
    for(i in 1:length(unique(intfactr))){
      points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col=ptcol[i], bg=ptcol[i],cex=cexsize)
      polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
      legpch[i] <- ppch[i]
      legcol[i] <- ptcol[i]
    }
  } else {
    for(i in 1:length(unique(intfactr))){
      points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col="black", bg=ptcol[i],cex=cexsize)
      polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
      legpch[i] <- ppch[i]
      legcol[i] <- ptcol[i]       
    }
  }
  ## Legend
  legend(x=legpos, legend=levels(f), pch=legpch, 
         pt.bg=legcol, col=legcol, bty="n", border=FALSE, pt.cex=legptsize, cex=legcexsize)
}   

## Axis legends for PCA output using prcomp() function
## Function 3 

PCAvarAxis <- function(PCA, decimal=1) {
  pcavar <- round((PCA$sdev^2)/sum((PCA$sdev^2)),3)*100   #Calculate % variance explained
  PC1var <- paste("Principal Component 1 (", pcavar[1], "%)", sep="")
  PC2var <- paste("Principal Component 2 (", pcavar[2], "%)", sep="")
  PC3var <- paste("Principal Component 3 (", pcavar[3], "%)", sep="")
  PC4var <- paste("Principal Component 4 (", pcavar[4], "%)", sep="") 
  PC5var <- paste("Principal Component 5 (", pcavar[5], "%)", sep="")     
  return(list(PC1=PC1var, PC2=PC2var, PC3=PC3var, PC4=PC4var, PC5=PC5var))
} 

# This last section was copied and pasted without any modification 
# check if it works for my data

#### Cleaner base plots tutorial 1 ####

explainPCAvar <- PCAvarAxis(vino_PCA)

# Plot 1: with ellipses 
ellipseplot(PCAscores[,1],                          # data for x-axis
            PCAscores[,2],                          # data for y-axis
            vint,                                   # factor with classes
            pcol=c("#66c2a5","#fc8d62","#8da0cb"),  # colors for plotting (must match # of factors)
            pbgcol=FALSE,                           # point borders black?
            cexsize=1.5,                            # size of points 
            ppch=c(21:23),                          # shape of points (must match # of factors)
            legpos="bottomright",                   # position of legend           
            legcexsize=1.5,                         # legend text size
            legptsize=1.5,                          # legend point size 
            axissize=1.5,                           # Set axis text size
            linewidth=1.5                           # Set axis line size
)                         
title(xlab=explainPCAvar[["PC1"]],    # % variance explained on PC1
      ylab=explainPCAvar[["PC2"]],    # % variance explained on PC2 
      main="Scores",                  # Title
      cex.lab=1.5,                    # size of label text
      cex.main=1.5                    # size of title text
)

## Plot 2
plot(PCAloadings[,1:2],   # x and y data
     pch=21,              # point shape
     bg="black",          # point color
     cex=1.5,             # point size
     # type="n",           # does not plot points
     axes=FALSE,          # does not print axes
     xlab="",             # removes x label
     ylab=""              # removes y label
)
pointLabel(PCAloadings[,1:2],             # set position of labels
           labels=rownames(PCAloadings),  # print labels
           cex=1.5                          # set size of label
) # pointLabel will try to position the text around the points
axis(1,                 # display x-axis
     cex.axis=1.5,      # set size of text
     lwd=1.5            # set size of axis line
)
axis(2,                 # display y-axis
     las=2,             # argument sets direction of text, 2 is perpendicular
     cex.axis=1.5,      # set size of text
     lwd=1.5            # set size of axis line
)
box(lwd=1.5             # line width of box surrounding plot
)
title(xlab=explainPCAvar[["PC1"]],    # % variance explained on PC1
      ylab=explainPCAvar[["PC2"]],    # % variance explained on PC2 
      main="Loading",                 # Title
      cex.lab=1.5,                    # size of label text
      cex.main=1.5                    # size of title text
)

