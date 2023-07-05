#### Metagenassist analyst

library(pheatmap)




metaheat <- read.csv("~/Grad_School/Maestria/Processed_data/METABOLISM.filtered.csv")
row.names(metaheat) <- metaheat[,1]
metaheat <- metaheat[,-1]

samplemeta <- read.csv("~/Grad_School/Maestria/Processed_data/metadata.csv")

metaheat <- metaheat %>%
  add_column(Plant_type = samplemeta$Plant_Type) %>%
  add_column(Treatment = samplemeta$Treatment)

for(i in 1:ncol(metaheat)){
  shapiro <- shapiro.test(metaheat[,i])
  normal <- ifelse(shapiro[["p.value"]]>0.05, "YES","NO")
  print(c(i,normal))
}

write.csv(metabolism,"~/Grad_School/Maestria/Processed_data/Metabolic_groups.csv")

### sulfur oxidizing
sulfur_anova <- aov(Sulfur.metabolizing ~ Plant_type * Treatment, 
                   data = metaheat) 

# Test homogeneity of variance
plot(sulfur_anova, 1) # check outliers
leveneTest(Sulfur.metabolizing ~ Plant_type * Treatment, 
           data = metaheat) # higher than 0.05 so it meets ANOVA assumption

# Anova
Anova(sulfur_anova) # Treatments and interactions are significant 
# Tukey HSD 
sulfur_Tukey <- TukeyHSD(sulfur_anova)

## Graphs
ggplot(metaheat, aes(x = Plant_type, y = Sulfur.metabolizing)) +
  geom_boxplot()+
  theme_minimal(base_size = 30)+
  scale_fill_manual(name= NULL, values = c("steelblue","burlywood3","lightpink"))+
  xlab("Treatment") + ylab("Ratio")


metabolism <- read.csv("~/Grad_School/Maestria/Processed_data/METABOLISM.filtered.csv")
row.names(metabolism) <- metabolism[,1]
metabolism <- metabolism[,-1]
metabolism <- as.matrix(metabolism)
sapply(metabolism, class)  
metabolism <- t(metabolism)

## Remove empty groups
metabolism <- metabolism[-c(92:99),]
metabolism <- metabolism[-c(27:90),]
metabolism <- metabolism[-26,]

## reorganize columns
Trtorder <- c("D1","D2","D3","D4","A2","B2","C2","A4","B4","C4",
              "A5","B5","C5","A7","B7","C7","A1","B1","C1","A3",
              "B3","C3","A6","B6","C6")
metabolism <- metabolism[,Trtorder] 

## base R
# no dendogram
heatmap(metabolism, Colv = NA) 
# heatmap(metabolism, scale = "column") not useful check scale function

coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

## ggplot all groups
meta_melt <- melt(metabolism)                                       

meta_melt$groups <- cut(meta_melt$value,               # Add group column
                       breaks = c(0, 10, 100, 1000, 5000, 7000, 10000,
                                  20000,30000,40000,50000,70000))

colorRampPalette(c("blue", "red"))( 11 )

ggp <- ggplot(meta_melt, aes(X2, X1)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill=groups))+
  scale_fill_manual(breaks = levels(meta_melt$groups),
                    values = c("#0000FF","#1900E5", "#3300CC",
                              "#4C00B2", "#660099", "#7F007F",
                            "#990065", "#B2004C" ,"#CC0032", "#E50019" ,
                            "#FF0000")) # 11 colors 
ggp

## Top 15  

top_meta <- metabolism[1:15,]
top_melt <- melt(top_meta)

# base R
heatmap(top_meta, Colv = NA, Rowv = NA, 
        col = colorRampPalette(c("oldlace", "darkslateblue"))(30))
legend(x = "bottomright", cex = 0.8, 
       legend= c(1:20),
       fill = colorRampPalette(brewer.pal(8, "Blues"))(20))

# ggplot 
 "oldlace"

top_melt$groups <- cut(top_melt$value,               # Add group column
                        breaks = c(400, 500, 600, 700, 800, 900, 1000,
                                   2000,3000,4000,5000,6000,7000,10000,
                                   15000,16000,18000,20000,22000,26000,
                                   30000,35000,40000,45000,50000,55000,
                                   70000))

ggplot(top_melt, aes(X2, X1)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill=groups),color = "black")+
  scale_fill_manual(breaks = levels(top_melt$groups),
                    values = colorRampPalette(c("rosybrown1", "darkslateblue"))(36))

### Using pheatmap

pheatmap(metabolism, cluster_cols = FALSE, cluster_rows = FALSE,
         scale = "row")

pheatmap(metabolism[1:7,], cluster_cols = FALSE, 
         cluster_rows = FALSE)

#### Data transformation for heatmap ####

## Log transformation 

metalog <- log10(metabolism)
metalog[is.infinite(metalog)] <- NA # remove infinite values

# base R
heatmap(metalog, Colv = NA, Rowv = NA, 
        col = colorRampPalette(c("oldlace", "darkslateblue"))(30))

heatmap(metalog[1:15,], Colv = NA, Rowv = NA, 
        col = colorRampPalette(c("oldlace", "darkslateblue"))(30))

# ggplot
metalog_melt <- melt(metalog)  

ggplot(metalog_melt, aes(X2, X1)) +                  
  geom_tile(aes(fill=value),color = "black")+
  scale_fill_gradient(low = "rosybrown1", high = "darkslateblue")

# pheatmap
pheatmap(metalog, cluster_cols = FALSE, cluster_rows = FALSE,
         scale = "row") # row scaling this might be the best one
# change colors and add column annotation


## Relative abundance transformation

# Retrieve max value per row

for(i in 1:nrow(metabolism)){
  nmax <- max(metabolism[i,])
  print(c(i,nmax))
}

rel_meta <- matrix(data = NA, nrow = 26, ncol = 25)

# for loop to relative abundance 
for(j in 1:25){
  for (i in 1:26) {
    rel_meta[i,j]= metabolism[i,j]/max(metabolism[i,])
  }
}

colnames(rel_meta) <- colnames(metabolism)
row.names(rel_meta) <- row.names(metabolism)

# base R: doesn't work
heatmap.2(rel_meta, Colv = NA, Rowv = NA, 
        col = colorRampPalette(c("rosybrown1", "darkslateblue"))(30))

heatmap(rel_meta[1:15,], Colv = NA, Rowv = NA, scale = "column",
        col = colorRampPalette(c("rosybrown1", "darkslateblue"))(30))

# ggplot
relog_melt <- melt(rel_meta)  

ggplot(relog_melt, aes(X2, X1)) +                  
  geom_tile(aes(fill=value))+
  scale_fill_gradient(low = "rosybrown1", high = "darkslateblue")

# pheatmap
pheatmap(rel_meta, cluster_cols = FALSE, cluster_rows = FALSE)

### Gplot heatmap2

breaks <- c(seq(0,30,length.out= 10),
            seq(31,100,length.out=70),
            seq(101,1000,length.out=100),
            seq(1001,3000,length.out=70),
            seq(3001,10000,length.out=180),
            seq(10001,25000,length.out=70),
            seq(25001,66000,length.out=180))

gradient3 <-  c(colorpanel(sum(breaks[-1]<=30),"snow2", "lightskyblue1"),
             colorpanel(sum(breaks[-1]>30 & breaks[-1]<=100), "lightskyblue1","plum"),
             colorpanel(sum(breaks[-1]>100 & breaks[-1]<=1000), "plum", "mediumorchid"),
             colorpanel(sum(breaks[-1]>1000 & breaks[-1]<=3000), "mediumorchid","magenta4"),
             colorpanel(sum(breaks[-1]>3000 & breaks[-1]<=10000), "magenta4","steelblue4"),
             colorpanel(sum(breaks[-1]>10000 & breaks[-1]<=25000), "steelblue4", "royalblue4"),
             colorpanel(sum(breaks[-1]>25000 & breaks[-1]<=66000),"royalblue4","black"))
             
             
breaks2 <- c(seq(0,5,length.out= 10),
            seq(6,1000,length.out=250),
            seq(1001,20000,length.out=150),
            seq(20001,66000,length.out=280))

gradient2 <-  c(colorpanel(sum(breaks2[-1]<=20),"snow2", "lightskyblue1"),
                colorpanel(sum(breaks2[-1]>20 & breaks2[-1]<=1000), "lightskyblue1","thistle"),
                colorpanel(sum(breaks2[-1]>1000 & breaks2[-1]<=12000), "thistle","magenta4"),
                colorpanel(sum(breaks2[-1]>12000 & breaks2[-1]<=20000),"magenta4","steelblue4"),
                colorpanel(sum(breaks2[-1]>20000 & breaks2[-1]<=66000),"steelblue4","royalblue4","black"))

# colors <- c(colorRampPalette(c("aliceblue","lightskyblue1"))(5),
#            colorRampPalette(c("lightskyblue1","lightpink"))(5),
#            colorRampPalette(c("lightpink","plum"))(10),
#            colorRampPalette(c("plum","mediumpurple1"))(5),
#            colorRampPalette(c("mediumpurple1","cornflowerblue"))(15),
#            colorRampPalette(c("cornflowerblue","magenta4"))(5),
#            colorRampPalette(c("navyblue","black"))(10))



heatmap.2(metabolism, dendrogram = "none", Rowv = F, Colv = F, 
          trace="none",margin = c(1,5),
          breaks = breaks, col=gradient3,
          main="my heatmap", scale = "none",
          sepcolor="grey", sepwidth=c(0.009,0.009),
          colsep= 1:25 ,rowsep=1:26, key = TRUE)

# c(4,7,10,13,16,19,21,25)

#### Final choice pheatmap ####

## Annotation column
cat_df = data.frame("Plant Type" = c(rep("No Plant",4),
                                   rep("36 ptxD",12),
                                   rep("Wild Type",9)))
rownames(cat_df) = colnames(metabolism)

mat_colors <-  c("burlywood3","steelblue","lightpink")
names(mat_colors) <- unique(cat_df$category)
ann_colors <- list(Plant.Type= c('No Plant'="burlywood3",
                                   '36 ptxD'="steelblue", 'Wild Type'= "lightpink"))

### Using Original data 

mat_breaks <- seq(min(metabolism), max(metabolism), length.out = 10)

# change to matrix

metabolism <- data.matrix(metabolism)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(metabolism, n = 21)

pheatmap(metabolism,cluster_rows=FALSE, cluster_cols=FALSE, scale = "none",
  color = colorRampPalette(c("snow2", "lightskyblue1", "plum3","magenta4","steelblue4", "black"))(20),
  breaks            = mat_breaks,
  drop_levels       = F,
  fontsize          = 10,
  main              = "Community metabolic profiles",
  gaps_col = c(4,7,10,13,16,19,22,25),
  labels_col = c ("NPLP", "NPC", "NPPhi","NPPM",
                  rep("36LP",3),rep("36C",3),
                  rep("36Phi",3),rep("36PM",3),
                  rep("WTLP",3),rep("WTC",3),rep("WTPM",3)))

