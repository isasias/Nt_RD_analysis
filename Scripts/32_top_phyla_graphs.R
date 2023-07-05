#### Top phyla abundance ####

# Subset top phyla
Top_phyla <- tax_glom(Psoil_rel,taxrank = "Phylum", NArm = FALSE)

Top_phyla <- prune_taxa(names(sort(taxa_sums(Top_phyla),TRUE)[1:5]), Top_phyla)

sorder <- c("A1","B1","C1","A2","B2","C2",
            "A3","B3","C3","A4","B4","C4",
            "A5","B5","C5",
            "A6","B6","C6","A7","B7","C7",
            "D1","D2","D3","D4")

# Top phyla barplot

Phyla_plot <-plot_composition(Top_phyla, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling2"),name = "Phylum", labels = c("Actinobacteria", "Proteobacteria", "Bacteroidetes","Chloroflexi","Acidobacteria"))+
  guides(x=guide_axis(angle = 35))+
  scale_x_discrete( labels=c(rep("36LP",3),rep("WTLP",3),  
                             rep("36C",3), rep("WTC",3),  
                             rep("36Phi",3), 
                             rep("36PM",3), rep("WTPM",3),
                             "NPLP","NPC","NPPhi","NPPM"))+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

### Subset actino and proteobacteria ###

Top_phyla <- subset_taxa(Psoil_rel, Phylum == "Actinobacteriota")
Act_phyla <- tax_glom(Top_phyla,taxrank = "Order", NArm = FALSE)
Act_phyla <- prune_taxa(names(sort(taxa_sums(Act_phyla),TRUE)[1:5]), Act_phyla)

plot_composition(Act_phyla, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling2"),name = "Phylum")+
  guides(x=guide_axis(angle = 35))+
  scale_x_discrete( labels=c(rep("36LP",3),rep("WTLP",3),  
                             rep("36C",3), rep("WTC",3),  
                             rep("36Phi",3), 
                             rep("36PM",3), rep("WTPM",3),
                             "NPLP","NPC","NPPhi","NPPM"))


Top_phyla <- subset_taxa(Psoil_rel, Phylum == "Proteobacteria")
Pro_phyla <- tax_glom(Top_phyla,taxrank = "Order", NArm = FALSE)
Pro_phyla <- prune_taxa(names(sort(taxa_sums(Pro_phyla),TRUE)[1:5]), Pro_phyla)

plot_composition(Pro_phyla, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling2"),name = "Phylum")+
  guides(x=guide_axis(angle = 35))+
  scale_x_discrete( labels=c(rep("36LP",3),rep("WTLP",3),  
                             rep("36C",3), rep("WTC",3),  
                             rep("36Phi",3), 
                             rep("36PM",3), rep("WTPM",3),
                             "NPLP","NPC","NPPhi","NPPM"))

### Subset top orders ###

Top_ordersP <- subset_taxa(Top_phyla, Order == "Burkholderiales")
Top_ordersP <- tax_glom(Top_ordersP,taxrank = "Genus", NArm = FALSE)
OrdersP <- prune_taxa(names(sort(taxa_sums(Top_ordersP),TRUE)[1:5]), Top_ordersP)

### Final Actinobacteria top orders graphs ###

Top_phyla <- subset_taxa(Psoil_rel, Phylum == "Actinobacteriota")
Top_ordersA <- subset_taxa(Top_phyla, Genus == c("Rubrobacter","Cellulomonas",
                                                 "Pseudarthrobacter","Agromyces"))
OrdersA <- tax_glom(Top_ordersA,taxrank = "Genus", NArm = FALSE)

plot_composition(OrdersA, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling2"),name = "Phylum")+
  guides(x=guide_axis(angle = 35))+
  scale_x_discrete( labels=c(rep("36LP",3),rep("WTLP",3),  
                             rep("36C",3), rep("WTC",3),  
                             rep("36Phi",3), 
                             rep("36PM",3), rep("WTPM",3),
                             "NPLP","NPC","NPPhi","NPPM"))

### Final Protebacteria top orders graphs ###

Top_phyla <- subset_taxa(Psoil_rel, Phylum == "Proteobacteria")
Top_ordersP <- subset_taxa(Top_phyla, Genus == c("Microvirga","Sphingomonas",
                                                 "MND1","Nitrosomonas"))
OrdersP <- tax_glom(Top_ordersP,taxrank = "Genus", NArm = FALSE)

plot_composition(OrdersP, plot.type = "barplot",sample.sort = sorder, x.label = "Sample_ID") +
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling2"),name = "Phylum")+
  guides(x=guide_axis(angle = 35))+
  scale_x_discrete( labels=c(rep("36LP",3),rep("WTLP",3),  
                             rep("36C",3), rep("WTC",3),  
                             rep("36Phi",3), 
                             rep("36PM",3), rep("WTPM",3),
                             "NPLP","NPC","NPPhi","NPPM"))
