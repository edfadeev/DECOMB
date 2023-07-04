########################################
#Metaproteome composition analysis
########################################
#load libraries
require(dplyr)
require(tidyr)
require(phyloseq)
require(ggplot2)
require(ggpubr)

#extra functions
source("scripts/extra_functions.R")

#########################################################
#Generate phyloseq object from the proteomics dataset
#!!! Run only once for Fig3-6 scripts!!!
#########################################################
source("scripts/Proteins2phyloseq.R")

#load metaproteome phyloseq object
metaP_merged<- readRDS("data/metaproteome/metaP_merged.rds")

###################
#Taxonomic compositions
###################
#conduct NSAF transformation and melt into df
prot_nsaf.long<- add_nsaf(metaP_merged, "prot_length") %>% 
                    psmelt(.)

#sum up each taxa
prot_nsaf.class.agg <- prot_nsaf.long %>% 
  mutate(Type = factor(Type, levels = c("Cellular","Exocellular"))) %>%
  group_by(Type, Sample_name, Treatment, Class, Order) %>% 
  summarize(Tot.abundance = sum(Abundance))

#calculate mean per fraction
prot_nsaf.frac.mean <- prot_nsaf.class.agg %>% 
  group_by(Treatment, Type, Class, Order) %>% 
  summarise(mean = mean(Tot.abundance), se = se(Tot.abundance))

#remove below 1% ra
taxa_orders <- sort(as.character(unique(prot_nsaf.frac.mean$Order[!prot_nsaf.frac.mean$mean<0.01])))

prot_nsaf.frac.mean$Order[prot_nsaf.frac.mean$mean<0.01] <- "Other taxa"
prot_nsaf.frac.mean$Order[is.na(prot_nsaf.frac.mean$Order)] <- "Other taxa"

prot_nsaf.frac.mean$Order <- factor(prot_nsaf.frac.mean$Order,
                                      levels=c(taxa_orders,"Other taxa"))


#plot
prot_tax_comp.p<- ggplot(prot_nsaf.frac.mean, 
       aes(x = Treatment, y = mean*100,
           fill = Order)) + 
  facet_wrap(.~Type) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Protein proportion (%) \n")+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Figure_4-MetaP_tax_comp.png", 
       plot = prot_tax_comp.p,
       units = "mm",
       width = 120, height = 90, 
       scale = 3,
       dpi = 300)


#summarize taxonomy on class level
prot_nsaf.frac.mean_Class <- prot_nsaf.frac.mean %>% 
  group_by(Treatment, Type, Class) %>% 
  summarize(Tot.abundance = sum(mean))


#check the different orders
top_orders <- prot_nsaf.long %>% 
  filter(Treatment!= "Inoculum",Order != "Other taxa") %>%
  mutate(Type = factor(Type, levels = c("Exocellular","Cellular"))) %>% 
  group_by(Treatment, Type, Order) %>% 
  summarize(Tot.abundance = sum(Abundance))
  
###################
#print session info and clean the workspace
###################
sessionInfo()
rm(list = ls())
gc()

