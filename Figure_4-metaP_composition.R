########################################
#Metaproteome composition analysis
########################################
#load libraries
require(dplyr)
require(tidyr)
require(phyloseq)
require(ggplot2)
require(ggpubr)

source("scripts/extra_functions.R")

#load metaproteome phyloseq object
metaP_runB<- readRDS("data/metaproteome/metaP_ps_runB.rds")

###################
#Taxonomic compositions
###################
#conduct NSAF transformation
metaP_agg_nsaf<- add_nsaf(metaP_runB, "prot_length")

#melt phyloseq into a dataframe for ploting
prot_nsaf.long <- psmelt(metaP_agg_nsaf) 

#sum up each taxa
prot_nsaf.class.agg <- prot_nsaf.long %>% group_by(Fraction, Sample_name, Group, Class, Order) %>% 
  summarize(Tot.abundance = sum(Abundance)) %>% 
  mutate(Fraction = factor(Fraction, levels =c("MP","EH","EL")))

#remove below 1% ra
taxa_classes <- sort(as.character(unique(prot_nsaf.class.agg$Class[!prot_nsaf.class.agg$Tot.abundance<0.01])))

prot_nsaf.class.agg$Class[prot_nsaf.class.agg$Tot.abundance<0.01] <- "Other taxa"
prot_nsaf.class.agg$Class[is.na(prot_nsaf.class.agg$Class)] <- "Other taxa"

prot_nsaf.class.agg$Class <- factor(prot_nsaf.class.agg$Class,
                                      levels=c(taxa_classes,"Other taxa"))

prot_tax_comp.p<- ggplot(prot_nsaf.class.agg, 
       aes(x = Sample_name, y = Tot.abundance,
           fill = Class)) + 
  facet_wrap(.~Fraction, scales = "free") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Protein proportions (>1%) \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        text=element_text(size=20),legend.position = "right", 
        axis.title.x = element_blank())


#save the plot
ggsave("./Figures/metaP_tax_order_comp.png", 
       plot = prot_tax_comp.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

###################
#Functional compositions
###################
#sum up each COG
prot_nsaf.COG.agg <- prot_nsaf.long %>% 
  mutate(COG20_CATEGORY_function = gsub("!!!.*","",COG20_CATEGORY_function)) %>% 
  mutate(COG20_CATEGORY_function = case_when(is.na(COG20_CATEGORY_function) ~ "Unk",
                                             TRUE ~ COG20_CATEGORY_function)) %>%
  group_by(Fraction, Sample_name, Group, COG20_CATEGORY_function) %>% 
  summarize(Tot.abundance = sum(Abundance)) %>% 
  mutate(Fraction = factor(Fraction, levels =c("MP","EH","EL")))

#remove below 1% ra
COG_categories <- sort(as.character(unique(prot_nsaf.COG.agg$COG20_CATEGORY_function[!prot_nsaf.COG.agg$Tot.abundance<0.02])))

prot_nsaf.COG.agg$COG20_CATEGORY_function[prot_nsaf.COG.agg$Tot.abundance<0.02] <- "Other"
prot_nsaf.COG.agg$COG20_CATEGORY_function[is.na(prot_nsaf.COG.agg$COG20_CATEGORY_function)] <- "Other"

prot_nsaf.COG.agg$COG20_CATEGORY_function <- factor(prot_nsaf.COG.agg$COG20_CATEGORY_function,
                                    levels=c(COG_categories,"Other"))

#plot
prot_COG.p<- ggplot(prot_nsaf.COG.agg, 
                         aes(x = Sample_name, y = Tot.abundance,
                             fill = COG20_CATEGORY_function)) + 
  facet_grid(.~Fraction, scales = "free") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("COG category ( >2% of all proteins) \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        text=element_text(size=20),legend.position = "right", 
        axis.title.x = element_blank())


#save the plot
ggsave("./Figures/metaP_COG20_comp.png", 
       plot = prot_COG.p,
       units = "cm",
       width = 50, height = 30, 
       #scale = 1,
       dpi = 300)



#Plot both composition figures together
ggarrange(prot_tax_comp.p,prot_COG.p,
          ncol = 2, nrow = 1, align = "hv")



###################
#Functional composition by taxa
###################
#sum up each COG
prot_nsaf.COG.agg.taxa <- prot_nsaf.long %>%
  mutate(COG20_CATEGORY_function = gsub("!!!.*","",COG20_CATEGORY_function)) %>% 
  mutate(COG20_CATEGORY_function = case_when(is.na(COG20_CATEGORY_function) ~ "Unk",
                                             TRUE ~ COG20_CATEGORY_function)) %>%
  group_by(Fraction,Sample_name,  Group, Class, COG20_CATEGORY_function) %>% 
  summarize(Tot.abundance = sum(Abundance)) %>% 
  mutate(Fraction = factor(Fraction, levels =c("MP","EH","EL")))

#merge rare classes (below 1% )
taxa_classes <- sort(as.character(unique(prot_nsaf.COG.agg.taxa$Class[!prot_nsaf.COG.agg.taxa$Tot.abundance<0.01])))

prot_nsaf.COG.agg.taxa$Class[prot_nsaf.COG.agg.taxa$Tot.abundance<0.01] <- "Other taxa"
prot_nsaf.COG.agg.taxa$Class[is.na(prot_nsaf.COG.agg.taxa$Class)] <- "Other taxa"
prot_nsaf.COG.agg.taxa$Class <- factor(prot_nsaf.COG.agg.taxa$Class,
                                    levels=c(taxa_classes,"Other taxa"))

#omit rare COG categories by changing to "Other"
prot_nsaf.COG.agg.taxa$COG20_CATEGORY_function <- factor(prot_nsaf.COG.agg.taxa$COG20_CATEGORY_function,
                                                         levels=c(COG_categories,"Other"))

#plot
prot_nsaf.COG.agg.taxa %>% 
  filter(Class %in% c("Alphaproteobacteria", "Gammaproteobacteria",
                      "Unknown_Synechococcales", "Other taxa")) %>% 
ggplot(aes(x = Sample_name, y = Tot.abundance,
                        fill = COG20_CATEGORY_function)) + 
  facet_grid(Class~Fraction, scales = "free_x") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Protein proportions (>1%) \n")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90),
        text=element_text(size=20),legend.position = "right", 
        axis.title.x = element_blank())


#save the plot
ggsave("./Figures/metaP_fun_comp_by_taxa.png", 
       plot = last_plot(),
       units = "cm",
       width = 50, height = 30, 
       #scale = 1,
       dpi = 300)


#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()

