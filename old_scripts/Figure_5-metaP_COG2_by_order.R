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
metaP_merged<- readRDS("data/metaproteome/metaP_runB_merged.rds")

###################
#Functional compositions
###################
#conduct NSAF transformation and melt into df
prot_nsaf.long<- add_nsaf(metaP_merged, "prot_length") %>% 
  psmelt(.)

#sum up each COG
prot_nsaf.COG.agg.taxa <- prot_nsaf.long %>%
  mutate(COG20_CATEGORY_accession = gsub("!!!.*","",COG20_CATEGORY_accession)) %>% 
  mutate(COG20_CATEGORY_accession = case_when(is.na(COG20_CATEGORY_accession) ~ "Unk",
                                              TRUE ~ COG20_CATEGORY_accession)) %>%
  group_by(Type,Sample_name, Treatment, Class, Order, COG20_CATEGORY_accession) %>% 
  filter(Abundance> 0) %>% 
  mutate(Type = factor(Type, levels = c("Exocellular","Cellular"))) %>% 
  summarize(Tot.abundance = sum(Abundance), n_prot = length(Abundance)) 

#mean of replicates
prot_nsaf.COG.agg.taxa.mean <-  prot_nsaf.COG.agg.taxa%>% 
  group_by(Type, Treatment, Order, COG20_CATEGORY_accession) %>% 
  summarize(Mean.abundance = mean(Tot.abundance),
            n_prot_mean = mean(n_prot),
            n_prot_se = se(n_prot)) 

#plot
prot_nsaf.COG.p <- prot_nsaf.COG.agg.taxa.mean %>% 
  filter(Order %in% c("Alteromonadales", "Oceanospirillales",
                      "Pelagibacterales")) %>% 
  ggplot(aes(x = Treatment, y = Mean.abundance*100,
             fill = COG20_CATEGORY_accession)) + 
  facet_grid(Order~Type, scales = "free_x") +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Protein proportions (%) \n")+
  theme_EF+
  theme(legend.position = "bottom")#+
#guides(fill=guide_legend(ncol=2))


#save the plot
ggsave("./Figures/Figure_5-metaP_COG20_by_order.pdf", 
       plot = prot_nsaf.COG.p,
       units = "mm",
       width = 120, height = 90, 
       scale = 3,
       dpi = 300)


prot_nsaf.COG.agg.taxa_sub <- prot_nsaf.long %>%
  mutate(COG20_CATEGORY_accession = gsub("!!!.*","",COG20_CATEGORY_accession),
         COG20_FUNCTION_function = gsub("!!!.*","",COG20_FUNCTION_function)) %>% 
  mutate(COG20_CATEGORY_accession = case_when(is.na(COG20_CATEGORY_accession) ~ "Unk",
                                              TRUE ~ COG20_CATEGORY_accession)) %>%
  group_by(Type,Sample_name, Treatment, Class, Order, COG20_FUNCTION_accession,COG20_FUNCTION_function) %>% 
  filter(Abundance> 0) %>% 
  mutate(Type = factor(Type, levels = c("Exocellular","Cellular"))) %>% 
  summarize(Tot.abundance = sum(Abundance), n_prot = length(Abundance)) %>% 
  filter(Treatment != "Inoculum")


prot_nsaf.COG.agg.taxa_mean <- prot_nsaf.COG.agg.taxa_sub%>% 
  group_by(Type, Treatment, Order, COG20_FUNCTION_accession,COG20_FUNCTION_function) %>% 
  summarize(Mean.abundance = mean(Tot.abundance),
            n_prot_mean = mean(n_prot),
            n_prot_se = se(n_prot)) 


COGs<- prot_nsaf.long %>% 
  select(COG20_CATEGORY_function, COG20_CATEGORY_accession) %>% 
  unique()


###################
#Functional composition by taxa
###################
#sum up each COG
prot_nsaf.COG.agg <- prot_nsaf.long %>% 
  mutate(COG20_CATEGORY_function = gsub("!!!.*","",COG20_CATEGORY_function)) %>% 
  mutate(COG20_CATEGORY_function = case_when(is.na(COG20_CATEGORY_function) ~ "Unk",
                                             TRUE ~ COG20_CATEGORY_function)) %>%
  group_by(Treatment, Type, Sample_name, COG20_CATEGORY_function) %>% 
  summarize(Tot.abundance = sum(Abundance)) 

#calculate mean per fraction
prot_nsaf.COG.frac.mean <- prot_nsaf.COG.agg %>% 
  group_by(Treatment, Type, COG20_CATEGORY_function) %>% 
  summarise(mean = mean(Tot.abundance), se = se(Tot.abundance))


#remove below 1% ra
COG_categories <- sort(as.character(unique(prot_nsaf.COG.frac.mean$COG20_CATEGORY_function[!prot_nsaf.COG.frac.mean$mean<0.01])))
prot_nsaf.COG.frac.mean$COG20_CATEGORY_function[prot_nsaf.COG.frac.mean$mean<0.01] <- "Other"
prot_nsaf.COG.frac.mean$COG20_CATEGORY_function[is.na(prot_nsaf.COG.frac.mean$COG20_CATEGORY_function)] <- "Other"
prot_nsaf.COG.frac.mean$COG20_CATEGORY_function <- factor(prot_nsaf.COG.frac.mean$COG20_CATEGORY_function,
                                                          levels=c(COG_categories,"Other"))

#plot
prot_COG.p<- ggplot(prot_nsaf.COG.frac.mean, 
                    aes(x = Treatment, y = mean*100,
                        fill = COG20_CATEGORY_function)) + 
  facet_wrap(.~Type) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = tol21rainbow)+ 
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Protein proportion (%) \n")+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Figure_S4-MetaP_COG20_comp.pdf", 
       plot = prot_COG.p,
       units = "mm",
       width = 120, height = 90, 
       scale = 3,
       dpi = 300)



#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()

