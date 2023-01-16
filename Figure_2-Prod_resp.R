require(dplyr)
require(ggplot2)

source("scripts/extra_functions.R")

#import the dataset
fish_raw <- read.table("data/microscopy/FiSH_results.txt",
                          h = TRUE, sep="\t", dec = ",", blank.lines.skip = TRUE) 

######################################
# Plot productivity and respiration by taxa
######################################
FiSH.mean <- fish_raw %>% 
            select(-c(Abundance..cells.mL.1.)) %>% 
            reshape2::melt(id.vars =c("Time", "Treatment","Bottle")) %>%
            mutate(value = as.numeric(value),
                   Method = case_when(gsub("_.*","",variable)=="P" ~ "Production",
                                      gsub("_.*","",variable)=="R" ~ "Respiration"),
                   Taxa = case_when(gsub("R_|P_","",variable)=="Bac" ~ "Total",
                                    gsub("R_|P_","",variable)=="Psu" ~ "Pseudo.",
                                    gsub("R_|P_","",variable)=="Vib" ~ "Vibrio.",
                                    gsub("R_|P_","",variable)=="Alt" ~ "Altero."),
                   Time = factor(Time, levels= c("0","21","154")),
                   Treatment = case_when(Treatment == "J" ~ "Cteno-OM",
                                                Treatment == "C" ~ "Control")) %>% 
            filter(!is.na(value),
                   Time != "0"
                   ) %>% 
            group_by(Time, Treatment, Method,Taxa) %>%
            summarise(across(where(is.numeric), list(mean = mean, se = se))) %>% 
  mutate(Taxa = factor(Taxa, levels= c("Total","Altero.","Pseudo.","Vibrio.")))
            
  
#plot production  
prod.p <-  FiSH.mean %>% 
  filter(Method =="Production") %>% 
  ggplot(aes(x = interaction(Taxa, Time), y = value_mean, colour = Taxa, group = Taxa))+
  geom_errorbar(aes(ymin= value_mean-value_se, ymax= value_mean+value_se), width = 0.3)+
  geom_point(size = 5, colour = "black")+
  geom_point(size = 4)+
  scale_colour_manual(values = tol21rainbow)+
  facet_grid( ~ Treatment, scales = "free")+
  theme_EF+
  scale_y_continuous(labels = scales::scientific_format(digits = 0), n.breaks=4)+
  theme(axis.text.x = element_text(angle = 90))

#save the plot
ggsave("./Figures/Figure_2-Bac_prod.pdf",
       plot = prod.p,
       units = "mm",
       width = 90, height = 90, 
       scale = 2,
       dpi = 300)

res.p <-  FiSH.mean %>% 
  filter(Method =="Respiration") %>% 
  ggplot(aes(x = interaction(Taxa, Time), y = value_mean, colour = Taxa, group = Taxa))+
  geom_point(size = 5, colour = "black")+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin= value_mean-value_se, ymax= value_mean+value_se), width = 0.3)+
  scale_colour_manual(values = tol21rainbow)+
  facet_grid( ~ Treatment, scales = "free")+
  theme_EF+
  scale_y_continuous(labels = scales::scientific_format(digits = 0), n.breaks=4)+
  theme(axis.text.x = element_text(angle = 90))

#save the plot
ggsave("./Figures/Figure_S2-Bac_res.pdf",
       plot = res.p,
       units = "mm",
       width = 90, height = 90, 
       scale = 2,
       dpi = 300)

#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
