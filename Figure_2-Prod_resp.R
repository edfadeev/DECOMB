require(dplyr)
require(ggplot2)

source("scripts/extra_functions.R")

#import the dataset
fish_raw <- read.table("data/microscopy/FiSH_results.txt",
                          h = TRUE, sep="\t", dec = ",", blank.lines.skip = TRUE) 

######################################
# Plot productivity and respiration by taxa
######################################
FiSH.p <- fish_raw %>% 
            select(-c(Abundance..cells.mL.1.)) %>% 
            reshape2::melt(id.vars =c("Time", "Treatment","Bottle")) %>%
            mutate(value = as.numeric(value),
                   Method = case_when(gsub("_.*","",variable)=="P" ~ "Production",
                                      gsub("_.*","",variable)=="R" ~ "Respiration"),
                   Taxa = case_when(gsub("R_|P_","",variable)=="Bac" ~ "Bact.",
                                    gsub("R_|P_","",variable)=="Psu" ~ "Pseudo.",
                                    gsub("R_|P_","",variable)=="Vib" ~ "Vibrio.",
                                    gsub("R_|P_","",variable)=="Alt" ~ "Altero."),
                   Time = factor(Time, levels= c("0","21","154"))) %>% 
            filter(!is.na(value),
                   #  !Time == "0"
                   ) %>% 
            group_by(Time, Treatment, Method,Taxa) %>%
            summarise(across(where(is.numeric), list(mean = mean, se = se))) %>% 
            #ggplot(aes(x = Taxa, y = value_mean, fill = Taxa))+
            ggplot(aes(x = Time, y = value_mean, colour = Treatment))+
                #geom_bar(size = 3, stat = "identity")+
                geom_point(size = 5, colour = "black")+
                geom_point(size = 4)+
                geom_errorbar(aes(ymin= value_mean-value_se, ymax= value_mean+value_se), width = 0.3)+
                #facet_grid(Method~Time+Treatment)+
                scale_colour_manual(values =c("Innoculum"="gray50",
                                              "J"="red",
                                              "C"="blue"))+
                facet_grid(Method~Taxa)+
                theme_bw()+
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      text = element_text(size = 20))

#save the plot
ggsave("./Figures/Figure_2-Prod_resp.pdf",
       plot = FiSH.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()