require(dplyr)
require(ggplot2)

source("scripts/extra_functions.R")

#import the dataset
ML_metadata <- read.table("data/ML_metadata.txt",
                         h = TRUE, sep="\t", dec = ",", blank.lines.skip = TRUE) %>% 
                  filter(!Treatment=="") %>% 
                  mutate(Treatment = case_when(Treatment == "J" ~ "Cteno-OM",
                                               Treatment == "C" ~ "Control"),
                         Nitrogen = NO3+NO2)

ML_metadata_long <- ML_metadata %>% 
  reshape2::melt(id.vars =c("Time", "Treatment", "Bottle"))

ML_metadata_mean <- ML_metadata_long %>% 
  group_by(Time, Treatment, variable) %>%
  mutate(variable = factor(variable, levels = c("Abundance","TDN","NH4",
                                                "NO3","NO2","PO43","Nitrogen",
                                                "DOC","DCAA","DFAA")),
         value = as.numeric(value)) %>% 
  summarise(mean = mean(value), se = se(value))
######################################
# Plot cell abundances, DON, ammonia and phosphate
######################################
metadata.p <- ML_metadata_mean %>% 
  filter(!is.na(mean),
         #Time < 100,
         variable %in% c("Abundance","TDN","NH4","PO43")) %>% 
  ggplot(aes(x = Time, y = mean, group = Treatment, colour = Treatment))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.2)+
  geom_point(size = 4, colour = "black")+
  geom_point(size = 3)+
  geom_line(linetype=2)+
  scale_y_continuous(n.breaks=4)+
  labs(x="Time (hours)", 
       y= expression(paste("Concentration (# ",mL^-1,")")))+
  scale_colour_manual(values =c(#"Innoculum"="gray50",
                              "Cteno-OM"="red",
                              "Control"="blue"))+
  facet_wrap(~variable, scales = "free_y")+
  theme_EF

#save the plot
ggsave("./Figures/Figure_1-metadata.pdf",
       plot = metadata.p,
       units = "mm",
       width = 90, height = 90, 
       scale = 2,
       dpi = 300)

######################################
# Plot cell abundances, DON, ammonia and phosphate
######################################
metadata_S.p <- ML_metadata_mean %>% 
  filter(!is.na(mean),
         #Time < 100,
         variable %in% c("Nitrogen",#"NO3","NO2",
                         "DOC","DCAA","DFAA")) %>% 
  ggplot(aes(x = Time, y = mean, group = Treatment, colour = Treatment))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.2)+
  geom_point(size = 4, colour = "black")+
  geom_point(size = 3)+
  geom_line(linetype=2)+
  scale_y_continuous(n.breaks=4)+
  labs(x="Time (hours)", 
       y= expression(paste("Concentration (# ",mL^-1,")")))+
  scale_colour_manual(values =c(#"Innoculum"="gray50",
    "Cteno-OM"="red",
    "Control"="blue"))+
  facet_wrap(~variable, scales = "free_y")+
  theme_EF

#save the plot
ggsave("./Figures/Figure_S1-metadata.pdf",
       plot = metadata_S.p,
       units = "mm",
       width = 90, height = 120, 
       scale = 2,
       dpi = 300)

######################################
# Statistics (work in progress)
######################################
ML_metadata_long %>% 
  filter(!is.na(value),
         Time < 100,
         variable == "DOC") %>% 
  group_by(Treatment) %>% 
  t_test(value ~ Time, paired = TRUE)






######################################
# Plot cell abundances with logarithmic axis
######################################
cell_abund.p <- ML_metadata %>% 
                  group_by(Time, Treatment) %>% 
                  summarise(across(where(is.numeric), list(mean = mean, sd = sd))) %>% 
                  ggplot(aes(x = Time, y = Abundance_mean, group = Treatment, colour = Treatment))+
                    geom_point(size = 5)+
                    geom_line()+
                    geom_errorbar(aes(ymin= Abundance_mean-Abundance_sd, ymax= Abundance_mean+Abundance_sd))+
                    scale_y_log10()+
                    theme_bw()



#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()

  
  
  
  