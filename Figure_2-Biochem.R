require(dplyr)
require(ggplot2)
require(ggpubr)

source("scripts/extra_functions.R")

######################################
# Calculate means 
######################################
#import the dataset
ML_metadata_means <- read.table("data/ML_metadata.txt", h = TRUE) %>% 
  mutate(DON = TDN - (NH4+NO3+NO2)) %>%  # calculate dissolved organic nitrogen
  reshape2::melt(id.vars =c("Time", "Treatment", "Bottle"))%>% 
  filter(variable %in% c("PO43","DON","NH4",
                         "DOC","DFAA","DCAA")) %>% 
  group_by(Time, Treatment, variable) %>%
  mutate(value = as.numeric(value)) %>% 
  summarise(mean = mean(value), se = se(value)) %>% 
  mutate(variable = factor(variable, levels =c("PO43","NH4",
                                                  "DOC","DON","DFAA","DCAA")))

ML_metadata_nitrogen<- read.table("data/ML_metadata.txt",
                                h = TRUE, sep="\t", dec = ",", blank.lines.skip = TRUE) %>% 
  filter(!Treatment=="") %>% 
  mutate(Treatment = case_when(Treatment == "J" ~ "Cteno-OM",
                               Treatment == "C" ~ "Control"),
         Nitrogen = NO3+NO2) %>% 
  reshape2::melt(id.vars =c("Time", "Treatment", "Bottle"))%>% 
  filter(variable %in% c("NO3","TDN","NH4",
                    "NO2","DFAA","DCAA")) %>% 
  group_by(Time,Bottle, Treatment, variable) %>%
  mutate(value = as.numeric(value)) %>% 
  summarise(mean = mean(value, na.rm=TRUE), se = se(value)) 
######################################
# Plot concentrations
######################################
#plot only cell abundances
meta.p <- ML_metadata_means %>% 
  filter(!is.na(mean),
         Time < 100) %>% 
  ggplot(aes(x = Time, y = mean, group = Treatment, colour = Treatment))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.2, colour = "black")+
  geom_point(size = 4, colour = "black")+
  geom_point(size = 3)+
  geom_line(linetype=2)+
  scale_y_continuous(n.breaks=4)+
  labs(x="Time (hours)", 
       y= expression(paste("Concentration (# ",mL^-1,")")))+
  scale_colour_manual(values =c(#"Innoculum"="gray50",
    "Cteno-OM"="red",
    "Control"="blue"))+
  facet_wrap(~variable, scales = "free_y", ncol = 3)+
  theme_EF+
  theme(legend.position = "bottom")

#save the plot
ggsave("./Figures/Figure_2-Biochem.pdf",
       plot = last_plot(),
       units = "mm",
       width = 90, height = 50, 
       scale = 3,
       dpi = 300)

#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()
