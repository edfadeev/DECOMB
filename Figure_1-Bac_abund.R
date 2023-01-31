require(dplyr)
require(ggplot2)
require(ggpubr)

source("scripts/extra_functions.R")

######################################
# Calculate mean cell abundance
######################################
#import the dataset
ML_metadata <- read.table("data/ML_metadata.txt",
                         h = TRUE, sep="\t", dec = ",", blank.lines.skip = TRUE) %>% 
                  filter(!Treatment=="") %>% 
                  mutate(Treatment = case_when(Treatment == "J" ~ "Cteno-OM",
                                               Treatment == "C" ~ "Control"),
                         Nitrogen = NO3+NO2)

Bac_abund_mean <- ML_metadata %>% 
  reshape2::melt(id.vars =c("Time", "Treatment", "Bottle"))%>% 
  filter(variable == "Abundance") %>% 
  group_by(Time, Treatment) %>%
  mutate(value = as.numeric(value)) %>% 
  summarise(mean = mean(value), se = se(value))

######################################
# Calculate mean productivity and respiration by taxa
######################################
FiSH.mean <- read.table("data/microscopy/FiSH_results.txt",
                       h = TRUE, sep="\t", dec = ",", blank.lines.skip = TRUE) %>% 
  select(-c(Abundance..cells.mL.1.)) %>% 
  reshape2::melt(id.vars =c("Time", "Treatment","Bottle")) %>%
  filter(!is.na(value),
         Time != "0") %>% 
  mutate(value = as.numeric(value),
         Method = case_when(gsub("_.*","",variable)=="P" ~ "Production",
                            gsub("_.*","",variable)=="R" ~ "Respiration"),
         Taxa = case_when(gsub("R_|P_","",variable)=="Bac" ~ "T",
                          gsub("R_|P_","",variable)=="Psu" ~ "P",
                          gsub("R_|P_","",variable)=="Vib" ~ "V",
                          gsub("R_|P_","",variable)=="Alt" ~ "A"),
         Treatment = case_when(Treatment == "J" ~ "Cteno-OM",
                               Treatment == "C" ~ "Control")) %>% 
  group_by(Time, Treatment, Method,Taxa) %>%
  summarise(across(where(is.numeric), list(mean = mean, se = se))) %>% 
  ungroup() %>% 
  mutate(Taxa = factor(Taxa, levels= c("T","A","P","V")),
         Time = as.numeric(Time))

######################################
# Plot cell abundances, production and respiration
######################################
#plot only cell abundances
Abund.p <- Bac_abund_mean %>% 
  filter(!is.na(mean),
         Time < 100) %>% 
  ggplot(aes(x = Time, y = mean, group = Treatment, colour = Treatment))+
  geom_errorbar(aes(ymin= mean-se, ymax= mean+se), width = 0.2)+
  geom_point(size = 4, colour = "black")+
  geom_point(size = 3)+
  geom_line(linetype=2)+
  scale_y_continuous(n.breaks=4)+
  labs(x="Time (hours)", 
       y= expression(paste("Concentration (Cells ",mL^-1,")")))+
  scale_colour_manual(values =c(#"Innoculum"="gray50",
    "Cteno-OM"="red",
    "Control"="blue"))+
  #facet_wrap(~variable, scales = "free_y")+
  theme_EF


#produce the FISH figures
fish_plots <- lapply(c("Production","Respiration"), function (i)
  FiSH.mean %>% 
    filter(Method == i, Time == 21) %>% 
    ggplot(aes(x = Taxa, y = value_mean, colour = Treatment, group = Taxa))+
    geom_errorbar(aes(ymin= value_mean-value_se, ymax= value_mean+value_se), width = 0.2)+
    geom_point(size = 5, colour = "black")+
    geom_point(size = 4)+
    #scale_colour_manual(values = tol21rainbow)+
    scale_colour_manual(values =c(#"Innoculum"="gray50",
      "Cteno-OM"="red",
      "Control"="blue"))+
    #facet_grid( . ~ Method, scales = "free")+
    labs(x="Taxa")+
    theme_EF+
    scale_y_continuous(labels = scales::scientific_format(digits = 0), n.breaks=4)+
    theme(axis.title.y = element_blank(),
          legend.position = "none")
)


#define font label size
font_labels <- list(size = 24, color = "black", face = "bold", family = NULL)

#generate combined figure 
ggarrange(
  Abund.p,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(fish_plots[[1]], fish_plots[[2]], 
            nrow = 2, ncol = 1, labels = c("B", "C"),
            font.label = font_labels), 
  nrow = 1, 
  labels = "A",       # Label of the line plot
  font.label = font_labels,
  common.legend = TRUE, legend = "bottom"
) 


#save the plot
ggsave("./Figures/Figure_1-Bacteria_cells.pdf",
       plot = last_plot(),
       units = "mm",
       width = 90, height = 50, 
       scale = 3,
       dpi = 300)


######################################
# Calculate growth rates
######################################
ML_metadata %>% 
  select(Treatment, Bottle, Time,Abundance) %>% 
  filter(Time %in%c(2,21)) %>% 
  group_by(Treatment, Bottle) %>% 
  summarize(generations = log2(Abundance[Time == 21]/ Abundance[Time == 2])) %>% 
  mutate(gen_time = 21/generations,
         growth_rate = 24/gen_time)


#print session info and clean the workspace
sessionInfo()
rm(list = ls())
gc()


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










Bac_abund <- ML_metadata_mean %>% 
              filter(variable =="Abundance") %>% 
              select(Treatment, Time,mean, se) %>% 
              mutate(Time = as.numeric(Time),
                     Method = "Abundance",
                     Taxa = "T") %>%
              dplyr::rename("value_mean" = "mean",
                            "value_se" ="se")

FiSH.mean_new <- FiSH.mean %>%  mutate(Time= as.numeric(Time))



test <- rbind(FiSH.mean_new,Bac_abund)






  
  
  
  