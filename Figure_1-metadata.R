require(dplyr)
require(ggplot2)


#import the dataset
ML_metadata <- read.table("/Users/eduardfadeev/ucloud/Projects/DE-COMB/01_Results/Metadata/for_R/ML_metadata.txt",
                         h = TRUE, sep="\t", dec = ",", blank.lines.skip = TRUE) %>% 
                  filter(!Treatment=="")

######################################
# Plot nutrients, DOC/DON, and cell abundances
######################################
metadata.p <- ML_metadata %>% 
  reshape2::melt(id.vars =c("Time", "Treatment", "Bottle")) %>%
  group_by(Time, Treatment, variable) %>%
  mutate(variable = factor(variable, levels = c("Abundance","TDN","NH4",
                                                   "NO3","NO2","PO43",
                                                   "DOC","DCAA","DFAA")),
          value = as.numeric(value)) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  filter(!is.na(mean),
         Time < 100
         ) %>% 
  ggplot(aes(x = Time, y = mean, group = Treatment, colour = Treatment))+
  geom_point(size = 4, colour = "black")+
  geom_point(size = 3)+
  geom_line(linetype=2)+
  geom_errorbar(aes(ymin= mean-sd, ymax= mean+sd))+
  facet_wrap(~variable, scales = "free_y")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20))

#save the plot
ggsave("./Figures/Figure_1-metadata.pdf",
       plot = metadata.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

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





  
  
  
  