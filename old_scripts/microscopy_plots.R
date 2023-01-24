require(tidyr)
require(dplyr)
require(reshape2)
require(ggplot2)
require(stringr)
require(ggpattern) #remotes::install_github("coolbutuseless/ggpattern")

#calculation factor
filtration.surface.area <- 176714437.5 #um
counting.surface.area <- 14560 # lxw = 140*104 um
calc.factor <- filtration.surface.area/counting.surface.area
volume <- 1 #mL

###################
#Plot total DAPI in different treatments 
###################

# read abundance data
abund_table_raw <- read.csv("F:/My Drive/DECOMB/Microscopy/for_R/mic_abundance.csv",
                        header = TRUE, sep = ";", dec = ",") 

abund_table_raw<- abund_table_raw %>% 
                  separate(sample_ID, c("Time","Sample"), sep =" ")
                  
#transpose the table and calculate cell abundance from DAPI 
abund_table_long <- abund_table_raw %>% 
  melt(id.vars = c("Time","Sample","grid")) %>% 
  mutate(FOV = gsub("X", "", variable)) %>% 
  filter(!is.na(value)) %>% 
  mutate(DAPI.conc = value*(176714437.5/grid),
         Type = case_when(Sample %in% c("C1","C2","C3") ~ "Control",
                          Sample %in% c("J1","J2","J3") ~ "Jelly-OM",
                                TRUE ~ Sample))
  
#summarize cell abundance per replicate
abund_table_long_rep_mean <- abund_table_long %>% 
                          group_by(Type, Sample, Time) %>% 
                          summarise(DAPI.mn = mean(DAPI.conc),
                                    DAPI.sd = sd(DAPI.conc),
                                    FOVs = n())

#plot cell abundances per replicate
abund_table_long_rep_mean %>% 
  filter(Type !="Ml") %>% 
  ggplot(aes(x=Time, y = DAPI.mn, group = Sample))+
  geom_point()+
  geom_line(aes(colour = Type))+
  geom_errorbar(aes(ymin = DAPI.mn - DAPI.sd, ymax = DAPI.mn + DAPI.sd),
                width = 0.2, position = position_dodge(0.1))+
  labs(x= "Time", y = "Cell abundance (cells mL-1)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom")


###################
#Respiration and production data of different taxa 
###################
# read counts data and adjust labeling
resp_prod_table_raw <- read.csv("F:/My Drive/DECOMB/Microscopy/for_R/resp_and_prod_data_complete.csv",
                            header = TRUE, sep = ";", dec = ",") %>% 
                       mutate(sample_name =str_replace(sample_name,
                                                       "(C1|C2|C3|J1|J2|J3)(T4|T7)(ALT|EUB|GV|PSU)(HPG|RSG)", 
                                                       "\\1_\\2_\\3_\\4")) %>% 
                       mutate(sample_name =str_replace(sample_name,
                                  "(T0)(ALT|EUB|GV|PSU)(HPG|RSG)", 
                                  "\\C0_\\1_\\2_\\3")) %>% 
                        separate(sample_name, into =c("Sample","Time","Probe","Method"), sep ="_")
  
#calculate concentration per FOV
resp_prod_FOV_conc <- resp_prod_table_raw %>% 
  mutate(DAPI.conc = (DAPI*calc.factor)/volume,
         Positive.conc = (DAPI_FITC_FISH*calc.factor)/volume,
         DAPI.FISH.conc = (DAPI_FISH*calc.factor)/volume,
         DAPI.FITC.conc = (DAPI_FITC*calc.factor)/volume,
         Type = case_when(Sample %in% c("C1","C2","C3") ~ "Control",
                          Sample %in% c("J1","J2","J3") ~ "Jelly-OM",
                          TRUE ~ Sample))

#calculate mean concentration for all FOVs in each sample
resp_prod_mean.conc <- resp_prod_FOV_conc %>% 
  group_by(Type, Sample, Time, Probe, Method) %>% 
  summarise(DAPI.mn = mean(DAPI.conc),
            Pos.mean = mean(Positive.conc),
            FITC.mean = mean(DAPI.FITC.conc),
            FISH.mean = mean(DAPI.FISH.conc))

###################
#Combined DAPI plot
###################
#summarize DAPI per treatment from manual counts
abund_table_long_mean <- abund_table_long_rep_mean %>% 
  group_by(Type, Time) %>% 
  summarise(DAPI.mean = mean(DAPI.mn),
            DAPI.se = sd(DAPI.mn)/sqrt(length(DAPI.mn)))


#summarize DAPI per treatment from automatic counts
resp_prod_DAPI_mean <- resp_prod_mean.conc %>% 
                        group_by(Type, Time) %>% 
                        summarise(DAPI.auto.mean = mean(DAPI.mn),
                        DAPI.auto.se = sd(DAPI.mn)/sqrt(length(DAPI.mn)))

#merge DAPi counts
DAPI_table_comb <- abund_table_long_mean %>% 
                    left_join(resp_prod_DAPI_mean, by = c("Type", "Time"))

#plot
DAPI_table_comb %>% 
  filter(Type !="Ml") %>% 
  ggplot(aes(x=Time, y = DAPI.mean, group = Type))+
  geom_point(size = 3)+
  geom_line(aes(colour = Type), size = 1)+
  geom_errorbar(aes(ymin = DAPI.mean - DAPI.se, ymax = DAPI.mean + DAPI.se),
                width = 0.1)+
  labs(x= "Time", y = "Cell abundance (cells mL-1)")+
  geom_point(aes(x=Time, y = DAPI.auto.mean, fill = Type), size = 4, shape = 23)+
  geom_errorbar(aes(ymin = DAPI.auto.mean - DAPI.auto.se, ymax = DAPI.auto.mean + DAPI.auto.se),
                width = 0.1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom")



###################
#Plot total abundance, respiration and production of each taxa 
###################
resp_prod_treat_mean <- resp_prod_mean.conc %>% 
  group_by(Type, Time, Probe, Method) %>% 
  summarise(Pos.conc = mean(Pos.mean),
            FITC.conc = mean(FITC.mean),
            FISH.conc = mean(FISH.mean)) %>% 
  melt()%>% 
  dplyr::rename(mean = value)


resp_prod_treat_se <- resp_prod_mean.conc %>% 
  group_by(Type, Time, Probe, Method) %>% 
  summarise(Pos.conc = sd(Pos.mean)/sqrt(length(Pos.mean)),
            FITC.conc = sd(FITC.mean)/sqrt(length(FITC.mean)),
            FISH.conc = sd(FISH.mean)/sqrt(length(FISH.mean))) %>% 
  melt() %>% 
  dplyr::rename(se = value)

#merge means and se
resp_prod_treat <- left_join(resp_prod_treat_mean, resp_prod_treat_se)

#plot means
resp_prod_treat %>% 
  filter(!Time %in% c("T0"),
         variable != "FITC.conc") %>% 
  ggplot(aes(x=Time, y= mean, fill = variable, group = variable))+
  geom_bar(position = position_dodge(), stat="identity")+
  #geom_col(position = "dodge")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.1,position=position_dodge(.9))+
  facet_grid(Method+Type~Probe)+
  scale_y_continuous(labels=scales::scientific_format())+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom")
  
  
  
  
  
  
  melt() %>% 
  #filter(variable %in% c("Pos.mean")) %>% 
  ggplot(aes(x=Time, y= value, fill = Probe))+
  geom_col(position = "dodge")+
  facet_grid(Method~Sample, scales = "free_y")+
  scale_y_continuous(labels=scales::scientific_format())+
  theme_bw()

#plot means
resp_prod_mean.conc %>% 
  melt() %>% 
  filter(variable %in% c("DAPI.mean")) %>% 
  ggplot(aes(x=Time, y= value, fill = Probe))+
  geom_col(position = "dodge")+
  facet_grid(Method~Sample, scales = "free_y")+
  scale_y_continuous(labels=scales::scientific_format())+
  theme_bw()


#plot EUB means
resp_prod_mean.conc %>% 
  melt() %>% 
  filter(variable %in% c("Pos.mean"),
         Probe == "EUB") %>% 
  ggplot(aes(x=interaction(Type,Time), y= value, fill = Type, colour = Type, group = Sample))+
  geom_point(size = 3)+
  #geom_col(position = "dodge")+
  facet_grid(Method~., scales = "free_y")+
  scale_y_continuous(labels=scales::scientific_format())+
  theme_bw()



#calculate proportion of active cells
resp_prod_mean.prop <- resp_prod_mean.conc %>% 
  mutate(Probe.trt.pos = FITC.mean*Pos.mean, # Probe and Treatment positive of all active cells
         perc.active.cells.ofDapi = (100/DAPI.mn)*FITC.mn) %>% # Treatment positive of all DAPI+