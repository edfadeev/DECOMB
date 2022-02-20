require(dplyr)
require(reshape2)
require(ggplot2)


###################
#Plot organic matter and innorganic nutrient concentrations
###################
biochem_raw <- read.csv("F:/My Drive/DECOMB/metadata/for_R/ML_OMandInorganics.csv")

#transpose the table and add type variable
biochem_long<- biochem_raw %>% 
                mutate(NO2_cal = NO2. - NO2.[Sample == "ASW"],
                       NO3_cal = NO3. - NO3.[Sample == "ASW"],
                       PO4_cal = PO43. - PO43.[Sample == "ASW"],
                       NH4_cal = NH4. - NH4.[Sample == "ASW"],
                       DOC_cal = DOC - DOC[Sample == "ASW"],
                       TDN_cal = TDN - TDN[Sample == "ASW"]) %>% 
                melt() %>% 
                mutate(Type = Sample) %>% 
                mutate(Type= case_when(Type %in% c("C1","C2","C3") ~ "Control",
                                       Type %in% c("J1","J2","J3") ~ "Jelly-OM",
                                       TRUE ~ Type)) %>% 
                filter(variable %in% c("NO2_cal","NO3_cal",
                                       "PO4_cal","NH4_cal",
                                       "DOC_cal","TDN_cal")) %>% 
                mutate(value = case_when(value<0 ~ 0,
                                         TRUE ~ value))#if value is negative after blanking force to 0
                      

#plot biochemical paprameters  
biochem_long %>% 
  filter(!Type %in% c("Inoculum","ASW")) %>%
  ggplot(aes(x= Time, y = value, group = interaction(Sample, variable)))+
  geom_point(aes(shape = Type))+
  geom_line(aes(colour = Type))+
  labs(x= "Time", y = "Concentration (umol L-1)")+
  facet_wrap(~variable, scales = "free_y")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom")


###################
#Plot amino acids concentrations
###################
AAs_raw <- read.csv("F:/My Drive/DECOMB/metadata/for_R/AA_compiled.csv")

#transpose the table and add type variable
AAs_long<- AAs_raw %>% 
  melt() %>% 
  filter(!is.na(value)) %>% 
  mutate(Type = Sample) %>% 
  mutate(Type= case_when(Type %in% c("C1","C2","C3") ~ "Control",
                         Type %in% c("J1","J2","J3") ~ "Jelly-OM",
                         TRUE ~ Type))

#plot biochemical parameters  
AAs_long %>% 
  filter(!Type %in% c("Inoculum","ASW")) %>%
  ggplot(aes(x= Time, y = value, group = interaction(Sample, Fraction, variable)))+
  geom_point(aes(shape = Fraction))+
  geom_line(aes(colour = Type))+
  labs(x= "Time", y = "Concentration (umol L-1)")+
  facet_wrap(~variable, scales = "free_y")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom")
