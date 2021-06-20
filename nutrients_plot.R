require(tidyr)
require(dplyr)
require(reshape2)

nutrients<- read.csv("data/DECOMB-Nut_DOC_TDN.csv", h = T) %>% 
  separate(SampleID, into = c("Experiment","Time","Type")) %>% 
  mutate(NO2_cal = NO2 - NO2[Type == "ASW"],
         NO3_cal = NO3 - NO3[Type == "ASW"],
         PO4_cal = PO4 - PO4[Type == "ASW"],
         NH4_cal = NH4 - NH4[Type == "ASW"],
         DOC_cal = DOC - DOC[Type == "ASW"],
         TDN_cal = TDN - TDN[Type == "ASW"],
         Sample = case_when(Type %in% c("C1","C2","C3")~ "Control", 
                            Type %in% c("J1","J2","J3")~ "Jelly"),
         Replicate = gsub("C|J","", Type)) %>% 
  select(Type, Time, Sample, Replicate, NO2_cal,NO3_cal,  PO4_cal, NH4_cal, DOC_cal, TDN_cal) %>% 
  melt()


ggplot(nutrients, aes(x= Time, y = value, colour = Sample, shape = Replicate, group = Type))+
  geom_point()+
  geom_line()+
  facet_wrap(~variable, scales = "free_y")+
  theme_bw()
