#setwd("T:/Health Behavior and Policy/Faculty/Kimmel/Common/Personnel work/Portia/Rwanda-model/R codes and data for MSM/Data")
setwd("C:/Users/exper/Desktop/Models/RW")
rm(list=ls())

files <- list.files(pattern="*.RData", full.names=TRUE)

data_list <- list()

for (file in files){
load(file)
data_list[[file]] <-  data.frame(N_total, S_total, N_F_total, S_F_total, N_M_total, S_M_total,
             N_F_young_total, S_F_young_total, N_M_young_total, S_M_young_total,
             N_total_all, S_total_all, N_F_total_all, S_F_total_all, N_M_total_all, S_M_total_all,
             NI_total, NI_F_total, NI_M_total,NI_F_young_total, NI_M_young_total,On_ART_grand_total,
             Diagnosed_total,Diagnosed_F_total,Diagnosed_M_total,
             On_ART_total, On_ART_total_F, On_ART_total_M, 
             Suppressed_total, Suppressed_total_F,  Suppressed_total_M,
             Suppressed_total_all, Suppressed_total_F_all,  Suppressed_total_M_all)}

combined_outcome <- bind_rows(data_list)
combined_outcome <- combined_outcome %>% 
   mutate(row=row_number())

combined_outcome$trial = ceiling(combined_outcome$row/324)

save(combined_outcome, file="outcomes.RData")