library(matrixStats)
library(gridExtra)
library(ggpubr)

# Assign a vector to store the year variables
year <- seq(from=2004, to=2020, by=1)
color=c("Projected outcomes"="blue","PHIA targets"="red","DHS targets"="purple")
color2=c("Projected outcomes"="blue","UNAIDS targets"="red")
color3=c("Projected outcomes"="blue","Rwanda HMIS targets"="red")

for (t in 1:trial_select){
  Prevalence_total[t,1,] <- (N_total[t,1,])/(N_total[t,1,]+S_total[t,1,])
  Prevalence_total_y[t,] <- Prevalence_total[t,1,][seq(1,length(Prevalence_total[t,1,]),12)]*100}
df_prev <- as.data.frame(Prevalence_total_y[,c(1:17)])
df_prev_long <- t(df_prev)
mean <- rowMeans(df_prev_long)
min <- rowMins(df_prev_long)
max <- rowMaxs(df_prev_long)
Prev <- Target_all$Prevalence*100
Prev_LB <- Target_all$Prevalence_LB*100
Prev_UB <- Target_all$Prevalence_UB*100
DHS <- Target_all$DHS_Prevalence*100
DHS_LB <- Target_all$DHS_Prevalence_LB*100
DHS_UB <- Target_all$DHS_Prevalence_UB*100
df_prev_new <- as.data.frame(cbind(year,mean,min,max,Prev,Prev_LB,Prev_UB,DHS,DHS_LB,DHS_UB))

ragg::agg_tiff("HIV prevalence_total_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_prev_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Prev,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Prev_LB,ymax=Prev_UB,color="PHIA targets"), size = 0.5)+
  geom_point(aes(y=DHS,color="DHS targets"), size = 2)+
  geom_errorbar(aes(ymin=DHS_LB,ymax=DHS_UB,color="DHS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV prevalence,total (15-49)",
       x="Year",
       y="Prevalence (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV prevalence female
for (t in 1:trial_select){
  Prevalence_F_total[t,1,] <- (N_F_total[t,1,])/(N_F_total[t,1,]+S_F_total[t,1,])
  Prevalence_F_total_y[t,] <- Prevalence_F_total[t,1,][seq(1,length(Prevalence_F_total[t,1,]),12)]*100}
df_prev <- as.data.frame(Prevalence_F_total_y[,c(1:17)])
df_prev_long <- t(df_prev)
mean <- rowMeans(df_prev_long)
min <- rowMins(df_prev_long)
max <- rowMaxs(df_prev_long)
Prev <- Target_all$Prevalence_F*100
Prev_LB <- Target_all$Prevalence_LB_F*100
Prev_UB <- Target_all$Prevalence_UB_F*100
DHS <- Target_all$DHS_Prevalence_F*100
DHS_LB <- Target_all$DHS_Prevalence_LB_F*100
DHS_UB <- Target_all$DHS_Prevalence_UB_F*100
df_prev_new <- as.data.frame(cbind(year,mean,min,max,Prev,Prev_LB,Prev_UB,DHS,DHS_LB,DHS_UB))

ragg::agg_tiff("HIV prevalence_female_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_prev_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Prev,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Prev_LB,ymax=Prev_UB,color="PHIA targets"), size = 0.5)+
  geom_point(aes(y=DHS,color="DHS targets"), size = 2)+
  geom_errorbar(aes(ymin=DHS_LB,ymax=DHS_UB,color="DHS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV prevalence,female (15-49)",
       x="Year",
       y="Prevalence (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV prevalence male
for (t in 1:trial_select){
  Prevalence_M_total[t,1,] <- (N_M_total[t,1,])/(N_M_total[t,1,]+S_M_total[t,1,])
  Prevalence_M_total_y[t,] <- Prevalence_M_total[t,1,][seq(1,length(Prevalence_M_total[t,1,]),12)]*100}
df_prev <- as.data.frame(Prevalence_M_total_y[,c(1:17)])
df_prev_long <- t(df_prev)
mean <- rowMeans(df_prev_long)
min <- rowMins(df_prev_long)
max <- rowMaxs(df_prev_long)
Prev <- Target_all$Prevalence_M*100
Prev_LB <- Target_all$Prevalence_LB_M*100
Prev_UB <- Target_all$Prevalence_UB_M*100
DHS <- Target_all$DHS_Prevalence_M*100
DHS_LB <- Target_all$DHS_Prevalence_LB_M*100
DHS_UB <- Target_all$DHS_Prevalence_UB_M*100
df_prev_new <- as.data.frame(cbind(year,mean,min,max,Prev,Prev_LB,Prev_UB,DHS,DHS_LB,DHS_UB))

ragg::agg_tiff("HIV prevalence_male_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_prev_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Prev,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Prev_LB,ymax=Prev_UB,color="PHIA targets"), size = 0.5)+
  geom_point(aes(y=DHS,color="DHS targets"), size = 2)+
  geom_errorbar(aes(ymin=DHS_LB,ymax=DHS_UB,color="DHS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV prevalence,male (15-49)",
       x="Year",
       y="Prevalence (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV prevalence young female
for (t in 1:trial_select){
  Prevalence_F_young_total[t,1,] <- (N_F_young_total[t,1,])/(N_F_young_total[t,1,]+S_F_young_total[t,1,])
  Prevalence_F_young_total_y[t,] <- Prevalence_F_young_total[t,1,][seq(1,length(Prevalence_F_young_total[t,1,]),12)]*100}
df_prev <- as.data.frame(Prevalence_F_young_total_y[,c(1:17)])
df_prev_long <- t(df_prev)
mean <- rowMeans(df_prev_long)
min <- rowMins(df_prev_long)
max <- rowMaxs(df_prev_long)
Prev <- Target_all$Prevalence_young_F*100
Prev_LB <- Target_all$Prevalence_LB_young_F*100
Prev_UB <- Target_all$Prevalence_UB_young_F*100
DHS <- Target_all$DHS_Prevalence_young_F*100
DHS_LB <- Target_all$DHS_Prevalence_LB_young_F*100
DHS_UB <- Target_all$DHS_Prevalence_UB_young_F*100
df_prev_new <- as.data.frame(cbind(year,mean,min,max,Prev,Prev_LB,Prev_UB,DHS,DHS_LB,DHS_UB))

ragg::agg_tiff("HIV prevalence_young female_15-24.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_prev_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Prev,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Prev_LB,ymax=Prev_UB,color="PHIA targets"), size = 0.5)+
  geom_point(aes(y=DHS,color="DHS targets"), size = 2)+
  geom_errorbar(aes(ymin=DHS_LB,ymax=DHS_UB,color="DHS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV prevalence,young female (15-24)",
       x="Year",
       y="Prevalence (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV prevalence young male
for (t in 1:trial_select){
  Prevalence_M_young_total[t,1,] <- (N_M_young_total[t,1,])/(N_M_young_total[t,1,]+S_M_young_total[t,1,])
  Prevalence_M_young_total_y[t,] <- Prevalence_M_young_total[t,1,][seq(1,length(Prevalence_M_young_total[t,1,]),12)]*100}
df_prev <- as.data.frame(Prevalence_M_young_total_y[,c(1:17)])
df_prev_long <- t(df_prev)
mean <- rowMeans(df_prev_long)
min <- rowMins(df_prev_long)
max <- rowMaxs(df_prev_long)
Prev <- Target_all$Prevalence_young_M*100
Prev_LB <- Target_all$Prevalence_LB_young_M*100
Prev_UB <- Target_all$Prevalence_UB_young_M*100
DHS <- Target_all$DHS_Prevalence_young_M*100
DHS_LB <- Target_all$DHS_Prevalence_LB_young_M*100
DHS_UB <- Target_all$DHS_Prevalence_UB_young_M*100
df_prev_new <- as.data.frame(cbind(year,mean,min,max,Prev,Prev_LB,Prev_UB,DHS,DHS_LB,DHS_UB))

ragg::agg_tiff("HIV prevalence_male_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_prev_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Prev,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Prev_LB,ymax=Prev_UB,color="PHIA targets"), size = 0.5)+
  geom_point(aes(y=DHS,color="DHS targets"), size = 2)+
  geom_errorbar(aes(ymin=DHS_LB,ymax=DHS_UB,color="DHS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV prevalence,young male (15-24)",
       x="Year",
       y="Prevalence (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV prevalence, total 15-64
for (t in 1:trial_select){
  Prevalence_total[t,1,] <- (N_total_all[t,1,])/(N_total_all[t,1,]+S_total_all[t,1,])
  Prevalence_total_y[t,] <- Prevalence_total[t,1,][seq(1,length(Prevalence_total[t,1,]),12)]*100}
df_prev <- as.data.frame(Prevalence_total_y[,c(1:17)])
df_prev_long <- t(df_prev)
mean <- rowMeans(df_prev_long)
min <- rowMins(df_prev_long)
max <- rowMaxs(df_prev_long)
Prev <- Target_new$Prevalence*100
Prev_LB <- Target_new$Prevalence_LB*100
Prev_UB <- Target_new$Prevalence_UB*100
df_prev_new <- as.data.frame(cbind(year,mean,min,max,Prev,Prev_LB,Prev_UB))

ragg::agg_tiff("HIV prevalence_total_15-64.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_prev_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Prev,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Prev_LB,ymax=Prev_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV prevalence,total (15-64)",
       x="Year",
       y="Prevalence (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV prevalence female 15-64
for (t in 1:trial_select){
  Prevalence_F_total[t,1,] <- (N_F_total_all[t,1,])/(N_F_total_all[t,1,]+S_F_total_all[t,1,])
  Prevalence_F_total_y[t,] <- Prevalence_F_total[t,1,][seq(1,length(Prevalence_F_total[t,1,]),12)]*100}
df_prev <- as.data.frame(Prevalence_F_total_y[,c(1:17)])
df_prev_long <- t(df_prev)
mean <- rowMeans(df_prev_long)
min <- rowMins(df_prev_long)
max <- rowMaxs(df_prev_long)
Prev <- Target_new$Prevalence_F*100
Prev_LB <- Target_new$Prevalence_LB_F*100
Prev_UB <- Target_new$Prevalence_UB_F*100
df_prev_new <- as.data.frame(cbind(year,mean,min,max,Prev,Prev_LB,Prev_UB))

ragg::agg_tiff("HIV prevalence_female_15-64.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_prev_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Prev,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Prev_LB,ymax=Prev_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV prevalence,female (15-64)",
       x="Year",
       y="Prevalence (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV prevalence male, 15-64
for (t in 1:trial_select){
  Prevalence_M_total[t,1,] <- (N_M_total_all[t,1,])/(N_M_total_all[t,1,]+S_F_total_all[t,1,])
  Prevalence_M_total_y[t,] <- Prevalence_M_total[t,1,][seq(1,length(Prevalence_M_total[t,1,]),12)]*100}
df_prev <- as.data.frame(Prevalence_M_total_y[,c(1:17)])
df_prev_long <- t(df_prev)
mean <- rowMeans(df_prev_long)
min <- rowMins(df_prev_long)
max <- rowMaxs(df_prev_long)
Prev <- Target_new$Prevalence_M*100
Prev_LB <- Target_new$Prevalence_LB_M*100
Prev_UB <- Target_new$Prevalence_UB_M*100
df_prev_new <- as.data.frame(cbind(year,mean,min,max,Prev,Prev_LB,Prev_UB))

ragg::agg_tiff("HIV prevalence_male_15-64.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_prev_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Prev,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Prev_LB,ymax=Prev_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV prevalence,male (15-64)",
       x="Year",
       y="Prevalence (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV incidence, total 15-44
for (t in 1:trial_select) {
  NI_total_y[t,] <- rollapply(NI_total[t,1,],12,sum,by=12)
  S_total_y[t,] <- S_total[t,1,][seq(1, length(S_total[t,1,]), 12)]
  N_total_y[t,] <- N_total[t,1,][seq(1, length(N_total[t,1,]), 12)]
  Incidence_total_y[t,] <- NI_total_y[t,]/S_total_y[t,]*1000}
df_inci <- as.data.frame(Incidence_total_y[,c(1:17)])
df_inci_long <- t(df_inci)
mean <- rowMeans(df_inci_long)
min <- rowMins(df_inci_long)
max <- rowMaxs(df_inci_long)
Inci <- Target_all$Incidence*1000
Inci_LB <- Target_all$Incidence_LB*1000
Inci_UB <- Target_all$Incidence_UB*1000
df_inci_new <- as.data.frame(cbind(year,mean,min,max,Inci,Inci_LB,Inci_UB))

ragg::agg_tiff("HIV incidence_total_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_inci_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Inci,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Inci_LB,ymax=Inci_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV incidence,total (15-49)",
       x="Year",
       y="Incidence (per 1,000 population)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV incidence female
for (t in 1:trial_select) {
  NI_F_total_y[t,] <- rollapply(NI_F_total[t,1,],12,sum,by=12)
  S_F_total_y[t,] <- S_F_total[t,1,][seq(1, length(S_F_total[t,1,]), 12)]
  N_F_total_y[t,] <- N_F_total[t,1,][seq(1, length(N_F_total[t,1,]), 12)]
  Incidence_F_total_y[t,] <- NI_F_total_y[t,]/S_F_total_y[t,]*1000}
df_inci <- as.data.frame(Incidence_F_total_y[,c(1:17)])
df_inci_long <- t(df_inci)
mean <- rowMeans(df_inci_long)
min <- rowMins(df_inci_long)
max <- rowMaxs(df_inci_long)
Inci <- Target_all$Incidence_F*1000
Inci_LB <- Target_all$Incidence_F_LB*1000
Inci_UB <- Target_all$Incidence_F_UB*1000
df_inci_new <- as.data.frame(cbind(year,mean,min,max,Inci,Inci_LB,Inci_UB))

ragg::agg_tiff("HIV incidence_female_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_inci_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Inci,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Inci_LB,ymax=Inci_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV incidence,female (15-49)",
       x="Year",
       y="Incidence (per 1,000 population)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV incidence male
for (t in 1:trial_select) {
  NI_M_total_y[t,] <- rollapply(NI_M_total[t,1,],12,sum,by=12)
  S_M_total_y[t,] <- S_M_total[t,1,][seq(1, length(S_M_total[t,1,]), 12)]
  N_M_total_y[t,] <- N_M_total[t,1,][seq(1, length(N_M_total[t,1,]), 12)]
  Incidence_M_total_y[t,] <- NI_M_total_y[t,]/S_M_total_y[t,]*1000}
df_inci <- as.data.frame(Incidence_M_total_y[,c(1:17)])
df_inci_long <- t(df_inci)
mean <- rowMeans(df_inci_long)
min <- rowMins(df_inci_long)
max <- rowMaxs(df_inci_long)
Inci <- Target_all$Incidence_M*1000
Inci_LB <- Target_all$Incidence_M_LB*1000
Inci_UB <- Target_all$Incidence_M_UB*1000
df_inci_new <- as.data.frame(cbind(year,mean,min,max,Inci,Inci_LB,Inci_UB))

ragg::agg_tiff("HIV incidence_male_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_inci_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Inci,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Inci_LB,ymax=Inci_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV incidence,male (15-49)",
       x="Year",
       y="Incidence (per 1,000 population)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV incidence young female
for (t in 1:trial_select) {
  NI_F_young_total_y[t,] <- rollapply(NI_F_young_total[t,1,],12,sum,by=12)
  S_F_young_total_y[t,] <- S_F_young_total[t,1,][seq(1, length(S_F_young_total[t,1,]), 12)]
  N_F_young_total_y[t,] <- N_F_young_total[t,1,][seq(1, length(N_F_young_total[t,1,]), 12)]
  Incidence_F_young_total_y[t,] <- NI_F_young_total_y[t,]/S_F_young_total_y[t,]*1000}
df_inci <- as.data.frame(Incidence_F_young_total_y[,c(1:17)])
df_inci_long <- t(df_inci)
mean <- rowMeans(df_inci_long)
min <- rowMins(df_inci_long)
max <- rowMaxs(df_inci_long)
Inci <- Target_all$Incidence_F_young*1000
Inci_LB <- Target_all$Incidence_F_young_LB*1000
Inci_UB <- Target_all$Incidence_F_young_UB*1000
df_inci_new <- as.data.frame(cbind(year,mean,min,max,Inci,Inci_LB,Inci_UB))

ragg::agg_tiff("HIV incidence_young_female_15-24.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_inci_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Inci,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Inci_LB,ymax=Inci_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV incidence,female (15-24)",
       x="Year",
       y="Incidence (per 1,000 population)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# HIV incidence young male
for (t in 1:trial_select) {
  NI_M_young_total_y[t,] <- rollapply(NI_M_young_total[t,1,],12,sum,by=12)
  S_M_young_total_y[t,] <- S_M_young_total[t,1,][seq(1, length(S_M_young_total[t,1,]), 12)]
  N_M_young_total_y[t,] <- N_M_young_total[t,1,][seq(1, length(N_M_young_total[t,1,]), 12)]
  Incidence_M_young_total_y[t,] <- NI_M_young_total_y[t,]/S_M_young_total_y[t,]*1000}
df_inci <- as.data.frame(Incidence_M_young_total_y[,c(1:17)])
df_inci_long <- t(df_inci)
mean <- rowMeans(df_inci_long)
min <- rowMins(df_inci_long)
max <- rowMaxs(df_inci_long)
Inci <- Target_all$Incidence_M*1000
Inci_LB <- Target_all$Incidence_M_LB*1000
Inci_UB <- Target_all$Incidence_M_UB*1000
df_inci_new <- as.data.frame(cbind(year,mean,min,max,Inci,Inci_LB,Inci_UB))

ragg::agg_tiff("HIV incidence_young_male_15-24.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_inci_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Inci,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Inci_LB,ymax=Inci_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="HIV incidence,male (15-24)",
       x="Year",
       y="Incidence (per 1,000 population)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Number on ART 
for (t in 1:trial_select){
  On_ART_grand_total_y[t,] <- On_ART_grand_total[t,1,][seq(1,length(On_ART_grand_total[t,1,]),12)]}
df_art <- as.data.frame(On_ART_grand_total_y[,c(1:17)])
df_art_long <- t(df_art)
mean <- rowMeans(df_art_long)
min <- rowMins(df_art_long)
max <- rowMaxs(df_art_long)
ART <- Target_all$N_ART
df_art_new <- as.data.frame(cbind(year,mean,min,max,ART))

ragg::agg_tiff("Number on ART_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_art_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=ART,color="Rwanda HMIS targets"), size = 2)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color3)+
  ylim(0,NA)+
  labs(title="Number of people on ART (15-49)",
       x="Year",
       y="Numbers")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent diagnosed, total
for (t in 1:trial_select) {
  Diagnosed_total_y[t,] <- Diagnosed_total[t,1,][seq(1, length(Diagnosed_total[t,1,]), 12)]
  Rate_diagnosed_total_y[t,] <- Diagnosed_total_y[t,]/(N_total_y[t,])*100}
df_diag <- as.data.frame(Rate_diagnosed_total_y[,c(1:17)])
df_diag_long <- t(df_diag)
mean <- rowMeans(df_diag_long)
min <- rowMins(df_diag_long)
max <- rowMaxs(df_diag_long)
Diag <- Target_all$Percent_diagnosed_total
Diag_LB <- Target_all$Percent_diagnosed_total_LB
Diag_UB <- Target_all$Percent_diagnosed_total_UB
df_diag_new <- as.data.frame(cbind(year,mean,min,max,Diag,Diag_LB,Diag_UB))

ragg::agg_tiff("Percent diagnosed_total_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_diag_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Diag,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Diag_LB,ymax=Diag_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent diagnosed,total (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent diagnosed, female
for (t in 1:trial_select) {
  Diagnosed_F_total_y[t,] <- Diagnosed_F_total[t,1,][seq(1, length(Diagnosed_F_total[t,1,]), 12)]
  Rate_diagnosed_F_total_y[t,] <- Diagnosed_F_total_y[t,]/(N_F_total_y[t,])*100}
df_diag <- as.data.frame(Rate_diagnosed_F_total_y[,c(1:17)])
df_diag_long <- t(df_diag)
mean <- rowMeans(df_diag_long)
min <- rowMins(df_diag_long)
max <- rowMaxs(df_diag_long)
Diag <- Target_all$Percent_diagnosed_Female
Diag_LB <- Target_all$Percent_diagnosed_Female_LB
Diag_UB <- Target_all$Percent_diagnosed_Female_UB
df_diag_new <- as.data.frame(cbind(year,mean,min,max,Diag,Diag_LB,Diag_UB))

ragg::agg_tiff("Percent diagnosed_female_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_diag_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Diag,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Diag_LB,ymax=Diag_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent diagnosed,female (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent diagnosed, male
for (t in 1:trial_select) {
  Diagnosed_M_total_y[t,] <- Diagnosed_M_total[t,1,][seq(1, length(Diagnosed_M_total[t,1,]), 12)]
  Rate_diagnosed_M_total_y[t,] <- Diagnosed_M_total_y[t,]/(N_M_total_y[t,])*100}
df_diag <- as.data.frame(Rate_diagnosed_M_total_y[,c(1:17)])
df_diag_long <- t(df_diag)
mean <- rowMeans(df_diag_long)
min <- rowMins(df_diag_long)
max <- rowMaxs(df_diag_long)
Diag <- Target_all$Percent_diagnosed_Male
Diag_LB <- Target_all$Percent_diagnosed_Male_LB
Diag_UB <- Target_all$Percent_diagnosed_Male_UB
df_diag_new <- as.data.frame(cbind(year,mean,min,max,Diag,Diag_LB,Diag_UB))

ragg::agg_tiff("Percent diagnosed_male_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_diag_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=Diag,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=Diag_LB,ymax=Diag_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent diagnosed,male (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent on ART among diagnosed
for (t in 1:trial_select){
  On_ART_p_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]/Diagnosed_total_y[t,]*100}
df_art <- as.data.frame( On_ART_p_y[,c(1:17)])
df_art_long <- t(df_art)
mean <- rowMeans(df_art_long)
min <- rowMins(df_art_long)
max <- rowMaxs(df_art_long)
ART <- Target_new$P_ART_All
ART_LB <- Target_new$P_ART_All_LB
ART_UB <- Target_new$P_ART_All_UB
df_art_new <- as.data.frame(cbind(year,mean,min,max,ART,ART_LB,ART_UB))

ragg::agg_tiff("Percent ART among diagnosed_total_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_art_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=ART,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=ART_LB,ymax=ART_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent on ART among diagnosed,total (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()


# Percent on ART female
for (t in 1:trial_select){
  On_ART_p_F_y[t,] <- On_ART_total_F[t,1,][seq(1,length(On_ART_total_F[t,1,]),12)]/Diagnosed_F_total_y[t,]*100}
df_art <- as.data.frame( On_ART_p_F_y[,c(1:17)])
df_art_long <- t(df_art)
mean <- rowMeans(df_art_long)
min <- rowMins(df_art_long)
max <- rowMaxs(df_art_long)
ART <- Target_new$P_ART_All_F
ART_LB <- Target_new$P_ART_All_F_LB
ART_UB <- Target_new$P_ART_All_F_UB
df_art_new <- as.data.frame(cbind(year,mean,min,max,ART,ART_LB,ART_UB))

ragg::agg_tiff("Percent ART among diagnosed_female_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_art_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=ART,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=ART_LB,ymax=ART_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent on ART among diagnosed,female (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent on ART male
for (t in 1:trial_select){
  On_ART_p_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]/Diagnosed_M_total_y[t,]*100}
df_art <- as.data.frame( On_ART_p_M_y[,c(1:17)])
df_art_long <- t(df_art)
mean <- rowMeans(df_art_long)
min <- rowMins(df_art_long)
max <- rowMaxs(df_art_long)
ART <- Target_new$P_ART_All_M
ART_LB <- Target_new$P_ART_All_M_LB
ART_UB <- Target_new$P_ART_All_M_UB
df_art_new <- as.data.frame(cbind(year,mean,min,max,ART,ART_LB,ART_UB))

ragg::agg_tiff("Percent ART among diagnosed_male_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_art_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=ART,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=ART_LB,ymax=ART_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent on ART among diagnosed,male (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent viral suppression (among ART)
for (t in 1:trial_select){
  On_ART_total_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
  Suppressed_total_y[t,] <- Suppressed_total[t,1,][seq(1,length(Suppressed_total[t,1,]),12)]
  P_suppressed_art_total_y[t,] <- Suppressed_total_y[t,]/On_ART_total_y[t,]*100}
df_vs_art <- as.data.frame(P_suppressed_art_total_y[,c(1:17)])
df_vs_art_long <- t(df_vs_art)
mean <- rowMeans(df_vs_art_long)
min <- rowMins(df_vs_art_long)
max <- rowMaxs(df_vs_art_long)
VS_ART <- Target_all$Percent_VS_on_ART_total
VS_ART_LB <- Target_all$Percent_VS_on_ART_total_LB
VS_ART_UB <- Target_all$Percent_VS_on_ART_total_UB
df_vs_art_new <- as.data.frame(cbind(year,mean,min,max,VS_ART,VS_ART_LB,VS_ART_UB))

ragg::agg_tiff("Percent VS among ART_total_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_vs_art_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=VS_ART,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=VS_ART_LB,ymax=VS_ART_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent suppressed among ART,total (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent viral suppression female (among ART)
for (t in 1:trial_select){
  On_ART_total_F_y[t,] <- On_ART_total_F[t,1,][seq(1,length(On_ART_total_F[t,1,]),12)]
  Suppressed_total_F_y[t,] <- Suppressed_total_F[t,1,][seq(1,length(Suppressed_total_F[t,1,]),12)]
  P_suppressed_art_F_total_y[t,] <- Suppressed_total_F_y[t,]/On_ART_total_F_y[t,]*100}
df_vs_art <- as.data.frame(P_suppressed_art_F_total_y[,c(1:17)])
df_vs_art_long <- t(df_vs_art)
mean <- rowMeans(df_vs_art_long)
min <- rowMins(df_vs_art_long)
max <- rowMaxs(df_vs_art_long)
VS_ART <- Target_all$Percent_VS_on_ART_Female
VS_ART_LB <- Target_all$Percent_VS_on_ART_Female_LB
VS_ART_UB <- Target_all$Percent_VS_on_ART_Female_UB
df_vs_art_new <- as.data.frame(cbind(year,mean,min,max,VS_ART,VS_ART_LB,VS_ART_UB))

ragg::agg_tiff("Percent VS among ART_female_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_vs_art_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=VS_ART,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=VS_ART_LB,ymax=VS_ART_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent suppressed among ART,female (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent viral suppression male (among ART)
for (t in 1:trial_select){
  On_ART_total_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]
  Suppressed_total_M_y[t,] <- Suppressed_total_M[t,1,][seq(1,length(Suppressed_total_M[t,1,]),12)]
  P_suppressed_art_M_total_y[t,] <- Suppressed_total_M_y[t,]/On_ART_total_M_y[t,]*100}
df_vs_art <- as.data.frame(P_suppressed_art_M_total_y[,c(1:17)])
df_vs_art_long <- t(df_vs_art)
mean <- rowMeans(df_vs_art_long)
min <- rowMins(df_vs_art_long)
max <- rowMaxs(df_vs_art_long)
VS_ART <- Target_all$Percent_VS_on_ART_Male
VS_ART_LB <- Target_all$Percent_VS_on_ART_Male_LB
VS_ART_UB <- Target_all$Percent_VS_on_ART_Male_UB
df_vs_art_new <- as.data.frame(cbind(year,mean,min,max,VS_ART,VS_ART_LB,VS_ART_UB))

ragg::agg_tiff("Percent VS among ART_male_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_vs_art_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=VS_ART,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=VS_ART_LB,ymax=VS_ART_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent suppressed among ART,male (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent viral suppression (among ALL)
for (t in 1:trial_select){
  Suppressed_total_y[t,] <- Suppressed_total[t,1,][seq(1,length(Suppressed_total[t,1,]),12)]
  P_suppressed_total_y[t,] <- Suppressed_total_y[t,]/N_total_y[t,]*100}
df_vs <- as.data.frame(P_suppressed_total_y[,c(1:17)])
df_vs_long <- t(df_vs)
mean <- rowMeans(df_vs_long)
min <- rowMins(df_vs_long)
max <- rowMaxs(df_vs_long)
VS <- Target_all$Percent_VS_total
VS_LB <- Target_all$Percent_VS_total_LB
VS_UB <- Target_all$Percent_VS_total_UB
df_vs_new <- as.data.frame(cbind(year,mean,min,max,VS,VS_LB,VS_UB))

ragg::agg_tiff("Percent VS among living with HIV_total_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_vs_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=VS,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=VS_LB,ymax=VS_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent suppressed among people living with HIV,total (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent viral suppression female
for (t in 1:trial_select){
  Suppressed_total_F_y[t,] <- Suppressed_total_F[t,1,][seq(1,length(Suppressed_total_F[t,1,]),12)]
  P_suppressed_F_total_y[t,] <- Suppressed_total_F_y[t,]/N_F_total_y[t,]*100}
df_vs <- as.data.frame(P_suppressed_F_total_y[,c(1:17)])
df_vs_long <- t(df_vs)
mean <- rowMeans(df_vs_long)
min <- rowMins(df_vs_long)
max <- rowMaxs(df_vs_long)
VS <- Target_all$Percent_VS_Female
VS_LB <- Target_all$Percent_VS_Female_LB
VS_UB <- Target_all$Percent_VS_Female_UB
df_vs_new <- as.data.frame(cbind(year,mean,min,max,VS,VS_LB,VS_UB))

ragg::agg_tiff("Percent VS among living with HIV_female_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_vs_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=VS,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=VS_LB,ymax=VS_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent suppressed among people living with HIV,female (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent viral suppression male
for (t in 1:trial_select){
  #On_ART_total_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]
  Suppressed_total_M_y[t,] <- Suppressed_total_M[t,1,][seq(1,length(Suppressed_total_M[t,1,]),12)]
  P_suppressed_M_total_y[t,] <- Suppressed_total_M_y[t,]/N_M_total_y[t,]*100}
df_vs <- as.data.frame(P_suppressed_M_total_y[,c(1:17)])
df_vs_long <- t(df_vs)
mean <- rowMeans(df_vs_long)
min <- rowMins(df_vs_long)
max <- rowMaxs(df_vs_long)
VS <- Target_all$Percent_VS_Male
VS_LB <- Target_all$Percent_VS_Male_LB
VS_UB <- Target_all$Percent_VS_Male_UB
df_vs_new <- as.data.frame(cbind(year,mean,min,max,VS,VS_LB,VS_UB))

ragg::agg_tiff("Percent VS among living with HIV_male_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_vs_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=VS,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=VS_LB,ymax=VS_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent suppressed among people living with HIV,male (15-49)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent viral suppression (among ALL)
for (t in 1:trial_select){
  #On_ART_total_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
  N_total_y[t,] <- N_total_all[t,1,][seq(1, length(N_total_all[t,1,]), 12)]
  Suppressed_total_y[t,] <- Suppressed_total_all[t,1,][seq(1,length(Suppressed_total_all[t,1,]),12)]
  P_suppressed_total_y[t,] <- Suppressed_total_y[t,]/N_total_y[t,]*100}
df_vs <- as.data.frame(P_suppressed_total_y[,c(1:17)])
df_vs_long <- t(df_vs)
mean <- rowMeans(df_vs_long)
min <- rowMins(df_vs_long)
max <- rowMaxs(df_vs_long)
VS <- Target_new$Percent_VS_total
VS_LB <- Target_new$Percent_VS_total_LB
VS_UB <- Target_new$Percent_VS_total_UB
df_vs_new <- as.data.frame(cbind(year,mean,min,max,VS,VS_LB,VS_UB))

ragg::agg_tiff("Percent VS among living with HIV_total_15-64.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_vs_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=VS,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=VS_LB,ymax=VS_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent suppressed among people living with HIV,total (15-64)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent viral suppression female
for (t in 1:trial_select){
  N_F_total_y[t,] <- N_F_total_all[t,1,][seq(1, length(N_F_total_all[t,1,]), 12)]
  Suppressed_total_F_y[t,] <- Suppressed_total_F[t,1,][seq(1,length(Suppressed_total_F[t,1,]),12)]
  P_suppressed_F_total_y[t,] <- Suppressed_total_F_y[t,]/N_F_total_y[t,]*100}
df_vs <- as.data.frame(P_suppressed_F_total_y[,c(1:17)])
df_vs_long <- t(df_vs)
mean <- rowMeans(df_vs_long)
min <- rowMins(df_vs_long)
max <- rowMaxs(df_vs_long)
VS <- Target_new$Percent_VS_Female
VS_LB <- Target_new$Percent_VS_Female_LB
VS_UB <- Target_new$Percent_VS_Female_UB
df_vs_new <- as.data.frame(cbind(year,mean,min,max,VS,VS_LB,VS_UB))

ragg::agg_tiff("Percent VS among living with HIV_female_15-64.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_vs_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=VS,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=VS_LB,ymax=VS_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent suppressed among people living with HIV,female (15-64)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Percent viral suppression male
for (t in 1:trial_select){
  #On_ART_total_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]
  N_M_total_y[t,] <- N_M_total_all[t,1,][seq(1, length(N_M_total_all[t,1,]), 12)]
  Suppressed_total_M_y[t,] <- Suppressed_total_M[t,1,][seq(1,length(Suppressed_total_M[t,1,]),12)]
  P_suppressed_M_total_y[t,] <- Suppressed_total_M_y[t,]/N_M_total_y[t,]*100}
df_vs <- as.data.frame(P_suppressed_M_total_y[,c(1:17)])
df_vs_long <- t(df_vs)
mean <- rowMeans(df_vs_long)
min <- rowMins(df_vs_long)
max <- rowMaxs(df_vs_long)
VS <- Target_new$Percent_VS_Male
VS_LB <- Target_new$Percent_VS_Male_LB
VS_UB <- Target_new$Percent_VS_Male_UB
df_vs_new <- as.data.frame(cbind(year,mean,min,max,VS,VS_LB,VS_UB))

ragg::agg_tiff("Percent VS among living with HIV_male_15-64.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_vs_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=VS,color="PHIA targets"), size = 2)+
  geom_errorbar(aes(ymin=VS_LB,ymax=VS_UB,color="PHIA targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color)+
  ylim(0,NA)+
  labs(title="Percent suppressed among people living with HIV,male (15-64)",
       x="Year",
       y="Percentage (%)")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# new infections
for (t in 1:trial_select) {
  NI_total_y[t,] <- rollapply(NI_total[t,1,],12,sum,by=12)}
df_ni <- as.data.frame(NI_total_y[,c(1:17)])
df_ni_long <- t(df_ni)
mean <- rowMeans(df_ni_long)
min <- rowMins(df_ni_long)
max <- rowMaxs(df_ni_long)
NI <- Target_all$N_NEW_ALL
NI_LB <- Target_all$N_NEW_ALL_LB
NI_UB <- Target_all$N_NEW_ALL_UB
df_ni_new <- as.data.frame(cbind(year,mean,min,max,NI,NI_LB,NI_UB))

ragg::agg_tiff("Corroboration_New infections_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_ni_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=NI,color="UNAIDS targets"), size = 2)+
  geom_errorbar(aes(ymin=NI_LB,ymax=NI_UB,color="UNAIDS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color2)+
  ylim(0,NA)+
  labs(title="Number of new infections,total (15-49)",
       x="Year",
       y="Numbers")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# new infections females
for (t in 1:trial_select) {
  NI_F_total_y[t,] <- rollapply(NI_F_total[t,1,],12,sum,by=12)}
df_ni <- as.data.frame(NI_F_total_y[,c(1:17)])
df_ni_long <- t(df_ni)
mean <- rowMeans(df_ni_long)
min <- rowMins(df_ni_long)
max <- rowMaxs(df_ni_long)
NI <- Target_all$N_NEW_Female
NI_LB <- Target_all$N_NEW_Female_LB
NI_UB <- Target_all$N_NEW_Female_UB
df_ni_new <- as.data.frame(cbind(year,mean,min,max,NI,NI_LB,NI_UB))

ragg::agg_tiff("Corroboration_New infections_female_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_ni_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=NI,color="UNAIDS targets"), size = 2)+
  geom_errorbar(aes(ymin=NI_LB,ymax=NI_UB,color="UNAIDS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color2)+
  ylim(0,NA)+
  labs(title="Number of new infections,female (15-49)",
       x="Year",
       y="Numbers")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()


# new infections males
for (t in 1:trial_select) {
  NI_M_total_y[t,] <- rollapply(NI_M_total[t,1,],12,sum,by=12)}
df_ni <- as.data.frame(NI_M_total_y[,c(1:17)])
df_ni_long <- t(df_ni)
mean <- rowMeans(df_ni_long)
min <- rowMins(df_ni_long)
max <- rowMaxs(df_ni_long)
NI <- Target_all$N_NEW_Male
NI_LB <- Target_all$N_NEW_Male_LB
NI_UB <- Target_all$N_NEW_Male_UB
df_ni_new <- as.data.frame(cbind(year,mean,min,max,NI,NI_LB,NI_UB))

ragg::agg_tiff("Corroboration_New infections_male_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_ni_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=NI,color="UNAIDS targets"), size = 2)+
  geom_errorbar(aes(ymin=NI_LB,ymax=NI_UB,color="UNAIDS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color2)+
  ylim(0,NA)+
  labs(title="Number of new infections,male (15-49)",
       x="Year",
       y="Numbers")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()


# new infections young (15-24)
for (t in 1:trial_select) {
  NI_young_total_y[t,] <- rollapply(NI_young_total[t,1,],12,sum,by=12)}
df_ni <- as.data.frame(NI_young_total_y[,c(1:17)])
df_ni_long <- t(df_ni)
mean <- rowMeans(df_ni_long)
min <- rowMins(df_ni_long)
max <- rowMaxs(df_ni_long)
NI <- Target_all$N_NEW_Young
NI_LB <- Target_all$N_NEW_Young_LB
NI_UB <- Target_all$N_NEW_Young_UB
df_ni_new <- as.data.frame(cbind(year,mean,min,max,NI,NI_LB,NI_UB))

ragg::agg_tiff("Corroboration_New infections_youth_15-24.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_ni_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=NI,color="UNAIDS targets"), size = 2)+
  geom_errorbar(aes(ymin=NI_LB,ymax=NI_UB,color="UNAIDS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color2)+
  ylim(0,NA)+
  labs(title="Number of new infections,youth (15-24)",
       x="Year",
       y="Numbers")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()


# Number with HIV
for (t in 1:trial_select) {
  N_total_y[t,] <- N_total[t,1,][seq(1,length(N_total[t,1,]),12)]}
df_n <- as.data.frame(N_total_y[,c(1:17)])
df_n_long <- t(df_n)
mean <- rowMeans(df_n_long)
min <- rowMins(df_n_long)
max <- rowMaxs(df_n_long)
N <- Target_all$N_HIV_ALL
N_LB <- Target_all$N_HIV_ALL_LB
N_UB <- Target_all$N_HIV_ALL_UB
df_n_new <- as.data.frame(cbind(year,mean,min,max,N,N_LB,N_UB))

ragg::agg_tiff("Corroboration_Living with HIV_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_n_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=N,color="UNAIDS targets"), size = 2)+
  geom_errorbar(aes(ymin=N_LB,ymax=N_UB,color="UNAIDS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color2)+
  ylim(0,NA)+
  labs(title="Number of people living with HIV,total (15-49)",
       x="Year",
       y="Numbers")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()


# Number with HIV females
for (t in 1:trial_select) {
  N_F_total_y[t,] <- N_F_total[t,1,][seq(1,length(N_F_total[t,1,]),12)]}
df_n <- as.data.frame(N_F_total_y[,c(1:17)])
df_n_long <- t(df_n)
mean <- rowMeans(df_n_long)
min <- rowMins(df_n_long)
max <- rowMaxs(df_n_long)
N <- Target_all$N_HIV_Female
N_LB <- Target_all$N_HIV_Female_LB
N_UB <- Target_all$N_HIV_Female_UB
df_n_new <- as.data.frame(cbind(year,mean,min,max,N,N_LB,N_UB))

ragg::agg_tiff("Corroboration_Living with HIV_female_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_n_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=N,color="UNAIDS targets"), size = 2)+
  geom_errorbar(aes(ymin=N_LB,ymax=N_UB,color="UNAIDS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color2)+
  ylim(0,NA)+
  labs(title="Number of people living with HIV,female (15-49)",
       x="Year",
       y="Numbers")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()


# Number with HIV males
for (t in 1:trial_select) {
  N_M_total_y[t,] <- N_M_total[t,1,][seq(1,length(N_M_total[t,1,]),12)]}
df_n <- as.data.frame(N_M_total_y[,c(1:17)])
df_n_long <- t(df_n)
mean <- rowMeans(df_n_long)
min <- rowMins(df_n_long)
max <- rowMaxs(df_n_long)
N <- Target_all$N_HIV_Male
N_LB <- Target_all$N_HIV_Male_LB
N_UB <- Target_all$N_HIV_Male_UB
df_n_new <- as.data.frame(cbind(year,mean,min,max,N,N_LB,N_UB))

ragg::agg_tiff("Corroboration_Living with HIV_male_15-49.tiff",units="in",height=5,width=6,res=300)
ggplot(data=df_n_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=N,color="UNAIDS targets"), size = 2)+
  geom_errorbar(aes(ymin=N_LB,ymax=N_UB,color="UNAIDS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color2)+
  ylim(0,NA)+
  labs(title="Number of people living with HIV,male (15-49)",
       x="Year",
       y="Numbers")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

# Number with HIV young (15-24)
for (t in 1:trial_select) {
  N_young_total_y[t,] <- N_young_total[t,1,][seq(1,length(N_young_total[t,1,]),12)]}
df_n <- as.data.frame(N_young_total_y[,c(1:17)])
df_n_long <- t(df_n)
mean <- rowMeans(df_n_long)
min <- rowMins(df_n_long)
max <- rowMaxs(df_n_long)
N <- Target_all$N_HIV_Young
N_LB <- Target_all$N_HIV_Young_LB
N_UB <- Target_all$N_HIV_Young_UB
df_n_new <- as.data.frame(cbind(year,mean,min,max,N,N_LB,N_UB))

ragg::agg_ragg::agg_tiff("Corroboration_Living with HIV_youth_15-24.tiff",units="in",height=5,width=6,res=300,units="in",height=5,width=6,res=300)
ggplot(data=df_n_new, aes(x=year))+
  geom_line(aes(y=mean, color="Projected outcomes"), size = 0.5)+
  geom_ribbon(aes(ymin=min, ymax=max, color="Projected outcomes"), fill="blue", alpha=0.1, linetype = "dotted")+
  geom_point(aes(y=N,color="UNAIDS targets"), size = 2)+
  geom_errorbar(aes(ymin=N_LB,ymax=N_UB,color="UNAIDS targets"), size = 0.5)+
  scale_color_manual(name="Legend", breaks=waiver(), values=color2)+
  ylim(0,NA)+
  labs(title="Number of people living with HIV,youth (15-24)",
       x="Year",
       y="Numbers")+
  theme_bw()+
  theme(plot.title = element_text(size = 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="bottom")
dev.off()

