# CODE FOR A DYNAMIC MODEL OF HIV TRANSMISSION FOR US #
# # START OF THE CODE # #

#Install additional packages (if not yet)
#Please restart R when error message occurs as the file is corrupted. 
# install.packages('EnvStats')
library('EnvStats')
library(zoo)
library(dplyr)
library(tidyr)
library(ggplot2)
# Clean global environment (clear previous data and results)
rm(list=ls())
set.seed(1)

# Specify the number of trials that alter parameter values 
trials <- 100
# Specify the number of sub-groups we are going to run
sgroup <-25
# Command specifing the number of monthly time periods in the model (432 monthly time periods from 2012 - 2020)
# Note that the generation shall be divisible by 12. Otherwise, errors would occur in calibration
generations <- 324
deathmulti<-1
initialvf<-0.25

# If error mesage "incomplete final line" occurs, open csv file with text edit and add a new empty line
# setwd("/gpfs_fs/home/panz/US model/")
#setwd("A:/Common/Personnel work/Rasnick/Rwanda Model/Rwanda Calibration Model 1.0")
#setwd("C:/Users/rrasn/Desktop/Rwanda Calibration Model 1.0")
#setwd("~/Dropbox/VCU_PhD_Year 3&4/GRA/Rwanda-Model/Full model/Data/")
# Please change the pathway to shared drive when using the VCU computer
# Calibration targets #
# Import calibration targets #
#files <- list.files(path="~/Dropbox/VCU_PhD_Year 3&4/Paper 1/R code/Calibration/", pattern = "csv")
#setwd("~/Dropbox/VCU_PhD_Year 3&4/Paper 1/R code/Calibration")
#setwd("~/Dropbox/VCU_PhD_Year 3&4/Paper 1/R code/Calibration")
setwd("C:/Users/exper/Desktop/Models/RW")
Target_all <- read.csv("Targets_00_15-54_new.csv")
Target_new <- read.csv("Targets_01.csv")

# Read parameters for HIV transmission probability
p_transmission <- read.csv(file='RW_p_transmission.csv')

# Read parameters for Condom effectiveness
p_condomEffect <- read.csv(file="RW_p_condomEffect.csv")

# Read parameters for Consistent condom use
p_condomUse <- read.csv(file="RW_p_condomUse.csv")

# Read parameters for viral suppression effectiveness
p_VSEffect <- read.csv(file="RW_p_VSEffect.csv")

# Read parameters for Number of sex acts
n_sexact <- read.csv(file="RW_n_sexact.csv")

# Read parameters for Initial population
Pop <- read.csv(file="RW_Pop.csv")
Pop_HIVPrev <- read.csv(file="RW_Pop_HIVPrev_2.csv")
Pop_CD4 <- read.csv(file="RW_Pop_CD4.csv")
Pop_comp_diag <- read.csv(file="RW_Pop_comp_diag.csv")
Pop_comp_link <- read.csv(file="RW_Pop_comp_link.csv")
# Pop_comp_ltfu <- read.csv(file="US_Pop_comp_ltfu_test.csv")
Pop_comp_artvs <- read.csv(file="RW_Pop_comp_artvs.csv")
Pop_comp_artvf <- read.csv(file="RW_Pop_comp_artvf.csv")
# Pop_2 <- read.csv(file="CM_Pop_v2.csv")

# Read parameters for population growth rate
p_popgrowth <- read.csv(file="RW_p_popgrowth.csv")

# Read parameters for natural history
p_DisProgress <- read.csv(file="RW_p_DisProgress.csv")

# Read parameters for Probability of HIV testing
p_diagnosis <- read.csv(file="RW_p_diagnosis.csv")

# Read parameters for Probability of linkage
p_link <- read.csv(file="RW_p_linkage.csv")

# Read parameters for Probability of LTFU
p_LTFU <- read.csv(file="RW_p_LTFU.csv")

# Read parameters for Probability of On ART and Suppressed
p_ARTVS <- read.csv(file="RW_p_ARTVS.csv")

# Read parameters for Probability of On ART and Not suppressed
p_ARTnotVS <- read.csv(file="RW_p_ARTnotVS.csv")

# Read parameters for Probability of Death
p_death <- read.csv(file="RW_p_death.csv")

# Read parameters for Probability of Re-engagement in care
p_reengage <- read.csv(file="RW_p_reengage.csv")

# Array to save the percent deviation
pd_top50 <- array(NA, dim=c(50,3))
pd_top50[,2] <- 999
pd_top50[,3] <- 0
colnames(pd_top50) <- c("trial","pd_total","N_within")


for ( t in 1:trials) {
# Set up empty variables to store the value for the parameter inputs 
# Note that the initial value of lambda shall set to 1 for separate rates otherwise error may occur (since we are multiplying the values)
lambda_inf <- array(1,dim=c(sgroup, generations))
lambda_t_r <- array(0,dim=c(sgroup, generations))

# Sets up variables to store the value for each compartment
# The initial value is set as null.
S <- array(NA,dim=c(sgroup, generations))
S_y <- array(NA,dim=c(sgroup, generations/12))
S_total <- array(0,dim=c(generations))
S_total_y <- array(NA,dim=c(generations/12))

S_F <- array(NA,dim=c(sgroup, generations))
S_F_y <- array(NA,dim=c(sgroup, generations/12))
S_F_total <- array(0,dim=c(generations))
S_F_total_y <- array(NA,dim=c(generations/12))

S_M <- array(NA,dim=c(sgroup, generations))
S_M_y <- array(NA,dim=c(sgroup, generations/12))
S_M_total <- array(0,dim=c(generations))
S_M_total_y <- array(NA,dim=c(generations/12))

S_F_young <- array(NA,dim=c(sgroup, generations))
S_F_young_y <- array(NA,dim=c(sgroup, generations/12))
S_F_young_total <- array(0,dim=c(generations))
S_F_young_total_y <- array(NA,dim=c(generations/12))

S_M_young <- array(NA,dim=c(sgroup, generations))
S_M_young_y <- array(NA,dim=c(sgroup, generations/12))
S_M_young_total <- array(0,dim=c(generations))
S_M_young_total_y <- array(NA,dim=c(generations/12))

D <- array(0,dim=c(sgroup, generations))
X1 <- array(NA,dim=c(sgroup, generations))
X2 <- array(NA,dim=c(sgroup, generations))
X3 <- array(NA,dim=c(sgroup, generations))
X4 <- array(NA,dim=c(sgroup, generations))
X5 <- array(NA,dim=c(sgroup, generations))
X6 <- array(NA,dim=c(sgroup, generations))
X7 <- array(NA,dim=c(sgroup, generations))
X8 <- array(NA,dim=c(sgroup, generations))
X9 <- array(NA,dim=c(sgroup, generations))
X10 <- array(NA,dim=c(sgroup, generations))
X11 <- array(NA,dim=c(sgroup, generations))
X12 <- array(NA,dim=c(sgroup, generations))
X13 <- array(NA,dim=c(sgroup, generations))
X14 <- array(NA,dim=c(sgroup, generations))
X15 <- array(NA,dim=c(sgroup, generations))
X16 <- array(NA,dim=c(sgroup, generations))
X17 <- array(NA,dim=c(sgroup, generations))
X18 <- array(NA,dim=c(sgroup, generations))
X19 <- array(NA,dim=c(sgroup, generations))
X20 <- array(NA,dim=c(sgroup, generations))
X21 <- array(NA,dim=c(sgroup, generations))
X22 <- array(NA,dim=c(sgroup, generations))
X23 <- array(NA,dim=c(sgroup, generations))
X24 <- array(NA,dim=c(sgroup, generations))


N_lr_u_f <- array(NA,dim=c(sgroup, generations))
N_hr_f <- array(NA,dim=c(sgroup, generations))
N_lr_r_f <- array(NA,dim=c(sgroup, generations))
N_u_m <- array(NA,dim=c(sgroup, generations))
N_r_m <- array(NA,dim=c(sgroup, generations))

# This command sets up variables to store the value of model parameters.
# Rate at which people enter the susceptible population compartment (population growth).
rho <- array(NA,dim=c(sgroup, generations))

# Rate of natural history HIV disease pregression.
delta_1 <- array(NA,dim=c(sgroup, generations))
delta_2 <- array(NA,dim=c(sgroup, generations))
delta_3 <- array(NA,dim=c(sgroup, generations))
delta_4 <- array(NA,dim=c(sgroup, generations))
delta_5 <- array(NA,dim=c(sgroup, generations))
delta_6 <- array(NA,dim=c(sgroup, generations))
m_delta <- array(NA,dim=c(sgroup, generations))

# Rate of getting diagnosed with HIV.
alpha_1 <- array(NA,dim=c(sgroup, generations))
alpha_2 <- array(NA,dim=c(sgroup, generations))
alpha_3 <- array(NA,dim=c(sgroup, generations))
alpha_4 <- array(NA,dim=c(sgroup, generations))

# Rate of being linked.
sigma_1 <- array(NA,dim=c(sgroup, generations))
sigma_2 <- array(NA,dim=c(sgroup, generations))
sigma_3 <- array(NA,dim=c(sgroup, generations))
sigma_4 <- array(NA,dim=c(sgroup, generations))

# Rate of getting lost from care.
gamma_1 <- array(NA,dim=c(sgroup, generations))
gamma_2 <- array(NA,dim=c(sgroup, generations))
gamma_3 <- array(NA,dim=c(sgroup, generations))
gamma_4 <- array(NA,dim=c(sgroup, generations))
gamma_5 <- array(NA,dim=c(sgroup, generations))
gamma_6 <- array(NA,dim=c(sgroup, generations))
gamma_7 <- array(NA,dim=c(sgroup, generations))
gamma_8 <- array(NA,dim=c(sgroup, generations))
gamma_9 <- array(NA,dim=c(sgroup, generations))
gamma_10 <- array(NA,dim=c(sgroup, generations))
gamma_11 <- array(NA,dim=c(sgroup, generations))
gamma_12 <- array(NA,dim=c(sgroup, generations))

# Rate of failure to maintain viral suppression.
psi_1 <- array(NA,dim=c(sgroup, generations))
psi_2 <- array(NA,dim=c(sgroup, generations))
psi_3 <- array(NA,dim=c(sgroup, generations))
psi_4 <- array(NA,dim=c(sgroup, generations))

# Rate of being on ART and suppressed.
theta_1 <- array(NA,dim=c(sgroup, generations))
theta_2 <- array(NA,dim=c(sgroup, generations))
theta_3 <- array(NA,dim=c(sgroup, generations))
theta_4 <- array(NA,dim=c(sgroup, generations))

# Rate of return to ART and being suppressed.
tau_1 <- array(NA,dim=c(sgroup, generations))
tau_2 <- array(NA,dim=c(sgroup, generations))
tau_3 <- array(NA,dim=c(sgroup, generations))
tau_4 <- array(NA,dim=c(sgroup, generations))

# Mortality rate.
mu_0 <- array(NA,dim=c(sgroup, generations))
mu_1 <- array(NA,dim=c(sgroup, generations))
mu_2 <- array(NA,dim=c(sgroup, generations))
mu_3 <- array(NA,dim=c(sgroup, generations))
mu_4 <- array(NA,dim=c(sgroup, generations))
mu_5 <- array(NA,dim=c(sgroup, generations))
mu_6 <- array(NA,dim=c(sgroup, generations))
mu_7 <- array(NA,dim=c(sgroup, generations))
mu_8 <- array(NA,dim=c(sgroup, generations))
mu_9 <- array(NA,dim=c(sgroup, generations))
mu_10 <- array(NA,dim=c(sgroup, generations))
mu_11 <- array(NA,dim=c(sgroup, generations))
mu_12 <- array(NA,dim=c(sgroup, generations))
mu_13 <- array(NA,dim=c(sgroup, generations))
mu_14 <- array(NA,dim=c(sgroup, generations))
mu_15 <- array(NA,dim=c(sgroup, generations))
mu_16 <- array(NA,dim=c(sgroup, generations))
m_mu <- array(NA,dim=c(sgroup, generations))

# Percentage of individuas consistently using a condom in sub-group r.
c_r <- array(NA,dim=c(sgroup, generations))
# Reduction in probability of HIV transmission when using a condom. 
epsilon <- array(NA,dim=c(sgroup, generations))
# Weight appleid for condom use in sub-group r.
omega_r <- array(NA,dim=c(sgroup, generations))
# Reduction in probability of HIV transmission when on ART and viral suppressed.
upsilon <- array(NA,dim=c(sgroup, generations))
# Weight applied for viral suppreessed individuals
kappa <- array(NA,dim=c(sgroup, generations))
# probability of HIV transmission per unprotected sex contact when not virally suppressed.
beta <- array(NA,dim=c(sgroup, generations))
# Average number of sexual acts per time period in sub-group r
n_r <- array(NA,dim=c(sgroup, generations))

# Total number of individuals in sub-group r
Total <- array(NA,dim=c(sgroup, generations))
initial_diag <- array(NA,dim=c(sgroup, generations))
initial_link <- array(NA,dim=c(sgroup, generations))
initial_ltfu <- array(NA,dim=c(sgroup, generations))
initial_artvs <- array(NA,dim=c(sgroup, generations))
initial_artvf <- array(NA,dim=c(sgroup, generations))
initial_cd4_unlink_1 <- array(NA,dim=c(sgroup, generations))
initial_cd4_unlink_2 <- array(NA,dim=c(sgroup, generations))
initial_cd4_unlink_3 <- array(NA,dim=c(sgroup, generations))
initial_cd4_unlink_4 <- array(NA,dim=c(sgroup, generations))
initial_cd4_link_1 <- array(NA,dim=c(sgroup, generations))
initial_cd4_link_2 <- array(NA,dim=c(sgroup, generations))
initial_cd4_link_3 <- array(NA,dim=c(sgroup, generations))
initial_cd4_link_4 <- array(NA,dim=c(sgroup, generations))
initial_cd4_artvs_1 <- array(NA,dim=c(sgroup, generations))
initial_cd4_artvs_2 <- array(NA,dim=c(sgroup, generations))
initial_cd4_artvs_3 <- array(NA,dim=c(sgroup, generations))
initial_cd4_artvs_4 <- array(NA,dim=c(sgroup, generations))
initial_cd4_artvf_1 <- array(NA,dim=c(sgroup, generations))
initial_cd4_artvf_2 <- array(NA,dim=c(sgroup, generations))
initial_cd4_artvf_3 <- array(NA,dim=c(sgroup, generations))
initial_cd4_artvf_4 <- array(NA,dim=c(sgroup, generations))

# Set up variables to store model outputs
# Set up variables to store the number of individuals infected.
N <- array(0,dim=c(sgroup, generations))
N_y <- array(0,dim=c(sgroup, generations/12))
N_total <- array(0,dim=c(generations))
N_total_y <- array(0,dim=c(generations/12))
N_F <- array(0,dim=c(sgroup, generations))
N_F_y <- array(0,dim=c(sgroup, generations/12))
N_F_total <- array(0,dim=c(generations))
N_F_total_y <- array(0,dim=c(generations/12))
N_M <- array(0,dim=c(sgroup, generations))
N_M_y <- array(0,dim=c(sgroup, generations/12))
N_M_total <- array(0,dim=c(generations))
N_M_total_y <- array(0,dim=c(generations/12))
N_young <- array(0,dim=c(sgroup, generations))
N_young_y <- array(0,dim=c(sgroup, generations/12))
N_young_total <- array(0,dim=c(generations))
N_young_total_y <- array(0,dim=c(generations/12))
N_F_young_total <- array(0,dim=c(generations))
N_F_young_total_y <- array(0,dim=c(generations/12))
N_M_young_total <- array(0,dim=c(generations))
N_M_young_total_y <- array(0,dim=c(generations/12))
N_adult <- array(0,dim=c(sgroup, generations))
N_adult_y <- array(0,dim=c(sgroup, generations/12))
N_adult_total <- array(0,dim=c(generations))
N_adult_total_y <- array(0,dim=c(generations/12))

S_F_total<- array(0,dim=c(generations))
S_M_total<- array(0,dim=c(generations))
S_young_total<- array(0,dim=c(generations))
S_F_young_total<- array(0,dim=c(generations))
S_M_young_total<- array(0,dim=c(generations))

#This command sets up variables to store the number of individuals virally not suppresed.
Z <- array(0,dim=c(sgroup, generations))
Z_total <- array(0,dim=c(generations))
Z_total_y <- array(0,dim=c(generations/12))
#This command sets up variables to store the number of individuals virally suppressed.
V <- array(0,dim=c(sgroup, generations))
V_total <- array(0,dim=c(generations))
V_total_y <- array(0,dim=c(generations/12))

# This command sets up variables to store the percent of individuals virally not suppressed.
p_r <- array(0,dim=c(sgroup, generations))

# This command sets up variables to store the percent of individuals virally suppressed.
p_r_v <- array(0,dim=c(sgroup, generations))

# This command sets up variables to store HIV treatment engagement status in the cohort. 
Undiagnosed <- array(NA,dim=c(sgroup, generations))
Undiagnosed_y <- array(NA,dim=c(sgroup, generations/12))
Undiagnosed_total <- array(0,dim=c(generations))
Undiagnosed_total_y <- array(NA,dim=c(generations/12))

Diagnosed <- array(NA,dim=c(sgroup, generations))
Diagnosed_y <- array(NA,dim=c(sgroup, generations/12))
Diagnosed_total <- array(0,dim=c(generations))
Diagnosed_total_y <- array(NA,dim=c(generations/12))
Diagnosed_F_total <- array(0,dim=c(generations))
Diagnosed_F_total_y <- array(NA,dim=c(generations/12))
Diagnosed_M_total <- array(0,dim=c(generations))
Diagnosed_M_total_y <- array(NA,dim=c(generations/12))
Diagnosed_F_young_total <- array(0,dim=c(generations))
Diagnosed_F_young_total_y <- array(NA,dim=c(generations/12))
Diagnosed_M_young_total <- array(0,dim=c(generations))
Diagnosed_M_young_total_y <- array(NA,dim=c(generations/12))
Diagnosed_total_py <- array(NA,dim=c(generations/12))

Rate_diagnosed <- array(NA,dim=c(sgroup, generations))
Rate_diagnosed_y <- array(NA,dim=c(sgroup, generations/12))
Rate_diagnosed_total <- array(0,dim=c(generations))
Rate_diagnosed_total_y <- array(NA,dim=c(generations/12))
Rate_diagnosed_F_total <- array(0,dim=c(generations))
Rate_diagnosed_F_total_y <- array(NA,dim=c(generations/12))
Rate_diagnosed_M_total <- array(0,dim=c(generations))
Rate_diagnosed_M_total_y <- array(NA,dim=c(generations/12))
Rate_diagnosed_F_young_total <- array(0,dim=c(generations))
Rate_diagnosed_F_young_total_y <- array(NA,dim=c(generations/12))
Rate_diagnosed_M_young_total <- array(0,dim=c(generations))
Rate_diagnosed_M_young_total_y <- array(NA,dim=c(generations/12))

New_Diagnosed <- array(NA,dim=c(sgroup, generations))
New_Diagnosed_y <- array(NA,dim=c(sgroup, generations/12))
New_Diagnosed_total <- array(0,dim=c(generations))
New_Diagnosed_total_y <- array(NA,dim=c(generations/12))

New_Diagnosed_MSM_total <- array(0,dim=c(generations))
New_Diagnosed_MSM_total_y <- array(NA,dim=c(generations/12))
New_Diagnosed_MSMIDU_total <- array(0,dim=c(generations))
New_Diagnosed_MSMIDU_total_y <- array(NA,dim=c(generations/12))
New_Diagnosed_IDUM_total <- array(0,dim=c(generations))
New_Diagnosed_IDUM_total_y <- array(NA,dim=c(generations/12))
New_Diagnosed_IDUF_total <- array(0,dim=c(generations))
New_Diagnosed_IDUF_total_y <- array(NA,dim=c(generations/12))
New_Diagnosed_HETM_total <- array(0,dim=c(generations))
New_Diagnosed_HETM_total_y <- array(NA,dim=c(generations/12))
New_Diagnosed_HETF_total <- array(0,dim=c(generations))
New_Diagnosed_HETF_total_y <- array(NA,dim=c(generations/12))

Rate_new_Diagnosed <- array(NA,dim=c(sgroup, generations))
Rate_new_Diagnosed_y <- array(NA,dim=c(sgroup, generations/12))
Rate_new_Diagnosed_total <- array(0,dim=c(generations))
Rate_new_Diagnosed_total_y <- array(NA,dim=c(generations/12))

Linked <- array(NA,dim=c(sgroup, generations))
Linked_y <- array(NA,dim=c(sgroup, generations/12))
Linked_total <- array(0,dim=c(generations))
Linked_total_y <- array(NA,dim=c(generations/12))

Lost <- array(NA,dim=c(sgroup, generations))
Lost_y <- array(NA,dim=c(sgroup, generations/12))
Lost_total <- array(0,dim=c(generations))
Lost_total_y <- array(NA,dim=c(generations/12))

Suppressed <- array(NA,dim=c(sgroup, generations))
Suppressed_y <- array(NA,dim=c(sgroup, generations/12))
Suppressed_total <- array(0,dim=c(generations))
Suppressed_total_y <- array(NA,dim=c(generations/12))
Suppressed_total_F <- array(0,dim=c(generations))
Suppressed_total_F_y <- array(NA,dim=c(generations/12))
Suppressed_total_M <- array(0,dim=c(generations))
Suppressed_total_M_y <- array(NA,dim=c(generations/12))
Suppressed_total_F_young <- array(0,dim=c(generations))
Suppressed_total_F_young_y <- array(NA,dim=c(generations/12))
Suppressed_total_M_young <- array(0,dim=c(generations))
Suppressed_total_M_young_y <- array(NA,dim=c(generations/12))
Target_ART <- array(NA,dim=c(generations/12))
P_Target_ART <- array(NA,dim=c(generations/12))

P_suppressed <- array(NA,dim=c(sgroup, generations))
P_suppressed_y <- array(NA,dim=c(sgroup, generations/12))
P_suppressed_total <- array(0,dim=c(generations))
P_suppressed_total_y <- array(NA,dim=c(generations/12))
P_suppressed_F_total <- array(0,dim=c(generations))
P_suppressed_F_total_y <- array(NA,dim=c(generations/12))
P_suppressed_M_total <- array(0,dim=c(generations))
P_suppressed_M_total_y <- array(NA,dim=c(generations/12))
P_suppressed_F_young_total <- array(0,dim=c(generations))
P_suppressed_F_young_total_y <- array(NA,dim=c(generations/12))
P_suppressed_M_young_total <- array(0,dim=c(generations))
P_suppressed_M_young_total_y <- array(NA,dim=c(generations/12))

P_suppressed_art <- array(NA,dim=c(sgroup, generations))
P_suppressed_art_y <- array(NA,dim=c(sgroup, generations/12))
P_suppressed_art_total <- array(0,dim=c(generations))
P_suppressed_art_total_y <- array(NA,dim=c(generations/12))
P_suppressed_art_F_total <- array(0,dim=c(generations))
P_suppressed_art_F_total_y <- array(NA,dim=c(generations/12))
P_suppressed_art_M_total <- array(0,dim=c(generations))
P_suppressed_art_M_total_y <- array(NA,dim=c(generations/12))
P_suppressed_art_F_young_total <- array(0,dim=c(generations))
P_suppressed_art_F_young_total_y <- array(NA,dim=c(generations/12))
P_suppressed_art_M_young_total <- array(0,dim=c(generations))
P_suppressed_art_M_young_total_y <- array(NA,dim=c(generations/12))


Suppressed_onART_y <- array(NA,dim=c(sgroup, generations/12))
Suppressed_onART_total_y <- array(NA,dim=c(generations/12))
Target_VSonART <- array(NA,dim=c(generations/12))
Target_VS_art_F<- array(NA,dim=c(generations/12))
Target_VS_art_M<- array(NA,dim=c(generations/12))
Target_VS_art_F_young<- array(NA,dim=c(generations/12))
Target_VS_art_M_young<- array(NA,dim=c(generations/12))

Target_VS <- array(NA,dim=c(generations/12))
Target_VS_F<- array(NA,dim=c(generations/12))
Target_VS_M<- array(NA,dim=c(generations/12))
Target_VS_F_young<- array(NA,dim=c(generations/12))
Target_VS_M_young<- array(NA,dim=c(generations/12))

Not_Suppressed <- array(NA,dim=c(sgroup, generations))
Not_Suppressed_y <- array(NA,dim=c(sgroup, generations/12))
Not_Suppressed_total <- array(0,dim=c(generations))
Not_Suppressed_total_y <- array(NA,dim=c(generations/12))

On_ART <- array(NA,dim=c(sgroup, generations))
On_ART_y <- array(NA,dim=c(sgroup, generations/12))
On_ART_p <- array(0,dim=c(generations))
On_ART_total <- array(0,dim=c(generations))
On_ART_total_y <- array(NA,dim=c(generations/12))
On_ART_total_F <- array(0,dim=c(generations))
On_ART_total_F_y <- array(NA,dim=c(generations/12))
On_ART_total_M <- array(0,dim=c(generations))
On_ART_total_M_y <- array(NA,dim=c(generations/12))
On_ART_total_F_young <- array(0,dim=c(generations))
On_ART_total_F_young_y <- array(NA,dim=c(generations/12))
On_ART_total_M_young <- array(0,dim=c(generations))
On_ART_total_M_young_y <- array(NA,dim=c(generations/12))
On_ART_grand_total<- array(0,dim=c(generations))
On_ART_grand_total_y <- array(NA,dim=c(generations/12))

On_ART_p_y <- array(NA,dim=c(generations/12))
On_ART_total_ny <- array(NA,dim=c(generations/12))
Target_onART_N <- array(NA,dim=c(generations/12))
Target_onART_P <- array(NA,dim=c(generations/12))

On_ART_p_F <- array(0,dim=c(generations))
On_ART_p_F_y <- array(NA,dim=c(generations/12))
On_ART_p_M <- array(0,dim=c(generations))
On_ART_p_M_y <- array(NA,dim=c(generations/12))

On_ART_p_F_young <- array(0,dim=c(generations))
On_ART_p_F_young_y <- array(NA,dim=c(generations/12))
On_ART_p_M_young <- array(0,dim=c(generations))
On_ART_p_M_young_y <- array(NA,dim=c(generations/12))

# This command sets up variables to store HIV incidence
Incidence <- array(NA,dim=c(sgroup, generations))
Incidence_y <- array(NA,dim=c(sgroup, generations/12))
Incidence_total <- array(0,dim=c(generations))
Incidence_total_y <- array(NA,dim=c(generations/12))
Target_incidence <- array(NA,dim=c(generations/12))
Incidence_F <- array(NA,dim=c(sgroup, generations))
Incidence_F_y <- array(NA,dim=c(sgroup, generations/12))
Incidence_F_total <- array(0,dim=c(generations))
Incidence_F_total_y <- array(NA,dim=c(generations/12))
Target_incidence_F <- array(NA,dim=c(generations/12))
Incidence_M <- array(NA,dim=c(sgroup, generations))
Incidence_M_y <- array(NA,dim=c(sgroup, generations/12))
Incidence_M_total <- array(0,dim=c(generations))
Incidence_M_total_y <- array(NA,dim=c(generations/12))
Target_incidence_M <- array(NA,dim=c(generations/12))
Incidence_F_young <- array(NA,dim=c(sgroup, generations))
Incidence_F_young_y <- array(NA,dim=c(sgroup, generations/12))
Incidence_F_young_total <- array(0,dim=c(generations))
Incidence_F_young_total_y <- array(NA,dim=c(generations/12))
Target_incidence_F_young  <- array(NA,dim=c(generations/12))
Incidence_M_young  <- array(NA,dim=c(sgroup, generations))
Incidence_M_young_y <- array(NA,dim=c(sgroup, generations/12))
Incidence_M_young_total <- array(0,dim=c(generations))
Incidence_M_young_total_y <- array(NA,dim=c(generations/12))
Target_incidence_M_young  <- array(NA,dim=c(generations/12))

# This command sets up variables to store HIV prevalence.
Prevalence <- array(NA,dim=c(sgroup, generations))
Prevalence_y <- array(NA,dim=c(sgroup, generations/12))
Prevalence_total <- array(0,dim=c(generations))
Prevalence_total_y <- array(NA,dim=c(generations/12))
Target_prev <- array(NA,dim=c(generations/12))
Target_prev_F <- array(NA,dim=c(generations/12))
Target_prev_M <- array(NA,dim=c(generations/12))
Target_prev_young <- array(NA,dim=c(generations/12))
Target_prev_F_young <- array(NA,dim=c(generations/12))
Target_prev_M_young <- array(NA,dim=c(generations/12))

DHS_Prevalence <- array(NA,dim=c(sgroup, generations))
DHS_Prevalence_y <- array(NA,dim=c(sgroup, generations/12))
DHS_Prevalence_total <- array(0,dim=c(generations))
DHS_Prevalence_total_y <- array(NA,dim=c(generations/12))
DHS_Target_prev <- array(NA,dim=c(generations/12))

DHS_Prevalence_F <- array(NA,dim=c(sgroup, generations))
DHS_Prevalence_y_F<- array(NA,dim=c(sgroup, generations/12))
DHS_Prevalence_F_total <- array(0,dim=c(generations))
DHS_Prevalence_F_total_y <- array(NA,dim=c(generations/12))
DHS_Target_prev_F <- array(NA,dim=c(generations/12))

DHS_Prevalence_M <- array(NA,dim=c(sgroup, generations))
DHS_Prevalence_y_M <- array(NA,dim=c(sgroup, generations/12))
DHS_Prevalence_M_total <- array(0,dim=c(generations))
DHS_Prevalence_M_total_y <- array(NA,dim=c(generations/12))
DHS_Target_prev_M <- array(NA,dim=c(generations/12))

DHS_Prevalence_F_young <- array(NA,dim=c(sgroup, generations))
DHS_Prevalence_y_F_young <- array(NA,dim=c(sgroup, generations/12))
DHS_Prevalence_F_total_young <- array(0,dim=c(generations))
DHS_Prevalence_F_total_y_young <- array(NA,dim=c(generations/12))
DHS_Target_prev_F_young <- array(NA,dim=c(generations/12))

DHS_Prevalence_M_young <- array(NA,dim=c(sgroup, generations))
DHS_Prevalence_y_M_young <- array(NA,dim=c(sgroup, generations/12))
DHS_Prevalence_M_total_young <- array(0,dim=c(generations))
DHS_Prevalence_M_total_y_young <- array(NA,dim=c(generations/12))
DHS_Target_prev_M_young <- array(NA,dim=c(generations/12))

Prevalence_F_total <- array(0,dim=c(generations))
Prevalence_F_total_y <- array(NA,dim=c(generations/12))
Target_prev_F <- array(NA,dim=c(generations/12))

Prevalence_M_total <- array(0,dim=c(generations))
Prevalence_M_total_y <- array(NA,dim=c(generations/12))
Target_prev_M <- array(NA,dim=c(generations/12))

Prevalence_young_total <- array(0,dim=c(generations))
Prevalence_young_total_y <- array(NA,dim=c(generations/12))
Target_prev_young <- array(NA,dim=c(generations/12))

Prevalence_F_young_total <- array(0,dim=c(generations))
Prevalence_F_young_total_y <- array(NA,dim=c(generations/12))
Target_prev_F_young <- array(NA,dim=c(generations/12))

Prevalence_M_young_total <- array(0,dim=c(generations))
Prevalence_M_young_total_y <- array(NA,dim=c(generations/12))
Target_prev_M_young <- array(NA,dim=c(generations/12))

# This command sets up variables to store the number of new infection.
NI <- array(NA,dim=c(sgroup, generations))
NI_y <- array(NA,dim=c(sgroup, generations/12))
NI_total <- array(0,dim=c(generations))
NI_total_y <- array(NA,dim=c(generations/12))
NI_F <- array(NA,dim=c(sgroup, generations))
NI_F_y <- array(NA,dim=c(sgroup, generations/12))
NI_F_total <- array(0,dim=c(generations))
NI_F_total_y <- array(NA,dim=c(generations/12))
NI_M <- array(NA,dim=c(sgroup, generations))
NI_M_y <- array(NA,dim=c(sgroup, generations/12))
NI_M_total <- array(0,dim=c(generations))
NI_M_total_y <- array(NA,dim=c(generations/12))
NI_young <- array(NA,dim=c(sgroup, generations))
NI_young_y <- array(NA,dim=c(sgroup, generations/12))
NI_young_total <- array(0,dim=c(generations))
NI_young_total_y <- array(NA,dim=c(generations/12))
NI_M_young <- array(NA,dim=c(sgroup, generations))
NI_M_young_y <- array(NA,dim=c(sgroup, generations/12))
NI_M_young_total <- array(0,dim=c(generations))
NI_M_young_total_y <- array(NA,dim=c(generations/12))
NI_F_young <- array(NA,dim=c(sgroup, generations))
NI_F_young_y <- array(NA,dim=c(sgroup, generations/12))
NI_F_young_total <- array(0,dim=c(generations))
NI_F_young_total_y <- array(NA,dim=c(generations/12))
NI_adult <- array(NA,dim=c(sgroup, generations))
NI_adult_y <- array(NA,dim=c(sgroup, generations/12))
NI_adult_total <- array(0,dim=c(generations))
NI_adult_total_y <- array(NA,dim=c(generations/12))

Target_N <- array(NA,dim=c(generations/12))
Target_N_F <- array(NA,dim=c(generations/12))
Target_N_M <- array(NA,dim=c(generations/12))
Target_N_young <- array(NA,dim=c(generations/12))
Target_N_adult <- array(NA,dim=c(generations/12))
Target_NI <- array(NA,dim=c(generations/12))
Target_NI_F <- array(NA,dim=c(generations/12))
Target_NI_M <- array(NA,dim=c(generations/12))
Target_NI_young <- array(NA,dim=c(generations/12))
Target_NI_adult <- array(NA,dim=c(generations/12))

# This command sets up variables to store total unprotected sex for heterosexual female.
N_het_f_u<- array(0,dim=c(sgroup, generations))
N_het_f_r<- array(0,dim=c(sgroup, generations))
# This command sets up variables to store total unprotected sex for heterosexual male.
N_het_m_u<- array(0,dim=c(sgroup, generations))
N_het_m_r<- array(0,dim=c(sgroup, generations))

Target_N_diag <- array(NA,dim=c(generations/12))
Target_R_diag <- array(NA,dim=c(generations/12))
Target_N_HIV <- array(NA,dim=c(generations/12))
Target_R_HIV_F <- array(NA,dim=c(generations/12))
Target_R_HIV_M <- array(NA,dim=c(generations/12))
Target_R_HIV_F_young <- array(NA,dim=c(generations/12))
Target_R_HIV_M_young <- array(NA,dim=c(generations/12))
Target_R_HIV <- array(NA,dim=c(generations/12))
Target_VS_art <- array(NA,dim=c(generations/12))

#Initialize the parameter inputs by random draw
# Initial population
# Susceptible population
  # Total population for each sub-population
  TotalUrbanWomen=sum(Pop[1:10,"Total"])
  TotalHighRisk = runif(1,min=11000, max=45000)
  TotalLowRiskUrban=TotalUrbanWomen-TotalHighRisk
  # High risk women
  Total[6,1]=TotalHighRisk*(157+470)/1338 #15-24
  Total[7,1]=TotalHighRisk*(342+288/2)/1338 #25-34
  Total[8,1]=TotalHighRisk*(288/2+81/5)/1338 #35-44
  Total[9,1]=TotalHighRisk*(81/5+81/5)/1338 #45-54
  Total[10,1]=TotalHighRisk*(81/5+81/5)/1338 #55-64
  # low risk urban women
  women=c(535205,460559,327619,256033,201099,209695,154203,124472,89338,67466)
  totalwomen=sum(women)
  Total[1,1]=TotalUrbanWomen*(sum(535205,460559)/totalwomen)-Total[6,1] #15-24
  Total[2,1]=TotalUrbanWomen*(sum(327619,256033)/totalwomen)-Total[7,1] #25-34
  Total[3,1]=TotalUrbanWomen*(sum(201099,209695)/totalwomen)-Total[8,1] #35-44
  Total[4,1]=TotalUrbanWomen*(sum(154203,124472)/totalwomen)-Total[9,1]#45-54
  Total[5,1]=TotalUrbanWomen*(sum(89338,67466)/totalwomen)-Total[10,1] #55-64
  
  for (s in 11:sgroup) {
    Total[s,1]=Pop[s,"Total"]
  }
  # HIV prevalence
  for (l in 1:25) {
    a=Pop_HIVPrev[l,"Alpha"]
    b=Pop_HIVPrev[l,"Beta"]
    Prevalence[l,1]=rbeta(1,a,b)
  }
  # Susceptible population
  for (s in 1:sgroup) {S[s,1]=Total[s,1]*(1-Prevalence[s,1])}

# By infectious compartment
# Diagnosis
  # By age
  multi_diag=rtri(1,0.15,0.8,0.3)
  
  a=Pop_comp_diag[1,"Alpha"]
  b=Pop_comp_diag[1,"Beta"]
  #mode=Pop_comp_diag[1,"Mode"]
  initial_diag[1,1]=rbeta(1,a,b)*multi_diag

  # Apply to all sub-population
  for (s in 1:sgroup) {initial_diag[s,1]=initial_diag[1,1]}

# Linkage
  # By age
  a=Pop_comp_link[1,"Alpha"]
  b=Pop_comp_link[1,"Beta"]
  initial_link[1,1]=rbeta(1,a,b)

  # Apply to all sub-population
  for (s in 1:sgroup) {initial_link[s,1]=initial_link[1,1]}

# LTFU (No one is LTFU)

# On ART and suppressed
  # By age
  a=Pop_comp_artvs[1,"Alpha"]
  b=Pop_comp_artvs[1,"Beta"]
  initial_artvs[1,1]=rbeta(1,a,b)

  # Apply to all sub-population
  for (s in 1:sgroup) {initial_artvs[s,1]=initial_artvs[1,1]}

# On ART and not suppressed
  # By age
  a=Pop_comp_artvf[1,"Alpha"]
  b=Pop_comp_artvf[1,"Beta"]
  initial_artvf[1,1]=rbeta(1,a,b)

  # Apply to all sub-population
  for (s in 1:sgroup) {initial_artvf[s,1]=initial_artvf[1,1]}

# Distribution by CD4
  # Not linked, CD4>500
  #a=Pop_CD4[1,"Alpha"]
  #b=Pop_CD4[1,"Beta"]
  initial_cd4_unlink_1[1,1]=Pop_CD4[1,"Mean"]
  # Not linked, CD4 350-500
  #a=Pop_CD4[2,"Alpha"]
  #b=Pop_CD4[2,"Beta"]
  initial_cd4_unlink_2[1,1]=Pop_CD4[2,"Mean"]
  # Not linked, CD4 200-350
  #a=Pop_CD4[3,"Alpha"]
  #b=Pop_CD4[3,"Beta"]
  initial_cd4_unlink_3[1,1]=Pop_CD4[3,"Mean"]
  # Not linked, CD4 <200
  initial_cd4_unlink_4[1,1]=1-initial_cd4_unlink_1[1,1]-initial_cd4_unlink_2[1,1]-initial_cd4_unlink_3[1,1]
  # Linked, CD4>500
  #a=Pop_CD4[4,"Alpha"]
  #b=Pop_CD4[4,"Beta"]
  initial_cd4_link_1[1,1]=Pop_CD4[4,"Mean"]
  # Linked, CD4 350-500
  #a=Pop_CD4[5,"Alpha"]
  #b=Pop_CD4[5,"Beta"]
  initial_cd4_link_2[1,1]=Pop_CD4[5,"Mean"]
  # Linked, CD4 200-350
  a=Pop_CD4[6,"Alpha"]
  b=Pop_CD4[6,"Beta"]
  initial_cd4_link_3[1,1]=Pop_CD4[6,"Mean"]
  # Linked, CD4 <200
  initial_cd4_link_4[1,1]=1-initial_cd4_link_1[1,1]-initial_cd4_link_2[1,1]-initial_cd4_link_3[1,1]
  # ARTVS, CD4>500
  #a=Pop_CD4[7,"Alpha"]
  #b=Pop_CD4[7,"Beta"]
  initial_cd4_artvs_1[1,1]=Pop_CD4[7,"Mean"]
  # ARTVS, CD4 350-500
  #a=Pop_CD4[8,"Alpha"]
  #b=Pop_CD4[8,"Beta"]
  initial_cd4_artvs_2[1,1]=Pop_CD4[8,"Mean"]
  # ARTVS, CD4 200-350
  #a=Pop_CD4[9,"Alpha"]
  #b=Pop_CD4[9,"Beta"]
  initial_cd4_artvs_3[1,1]=Pop_CD4[9,"Mean"]
  # ARTVS, CD4 <200
  initial_cd4_artvs_4[1,1]=1-initial_cd4_artvs_1[1,1]-initial_cd4_artvs_2[1,1]-initial_cd4_artvs_3[1,1]
  # ARTVF, CD4>500
  # a=Pop_CD4[7,"Alpha"]
  #  b=Pop_CD4[7,"Beta"]
  initial_cd4_artvf_1[1,1]=Pop_CD4[10,"Mean"]
  # ARTVF, CD4 350-500
  # a=Pop_CD4[8,"Alpha"]
  #  b=Pop_CD4[8,"Beta"]
  initial_cd4_artvf_2[1,1]=Pop_CD4[11,"Mean"]
  # ARTVF, CD4 200-350
  # a=Pop_CD4[9,"Alpha"]
  #  b=Pop_CD4[9,"Beta"]
  initial_cd4_artvf_3[1,1]=Pop_CD4[12,"Mean"]
  # ARTVF, CD4 <200
  initial_cd4_artvf_4[1,1]=1-initial_cd4_artvf_1[1,1]-initial_cd4_artvf_2[1,1]-initial_cd4_artvf_3[1,1]
  
  # Apply to all sub-population
  for (s in 1:sgroup){
    initial_cd4_unlink_1[s,1]=initial_cd4_unlink_1[1,1]
    initial_cd4_unlink_2[s,1]=initial_cd4_unlink_2[1,1]
    initial_cd4_unlink_3[s,1]=initial_cd4_unlink_3[1,1]
    initial_cd4_unlink_4[s,1]=initial_cd4_unlink_4[1,1]
    initial_cd4_link_1[s,1]=initial_cd4_link_1[1,1]
    initial_cd4_link_2[s,1]=initial_cd4_link_2[1,1]
    initial_cd4_link_3[s,1]=initial_cd4_link_3[1,1]
    initial_cd4_link_4[s,1]=initial_cd4_link_4[1,1]
    initial_cd4_artvs_1[s,1]=initial_cd4_artvs_1[1,1]
    initial_cd4_artvs_2[s,1]=initial_cd4_artvs_2[1,1]
    initial_cd4_artvs_3[s,1]=initial_cd4_artvs_3[1,1]
    initial_cd4_artvs_4[s,1]=initial_cd4_artvs_4[1,1]
    initial_cd4_artvf_1[s,1]=initial_cd4_artvf_1[1,1]
    initial_cd4_artvf_2[s,1]=initial_cd4_artvf_2[1,1]
    initial_cd4_artvf_3[s,1]=initial_cd4_artvf_3[1,1]
    initial_cd4_artvf_4[s,1]=initial_cd4_artvf_4[1,1]
  }


# Infectious compartment
  for (s in 1:sgroup){
    X1[s,1]=Total[s,1]*Prevalence[s,1]*(1-initial_diag[s,1])*initial_cd4_unlink_1[s,1]
    X2[s,1]=Total[s,1]*Prevalence[s,1]*(1-initial_diag[s,1])*initial_cd4_unlink_2[s,1]
    X3[s,1]=Total[s,1]*Prevalence[s,1]*(1-initial_diag[s,1])*initial_cd4_unlink_3[s,1]
    X4[s,1]=Total[s,1]*Prevalence[s,1]*(1-initial_diag[s,1])*initial_cd4_unlink_4[s,1]
    X5[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*(1-initial_link[s,1])*initial_cd4_link_1[s,1]
    X6[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*(1-initial_link[s,1])*initial_cd4_link_2[s,1]
    X7[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*(1-initial_link[s,1])*initial_cd4_link_3[s,1]
    X8[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*(1-initial_link[s,1])*initial_cd4_link_4[s,1]
    X9[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*(1-initial_artvs[s,1])*initial_cd4_link_1[s,1]
    X10[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*(1-initial_artvs[s,1])*initial_cd4_link_2[s,1]
    X11[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*(1-initial_artvs[s,1])*initial_cd4_link_3[s,1]
    X12[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*(1-initial_artvs[s,1])*initial_cd4_link_4[s,1]
    X13[s,1]=0
    X14[s,1]=0
    X15[s,1]=0
    X16[s,1]=0
    X17[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*initial_artvs[s,1]*(1-initial_artvf[s,1])*initial_cd4_artvs_1[s,1]
    X18[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*initial_artvs[s,1]*(1-initial_artvf[s,1])*initial_cd4_artvs_2[s,1]
    X19[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*initial_artvs[s,1]*(1-initial_artvf[s,1])*initial_cd4_artvs_3[s,1]
    X20[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*initial_artvs[s,1]*(1-initial_artvf[s,1])*initial_cd4_artvs_4[s,1]
    X21[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*initial_artvs[s,1]*initial_artvf[s,1]*initial_cd4_artvf_1[s,1]
    X22[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*initial_artvs[s,1]*initial_artvf[s,1]*initial_cd4_artvf_2[s,1]
    X23[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*initial_artvs[s,1]*initial_artvf[s,1]*initial_cd4_artvf_3[s,1]
    X24[s,1]=Total[s,1]*Prevalence[s,1]*initial_diag[s,1]*initial_link[s,1]*initial_artvs[s,1]*initial_artvf[s,1]*initial_cd4_artvf_4[s,1]
  }


# Probability of HIV transmission (heterosexual, converted to monthly probability) 
  # FSW
  min=p_transmission[1,"Min"]
  max=p_transmission[1,"Max"]
  mode=p_transmission[1,"Mode"]
  # random draw
  beta[6,1]=rtri(1,min,max,mode)
  # Heterosexual female
  min=p_transmission[2,"Min"]
  max=p_transmission[2,"Max"]
  mode=p_transmission[2,"Mode"]
  # random draw
  beta[1,1]=rtri(1,min,max,mode)
  # Heterosexual male
  min=p_transmission[3,"Min"]
  max=p_transmission[3,"Max"]
  mode=p_transmission[3,"Mode"]
  # random draw
  beta[21,1]=rtri(1,min,max,mode)
  # Apply the random draw value to all sub-population
  for (s in 1:sgroup) {
    if (s>=1&s<=5||s>=11&s<=15) {beta[s,1]=beta[1,1]}
    if (s>=6&s<=10){beta[s,1]=beta[6,1]}
    if (s>=16&s<=25){beta[s,1]=beta[21,1]}
    for (g in 2:generations) {beta[s,g]=beta[s,1]}}

# Condom effectiveness
  # HET
  min=p_condomEffect[1,"Min"]
  max=p_condomEffect[1,"Max"]
  mode=p_condomEffect[1,"Mode"]
  # random draw
  epsilon[1,1]=rtri(1,min,max,mode)
  # Apply the random draw value to all sub-population
  for (s in 1:sgroup) {epsilon[s,1]=epsilon[1,1]
  for (g in 2:generations) {
    epsilon[s,g]=epsilon[s,1]}}

# Proportion of consistent condom use
  # Low-risk urban women
  for (l in 1:5){
    a=p_condomUse[l,"Alpha"]
    b=p_condomUse[l,"Beta"]
    c_r[l,1]=rbeta(1,a,b)

    a=p_condomUse[l+5,"Alpha"]
    b=p_condomUse[l+5,"Beta"]
    c_r[l,73]=rbeta(1,a,b)

    a=p_condomUse[l+10,"Alpha"]
    b=p_condomUse[l+10,"Beta"]
    c_r[l,133]=rbeta(1,a,b)

    # Low-risk rural women
    a=p_condomUse[l+15,"Alpha"]
    b=p_condomUse[l+15,"Beta"]
    c_r[l+10,1]=rbeta(1,a,b)
    a=p_condomUse[l+20,"Alpha"]
    b=p_condomUse[l+20,"Beta"]
    c_r[l+10,73]=rbeta(1,a,b)
    a=p_condomUse[l+25,"Alpha"]
    b=p_condomUse[l+25,"Beta"]
    c_r[l+10,133]=rbeta(1,a,b)
    
    # Urban men
    a=p_condomUse[l+30,"Alpha"]
    b=p_condomUse[l+30,"Beta"]
    c_r[l+15,1]=rbeta(1,a,b)
    a=p_condomUse[l+35,"Alpha"]
    b=p_condomUse[l+35,"Beta"]
    c_r[l+15,73]=rbeta(1,a,b)
    a=p_condomUse[l+40,"Alpha"]
    b=p_condomUse[l+40,"Beta"]
    c_r[l+15,133]=rbeta(1,a,b)

    # Rural men
    a=p_condomUse[l+45,"Alpha"]
    b=p_condomUse[l+45,"Beta"]
    c_r[l+20,1]=rbeta(1,a,b)
    
    a=p_condomUse[l+50,"Alpha"]
    b=p_condomUse[l+50,"Beta"]
    c_r[l+20,73]=rbeta(1,a,b)
    
    a=p_condomUse[l+55,"Alpha"]
    b=p_condomUse[l+55,"Beta"]
    c_r[l+20,133]=rbeta(1,a,b)

    
    # High-risk women
    a=p_condomUse[l+60,"Alpha"]
    b=p_condomUse[l+60,"Beta"]
    c_r[l+5,1]=rbeta(1,a,b)
    c_r[l+5,73]=c_r[l+5,1]
    c_r[l+5,133]=c_r[l+5,1]
  }
  
  #Apply to all generations
  for (s in 1:sgroup) {
    for (g in 2:generations) {
      if (g>=1&g<=72) {
        c_r[s,g]=c_r[s,1]
        #c_r_ipv[s,g]=c_r_ipv[s,1]
      }
      if (g>=73&g<=132) {
        c_r[s,g]=c_r[s,73]
        #c_r_ipv[s,g]=c_r_ipv[s,49]
      }
      if (g>=133&g<=generations) {
        c_r[s,g]=c_r[s,133]
      }}}

# Effectiveness of viral suppression
  # Sexual
  min=p_VSEffect[1,"Min"]
  max=p_VSEffect[1,"Max"]
  mode=p_VSEffect[1,"Mode"]
  # random draw
  upsilon[1,1]=rtri(1,min,max,mode)
  
  # Apply the random draw value to all sub-population
  for (s in 1:sgroup) {
    upsilon[s,1]=upsilon[1,1]
    for (g in 2:generations) {
      upsilon[s,g]=upsilon[s,1]
    }}

# Number of sex acts
  # Low risk
  min=n_sexact[1,"Min"]
  max=n_sexact[1,"Max"]
  mode=n_sexact[1,"Mode"] 
  n_r[1,1]=rtri(1,min,max,mode)
  n_r[2,1]=n_r[1,1]
  # High risk
  min=n_sexact[2,"Min"]
  max=n_sexact[2,"Max"]
  mode=n_sexact[2,"Mode"] 
  n_r[6,1]=rtri(1,min,max,mode)*runif(n=1,min=1,max=5)

  #Youth multiplier
  min=0.2
  max=1.0
  mode=0.5
  youth_multi_f=rtri(1,min,max,mode)
  youth_multi_m=rtri(1,min,max,mode)
  
  for (s in 1:sgroup) {
    # Assign value to all other sub-populations 
      if (s==1||s==11){n_r[s,1]=n_r[2,1]*youth_multi_f}
      if (s==16||s==21){n_r[s,1]=n_r[2,1]*youth_multi_m}
    if (s>=2&s<=5||s>=12&s<=15||s>=17&s<=20||s>=22&s<=25) {n_r[s,1]=n_r[2,1]}
    if (s>=6&s<=10) {n_r[s,1]=n_r[6,1]}
    #Apply to all generations
    for (g in 2:generations) {n_r[s,g]=n_r[s,1]}
  }

# Parameter inputs for infectious stages
# Population growth
# for each trial, draw randomly from the input distribution (monthly)
  for(s in 1:sgroup){
    u=p_popgrowth[s,"Mean"]
    sd=p_popgrowth[s,"SD"]
    # Random draw (Note that rho can be either positive or negative)
    rho[s,1] <- rnorm(1,mean=u,sd=sd)
  }
  
  # For sub-population with same characters, we assume that the parameters are the same
  # For example, for all individuals that are white and aged 14-25 would share the same population growth rate
  for (s in 1:sgroup) {
    # We assume that the population growth does not change overtime
    for (g in 2:generations){rho[s,g]=rho[s,1]}
  }


# Natural history
# for each trial, draw randomly from the input distribution (monthly)
  # Assume that the natural history does not differ by engagement in care
  m_delta[1,1]=1
  
  # Disease progression from CD4>500 to CD4>350-500
  a=p_DisProgress[1,"Alpha"]
  b=p_DisProgress[1,"Beta"]
  # Random draw
  delta_1[1,1] <- rbeta(1,a,b)
  delta_1[1,1] <- 1-exp(log(1-delta_1[1,1])/12)

  a=p_DisProgress[4,"Alpha"]
  b=p_DisProgress[4,"Beta"]
  # Random draw
  delta_1[16,1] <- rbeta(1,a,b)
  delta_1[16,1] <- 1-exp(log(1-delta_1[16,1])/12)
  
  # Step length is set to 1 since sub-populations would have same natural history with the same CD4 level
  for (s in 1:sgroup) {
    if (s>=1 & s<=15) {
      delta_1[s,1] <- delta_1[1,1]
      delta_4[s,1] <- delta_1[1,1]*m_delta[1,1]}
    if (s>=16 & s<=25) {
      delta_1[s,1] <- delta_1[16,1]
      delta_4[s,1] <- delta_1[16,1]*m_delta[1,1]
    }
    # We assume that the natural history does not change overtime
    for (g in 2:generations){
      delta_1[s,g]=delta_1[s,1]
      delta_4[s,g]=delta_4[s,1]}
  }
  
  # Disease progression from CD4>350-500 to CD4>200-350
  a=p_DisProgress[2,"Alpha"]
  b=p_DisProgress[2,"Beta"]
  # Random draw
  delta_2[1,1] <- rbeta(1,a,b)
  delta_2[1,1] <- 1-exp(log(1-delta_2[1,1])/12)

  a=p_DisProgress[5,"Alpha"]
  b=p_DisProgress[5,"Beta"]
  # Random draw
  delta_2[16,1] <- rbeta(1,a,b)
  delta_2[16,1] <- 1-exp(log(1-delta_2[16,1])/12)

  for (s in 1:sgroup) {
    # Step length is set to 1 since sub-populations would have same natural history with the same CD4 level
    if (s>=1 & s<=15) {
      delta_2[s,1] <- delta_2[1,1]
      delta_5[s,1] <- delta_2[1,1]*m_delta[1,1]}
    if (s>=16 & s<=25) {
      delta_2[s,1] <- delta_2[16,1]
      delta_5[s,1] <- delta_2[16,1]*m_delta[1,1]}
    # We assume that the natural history does not change overtime
    for (g in 2:generations){
      delta_2[s,g]=delta_2[s,1]
      delta_5[s,g]=delta_5[s,1]}
  }
  
  # Disease progression from CD4>200-350 to CD4<200
  a=p_DisProgress[3,"Alpha"]
  b=p_DisProgress[3,"Beta"]
  # Random draw
  delta_3[1,1] <- rbeta(1,a,b)
  delta_3[1,1] <- 1-exp(log(1-delta_3[1,1])/12)

  a=p_DisProgress[6,"Alpha"]
  b=p_DisProgress[6,"Beta"]
  # Random draw
  delta_3[16,1] <- rbeta(1,a,b)
  delta_3[16,1] <- 1-exp(log(1-delta_3[16,1])/12)

  for (s in 1:sgroup) {
    # Step length is set to 1 since sub-populations would have same natural history with the same CD4 level
    if (s>=1 & s<=15) {
      delta_3[s,1] <- delta_3[1,1]
      delta_6[s,1] <- delta_3[1,1]*m_delta[1,1]}
    if (s>=16 & s<=25) {
      delta_3[s,1] <- delta_3[16,1]
      delta_6[s,1] <- delta_3[16,1]*m_delta[1,1]}
    # We assume that the natural history does not change overtime
    for (g in 2:generations){
      delta_3[s,g]=delta_3[s,1]
      delta_6[s,g]=delta_6[s,1]}
  }


# Probability of HIV testing
# For each trial, draw randomly from the distribution (monthly)
  # CD4 multiplier
  min=p_diagnosis[62,"Min"]
  max=p_diagnosis[62,"Max"]
  mode=p_diagnosis[62,"Mode"]
  cd4_multi=rtri(1,min,max,mode)
  
  correct_f=rtri(1,0.15,0.8,0.2)
  
  correct_m=rtri(1,0.2,0.8,0.3)
  
  # Low-risk urban women
  for (l in 1:5){
    a=p_diagnosis[l,"Alpha"]
    b=p_diagnosis[l,"Beta"]
    alpha_1[l,1]=rbeta(1,a,b)
    alpha_1[l,1]=-log(1-alpha_1[l,1])
    alpha_2[l,1]=alpha_1[l,1]
    alpha_3[l,1]=alpha_1[l,1]
    alpha_4[l,1]=alpha_1[l,1]*cd4_multi
    a=p_diagnosis[l+5,"Alpha"]
    b=p_diagnosis[l+5,"Beta"]
    alpha_1[l,73]=rbeta(1,a,b)
    alpha_1[l,73]=-log(1-alpha_1[l,73])
    alpha_2[l,73]=alpha_1[l,73]
    alpha_3[l,73]=alpha_1[l,73]
    alpha_4[l,73]=alpha_1[l,73]*cd4_multi
    a=p_diagnosis[l+10,"Alpha"]
    b=p_diagnosis[l+10,"Beta"]
    alpha_1[l,133]=rbeta(1,a,b)
    alpha_1[l,133]=-log(1-alpha_1[l,133])
    alpha_2[l,133]=alpha_1[l,133]
    alpha_3[l,133]=alpha_1[l,133]
    alpha_4[l,133]=alpha_1[l,133]*cd4_multi

    # Low-risk rural women
    a=p_diagnosis[l+15,"Alpha"]
    b=p_diagnosis[l+15,"Beta"]
    alpha_1[l+10,1]=rbeta(1,a,b)
    alpha_1[l+10,1]=-log(1-alpha_1[l+10,1])
    alpha_2[l+10,1]=alpha_1[l+10,1]
    alpha_3[l+10,1]=alpha_1[l+10,1]
    alpha_4[l+10,1]=alpha_1[l+10,1]*cd4_multi
    a=p_diagnosis[l+20,"Alpha"]
    b=p_diagnosis[l+20,"Beta"]
    alpha_1[l+10,73]=rbeta(1,a,b)
    alpha_1[l+10,73]=-log(1-alpha_1[l+10,73])
    alpha_2[l+10,73]=alpha_1[l+10,73]
    alpha_3[l+10,73]=alpha_1[l+10,73]
    alpha_4[l+10,73]=alpha_1[l+10,73]*cd4_multi
    a=p_diagnosis[l+25,"Alpha"]
    b=p_diagnosis[l+25,"Beta"]
    alpha_1[l+10,133]=rbeta(1,a,b)
    alpha_1[l+10,133]=-log(1-alpha_1[l+10,133])
    alpha_2[l+10,133]=alpha_1[l+10,133]
    alpha_3[l+10,133]=alpha_1[l+10,133]
    alpha_4[l+10,133]=alpha_1[l+10,133]*cd4_multi

    # Urban men
    a=p_diagnosis[l+30,"Alpha"]
    b=p_diagnosis[l+30,"Beta"]
    alpha_1[l+15,1]=rbeta(1,a,b)
    alpha_1[l+15,1]=-log(1-alpha_1[l+15,1])
    alpha_2[l+15,1]=alpha_1[l+15,1]
    alpha_3[l+15,1]=alpha_1[l+15,1]
    alpha_4[l+15,1]=alpha_1[l+15,1]*cd4_multi
    a=p_diagnosis[l+35,"Alpha"]
    b=p_diagnosis[l+35,"Beta"]
    alpha_1[l+15,73]=rbeta(1,a,b)
    alpha_1[l+15,73]=-log(1-alpha_1[l+15,73])
    alpha_2[l+15,73]=alpha_1[l+15,73]
    alpha_3[l+15,73]=alpha_1[l+15,73]
    alpha_4[l+15,73]=alpha_1[l+15,73]*cd4_multi
    a=p_diagnosis[l+40,"Alpha"]
    b=p_diagnosis[l+40,"Beta"]
    alpha_1[l+15,133]=rbeta(1,a,b)
    alpha_1[l+15,133]=-log(1-alpha_1[l+15,133])
    alpha_2[l+15,133]=alpha_1[l+15,133]
    alpha_3[l+15,133]=alpha_1[l+15,133]
    alpha_4[l+15,133]=alpha_1[l+15,133]*cd4_multi
    
    # Rural men
    a=p_diagnosis[l+45,"Alpha"]
    b=p_diagnosis[l+45,"Beta"]
    alpha_1[l+20,1]=rbeta(1,a,b)
    alpha_1[l+20,1]=-log(1-alpha_1[l+20,1])
    alpha_2[l+20,1]=alpha_1[l+20,1]
    alpha_3[l+20,1]=alpha_1[l+20,1]
    alpha_4[l+20,1]=alpha_1[l+20,1]*cd4_multi
    a=p_diagnosis[l+50,"Alpha"]
    b=p_diagnosis[l+50,"Beta"]
    alpha_1[l+20,73]=rbeta(1,a,b)
    alpha_1[l+20,73]=-log(1-alpha_1[l+20,73])
    alpha_2[l+20,73]=alpha_1[l+20,73]
    alpha_3[l+20,73]=alpha_1[l+20,73]
    alpha_4[l+20,73]=alpha_1[l+20,73]*cd4_multi
    a=p_diagnosis[l+55,"Alpha"]
    b=p_diagnosis[l+55,"Beta"]
    alpha_1[l+20,133]=rbeta(1,a,b)
    alpha_1[l+20,133]=-log(1-alpha_1[l+20,133])
    alpha_2[l+20,133]=alpha_1[l+20,133]
    alpha_3[l+20,133]=alpha_1[l+20,133]
    alpha_4[l+20,133]=alpha_1[l+20,133]*cd4_multi

    # High-risk women
    a=p_diagnosis[61,"Alpha"]
    b=p_diagnosis[61,"Beta"]
    alpha_1[6,1]=rbeta(1,a,b)
    alpha_1[6,1]=-log(1-alpha_1[6,1])
    alpha_2[6,1]=alpha_1[6,1]
    alpha_3[6,1]=alpha_1[6,1]
    alpha_4[6,1]=alpha_1[6,1]*cd4_multi
    for (s in 1:sgroup){
      if (s>=6&s<=10){
        alpha_1[s,1]=alpha_1[6,1]
        alpha_2[s,1]=alpha_2[6,1]
        alpha_3[s,1]=alpha_3[6,1]
        alpha_4[s,1]=alpha_4[6,1]
        alpha_1[s,73]=alpha_1[6,1]
        alpha_2[s,73]=alpha_2[6,1]
        alpha_3[s,73]=alpha_3[6,1]
        alpha_4[s,73]=alpha_4[6,1]
        alpha_1[s,133]=alpha_1[6,1]
        alpha_2[s,133]=alpha_2[6,1]
        alpha_3[s,133]=alpha_3[6,1]
        alpha_4[s,133]=alpha_4[6,1]}}}
  
  # Adjust for youth
  for (s in 1:sgroup) {
    if (s<=15){
      alpha_1[s,1]=alpha_1[s,1]*correct_f
      alpha_2[s,1]=alpha_2[s,1]*correct_f
      alpha_3[s,1]=alpha_3[s,1]*correct_f
      alpha_4[s,1]=alpha_4[s,1]*correct_f
      alpha_1[s,73]=alpha_1[s,73]*correct_f
      alpha_2[s,73]=alpha_2[s,73]*correct_f
      alpha_3[s,73]=alpha_3[s,73]*correct_f
      alpha_4[s,73]=alpha_4[s,73]*correct_f
      alpha_1[s,133]=alpha_1[s,133]*correct_f
      alpha_2[s,133]=alpha_2[s,133]*correct_f
      alpha_3[s,133]=alpha_3[s,133]*correct_f
      alpha_4[s,133]=alpha_4[s,133]*correct_f
    }
    if (s>15){
      alpha_1[s,1]=alpha_1[s,1]*correct_m
      alpha_2[s,1]=alpha_2[s,1]*correct_m
      alpha_3[s,1]=alpha_3[s,1]*correct_m
      alpha_4[s,1]=alpha_4[s,1]*correct_m
      alpha_1[s,73]=alpha_1[s,73]*correct_m
      alpha_2[s,73]=alpha_2[s,73]*correct_m
      alpha_3[s,73]=alpha_3[s,73]*correct_m
      alpha_4[s,73]=alpha_4[s,73]*correct_m
      alpha_1[s,133]=alpha_1[s,133]*correct_m
      alpha_2[s,133]=alpha_2[s,133]*correct_m
      alpha_3[s,133]=alpha_3[s,133]*correct_m
      alpha_4[s,133]=alpha_4[s,133]*correct_m
    }}
  
  for (s in 1:sgroup){
    # Convert to probability
    alpha_1[s,1]=1-exp(-alpha_1[s,1])
    alpha_2[s,1]=1-exp(-alpha_2[s,1])
    alpha_3[s,1]=1-exp(-alpha_3[s,1])
    alpha_4[s,1]=1-exp(-alpha_4[s,1])
    alpha_1[s,73]=1-exp(-alpha_1[s,73])
    alpha_2[s,73]=1-exp(-alpha_2[s,73])
    alpha_3[s,73]=1-exp(-alpha_3[s,73])
    alpha_4[s,73]=1-exp(-alpha_4[s,73])
    alpha_1[s,133]=1-exp(-alpha_1[s,133])
    alpha_2[s,133]=1-exp(-alpha_2[s,133])
    alpha_3[s,133]=1-exp(-alpha_3[s,133])
    alpha_4[s,133]=1-exp(-alpha_4[s,133])}
  
  #Apply to all generations
  for (s in 1:sgroup) {
    for (g in 2:generations) {
      if (g>=1&g<=72) {
        alpha_1[s,g]=alpha_1[s,1]
        alpha_2[s,g]=alpha_2[s,1]
        alpha_3[s,g]=alpha_3[s,1]
        alpha_4[s,g]=alpha_4[s,1]}
      if (g>=73&g<=132) {
        alpha_1[s,g]=alpha_1[s,73]
        alpha_2[s,g]=alpha_2[s,73]
        alpha_3[s,g]=alpha_3[s,73]
        alpha_4[s,g]=alpha_4[s,73]}
      if (g>=133&g<=generations) {
        alpha_1[s,g]=alpha_1[s,133]
        alpha_2[s,g]=alpha_2[s,133]
        alpha_3[s,g]=alpha_3[s,133]
        alpha_4[s,g]=alpha_4[s,133]}}}


# Probability of HIV linkage
# For each trial, draw randomly from the distribution (monthly)
  min=p_link[3,"Min"]
  max=p_link[3,"Max"]
  mode=p_link[3,"Mode"]
  cd4_multi=rtri(1,min,max,mode)
  
  # High-risk women
  min=p_link[1,"Min"]
  max=p_link[1,"Max"]
  mode=p_link[1,"Mode"]
  sigma_1[6,1]=rtri(1,min,max,mode)
  sigma_1[6,1]=-log(1-sigma_1[6,1])
  sigma_2[6,1]=sigma_1[6,1]
  sigma_3[6,1]=sigma_1[6,1]
  sigma_4[6,1]=sigma_1[6,1]*cd4_multi
  # Low-risk
  min=p_link[2,"Min"]
  max=p_link[2,"Max"]
  mode=p_link[2,"Mode"]
  sigma_1[1,1]=rtri(1,min,max,mode)
  sigma_1[1,1]=-log(1-sigma_1[1,1])
  sigma_2[1,1]=sigma_1[1,1]
  sigma_3[1,1]=sigma_1[1,1]
  sigma_4[1,1]=sigma_1[1,1]*cd4_multi
  
  # Convert to probabilities
  sigma_1[1,1]=1-exp(-sigma_1[1,1])
  sigma_2[1,1]=1-exp(-sigma_2[1,1])
  sigma_3[1,1]=1-exp(-sigma_3[1,1])
  sigma_4[1,1]=1-exp(-sigma_4[1,1])
  sigma_1[6,1]=1-exp(-sigma_1[6,1])
  sigma_2[6,1]=1-exp(-sigma_2[6,1])
  sigma_3[6,1]=1-exp(-sigma_3[6,1])
  sigma_4[6,1]=1-exp(-sigma_4[6,1])
  
  for (s in 1:sgroup){
    if (s>=6&s<=10){
      sigma_1[s,1]=sigma_1[6,1]
      sigma_2[s,1]=sigma_2[6,1]
      sigma_3[s,1]=sigma_3[6,1]
      sigma_4[s,1]=sigma_4[6,1]
    }else{sigma_1[s,1]=sigma_1[1,1]
    sigma_2[s,1]=sigma_2[1,1]
    sigma_3[s,1]=sigma_3[1,1]
    sigma_4[s,1]=sigma_4[1,1]}
    # We assume that the probability of diagnosis will change every two years (24 months)
    for (g in 2:generations){
      sigma_1[s,g]=sigma_1[s,1]
      sigma_2[s,g]=sigma_2[s,1]
      sigma_3[s,g]=sigma_3[s,1]
      sigma_4[s,g]=sigma_4[s,1]}}

# Probability of LTFU
# For each trial, draw randomly from the distribution 
  a=p_LTFU[1,"Alpha"]
  b=p_LTFU[1,"Beta"]
  #mode=p_LTFU[1,"Mode"]
  gamma_1[1,1]=rbeta(1,a,b)
  gamma_1[1,1]=1-exp(log(1-gamma_1[1,1])/12)
  
  a=p_LTFU[2,"Alpha"]
  b=p_LTFU[2,"Beta"]
  #mode=p_LTFU[2,"Mode"]
  gamma_2[1,1]=rbeta(1,a,b)
  gamma_2[1,1]=1-exp(log(1-gamma_2[1,1])/12)
  
  a=p_LTFU[3,"Alpha"]
  b=p_LTFU[3,"Beta"]
  #mode=p_LTFU[3,"Mode"]
  gamma_3[1,1]=rbeta(1,a,b)
  gamma_3[1,1]=1-exp(log(1-gamma_3[1,1])/12)
  
  a=p_LTFU[4,"Alpha"]
  b=p_LTFU[4,"Beta"]
  #mode=p_LTFU[4,"Mode"]
  gamma_4[1,1]=rbeta(1,a,b)
  gamma_4[1,1]=1-exp(log(1-gamma_4[1,1])/12)
  
  a=p_LTFU[5,"Alpha"]
  b=p_LTFU[5,"Beta"]
  #mode=p_LTFU[5,"Mode"]
  gamma_5[1,1]=rbeta(1,a,b)
  gamma_5[1,1]=1-exp(log(1-gamma_5[1,1])/12)
  gamma_9[1,1]=gamma_5[1,1]
  
  a=p_LTFU[6,"Alpha"]
  b=p_LTFU[6,"Beta"]
  #mode=p_LTFU[6,"Mode"]
  gamma_6[1,1]=rbeta(1,a,b)
  gamma_6[1,1]=1-exp(log(1-gamma_6[1,1])/12)
  gamma_10[1,1]=gamma_6[1,1]
  
  a=p_LTFU[7,"Alpha"]
  b=p_LTFU[7,"Beta"]
  #mode=p_LTFU[7,"Mode"]
  gamma_7[1,1]=rbeta(1,a,b)
  gamma_7[1,1]=1-exp(log(1-gamma_7[1,1])/12)
  gamma_11[1,1]=gamma_7[1,1]
  
  a=p_LTFU[8,"Alpha"]
  b=p_LTFU[8,"Beta"]
  #mode=p_LTFU[8,"Mode"]
  gamma_8[1,1]=rbeta(1,a,b)
  gamma_8[1,1]=1-exp(log(1-gamma_8[1,1])/12)
  gamma_12[1,1]=gamma_8[1,1]
  
  for (s in 1:sgroup){
    gamma_1[s,1]=gamma_1[1,1]
    gamma_2[s,1]=gamma_2[1,1]
    gamma_3[s,1]=gamma_3[1,1]
    gamma_4[s,1]=gamma_4[1,1]
    gamma_5[s,1]=gamma_5[1,1]
    gamma_6[s,1]=gamma_6[1,1]
    gamma_7[s,1]=gamma_7[1,1]
    gamma_8[s,1]=gamma_8[1,1]
    gamma_9[s,1]=gamma_9[1,1]
    gamma_10[s,1]=gamma_10[1,1]
    gamma_11[s,1]=gamma_11[1,1]
    gamma_12[s,1]=gamma_12[1,1]
    # We assume that the probability of diagnosis will be constant
    for (g in 2:generations){
      gamma_1[s,g]=gamma_1[s,1]
      gamma_2[s,g]=gamma_2[s,1]
      gamma_3[s,g]=gamma_3[s,1]
      gamma_4[s,g]=gamma_4[s,1]
      gamma_5[s,g]=gamma_5[s,1]
      gamma_6[s,g]=gamma_6[s,1]
      gamma_7[s,g]=gamma_7[s,1]
      gamma_8[s,g]=gamma_8[s,1]
      gamma_9[s,g]=gamma_9[s,1]
      gamma_10[s,g]=gamma_10[s,1]
      gamma_11[s,g]=gamma_11[s,1]
      gamma_12[s,g]=gamma_12[s,1]
    }
  }   

# Probability of On ART and Suppressed
# For each trial, draw randomly from the distribution 
  
  
  a=p_ARTVS[1,"Alpha"]
  b=p_ARTVS[1,"Beta"]
  #mode=p_ARTVS[1,"Mode"]
  theta_1[1,1]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[2,"Alpha"]
  b=p_ARTVS[2,"Beta"]
  #mode=p_ARTVS[2,"Mode"]
  theta_2[1,1]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[3,"Alpha"]
  b=p_ARTVS[3,"Beta"]
  #mode=p_ARTVS[3,"Mode"]
  theta_3[1,1]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[4,"Alpha"]
  b=p_ARTVS[4,"Beta"]
  #mode=p_ARTVS[4,"Mode"]
  theta_4[1,1]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[5,"Alpha"]
  b=p_ARTVS[5,"Beta"]
  #mode=p_ARTVS[5,"Mode"]
  theta_1[1,37]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[6,"Alpha"]
  b=p_ARTVS[6,"Beta"]
  #mode=p_ARTVS[6,"Mode"]
  theta_2[1,37]=1-exp(log(1-rbeta(1,a,b))/12)
  
  
  a=p_ARTVS[7,"Alpha"]
  b=p_ARTVS[7,"Beta"]
  #mode=p_ARTVS[7,"Mode"]
  theta_3[1,37]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[8,"Alpha"]
  b=p_ARTVS[8,"Beta"]
  #mode=p_ARTVS[8,"Mode"]
  theta_4[1,37]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[9,"Alpha"]
  b=p_ARTVS[9,"Beta"]
  #mode=p_ARTVS[9,"Mode"]
  theta_1[1,97]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[10,"Alpha"]
  b=p_ARTVS[10,"Beta"]
  #mode=p_ARTVS[10,"Mode"]
  theta_2[1,97]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[11,"Alpha"]
  b=p_ARTVS[11,"Beta"]
  #mode=p_ARTVS[11,"Mode"]
  theta_3[1,97]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[12,"Alpha"]
  b=p_ARTVS[12,"Beta"]
  #mode=p_ARTVS[12,"Mode"]
  theta_4[1,97]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[13,"Alpha"]
  b=p_ARTVS[13,"Beta"]
  #mode=p_ARTVS[13,"Mode"]
  theta_1[1,145]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[14,"Alpha"]
  b=p_ARTVS[14,"Beta"]
  #mode=p_ARTVS[14,"Mode"]
  theta_2[1,145]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[15,"Alpha"]
  b=p_ARTVS[15,"Beta"]
  #mode=p_ARTVS[15,"Mode"]
  theta_3[1,145]=1-exp(log(1-rbeta(1,a,b))/12)
  
  a=p_ARTVS[16,"Alpha"]
  b=p_ARTVS[16,"Beta"]
  #mode=p_ARTVS[16,"Mode"]
  theta_4[1,145]=1-exp(log(1-rbeta(1,a,b))/12)
  
  for (s in 1:sgroup) {
    theta_1[s,1]=theta_1[1,1]
    theta_2[s,1]=theta_2[1,1]
    theta_3[s,1]=theta_3[1,1]
    theta_4[s,1]=theta_4[1,1]
    theta_1[s,37]=theta_1[1,37]
    theta_2[s,37]=theta_2[1,37]
    theta_3[s,37]=theta_3[1,37]
    theta_4[s,37]=theta_4[1,37]
    theta_1[s,97]=theta_1[1,97]
    theta_2[s,97]=theta_2[1,97]
    theta_3[s,97]=theta_3[1,97]
    theta_4[s,97]=theta_4[1,97]
    theta_1[s,145]=theta_1[1,145]
    theta_2[s,145]=theta_2[1,145]
    theta_3[s,145]=theta_3[1,145]
    theta_4[s,145]=theta_4[1,145]
    if (s>=6&s<11){
      theta_1[s,1]=theta_1[1,1]#*runif(1,min=0.2,max=1)
      theta_2[s,1]=theta_2[1,1]#*runif(1,min=0.2,max=1)
      theta_3[s,1]=theta_3[1,1]#*runif(1,min=0.2,max=1)
      theta_4[s,1]=theta_4[1,1]#*runif(1,min=0.2,max=1)
      theta_1[s,37]=theta_1[1,37]#*runif(1,min=0.2,max=1)
      theta_2[s,37]=theta_2[1,37]#*runif(1,min=0.2,max=1)
      theta_3[s,37]=theta_3[1,37]#*runif(1,min=0.2,max=1)
      theta_4[s,37]=theta_4[1,37]#*runif(1,min=0.2,max=1)
      theta_1[s,97]=theta_1[1,97]#*runif(1,min=0.2,max=1)
      theta_2[s,97]=theta_2[1,97]#*runif(1,min=0.2,max=1)
      theta_3[s,97]=theta_3[1,97]#*runif(1,min=0.2,max=1)
      theta_4[s,97]=theta_4[1,97]#*runif(1,min=0.2,max=1)
      theta_1[s,145]=theta_1[1,145]#*runif(1,min=0.2,max=1)
      theta_2[s,145]=theta_2[1,145]#*runif(1,min=0.2,max=1)
      theta_3[s,145]=theta_3[1,145]#*runif(1,min=0.2,max=1)
      theta_4[s,145]=theta_4[1,145]#*runif(1,min=0.2,max=1)
    }
    # # We assume that the probability of diagnosis will change every years (12 months)
    for (g in 2:generations){
      if (g<=36) {
        theta_1[s,g]=theta_1[s,1]
        theta_2[s,g]=theta_2[s,1]
        theta_3[s,g]=theta_3[s,1]
        theta_4[s,g]=theta_4[s,1]}
      if (g>=37&g<=96) {
        theta_1[s,g]=theta_1[s,37]
        theta_2[s,g]=theta_2[s,37]
        theta_3[s,g]=theta_3[s,37]
        theta_4[s,g]=theta_4[s,37]}
      if (g>=97&g<=144) {
        theta_1[s,g]=theta_1[s,97]
        theta_2[s,g]=theta_2[s,97]
        theta_3[s,g]=theta_3[s,97]
        theta_4[s,g]=theta_4[s,97]}
      if (g>=145&g<=generations) {
        theta_1[s,g]=theta_1[s,145]
        theta_2[s,g]=theta_2[s,145]
        theta_3[s,g]=theta_3[s,145]
        theta_4[s,g]=theta_4[s,145]}}
  }   
  
    for (s in 1:sgroup){
      for (g in 1:generations){
        theta_1[s,g]=-log(1-theta_1[s,g])
        theta_2[s,g]=-log(1-theta_2[s,g])
        theta_3[s,g]=-log(1-theta_3[s,g])
        theta_4[s,g]=-log(1-theta_4[s,g])}}
  
  correct_female=rtri(1,1,2,1.5)
  correct_old=rtri(1,1,2,1.5)
    for (s in 1:sgroup){
      for (g in 1:generations){
        if (s>=2&s<=5||s>=7&s<=10||s>=12&s<=15){
          theta_1[s,g]=theta_1[s,g]*correct_female*correct_old
          theta_2[s,g]=theta_2[s,g]*correct_female*correct_old
          theta_3[s,g]=theta_3[s,g]*correct_female*correct_old
          theta_4[s,g]=theta_4[s,g]*correct_female*correct_old}
        if (s>=17&s<=20||s>=22&s<=25){
          theta_1[s,g]=theta_1[s,g]*correct_old
          theta_2[s,g]=theta_2[s,g]*correct_old
          theta_3[s,g]=theta_3[s,g]*correct_old
          theta_4[s,g]=theta_4[s,g]*correct_old}
      }
  }

  for (s in 1:sgroup){
    for (g in 1:generations){
      theta_1[s,g]=1-exp(-theta_1[s,g])
      theta_2[s,g]=1-exp(-theta_2[s,g])
      theta_3[s,g]=1-exp(-theta_3[s,g])
      theta_4[s,g]=1-exp(-theta_4[s,g])}}
  
# Probability of On ART and Not suppressed
# For each trial, draw randomly from the distribution 
  # Extract the baseline value
  a=p_ARTnotVS[1,"Alpha"]
  b=p_ARTnotVS[1,"Beta"]
  #mode=p_ARTnotVS[1,"Mode"]
  # Random draw
  psi_1[1,1]=rbeta(1,a,b)/12

  a=p_ARTnotVS[2,"Alpha"]
  b=p_ARTnotVS[2,"Beta"]
  #mode=p_ARTnotVS[2,"Mode"]
  # Random draw
  psi_2[1,1]=rbeta(1,a,b)/12

  a=p_ARTnotVS[3,"Alpha"]
  b=p_ARTnotVS[3,"Beta"]
  #mode=p_ARTnotVS[3,"Mode"]
  # Random draw
  psi_3[1,1]=rbeta(1,a,b)/12

  a=p_ARTnotVS[4,"Alpha"]
  b=p_ARTnotVS[4,"Beta"]
  #mode=p_ARTnotVS[4,"Mode"]
  # Random draw
  psi_4[1,1]=rbeta(1,a,b)/12
  # p_ARTnotVS_trial[1]=psi_1[1,1]
  
  min=0.2
  max=1.0
  mode=0.5
  # Random draw
  multi_female=rtri(1,min,max,mode)
  # p_ARTnotVS_trial[1]=psi_1[1,1]
  
  min=1.0
  max=3.0
  mode=2.5
  # Random draw
  multi_male=rtri(1,min,max,mode)
  # p_ARTnotVS_trial[1]=psi_1[1,1]
  
  for (s in 1:sgroup) {
    if (s>=1&s<=15){
    # Apply baseline value to all sub-population
    psi_1[s,1]=psi_1[1,1]*multi_female
    psi_2[s,1]=psi_2[1,1]*multi_female
    psi_3[s,1]=psi_3[1,1]*multi_female
    psi_4[s,1]=psi_4[1,1]*multi_female}
    if (s>=16&s<=25){
      # Apply baseline value to all sub-population
      psi_1[s,1]=psi_1[1,1]*multi_male
      psi_2[s,1]=psi_2[1,1]*multi_male
      psi_3[s,1]=psi_3[1,1]*multi_male
      psi_4[s,1]=psi_4[1,1]*multi_male}}
    
    for (s in 1:sgroup) {
      # Convert rate to probability
        psi_1[s,1]=1-exp(-psi_1[s,1])
        psi_2[s,1]=1-exp(-psi_2[s,1])
        psi_3[s,1]=1-exp(-psi_3[s,1])
        psi_4[s,1]=1-exp(-psi_4[s,1])
    
    
    # We assume that the probability is constant overtime
    for (g in 2:generations){
      psi_1[s,g]=psi_1[s,1]
      psi_2[s,g]=psi_2[s,1]
      psi_3[s,g]=psi_3[s,1]
      psi_4[s,g]=psi_4[s,1]}}

# Probability of death
  # Multiplier 
  min=p_death[41,"Min"]
  max=p_death[41,"Max"]
  mode=p_death[41,"Mode"]
  pre_cd4_1=1
  
  min=p_death[42,"Min"]
  max=p_death[42,"Max"]
  mode=p_death[42,"Mode"]
  pre_cd4_2=1
  
  min=p_death[43,"Min"]
  max=p_death[43,"Max"]
  mode=p_death[43,"Mode"]
  pre_cd4_3=1
  
  min=p_death[44,"Min"]
  max=p_death[44,"Max"]
  mode=p_death[44,"Mode"]
  pre_cd4_4=rtri(1,min,max,mode)
  
  # On ART
  for (l in 1:5){
    # Women
    min=p_death[l,"Min"]
    max=p_death[l,"Max"]
    mode=p_death[l,"Mode"]
    mu_9[l,1]=rtri(1,min,max,mode)
    mu_1[l,1]=mu_9[l,1]*pre_cd4_1
    mu_5[l,1]=mu_1[l,1]
    mu_13[l,1]=mu_1[l,1]*deathmulti
    mu_1[l+5,1]=mu_1[l,1]
    mu_5[l+5,1]=mu_5[l,1]
    mu_9[l+5,1]=mu_9[l,1]
    mu_13[l+5,1]=mu_13[l,1]
    mu_1[l+10,1]=mu_1[l,1]
    mu_5[l+10,1]=mu_5[l,1]
    mu_9[l+10,1]=mu_9[l,1]
    mu_13[l+10,1]=mu_13[l,1]
    # 
    min=p_death[l+5,"Min"]
    max=p_death[l+5,"Max"]
    mode=p_death[l+5,"Mode"]
    mu_10[l,1]=rtri(1,min,max,mode)
    mu_2[l,1]=mu_10[l,1]*pre_cd4_2
    mu_6[l,1]=mu_2[l,1]
    mu_14[l,1]=mu_2[l,1]*deathmulti
    mu_2[l+5,1]=mu_2[l,1]
    mu_6[l+5,1]=mu_6[l,1]
    mu_10[l+5,1]=mu_10[l,1]
    mu_14[l+5,1]=mu_14[l,1]
    mu_2[l+10,1]=mu_2[l,1]
    mu_6[l+10,1]=mu_6[l,1]
    mu_10[l+10,1]=mu_10[l,1]
    mu_14[l+10,1]=mu_14[l,1]
    # 
    min=p_death[l+10,"Min"]
    max=p_death[l+10,"Max"]
    mode=p_death[l+10,"Mode"]
    mu_11[l,1]=rtri(1,min,max,mode)
    mu_3[l,1]=mu_11[l,1]*pre_cd4_3
    mu_7[l,1]=mu_3[l,1]
    mu_15[l,1]=mu_3[l,1]*deathmulti
    mu_3[l+5,1]=mu_3[l,1]
    mu_7[l+5,1]=mu_7[l,1]
    mu_11[l+5,1]=mu_11[l,1]
    mu_15[l+5,1]=mu_15[l,1]
    mu_3[l+10,1]=mu_3[l,1]
    mu_7[l+10,1]=mu_7[l,1]
    mu_11[l+10,1]=mu_11[l,1]
    mu_15[l+10,1]=mu_15[l,1]
    # mu_3[l+15,1]=mu_3[l,1]
    # mu_7[l+15,1]=mu_7[l,1]
    # mu_11[l+15,1]=mu_11[l,1]
    # mu_15[l+15,1]=mu_15[l,1]
    # 
    min=p_death[l+15,"Min"]
    max=p_death[l+15,"Max"]
    mode=p_death[l+15,"Mode"]
    mu_12[l,1]=rtri(1,min,max,mode)
    mu_4[l,1]=mu_12[l,1]*pre_cd4_4
    mu_8[l,1]=mu_4[l,1]
    mu_16[l,1]=mu_4[l,1]*deathmulti
    mu_4[l+5,1]=mu_4[l,1]
    mu_8[l+5,1]=mu_8[l,1]
    mu_12[l+5,1]=mu_12[l,1]
    mu_16[l+5,1]=mu_16[l,1]
    mu_4[l+10,1]=mu_4[l,1]
    mu_8[l+10,1]=mu_8[l,1]
    mu_12[l+10,1]=mu_12[l,1]
    mu_16[l+10,1]=mu_16[l,1]
    # mu_4[l+15,1]=mu_4[l,1]
    # mu_8[l+15,1]=mu_8[l,1]
    # mu_12[l+15,1]=mu_12[l,1]
    # mu_16[l+15,1]=mu_16[l,1]
    
    #Men
    min=p_death[l+20,"Min"]
    max=p_death[l+20,"Max"]
    mode=p_death[l+20,"Mode"]
    mu_9[l+15,1]=rtri(1,min,max,mode)
    mu_1[l+15,1]=mu_9[l+15,1]*pre_cd4_1
    mu_5[l+15,1]=mu_1[l+15,1]
    mu_13[l+15,1]=mu_1[l+15,1]*deathmulti
    mu_1[l+20,1]=mu_1[l+15,1]
    mu_5[l+20,1]=mu_5[l+15,1]
    mu_9[l+20,1]=mu_9[l+15,1]
    mu_13[l+20,1]=mu_13[l+15,1]
    
    min=p_death[l+25,"Min"]
    max=p_death[l+25,"Max"]
    mode=p_death[l+25,"Mode"]
    mu_10[l+15,1]=rtri(1,min,max,mode)
    mu_2[l+15,1]=mu_10[l+15,1]*pre_cd4_2
    mu_6[l+15,1]=mu_2[l+15,1]
    mu_14[l+15,1]=mu_2[l+15,1]*deathmulti
    mu_2[l+20,1]=mu_2[l+15,1]
    mu_6[l+20,1]=mu_6[l+15,1]
    mu_10[l+20,1]=mu_10[l+15,1]
    mu_14[l+20,1]=mu_14[l+15,1]

    min=p_death[l+30,"Min"]
    max=p_death[l+30,"Max"]
    mode=p_death[l+30,"Mode"]
    mu_11[l+15,1]=rtri(1,min,max,mode)
    mu_3[l+15,1]=mu_11[l+15,1]*pre_cd4_3
    mu_7[l+15,1]=mu_3[l+15,1]
    mu_15[l+15,1]=mu_3[l+15,1]*deathmulti
    mu_3[l+20,1]=mu_3[l+15,1]
    mu_7[l+20,1]=mu_7[l+15,1]
    mu_11[l+20,1]=mu_11[l+15,1]
    mu_15[l+20,1]=mu_15[l+15,1]
    
    min=p_death[l+35,"Min"]
    max=p_death[l+35,"Max"]
    mode=p_death[l+35,"Mode"]
    mu_12[l+15,1]=rtri(1,min,max,mode)
    mu_4[l+15,1]=mu_12[l+15,1]*pre_cd4_4
    mu_8[l+15,1]=mu_4[l+15,1]
    mu_16[l+15,1]=mu_4[l+15,1]*deathmulti
    mu_4[l+20,1]=mu_4[l+15,1]
    mu_8[l+20,1]=mu_8[l+15,1]
    mu_12[l+20,1]=mu_12[l+15,1]
    mu_16[l+20,1]=mu_16[l+15,1]
  }
  
  for (s in 1:sgroup) {
    # Convert to probabilities
    mu_1[s,1]=1-exp(-mu_1[s,1]/12)
    mu_2[s,1]=1-exp(-mu_2[s,1]/12)
    mu_3[s,1]=1-exp(-mu_3[s,1]/12)
    mu_4[s,1]=1-exp(-mu_4[s,1]/12)
    mu_5[s,1]=1-exp(-mu_5[s,1]/12)
    mu_6[s,1]=1-exp(-mu_6[s,1]/12)
    mu_7[s,1]=1-exp(-mu_7[s,1]/12)
    mu_8[s,1]=1-exp(-mu_8[s,1]/12)
    mu_9[s,1]=1-exp(-mu_9[s,1]/12)
    mu_10[s,1]=1-exp(-mu_10[s,1]/12)
    mu_11[s,1]=1-exp(-mu_11[s,1]/12)
    mu_12[s,1]=1-exp(-mu_12[s,1]/12)
    mu_13[s,1]=1-exp(-mu_13[s,1]/12)
    mu_14[s,1]=1-exp(-mu_14[s,1]/12)
    mu_15[s,1]=1-exp(-mu_15[s,1]/12)
    mu_16[s,1]=1-exp(-mu_16[s,1]/12)}
  
  # We assume that death does not change overtime
  
  for (g in 1:generations) {
    for (s in 1:sgroup){
      #undiagnosed, diagnosed but not linked to care, LTFU
      mu_1[s,g]=mu_1[s,1]
      mu_2[s,g]=mu_2[s,1]
      mu_3[s,g]=mu_3[s,1]
      mu_4[s,g]=mu_4[s,1]
      # linked to care
      mu_5[s,g]=mu_5[s,1]
      mu_6[s,g]=mu_6[s,1]
      mu_7[s,g]=mu_7[s,1]
      mu_8[s,g]=mu_8[s,1]
      # On ART & Suppressed
      mu_9[s,g]=mu_9[s,1]
      mu_10[s,g]=mu_10[s,1]
      mu_11[s,g]=mu_11[s,1]
      mu_12[s,g]=mu_12[s,1]
      # On ART & Not suppressed
      mu_13[s,g]=mu_13[s,1]
      mu_14[s,g]=mu_14[s,1]
      mu_15[s,g]=mu_15[s,1]
      mu_16[s,g]=mu_16[s,1]}
  }


# Probability of re-engagement in care
# For each trial, draw randomly from the distribution 
  min=p_reengage[1,"Min"]
  max=p_reengage[1,"Max"]
  mode=p_reengage[1,"Mode"]
  # Random draw
  tau_1[1,1]=1-exp(-rtri(1,min,max,mode)/12)
  tau_2[1,1]=1-exp(-rtri(1,min,max,mode)/12)
  tau_3[1,1]=1-exp(-rtri(1,min,max,mode)/12)
  tau_4[1,1]=1-exp(-rtri(1,min,max,mode)/12)
  # p_reengage_trial[1]=tau_1[1,1]
  # We assume that it is the same for all sub-population and overtime
  for (s in 1:sgroup) {
    tau_1[s,1]=tau_1[1,1]
    tau_2[s,1]=tau_2[1,1]
    tau_3[s,1]=tau_3[1,1]
    tau_4[s,1]=tau_4[1,1]
    for (g in 2:generations) {
      tau_1[s,g]=tau_1[s,1]
      tau_2[s,g]=tau_2[s,1]
      tau_3[s,g]=tau_3[s,1]
      tau_4[s,g]=tau_4[s,1]
    }
  }


# Initalize baseline population characteristics
  #Initalize parameters for each sub-group
  for (s in 1:(sgroup)) {
    # Calculate the initial HIV prevalence 
    
    # Total HIV infected individuals
    N[s,1] <- X1[s,1]+X2[s,1]+X3[s,1]+X4[s,1]+X5[s,1]+X6[s,1]+X7[s,1]+X8[s,1]+X9[s,1]+X10[s,1]+X11[s,1]+X12[s,1]+X13[s,1]+X14[s,1]+X15[s,1]+X16[s,1]+X17[s,1]+X18[s,1]+X19[s,1]+X20[s,1]+X21[s,1]+X22[s,1]+X23[s,1]+X24[s,1]
    
    # Total population of individuals who are not viral suppressed
    Z[s,1] <- X1[s,1]+X2[s,1]+X3[s,1]+X4[s,1]+X5[s,1]+X6[s,1]+X7[s,1]+X8[s,1]+X9[s,1]+X10[s,1]+X11[s,1]+X12[s,1]+X13[s,1]+X14[s,1]+X15[s,1]+X16[s,1]+X21[s,1]+X22[s,1]+X23[s,1]+X24[s,1]
    
    # Total population of individuals who are viral suppressed
    V[s,1] <- X17[s,1]+X18[s,1]+X19[s,1]+X20[s,1]
    
    # Total population of individuals infected but undiagnosed 
    # Undiagnosed[s,1] <- X1[s,1]+X2[s,1]+X3[s,1]+X4[s,1]
    
    # Total population of individuals diagnosed with HIV
    Diagnosed[s,1] <- X5[s,1]+X6[s,1]+X7[s,1]+X8[s,1]+X9[s,1]+X10[s,1]+X11[s,1]+X12[s,1]+X13[s,1]+X14[s,1]+X15[s,1]+X16[s,1]+X17[s,1]+X18[s,1]+X19[s,1]+X20[s,1]+X21[s,1]+X22[s,1]+X23[s,1]+X24[s,1]
    
    # Total population of individuals newly diagnosed
    New_Diagnosed[s,1] <- X1[s,1]*alpha_1[s,g]+X2[s,1]*alpha_2[s,g]+X3[s,1]*alpha_3[s,g]+X4[s,1]*alpha_4[s,g]
    
    # Total population of inidivudlas linked to care
    # Linked[s,1] <- X9[s,1]+X10[s,1]+X11[s,1]+X12[s,1]
    
    # Total population of individuals lost
    # Lost[s,1] <-  X13[s,1]+X14[s,1]+X15[s,1]+X16[s,1]
    
    # Total population of individuals on ART and suppressed
    Suppressed[s,1] <- X17[s,1]+X18[s,1]+X19[s,1]+X20[s,1]
    
    # Total population of individuals on ART and not suppressed
    # Not_Suppressed[s,1] <- X21[s,1]+X22[s,1]+X23[s,1]+X24[s,1]
    
    # Total population of individuals on ART
    On_ART[s,1] <- X17[s,1]+X18[s,1]+X19[s,1]+X20[s,1]+X21[s,1]+X22[s,1]+X23[s,1]+X24[s,1]
    
    if (N[s,1]==0) {
      
      p_r[s,1] <- 1
      p_r_v[s,1] <- 0
      
    } else {
      # Prevalence of individuals who are viral supprressed among people living with HIV
      p_r_v[s,1] <- V[s,1]/N[s,1]  
      
      # percentage of individuals who are not viral supprressed among people living with HIV
      p_r[s,1] <- 1-p_r_v[s,1]
    }
    
    # initial value for force of infection and new infection (set to null)
    lambda_t_r[s,1] <- NA
    NI[s,1] <- NA
    N_total[1]<- NA
    #lambda_inf[s,1] <- 0
    # lambda_t_r[s,1] <- 0
    # NI[s,1] <- 0
    # N_total[1,1]<- 0
  }

# Start of the loop for each trial
  # Start the loop by generations
  # The loop of generations goes before loops of sub-groups since operations in time period t+1 depends on the outcomes from all sub-groups in time period t.
  # The generations starts from time 2 so that the equations in R are in line with the difference equation. 
  for (g in 2:(generations))
  {
    for (s in 1:sgroup) {
      #for (j in 1:sgroup) {
      
      # total population of low risk urban women from last generation, who will mix with urban men 
      N_lr_u_f[s,g-1] <- S[1,g-1]+S[2,g-1]+S[3,g-1]+S[4,g-1]+S[5,g-1]+N[1,g-1]+N[2,g-1]+N[3,g-1]+N[4,g-1]+N[5,g-1]
      # total population of high risk urban women from last generation, who will mix with urban men 
      N_hr_f[s,g-1] <- S[6,g-1]+S[7,g-1]+S[8,g-1]+S[9,g-1]+S[10,g-1]+N[6,g-1]+N[7,g-1]+N[8,g-1]+N[9,g-1]+N[10,g-1]
      # total population of low risk rural women from last generation, who will mix with rural men 
      N_lr_r_f[s,g-1] <- S[11,g-1]+S[12,g-1]+S[13,g-1]+S[14,g-1]+S[15,g-1]+N[11,g-1]+N[12,g-1]+N[13,g-1]+N[14,g-1]+N[15,g-1]
      # total population of urban men from last generation, who will mix with high risk women and low risk urban women
      N_u_m[s,g-1] <- S[16,g-1]+S[17,g-1]+S[18,g-1]+S[19,g-1]+S[20,g-1]+N[16,g-1]+N[17,g-1]+N[18,g-1]+N[19,g-1]+N[20,g-1]
      # total population of rural men from last generation, who will mix with high risk women and low risk rural women
      N_r_m[s,g-1] <- S[21,g-1]+S[22,g-1]+S[23,g-1]+S[24,g-1]+S[25,g-1]+N[21,g-1]+N[22,g-1]+N[23,g-1]+N[24,g-1]+N[25,g-1]
      #}
      # Growth rate in this generation
      # rho_2[s,g-1] <- rho[s,g]*(1+N[s,g]/S[s,g])
      
      # low risk urban women mixed with urban men
      if (s>=1&s<=5) 
      {
        for (j in 16:20) 
        {
          # force of infection low risk urban women
          #lambda_t_r[s,g] <-  lambda_t_r[s,g]+((N[j,g-1]+S[j,g-1])/N_u_m[s,g-1])*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*(1-(1-(omega_r[s,g]*beta[s,g]*(kappa[s,g]*p_r_v[j,g-1]+p_r[j,g-1])))^n_r[s,g])
          lambda_inf[s,g] <-  lambda_inf[s,g]*(1-((N[j,g-1]+S[j,g-1])/N_u_m[s,g-1])*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*beta[s,g]*((1-upsilon[s,g])*p_r_v[j,g-1]+p_r[j,g-1]))^((1-epsilon[s,g]*c_r[s,g])*n_r[s,g])
          
        } 
      }
      
      # high risk women mixed with urban men
      if (s>=6&s<=10) 
      {
        for (j in 16:20) 
        {
          # force of infection (high risk urban women)
          #lambda_t_r[s,g] <-  lambda_t_r[s,g]+((N[j,g-1]+S[j,g-1])/N_u_m[s,g-1])*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*(1-(1-(omega_r[s,g]*beta[s,g]*(kappa[s,g]*p_r_v[j,g-1]+p_r[j,g-1])))^n_r[s,g])
          lambda_inf[s,g] <-  lambda_inf[s,g]*(1-((N[j,g-1]+S[j,g-1])/N_u_m[s,g-1])*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*beta[s,g]*((1-upsilon[s,g])*p_r_v[j,g-1]+p_r[j,g-1]))^((1-epsilon[s,g]*c_r[s,g])*n_r[s,g])
          
        }
      }
      
      # low risk rural women mixed with rural men
      if (s>=11&s<=15) 
      {
        for (j in 21:25) {
          # force of infection (low risk rural women)
          #lambda_t_r[s,g] <-  lambda_t_r[s,g]+((N[j,g-1]+S[j,g-1])/N_u_m[s,g-1])*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*(1-(1-(omega_r[s,g]*beta[s,g]*(kappa[s,g]*p_r_v[j,g-1]+p_r[j,g-1])))^n_r[s,g])
          lambda_inf[s,g] <-  lambda_inf[s,g]*(1-((N[j,g-1]+S[j,g-1])/N_r_m[s,g-1])*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*beta[s,g]*((1-upsilon[s,g])*p_r_v[j,g-1]+p_r[j,g-1]))^((1-epsilon[s,g]*c_r[s,g])*n_r[s,g])
          
        }
      }
      
      # urban men mixed with low risk urban women and high risk women
      if (s>=16&s<=20) 
      {
        for (j in 1:10) 
        {
          # force of infection (urban men)
          #lambda_t_r[s,g] <-  lambda_t_r[s,g]+((N[j,g-1]+S[j,g-1])/N_u_m[s,g-1])*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*(1-(1-(omega_r[s,g]*beta[s,g]*(kappa[s,g]*p_r_v[j,g-1]+p_r[j,g-1])))^n_r[s,g])
          lambda_inf[s,g] <-  lambda_inf[s,g]*(1-((N[j,g-1]+S[j,g-1])/(N_lr_u_f[s,g-1]+N_hr_f[s,g-1]))*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*beta[s,g]*((1-upsilon[s,g])*p_r_v[j,g-1]+p_r[j,g-1]))^((1-epsilon[s,g]*c_r[s,g])*n_r[s,g])
          
        }
      }
      
      # rural men mixed with low risk rural women
      if (s>=21&s<=25) 
      {
        for (j in 11:15) 
        {
          # force of infection (rural men)
          #lambda_t_r[s,g] <-  lambda_t_r[s,g]+((N[j,g-1]+S[j,g-1])/N_u_m[s,g-1])*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*(1-(1-(omega_r[s,g]*beta[s,g]*(kappa[s,g]*p_r_v[j,g-1]+p_r[j,g-1])))^n_r[s,g])
          lambda_inf[s,g] <-  lambda_inf[s,g]*(1-((N[j,g-1]+S[j,g-1])/N_lr_r_f[s,g-1])*(N[j,g-1]/(N[j,g-1]+S[j,g-1]))*beta[s,g]*((1-upsilon[s,g])*p_r_v[j,g-1]+p_r[j,g-1]))^((1-epsilon[s,g]*c_r[s,g])*n_r[s,g])
          
        }
      }
      
      # Calculate new infections
      #NI[s,g] <- lambda_t_r[s,g]*S[s,g-1]
      lambda_t_r[s,g]<-(1-lambda_inf[s,g])
      NI[s,g] <- (1-lambda_inf[s,g])*S[s,g-1]
      
      # Difference equation for each compartments
      # Difference equation for susceptable population
      S[s,g] <- S[s,g-1]+rho[s,g]*S[s,g-1]-lambda_t_r[s,g]*S[s,g-1]
      
      # Difference equations for compartments X1 - X4 (undiagnosed individuals)
      
      X1[s,g] <- X1[s,g-1]+lambda_t_r[s,g]*S[s,g-1]-delta_1[s,g]*X1[s,g-1]-alpha_1[s,g]*X1[s,g-1]-mu_1[s,g]*X1[s,g-1]
      
      X2[s,g] <- X2[s,g-1]+delta_1[s,g]*X1[s,g-1]-delta_2[s,g]*X2[s,g-1]-alpha_2[s,g]*X2[s,g-1]-mu_2[s,g]*X2[s,g-1]
      
      X3[s,g] <- X3[s,g-1]+delta_2[s,g]*X2[s,g-1]-delta_3[s,g]*X3[s,g-1]-alpha_3[s,g]*X3[s,g-1]-mu_3[s,g]*X3[s,g-1]
      
      X4[s,g] <- X4[s,g-1]+delta_3[s,g]*X3[s,g-1]-alpha_4[s,g]*X4[s,g-1]-mu_4[s,g]*X4[s,g-1]
      
      
      
      # Difference equations for compartments X5 - X8 (diagnosed individuals)
      
      X5[s,g] <- X5[s,g-1]+alpha_1[s,g]*X1[s,g-1]-delta_1[s,g]*X5[s,g-1]-sigma_1[s,g]*X5[s,g-1]-mu_1[s,g]*X5[s,g-1]
      
      X6[s,g] <- X6[s,g-1]+alpha_2[s,g]*X2[s,g-1]+delta_1[s,g]*X5[s,g-1]-delta_2[s,g]*X6[s,g-1]-sigma_2[s,g]*X6[s,g-1]-mu_2[s,g]*X6[s,g-1]
      
      X7[s,g] <- X7[s,g-1]+alpha_3[s,g]*X3[s,g-1]+delta_2[s,g]*X6[s,g-1]-delta_3[s,g]*X7[s,g-1]-sigma_3[s,g]*X7[s,g-1]-mu_3[s,g]*X7[s,g-1]
      
      X8[s,g] <- X8[s,g-1]+alpha_4[s,g]*X4[s,g-1]+delta_3[s,g]*X7[s,g-1]-sigma_4[s,g]*X8[s,g-1]-mu_4[s,g]*X8[s,g-1]
      
      
      
      # Difference equations for compartments X9 - X12 (individuals linked to care)
      
      X9[s,g] <- X9[s,g-1]+sigma_1[s,g]*X5[s,g-1]-delta_4[s,g]*X9[s,g-1]-theta_1[s,g]*X9[s,g-1]-gamma_1[s,g]*X9[s,g-1]-mu_5[s,g]*X9[s,g-1]
      
      X10[s,g] <- X10[s,g-1]+sigma_2[s,g]*X6[s,g-1]+delta_4[s,g]*X9[s,g-1]-delta_5[s,g]*X10[s,g-1]-theta_2[s,g]*X10[s,g-1]-gamma_2[s,g]*X10[s,g-1]-mu_6[s,g]*X10[s,g-1]
      
      X11[s,g] <- X11[s,g-1]+sigma_3[s,g]*X7[s,g-1]+delta_5[s,g]*X10[s,g-1]-delta_6[s,g]*X11[s,g-1]-theta_3[s,g]*X11[s,g-1]-gamma_3[s,g]*X11[s,g-1]-mu_7[s,g]*X11[s,g-1]
      
      X12[s,g] <- X12[s,g-1]+sigma_4[s,g]*X8[s,g-1]+delta_6[s,g]*X11[s,g-1]-theta_4[s,g]*X12[s,g-1]-gamma_4[s,g]*X12[s,g-1]-mu_8[s,g]*X12[s,g-1]
      
      
      
      # Difference equations for compartments X13 - X16 (individuals lost from care)
      
      X13[s,g] <- X13[s,g-1]+gamma_1[s,g]*X9[s,g-1]+gamma_5[s,g]*X17[s,g-1]+gamma_9[s,g]*X21[s,g-1]-delta_1[s,g]*X13[s,g-1]-tau_1[s,g]*X13[s,g-1]-mu_1[s,g]*X13[s,g-1]
      
      X14[s,g] <- X14[s,g-1]+gamma_2[s,g]*X10[s,g-1]+gamma_6[s,g]*X18[s,g-1]+gamma_10[s,g]*X22[s,g-1]+delta_1[s,g]*X13[s,g-1]-tau_2[s,g]*X14[s,g-1]-delta_2[s,g]*X14[s,g-1]-mu_2[s,g]*X14[s,g-1]
      
      X15[s,g] <- X15[s,g-1]+gamma_3[s,g]*X11[s,g-1]+gamma_7[s,g]*X19[s,g-1]+gamma_11[s,g]*X23[s,g-1]+delta_2[s,g]*X14[s,g-1]-tau_3[s,g]*X15[s,g-1]-delta_3[s,g]*X15[s,g-1]-mu_3[s,g]*X15[s,g-1]
      
      X16[s,g] <- X16[s,g-1]+gamma_4[s,g]*X12[s,g-1]+gamma_8[s,g]*X20[s,g-1]+gamma_12[s,g]*X24[s,g-1]+delta_3[s,g]*X15[s,g-1]-tau_4[s,g]*X16[s,g-1]-mu_4[s,g]*X16[s,g-1]
      
      
      
      # Difference equations for compartments X17 - X20 (individuals who are viral suppressed)
      
      X17[s,g] <- X17[s,g-1]+theta_1[s,g]*X9[s,g-1]+tau_1[s,g]*X13[s,g-1]-gamma_5[s,g]*X17[s,g-1]-psi_1[s,g]*X17[s,g-1]-mu_9[s,g]*X17[s,g-1]
      
      X18[s,g] <- X18[s,g-1]+theta_2[s,g]*X10[s,g-1]+tau_2[s,g]*X14[s,g-1]-gamma_6[s,g]*X18[s,g-1]-psi_2[s,g]*X18[s,g-1]-mu_10[s,g]*X18[s,g-1]
      
      X19[s,g] <- X19[s,g-1]+theta_3[s,g]*X11[s,g-1]+tau_3[s,g]*X15[s,g-1]-gamma_7[s,g]*X19[s,g-1]-psi_3[s,g]*X19[s,g-1]-mu_11[s,g]*X19[s,g-1]
      
      X20[s,g] <- X20[s,g-1]+theta_4[s,g]*X12[s,g-1]+tau_4[s,g]*X16[s,g-1]-gamma_8[s,g]*X20[s,g-1]-psi_4[s,g]*X20[s,g-1]-mu_12[s,g]*X20[s,g-1]
      
      
      # Difference equations compartments X21 - X24 (individuals who are not viral suppressed)
      
      X21[s,g] <- X21[s,g-1]+psi_1[s,g]*X17[s,g-1]-gamma_9[s,g]*X21[s,g-1]-mu_13[s,g]*X21[s,g-1]
      
      X22[s,g] <- X22[s,g-1]+psi_2[s,g]*X18[s,g-1]-gamma_10[s,g]*X22[s,g-1]-mu_14[s,g]*X22[s,g-1]
      
      X23[s,g] <- X23[s,g-1]+psi_3[s,g]*X19[s,g-1]-gamma_11[s,g]*X23[s,g-1]-mu_15[s,g]*X23[s,g-1]
      
      X24[s,g] <- X24[s,g-1]+psi_4[s,g]*X20[s,g-1]-gamma_12[s,g]*X24[s,g-1]-mu_16[s,g]*X24[s,g-1]
      
      # Difference equations for compartment D (individuals who have died)
      
      D[s,g] <- D[s,g-1]+mu_1[s,g]*(X1[s,g-1]+X5[s,g-1]+X13[s,g-1])+mu_2[s,g]*(X2[s,g-1]+X6[s,g-1]+X14[s,g-1])+mu_3[s,g]*(X3[s,g-1]+X7[s,g-1]+X15[s,g-1])+ mu_4[s,g]*(X4[s,g-1]+X8[s,g-1]+X16[s,g-1])+mu_5[s,g]*X9[s,g-1]+mu_6[s,g]*X10[s,g-1]+mu_7[s,g]*X11[s,g-1]+mu_8[s,g]*X12[s,g-1]+mu_9[s,g]*X17[s,g-1]+mu_10[s,g]*X18[s,g-1]+mu_11[s,g]*X19[s,g-1]+mu_12[s,g]*X20[s,g-1]+mu_13[s,g]*X21[s,g-1]+mu_14[s,g]*X22[s,g-1]+mu_15[s,g]*X23[s,g-1]+mu_16[s,g]*X24[s,g-1]
      
      
      # Total population of HIV infected individuals
      N[s,g] <- X1[s,g]+X2[s,g]+X3[s,g]+X4[s,g]+X5[s,g]+X6[s,g]+X7[s,g]+X8[s,g]+X9[s,g]+X10[s,g]+X11[s,g]+X12[s,g]+X13[s,g]+X14[s,g]+X15[s,g]+X16[s,g]+X17[s,g]+X18[s,g]+X19[s,g]+X20[s,g]+X21[s,g]+X22[s,g]+X23[s,g]+X24[s,g]
      
      # Total population of individuals diagnosed with HIV
      Diagnosed[s,g] <- X5[s,g]+X6[s,g]+X7[s,g]+X8[s,g]+X9[s,g]+X10[s,g]+X11[s,g]+X12[s,g]+X13[s,g]+X14[s,g]+X15[s,g]+X16[s,g]+X17[s,g]+X18[s,g]+X19[s,g]+X20[s,g]+X21[s,g]+X22[s,g]+X23[s,g]+X24[s,g]#+
      #X1[s,g]*alpha_1[s,g]+X2[s,g]*alpha_2[s,g]+X3[s,g]*alpha_3[s,g]+X4[s,g]*alpha_4[s,g]
      
      # Total population of individuals newly diagnosed
      New_Diagnosed[s,g] <- X1[s,g]*alpha_1[s,g]+X2[s,g]*alpha_2[s,g]+X3[s,g]*alpha_3[s,g]+X4[s,g]*alpha_4[s,g]
      
      # Total population of individuals on ART
      On_ART[s,g] <- X17[s,g]+X18[s,g]+X19[s,g]+X20[s,g]+X21[s,g]+X22[s,g]+X23[s,g]+X24[s,g]
      
      # Total population of individuals who are not viral suppressed
      Z[s,g] <- X1[s,g]+X2[s,g]+X3[s,g]+X4[s,g]+X5[s,g]+X6[s,g]+X7[s,g]+X8[s,g]+X9[s,g]+X10[s,g]+X11[s,g]+X12[s,g]+X13[s,g]+X14[s,g]+X15[s,g]+X16[s,g]+X21[s,g]+X22[s,g]+X23[s,g]+X24[s,g]
      
      # Total population of individuals who are viral suppressed
      V[s,g] <- X17[s,g]+X18[s,g]+X19[s,g]+X20[s,g]
      
      if (N[s,g]==0) {
        
        p_r[s,g] <- 1
        p_r_v[s,g] <- 0
        
      } else {
        
        # Prevalence of HIV in the population based on individuals who are not viral supprressed
        p_r[s,g] <- Z[s,g]/N[s,g]
        # Prevalence of HIV in the population based on individuals who are viral supprressed
        p_r_v[s,g] <- V[s,g]/N[s,g]  
      }
    }
  }


# Model calibration_targets&GOF
# Generate predicted outcomes
  for (g in (1:generations)) {
    N_F_total[g] <- 0
    N_M_total[g] <- 0
    N_young_total[g] <- N[1,g]+N[6,g]+N[11,g]+N[16,g]+N[21,g]
    N_F_young_total[g] <- N[1,g]+N[6,g]+N[11,g]
    N_M_young_total[g] <- N[16,g]+N[21,g]
    N_adult_total[g] <- sum(N[1:3,g]+N[6:8,g]+N[11:13,g]+N[16:18,g]+N[21:23,g])
    NI_total[g] <- 0
    NI_F_total[g] <- 0
    NI_M_total[g] <- 0
    NI_young_total[g] <-  NI[1,g]+NI[6,g]+NI[11,g]+NI[16,g]+NI[21,g]
    NI_F_young_total[g] <- NI[1,g]+NI[6,g]+NI[11,g]
    NI_M_young_total[g] <- NI[16,g]+NI[21,g]
    NI_adult_total[g] <- sum(NI[1:3,g]+NI[6:8,g]+NI[11:13,g]+NI[16:18,g]+NI[21:23,g])
    N_total[g] <- 0
    Z_total[g] <- 0
    V_total[g] <- 0
    S_total[g] <- 0
    S_F_total[g] <- 0
    S_M_total[g] <- 0
    S_F_young_total[g] <- S[1,g]+S[6,g]+S[11,g]
    S_M_young_total[g] <- S[16,g]+S[21,g]
    
    Diagnosed_total[g] <- 0
    Diagnosed_F_total[g] <- 0
    Diagnosed_M_total[g] <- 0
    Diagnosed_F_young_total[g] <- Diagnosed[1,g]+Diagnosed[6,g]+Diagnosed[11,g]
    Diagnosed_M_young_total[g] <- Diagnosed[16,g]+Diagnosed[21,g]
    
    
    New_Diagnosed_total[g] <- 0
    
    Suppressed_total[g] <- 0
    Suppressed_total_F[g] <- 0
    Suppressed_total_M[g] <- 0
    Suppressed_total_F_young[g] <- V[1,g]+V[6,g]+V[11,g]
    Suppressed_total_M_young[g] <- V[16,g]+V[21,g]
    
    On_ART_total[g] <- 0
    On_ART_total_F[g] <- 0
    On_ART_total_M[g] <- 0
    On_ART_total_F_young[g] <- On_ART[1,g]+On_ART[6,g]+On_ART[11,g]
    On_ART_total_M_young[g] <- On_ART[16,g]+On_ART[21,g]
    
    
    Linked_total[g] <- 0
    for(s in c(1:3, 6:8, 11:13)){ # for female subpopulations
      N_F_total[g] <- N_F_total[g]+N[s,g]
      NI_F_total[g] <-  NI_F_total[g]+NI[s,g]
      S_F_total[g] = S_F_total[g]+S[s,g]
      On_ART_total_F[g]=On_ART_total_F[g]+On_ART[s,g]
      Suppressed_total_F[g]=Suppressed_total_F[g]+V[s,g]
      Diagnosed_F_total[g] = Diagnosed_F_total[g]+Diagnosed[s,g]
    }
    for(s in c(16:18, 21:23)){ # for male subpopulations
      N_M_total[g] <- N_M_total[g]+N[s,g]
      NI_M_total[g] <- NI_M_total[g]+NI[s,g]
      S_M_total[g]=S_M_total[g]+S[s,g]
      On_ART_total_M[g]=On_ART_total_M[g]+On_ART[s,g]
      Suppressed_total_M [g]=Suppressed_total_M[g]+V[s,g] 
      Diagnosed_M_total[g] = Diagnosed_M_total[g]+Diagnosed[s,g]
    }
    
    for (s in c(1:3, 6:8, 11:13, 16:18, 21:23)){
      # Total population    
      
      if(is.na(NI[s,g])){
        NI[s,g]=0
      }
      NI_total[g]=NI_total[g]+NI[s,g]
      N_total[g]=N_total[g]+N[s,g]
      S_total[g]=S_total[g]+S[s,g]
      Diagnosed_total[g]=Diagnosed_total[g]+Diagnosed[s,g]
      New_Diagnosed_total[g]=New_Diagnosed_total[g]+New_Diagnosed[s,g]
      Suppressed_total[g]=Suppressed_total[g]+V[s,g]
      Linked_total[g]=Linked_total[g]+Linked[s,g]
      On_ART_total[g]=On_ART_total[g]+On_ART[s,g]
    }
    for (s in c(1:3, 6:8, 11:13, 16:18, 21:23)){
      On_ART_grand_total[g]=On_ART_grand_total[g]+On_ART[s,g]
    }
  }

N_F_total_all <- array(0, dim=c(generations))
N_M_total_all<- array(0, dim=c(generations))
N_adult_total_all<- array(0, dim=c(generations))
S_total_all<- array(0, dim=c(generations))
S_F_total_all<- array(0, dim=c(generations))
S_M_total_all<- array(0, dim=c(generations))
Suppressed_total_all<- array(0, dim=c(generations))
Suppressed_total_F_all<- array(0, dim=c(generations))
Suppressed_total_M_all<- array(0, dim=c(generations))
N_total_all<- array(0, dim=c(generations))
 
  for (g in (1:generations)) {
    N_F_total_all[g] <- 0
    N_M_total_all[g] <- 0
    N_total_all[g] <- 0
    S_total_all[g] <- 0
    S_F_total_all[g] <- 0
    S_M_total_all[g] <- 0
    
    Suppressed_total_all[g] <- 0
    Suppressed_total_F_all[g] <- 0
    Suppressed_total_M_all[g] <- 0
    
    for(s in c(1:15)){ # for female subpopulations
      N_F_total_all[g] <- N_F_total_all[g]+N[s,g]
      S_F_total_all[g] = S_F_total_all[g]+S[s,g]
      Suppressed_total_F_all[g]=Suppressed_total_F_all[g]+V[s,g]
    }
    for(s in c(16:20, 21:25)){ # for male subpopulations
      N_M_total_all[g] <- N_M_total_all[g]+N[s,g]
      S_M_total_all[g]=S_M_total_all[g]+S[s,g]
      Suppressed_total_M_all [g]=Suppressed_total_M_all[g]+V[s,g] 
    }
    
    for (s in c(1:25)){
      # Total population    
      N_total_all[g]=N_total_all[g]+N[s,g]
      S_total_all[g]=S_total_all[g]+S[s,g]
      Suppressed_total_all[g]=Suppressed_total_all[g]+V[s,g]
    }
  }


# Save the generated parameters
# save.image("Data_trial500.RData")
  # number on ART
  pd_n_art_total<- 0
  
  # percent diagnosed
  pd_p_diagnosed<-0
  pd_p_diagnosed_F<-0
  pd_p_diagnosed_M<-0
  pd_p_diagnosed_F_young<-0
  pd_p_diagnosed_M_young<-0  
  
  # Percent on Treatment
  pd_art_p<- 0
  pd_art_p_F <- 0
  pd_art_p_M<- 0
  pd_art_p_F_young<- 0
  pd_art_p_M_young<- 0
  
  # Percent virally suppressed
  pd_vs_total <- 0
  pd_vs_total_F <- 0
  pd_vs_total_M <- 0
  pd_vs_total_M_young <- 0
  pd_vs_total_F_young <- 0
  pd_vs_total_all <- 0
  pd_vs_total_F_all<- 0
  pd_vs_total_M_all <- 0
  
  
  # Percent virally suppressed on ART
  pd_vs_ART_total <- 0
  pd_vs_ART_total_F <- 0
  pd_vs_ART_total_M <- 0
  pd_vs_ART_total_M_young <- 0
  pd_vs_ART_total_F_young <- 0
  
  
  # HIV Prevalence
  pd_prev_total <- 0
  pd_prev_total_F <- 0
  pd_prev_total_M <- 0
  pd_prev_total_young <- 0
  pd_prev_total_F_young <- 0
  pd_prev_total_M_young <- 0
  pd_prev_total_all <- 0
  pd_prev_total_F_all <- 0
  pd_prev_total_M_all <- 0
  
  # HIV Incidence
  pd_incidence_total <- 0
  pd_incidence_F_total <- 0
  pd_incidence_M_total <- 0
  pd_incidence_F_young_total <- 0
  pd_incidence_M_young_total <- 0
  
  # Number within bounds 
  N_within <- 0
  
    Prevalence_total[] <- (N_total[])/(N_total[]+S_total[])
    Prevalence_total_y[] <- Prevalence_total[][seq(1,length(Prevalence_total[]),12)]
  for (y in (1:length(Prevalence_total_y[]))) {
    # Prevalence
    Target_prev[y]<- Target_all[y,"Prevalence"]
    LB<- Target_all[y,"Prevalence_LB"]
    UB<- Target_all[y,"Prevalence_UB"]
    if (is.na(Target_prev[y])||is.na(Prevalence_total_y[y]))
    {
      pd_prev_total <- pd_prev_total
    } else {
      pd_prev_total <- pd_prev_total+abs(Prevalence_total_y[y]-Target_prev[y])/(Target_prev[y])
      if (Prevalence_total_y[y]>=LB&Prevalence_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
    Prevalence_F_total[] <- (N_F_total[])/(N_F_total[]+S_F_total[])
    Prevalence_F_total_y[] <- Prevalence_F_total[][seq(1,length(Prevalence_F_total[]),12)]
  for (y in (1:length(Prevalence_total_y[]))) {
    # Prevalence Female 
    Target_prev_F[y]<- Target_all[y,"Prevalence_F"]
    LB<- Target_all[y,"Prevalence_LB_F"]
    UB<- Target_all[y,"Prevalence_UB_F"]
    if (is.na(Target_prev_F[y])||is.na(Prevalence_F_total_y[y]))
    {
      pd_prev_total_F <- pd_prev_total_F
    } else {
      pd_prev_total_F <- pd_prev_total_F+abs(Prevalence_F_total_y[y]-Target_prev_F[y])/(Target_prev_F[y])
      if (Prevalence_F_total_y[y]>=LB&Prevalence_F_total_y[y]<=UB){N_within = N_within+1}
      }
  }
  
    Prevalence_M_total[] <- (N_M_total[])/(N_M_total[]+S_M_total[])
    Prevalence_M_total_y[] <- Prevalence_M_total[][seq(1,length(Prevalence_M_total[]),12)]
  for (y in (1:length(Prevalence_total_y[]))) {
    # Prevalence male
    Target_prev_M[y]<- Target_all[y,"Prevalence_M"]
    LB<- Target_all[y,"Prevalence_LB_M"]
    UB<- Target_all[y,"Prevalence_UB_M"]
    if (is.na(Target_prev_M[y])||is.na(Prevalence_M_total_y[y]))
    {
      pd_prev_total_M <- pd_prev_total_M
    } else {
      pd_prev_total_M <- pd_prev_total_M+abs(Prevalence_M_total_y[y]-Target_prev_M[y])/(Target_prev_M[y])
      if (Prevalence_M_total_y[y]>=LB&Prevalence_M_total_y[y]<=UB){N_within = N_within+1}
    }
  }
    
    Prevalence_F_young_total[] <- (N_F_young_total[])/(N_F_young_total[]+S_F_young_total[])
    Prevalence_F_young_total_y[] <- Prevalence_F_young_total[][seq(1,length(Prevalence_F_young_total[]),12)]
    for (y in (1:length(Prevalence_total_y[]))) {
      # Prevalence Female 
      Target_prev_F_young[y]<- Target_all[y,"Prevalence_young_F"]
      LB<- Target_all[y,"Prevalence_LB_young_F"]
      UB<- Target_all[y,"Prevalence_UB_young_F"]
      if (is.na(Target_prev_F_young[y])||is.na(Prevalence_F_young_total_y[y]))
      {
        pd_prev_total_F_young <- pd_prev_total_F_young
      } else {
        pd_prev_total_F_young <- pd_prev_total_F_young+abs(Prevalence_F_young_total_y[y]-Target_prev_F_young[y])/(Target_prev_F_young[y])
        if (Prevalence_F_young_total_y[y]>=LB&Prevalence_F_young_total_y[y]<=UB){N_within = N_within+1}
      }
    }
    
    Prevalence_M_young_total[] <- (N_M_young_total[])/(N_M_young_total[]+S_M_young_total[])
    Prevalence_M_young_total_y[] <- Prevalence_M_young_total[][seq(1,length(Prevalence_M_young_total[]),12)]
    for (y in (1:length(Prevalence_total_y[]))) {
      # Prevalence male
      Target_prev_M_young[y]<- Target_all[y,"Prevalence_young_M"]
      LB<- Target_all[y,"Prevalence_LB_young_M"]
      UB<- Target_all[y,"Prevalence_UB_young_M"]
      if (is.na(Target_prev_M_young[y])||is.na(Prevalence_M_young_total_y[y]))
      {
        pd_prev_total_M_young <- pd_prev_total_M_young
      } else {
        pd_prev_total_M_young <- pd_prev_total_M_young+abs(Prevalence_M_young_total_y[y]-Target_prev_M_young[y])/(Target_prev_M_young[y])
        if (Prevalence_M_young_total_y[y]>=LB&Prevalence_M_young_total_y[y]<=UB){N_within = N_within+1}
      }
    }
  
    Prevalence_total[] <- (N_total_all[])/(N_total_all[]+S_total_all[])
    Prevalence_total_y[] <- Prevalence_total[][seq(1,length(Prevalence_total[]),12)]
  for (y in (1:length(Prevalence_total_y[]))) {
    # Prevalence
    Target_prev[y]<- Target_new[y,"Prevalence"]
    LB<- Target_new[y,"Prevalence_LB"]
    UB<- Target_new[y,"Prevalence_UB"]
    if (is.na(Target_prev[y])||is.na(Prevalence_total_y[y]))
    {
      pd_prev_total_all <- pd_prev_total_all
    } else {
      pd_prev_total_all <- pd_prev_total_all+abs(Prevalence_total_y[y]-Target_prev[y])/(Target_prev[y])
      if (Prevalence_total_y[y]>=LB&Prevalence_total_y[y]<=UB){N_within = N_within+1}
    }
  }
    
    Prevalence_F_total[] <- (N_F_total_all[])/(N_F_total_all[]+S_F_total_all[])
    Prevalence_F_total_y[] <- Prevalence_F_total[][seq(1,length(Prevalence_F_total[]),12)]
  for (y in (1:length(Prevalence_total_y[]))) {
    # Prevalence Female 
    Target_prev_F[y]<- Target_new[y,"Prevalence_F"]
    LB<- Target_new[y,"Prevalence_LB_F"]
    UB<- Target_new[y,"Prevalence_UB_F"]
    if (is.na(Target_prev_F[y])||is.na(Prevalence_F_total_y[y]))
    {
      pd_prev_total_F_all <- pd_prev_total_F_all
    } else {
      pd_prev_total_F_all <- pd_prev_total_F_all+abs(Prevalence_F_total_y[y]-Target_prev_F[y])/(Target_prev_F[y])
      if (Prevalence_F_total_y[y]>=LB&Prevalence_F_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
      Prevalence_M_total[] <- (N_M_total_all[])/(N_M_total_all[]+S_M_total_all[])
      Prevalence_M_total_y[] <- Prevalence_M_total[][seq(1,length(Prevalence_M_total[]),12)]
  for (y in (1:length(Prevalence_total_y[]))) {
    # Prevalence male
    Target_prev_M[y]<- Target_new[y,"Prevalence_M"]
    LB<- Target_new[y,"Prevalence_LB_M"]
    UB<- Target_new[y,"Prevalence_UB_M"]
    if (is.na(Target_prev_M[y])||is.na(Prevalence_M_total_y[y]))
    {
      pd_prev_total_M_all <- pd_prev_total_M_all
    } else {
      pd_prev_total_M_all <- pd_prev_total_M_all+abs(Prevalence_M_total_y[y]-Target_prev_M[y])/(Target_prev_M[y])
      if (Prevalence_M_total_y[y]>=LB&Prevalence_M_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  #Incidence
    NI_total_y[] <- rollapply(NI_total[],12,sum,by=12)
    S_total_y[] <- S_total[][seq(1, length(S_total[]), 12)]
    N_total_y[] <- N_total[][seq(1, length(N_total[]), 12)]
    Incidence_total_y[] <- NI_total_y[]/S_total_y[]
  for (y in (1:length(Incidence_total_y[]))) {
    Target_incidence[y]<- Target_all[y,"Incidence"]
    LB<- Target_all[y,"Incidence_LB"]
    UB<- Target_all[y,"Incidence_UB"]
    if (is.na(Target_incidence[y])||is.na(Incidence_total_y[y]))
    {
      pd_incidence_total <- pd_incidence_total
    } else {
      pd_incidence_total <- pd_incidence_total+abs(Incidence_total_y[y]-Target_incidence[y])/ (Target_incidence[y])
      if (Incidence_total_y[y]>=LB&Incidence_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
    # HIV incidence female
      NI_F_total_y[] <- rollapply(NI_F_total[],12,sum,by=12)
      S_F_total_y[] <- S_F_total[][seq(1, length(S_F_total[]), 12)]
      N_F_total_y[] <- N_F_total[][seq(1, length(N_F_total[]), 12)]
      Incidence_F_total_y[] <- NI_F_total_y[]/S_F_total_y[]
  for (y in (1:length(Incidence_F_total_y[]))) {
    #Incidence female
    Target_incidence_F[y]<- Target_all[y,"Incidence_F"]
    LB<- Target_all[y,"Incidence_F_LB"]
    UB<- Target_all[y,"Incidence_F_UB"]
    if (is.na(Target_incidence_F[y])||is.na(Incidence_F_total_y[y]))
    {
      pd_incidence_F_total <- pd_incidence_F_total
    } else {
      pd_incidence_F_total <- pd_incidence_F_total+abs(Incidence_F_total_y[y]-Target_incidence_F[y])/ (Target_incidence_F[y])
      if (Incidence_total_y[y]>=LB&Incidence_total_y[y]<=UB){N_within = N_within+1}
    }
  }
      
  # HIV incidence male
  NI_M_total_y[] <- rollapply(NI_M_total[],12,sum,by=12)
  S_M_total_y[] <- S_M_total[][seq(1, length(S_M_total[]), 12)]
  N_M_total_y[] <- N_M_total[][seq(1, length(N_M_total[]), 12)]
  Incidence_M_total_y[] <- NI_M_total_y[]/S_M_total_y[]
  for (y in (1:length(Incidence_M_total_y[]))) {
    #Incidence male
    Target_incidence_M[y]<- Target_all[y,"Incidence_M"]
    LB<- Target_all[y,"Incidence_M_LB"]
    UB<- Target_all[y,"Incidence_M_UB"]
    if (is.na(Target_incidence_M[y])||is.na(Incidence_M_total_y[y]))
    {
      pd_incidence_M_total <- pd_incidence_M_total
    } else {
      pd_incidence_M_total <- pd_incidence_M_total+abs(Incidence_M_total_y[y]-Target_incidence_M[y])/ (Target_incidence_M[y])
      if (Incidence_total_y[y]>=LB&Incidence_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  # HIV incidence young female
  NI_F_young_total_y[] <- rollapply(NI_F_young_total[],12,sum,by=12)
  S_F_young_total_y[] <- S_F_young_total[][seq(1, length(S_F_young_total[]), 12)]
  N_F_young_total_y[] <- N_F_young_total[][seq(1, length(N_F_young_total[]), 12)]
  Incidence_F_young_total_y[] <- NI_F_young_total_y[]/S_F_young_total_y[]
  for (y in (1:length(Incidence_F_young_total_y[]))) {
    #Incidence female
    Target_incidence_F_young[y]<- Target_all[y,"Incidence_F_young"]
    LB<- Target_all[y,"Incidence_F_young_LB"]
    UB<- Target_all[y,"Incidence_F_young_UB"]
    if (is.na(Target_incidence_F_young[y])||is.na(Incidence_F_young_total_y[y]))
    {
      pd_incidence_F_young_total <- pd_incidence_F_young_total
    } else {
      pd_incidence_F_young_total <- pd_incidence_F_young_total+abs(Incidence_F_young_total_y[y]-Target_incidence_F_young[y])/ (Target_incidence_F_young[y])
      if (Incidence_total_y[y]>=LB&Incidence_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  # HIV incidence young male
  NI_M_young_total_y[] <- rollapply(NI_M_young_total[],12,sum,by=12)
  S_M_young_total_y[] <- S_M_young_total[][seq(1, length(S_M_young_total[]), 12)]
  N_M_young_total_y[] <- N_M_young_total[][seq(1, length(N_M_young_total[]), 12)]
  Incidence_M_young_total_y[] <- NI_M_young_total_y[]/S_M_young_total_y[]
  for (y in (1:length(Incidence_M_young_total_y[]))) {
    #Incidence male
    Target_incidence_M_young[y]<- Target_all[y,"Incidence_M_young"]
    LB<- Target_all[y,"Incidence_M_young_LB"]
    UB<- Target_all[y,"Incidence_M_young_UB"]
    if (is.na(Target_incidence_M_young[y])||is.na(Incidence_M_young_total_y[y]))
    {
      pd_incidence_M_young_total <- pd_incidence_M_young_total
    } else {
      pd_incidence_M_young_total <- pd_incidence_M_young_total+abs(Incidence_M_young_total_y[y]-Target_incidence_M_young[y])/ (Target_incidence_M_young[y])
      if (Incidence_total_y[y]>=LB&Incidence_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  On_ART_grand_total_y[] <- On_ART_grand_total[][seq(1,length(On_ART_grand_total[]),12)]
  for (y in (1:length(P_suppressed_total_y[]))) {
    # number on ART
    Target_ART[y]<- Target_all[y,"N_ART"]
    if (is.na(Target_ART[y])||is.na(On_ART_grand_total_y[y]))
    {
      pd_n_art_total <- pd_n_art_total
    } else {
      pd_n_art_total <- pd_n_art_total+abs(On_ART_grand_total_y[y]-Target_ART[y])/ (Target_ART[y])
    }
  }
  
  # Percent diag
  Diagnosed_total_y[] <- Diagnosed_total[][seq(1, length(Diagnosed_total[]), 12)]
  Rate_diagnosed_total_y[] <- Diagnosed_total_y[]/N_total_y[]*100
  for (y in (1:length(P_suppressed_total_y[]))) {
    # percent on ART
    Target_R_HIV[y]<- Target_all[y,"Percent_diagnosed_total"]
    LB<- Target_all[y,"Percent_diagnosed_total_LB"]
    UB<- Target_all[y,"Percent_diagnosed_total_UB"]
    if (is.na(Target_R_HIV[y])||is.na(Rate_diagnosed_total_y[y]))
    {
      pd_p_diagnosed <- pd_p_diagnosed
    } else {
      pd_p_diagnosed <- pd_p_diagnosed+abs(Rate_diagnosed_total_y[y]-Target_R_HIV[y])/ (Target_R_HIV[y])
      if (Rate_diagnosed_total_y[y]>=LB&Rate_diagnosed_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  
  
  # Percent on ART female
  Diagnosed_F_total_y[] <- Diagnosed_F_total[][seq(1, length(Diagnosed_F_total[]), 12)]
  Rate_diagnosed_F_total_y[] <- Diagnosed_F_total_y[]/N_F_total_y[]*100
  for (y in (1:length(P_suppressed_total_y[]))) {
    # percent on ART female
    Target_R_HIV[y]<- Target_all[y,"Percent_diagnosed_Female"]
    LB<- Target_all[y,"Percent_diagnosed_Female_LB"]
    UB<- Target_all[y,"Percent_diagnosed_Female_UB"]
    if (is.na(Target_R_HIV[y])||is.na(Rate_diagnosed_F_total_y[y]))
    {
      pd_p_diagnosed_F <- pd_p_diagnosed_F
    } else {
      pd_p_diagnosed_F <- pd_p_diagnosed_F+abs(Rate_diagnosed_F_total_y[y]-Target_R_HIV[y])/ (Target_R_HIV[y])
      if (Rate_diagnosed_F_total_y[y]>=LB&Rate_diagnosed_F_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  # Percent on ART male
  Diagnosed_M_total_y[] <- Diagnosed_M_total[][seq(1, length(Diagnosed_M_total[]), 12)]
  Rate_diagnosed_M_total_y[] <- Diagnosed_M_total_y[]/N_M_total_y[]*100
  for (y in (1:length(P_suppressed_total_y[]))) {
    # percent on ART male
    Target_R_HIV[y]<- Target_all[y,"Percent_diagnosed_Male"]
    LB<- Target_all[y,"Percent_diagnosed_Male_LB"]
    UB<- Target_all[y,"Percent_diagnosed_Male_UB"]
    if (is.na(Target_R_HIV[y])||is.na(Rate_diagnosed_M_total_y[y]))
    {
      pd_p_diagnosed_M <- pd_p_diagnosed_M
    } else {
      pd_p_diagnosed_M <- pd_p_diagnosed_M+abs(Rate_diagnosed_M_total_y[y]-Target_R_HIV[y])/ (Target_R_HIV[y])
      if (Rate_diagnosed_M_total_y[y]>=LB&Rate_diagnosed_M_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  # Percent diag young female
  Diagnosed_F_young_total_y[] <- Diagnosed_F_young_total[][seq(1, length(Diagnosed_F_young_total[]), 12)]
  Rate_diagnosed_F_young_total_y[] <- Diagnosed_F_young_total_y[]/N_F_young_total_y[]*100
  for (y in (1:length(P_suppressed_total_y[]))) {
    # percent on ART female
    Target_R_HIV[y]<- Target_all[y,"Percent_diagnosed_Female_YOUNG"]
    LB<- Target_all[y,"Percent_diagnosed_Female_YOUNG_LB"]
    UB<- Target_all[y,"Percent_diagnosed_Female_YOUNG_UB"]
    if (is.na(Target_R_HIV[y])||is.na(Rate_diagnosed_F_young_total_y[y]))
    {
      pd_p_diagnosed_F_young <- pd_p_diagnosed_F_young
    } else {
      pd_p_diagnosed_F_young <- pd_p_diagnosed_F_young+abs(Rate_diagnosed_F_young_total_y[y]-Target_R_HIV[y])/ (Target_R_HIV[y])
      if (Rate_diagnosed_F_young_total_y[y]>=LB&Rate_diagnosed_F_young_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  # Percent on ART 
  Diagnosed_total_y[] <- Diagnosed_total[][seq(1, length(Diagnosed_total[]), 12)]
  On_ART_p_y[] <- On_ART_total[][seq(1,length(On_ART_total[]),12)]/Diagnosed_total_y[]*100
  for (y in (1:length(P_suppressed_total_y[]))) {
    # percent on ART
    P_Target_ART[y]<- Target_new[y,"P_ART_All"]
    LB<- Target_new[y,"P_ART_All_LB"]
    UB<- Target_new[y,"P_ART_All_UB"]
    if (is.na(P_Target_ART[y])||is.na(On_ART_p_y[y]))
    {
      pd_art_p <- pd_art_p
    } else {
      pd_art_p <- pd_art_p+abs(On_ART_p_y[y]-P_Target_ART[y])/ (P_Target_ART[y])
      if (On_ART_p_y[y]>=LB&On_ART_p_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  
  
  # Percent on ART female
  Diagnosed_F_total_y[] <- Diagnosed_F_total[][seq(1, length(Diagnosed_F_total[]), 12)]
   On_ART_p_F_y[] <- On_ART_total_F[][seq(1,length(On_ART_total_F[]),12)]/Diagnosed_F_total_y[]*100
  for (y in (1:length(P_suppressed_total_y[]))) {
    # percent on ART female
    P_Target_ART[y]<- Target_new[y,"P_ART_All_F"]
    LB<- Target_new[y,"P_ART_All_F_LB"]
    UB<- Target_new[y,"P_ART_All_F_UB"]
    if (is.na(P_Target_ART[y])||is.na(On_ART_p_F_y[y]))
    {
      pd_art_p_F <- pd_art_p_F
    } else {
      pd_art_p_F <- pd_art_p_F+abs(On_ART_p_F_y[y]-P_Target_ART[y])/ (P_Target_ART[y])
      if (On_ART_p_F_y[y]>=LB&On_ART_p_F_y[y]<=UB){N_within = N_within+1}
    }
  }
  
   # Percent on ART male
   Diagnosed_M_total_y[] <- Diagnosed_M_total[][seq(1, length(Diagnosed_M_total[]), 12)]
   On_ART_p_M_y[] <- On_ART_total_M[][seq(1,length(On_ART_total_M[]),12)]/Diagnosed_M_total_y[]*100
  for (y in (1:length(P_suppressed_total_y[]))) {
    # percent on ART male
    P_Target_ART[y]<- Target_new[y,"P_ART_All_M"]
    LB<- Target_new[y,"P_ART_All_M_LB"]
    UB<- Target_new[y,"P_ART_All_M_UB"]
    if (is.na(P_Target_ART[y])||is.na(On_ART_p_M_y[y]))
    {
      pd_art_p_M <- pd_art_p_M
    } else {
      pd_art_p_M <- pd_art_p_M+abs(On_ART_p_M_y[y]-P_Target_ART[y])/ (P_Target_ART[y])
      if (On_ART_p_M_y[y]>=LB&On_ART_p_M_y[y]<=UB){N_within = N_within+1}
    }
  }
  
   # Percent on ART female
   Diagnosed_F_young_total_y[] <- Diagnosed_F_young_total[][seq(1, length(Diagnosed_F_young_total[]), 12)]
   On_ART_p_F_young_y[] <- On_ART_total_F_young[][seq(1,length(On_ART_total_F_young[]),12)]/Diagnosed_F_young_total_y[]*100
   for (y in (1:length(P_suppressed_total_y[]))) {
     # percent on ART female
     P_Target_ART[y]<- Target_new[y,"P_ART_All_F_young"]
     LB<- Target_new[y,"P_ART_All_F_young_LB"]
     UB<- Target_new[y,"P_ART_All_F_young_UB"]
     if (is.na(P_Target_ART[y])||is.na(On_ART_p_F_young_y[y]))
     {
       pd_art_p_F_young <- pd_art_p_F_young
     } else {
       pd_art_p_F_young <- pd_art_p_F_young+abs(On_ART_p_F_young_y[y]-P_Target_ART[y])/ (P_Target_ART[y])
       if (On_ART_p_F_young_y[y]>=LB&On_ART_p_F_young_y[y]<=UB){N_within = N_within+1}
     }
   }
   
   # Percent viral suppression (among ART)
   On_ART_total_y[] <- On_ART_total[][seq(1,length(On_ART_total[]),12)]
   Suppressed_total_y[] <- Suppressed_total[][seq(1,length(Suppressed_total[]),12)]
   P_suppressed_art_total_y[] <- Suppressed_total_y[]/On_ART_total_y[]*100
  for (y in (1:length(P_suppressed_art_total_y[]))) {
    # % VS among linked
    Target_VS_art[y]<- Target_all[y,"Percent_VS_on_ART_total"]
    LB<- Target_all[y,"Percent_VS_on_ART_total_LB"]
    UB<- Target_all[y,"Percent_VS_on_ART_total_UB"]
    if (is.na(Target_VS_art[y])||is.na(P_suppressed_art_total_y[y]))
    {
      pd_vs_ART_total <- pd_vs_ART_total
    } else {
      pd_vs_ART_total <- pd_vs_ART_total+abs(P_suppressed_art_total_y[y]-Target_VS_art[y])/ (Target_VS_art[y])
      if (P_suppressed_art_total_y[y]>=LB&P_suppressed_art_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
   # Percent viral suppression female (among ART)
   On_ART_total_F_y[] <- On_ART_total_F[][seq(1,length(On_ART_total_F[]),12)]
   Suppressed_total_F_y[] <- Suppressed_total_F[][seq(1,length(Suppressed_total_F[]),12)]
   P_suppressed_art_F_total_y[] <- Suppressed_total_F_y[]/On_ART_total_F_y[]*100
  for (y in (1:length(P_suppressed_art_total_y[]))) {
    # % VS among linked females
    Target_VS_art_F[y]<- Target_all[y,"Percent_VS_on_ART_Female"]
    LB<- Target_all[y,"Percent_VS_on_ART_Female_LB"]
    UB<- Target_all[y,"Percent_VS_on_ART_Female_UB"]
    if (is.na(Target_VS_art_F[y])||is.na(P_suppressed_art_F_total_y[y]))
    {
      pd_vs_ART_total_F <- pd_vs_ART_total_F
    } else {
      pd_vs_ART_total_F <- pd_vs_ART_total_F+abs(P_suppressed_art_F_total_y[y]-Target_VS_art_F[y])/ (Target_VS_art_F[y])
      if (P_suppressed_art_F_total_y[y]>=LB&P_suppressed_art_F_total_y[y]<=UB){N_within = N_within+1}
    }
  }
   
   # Percent viral suppression male (among ART)
   On_ART_total_M_y[] <- On_ART_total_M[][seq(1,length(On_ART_total_M[]),12)]
   Suppressed_total_M_y[] <- Suppressed_total_M[][seq(1,length(Suppressed_total_M[]),12)]
   P_suppressed_art_M_total_y[] <- Suppressed_total_M_y[]/On_ART_total_M_y[]*100
  for (y in (1:length(P_suppressed_art_total_y[]))) {
    # % VS among linked males
    Target_VS_art_M[y]<- Target_all[y,"Percent_VS_on_ART_Male"]
    LB<- Target_all[y,"Percent_VS_on_ART_Male_LB"]
    UB<- Target_all[y,"Percent_VS_on_ART_Male_UB"]
    if (is.na(Target_VS_art_M[y])||is.na(P_suppressed_art_M_total_y[y]))
    {
      pd_vs_ART_total_M <- pd_vs_ART_total_M
    } else {
      pd_vs_ART_total_M <- pd_vs_ART_total_M+abs(P_suppressed_art_M_total_y[y]-Target_VS_art_M[y])/ (Target_VS_art_M[y])
      if (P_suppressed_art_M_total_y[y]>=LB&P_suppressed_art_M_total_y[y]<=UB){N_within = N_within+1}
    }
  }
   
   # Percent viral suppression female (among ART)
   On_ART_total_F_young_y[] <- On_ART_total_F_young[][seq(1,length(On_ART_total_F_young[]),12)]
   Suppressed_total_F_young_y[] <- Suppressed_total_F_young[][seq(1,length(Suppressed_total_F_young[]),12)]
   P_suppressed_art_F_young_total_y[] <- Suppressed_total_F_young_y[]/On_ART_total_F_young_y[]*100
   for (y in (1:length(P_suppressed_art_total_y[]))) {
     # % VS among linked females
     Target_VS_art_F_young[y]<- Target_all[y,"Percent_VS_on_ART_Female_YOUNG"]
     LB<- Target_all[y,"Percent_VS_on_ART_Female_YOUNG_LB"]
     UB<- Target_all[y,"Percent_VS_on_ART_Female_YOUNG_UB"]
     if (is.na(Target_VS_art_F_young[y])||is.na(P_suppressed_art_F_young_total_y[y]))
     {
       pd_vs_ART_total_F_young <- pd_vs_ART_total_F_young
     } else {
       pd_vs_ART_total_F_young <- pd_vs_ART_total_F_young+abs(P_suppressed_art_F_young_total_y[y]-Target_VS_art_F_young[y])/ (Target_VS_art_F_young[y])
       if (P_suppressed_art_F_young_total_y[y]>=LB&P_suppressed_art_F_young_total_y[y]<=UB){N_within = N_within+1}
     }
   }
   
   # Percent viral suppression (among ALL)
   Suppressed_total_y[] <- Suppressed_total[][seq(1,length(Suppressed_total[]),12)]
   P_suppressed_total_y[] <- Suppressed_total_y[]/N_total_y[]*100
  for (y in (1:length(P_suppressed_art_total_y[]))) {
    # % VS linked
    Target_VS[y]<- Target_all[y,"Percent_VS_total"]
    LB<- Target_all[y,"Percent_VS_total_LB"]
    UB<- Target_all[y,"Percent_VS_total_UB"]
    if (is.na(Target_VS[y])||is.na(P_suppressed_total_y[y]))
    {
      pd_vs_total <- pd_vs_total
    } else {
      pd_vs_total <- pd_vs_total+abs(P_suppressed_total_y[y]-Target_VS[y])/ (Target_VS[y])
      if (P_suppressed_total_y[y]>=LB&P_suppressed_total_y[y]<=UB){N_within = N_within+1}
    }
  }
   
   # Percent viral suppression female
   Suppressed_total_F_y[] <- Suppressed_total_F[][seq(1,length(Suppressed_total_F[]),12)]
   P_suppressed_F_total_y[] <- Suppressed_total_F_y[]/N_F_total_y[]*100
  for (y in (1:length(P_suppressed_art_total_y[]))) {
    # % VS among females
    Target_VS_F[y]<- Target_all[y,"Percent_VS_Female"]
    LB<- Target_all[y,"Percent_VS_Female_LB"]
    UB<- Target_all[y,"Percent_VS_Female_UB"]
    if (is.na(Target_VS_F[y])||is.na(P_suppressed_F_total_y[y]))
    {
      pd_vs_total_F <- pd_vs_total_F
    } else {
      pd_vs_total_F <- pd_vs_total_F+abs(P_suppressed_F_total_y[y]-Target_VS_F[y])/ (Target_VS_F[y])
      if (P_suppressed_F_total_y[y]>=LB&P_suppressed_F_total_y[y]<=UB){N_within = N_within+1}
    }
  }
   
   # Percent viral suppression male
   Suppressed_total_M_y[] <- Suppressed_total_M[][seq(1,length(Suppressed_total_M[]),12)]
   P_suppressed_M_total_y[] <- Suppressed_total_M_y[]/N_M_total_y[]*100
  for (y in (1:length(P_suppressed_art_total_y[]))) {
    # % VS among males
    Target_VS_M[y]<- Target_all[y,"Percent_VS_Male"]
    LB<- Target_all[y,"Percent_VS_Male_LB"]
    UB<- Target_all[y,"Percent_VS_Male_UB"]
    if (is.na(Target_VS_M[y])||is.na(P_suppressed_M_total_y[y]))
    {
      pd_vs_total_M <- pd_vs_total_M
    } else {
      pd_vs_total_M <- pd_vs_total_M+abs(P_suppressed_M_total_y[y]-Target_VS_M[y])/ (Target_VS_M[y])
      if (P_suppressed_M_total_y[y]>=LB&P_suppressed_M_total_y[y]<=UB){N_within = N_within+1}
    }
  }
   
   Suppressed_total_F_young_y[] <- Suppressed_total_F_young[][seq(1,length(Suppressed_total_F_young[]),12)]
   P_suppressed_F_young_total_y[] <- Suppressed_total_F_young_y[]/N_F_young_total_y[]*100
   for (y in (1:length(P_suppressed_art_total_y[]))) {
     # % VS among females
     Target_VS_F_young[y]<- Target_all[y,"Percent_VS_Female_YOUNG"]
     LB<- Target_all[y,"Percent_VS_Female_YOUNG_LB"]
     UB<- Target_all[y,"Percent_VS_Female_YOUNG_UB"]
     if (is.na(Target_VS_F_young[y])||is.na(P_suppressed_F_young_total_y[y]))
     {
       pd_vs_total_F_young <- pd_vs_total_F_young
     } else {
       pd_vs_total_F_young <- pd_vs_total_F_young+abs(P_suppressed_F_young_total_y[y]-Target_VS_F_young[y])/ (Target_VS_F_young[y])
       if (P_suppressed_F_young_total_y[y]>=LB&P_suppressed_F_young_total_y[y]<=UB){N_within = N_within+1}
     }
   }
   
   # Percent viral suppression (among ALL)
   N_total_y[] <- N_total_all[][seq(1, length(N_total_all[]), 12)]
   Suppressed_total_y[] <- Suppressed_total_all[][seq(1,length(Suppressed_total_all[]),12)]
   P_suppressed_total_y[] <- Suppressed_total_y[]/N_total_y[]*100
  for (y in (1:length(P_suppressed_art_total_y[]))) {
    # % VS linked
    Target_VS[y]<- Target_new[y,"Percent_VS_total"]
    LB<- Target_new[y,"Percent_VS_total_LB"]
    UB<- Target_new[y,"Percent_VS_total_UB"]
    if (is.na(Target_VS[y])||is.na(P_suppressed_total_y[y]))
    {
      pd_vs_total_all <- pd_vs_total_all
    } else {
      pd_vs_total_all <- pd_vs_total_all+abs(P_suppressed_total_y[y]-Target_VS[y])/ (Target_VS[y])
      if (P_suppressed_total_y[y]>=LB&P_suppressed_total_y[y]<=UB){N_within = N_within+1}
    }
  }
   
   N_F_total_y[] <- N_F_total_all[][seq(1, length(N_F_total_all[]), 12)]
   Suppressed_total_F_y[] <- Suppressed_total_F_all[][seq(1,length(Suppressed_total_F_all[]),12)]
   P_suppressed_F_total_y[] <- Suppressed_total_F_y[]/N_F_total_y[]*100
  for (y in (1:length(P_suppressed_art_total_y[]))) {
    # % VS among females
    Target_VS_F[y]<- Target_new[y,"Percent_VS_Female"]
    LB<- Target_new[y,"Percent_VS_Female_LB"]
    UB<- Target_new[y,"Percent_VS_Female_UB"]
    if (is.na(Target_VS_F[y])||is.na(P_suppressed_F_total_y[y]))
    {
      pd_vs_total_F_all <- pd_vs_total_F_all
    } else {
      pd_vs_total_F_all <- pd_vs_total_F_all+abs(P_suppressed_F_total_y[y]-Target_VS_F[y])/ (Target_VS_F[y])
      if (P_suppressed_F_total_y[y]>=LB&P_suppressed_F_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  N_M_total_y[] <- N_M_total_all[][seq(1, length(N_M_total_all[]), 12)]
  Suppressed_total_M_y[] <- Suppressed_total_M_all[][seq(1,length(Suppressed_total_M_all[]),12)]
  P_suppressed_M_total_y[] <- Suppressed_total_M_y[]/N_M_total_y[]*100
  for (y in (1:length(P_suppressed_art_total_y[]))) {
    # % VS among males
    Target_VS_M[y]<- Target_new[y,"Percent_VS_Male"]
    LB<- Target_new[y,"Percent_VS_Male_LB"]
    UB<- Target_new[y,"Percent_VS_Male_UB"]
    if (is.na(Target_VS_M[y])||is.na(P_suppressed_M_total_y[y]))
    {
      pd_vs_total_M_all <- pd_vs_total_M_all
    } else {
      pd_vs_total_M_all <- pd_vs_total_M_all+abs(P_suppressed_M_total_y[y]-Target_VS_M[y])/ (Target_VS_M[y])
      if (P_suppressed_M_total_y[y]>=LB&P_suppressed_M_total_y[y]<=UB){N_within = N_within+1}
    }
  }
  
  pd_total= pd_n_art_total+
    pd_art_p + pd_art_p_F + pd_art_p_M +
    pd_vs_total + pd_vs_total_F + pd_vs_total_M +
    pd_vs_ART_total + pd_vs_ART_total_F + pd_vs_ART_total_M+
    pd_prev_total + pd_prev_total_F + pd_prev_total_M+
    pd_incidence_total+pd_incidence_F_total+pd_incidence_M_total+
    pd_prev_total_all + pd_prev_total_F_all + pd_prev_total_M_all+
    pd_vs_total_all + pd_vs_total_F_all + pd_vs_total_M_all+
    pd_prev_total_F_young+pd_prev_total_M_young+
    pd_incidence_F_young_total+pd_incidence_M_young_total+
    #pd_art_p_F_young+pd_vs_total_F_young+pd_vs_ART_total_F_young+
    pd_p_diagnosed+pd_p_diagnosed_F+pd_p_diagnosed_M #+pd_p_diagnosed_F_young
  
    pd_total_new <- c(t,pd_total,N_within)
    
    # for (i in 1:50){
    #   if (i==1&N_within>pd_top50[i,3]){
    #       pd_top50 <- rbind(pd_total_new,pd_top50[(1:50), ])
    #       pd_top50 <- head(pd_top50,-1)
    #     break}
    #   if (i==1&N_within==pd_top50[i,3]&pd_total<pd_top50[i,2]){
    #     pd_top50 <- rbind(pd_total_new,pd_top50[(1:50), ])
    #     pd_top50 <- head(pd_top50,-1)
    #     break}
    #   if (i>1&N_within>pd_top50[i,3]){
    #     pd_top50 <- rbind(pd_top50[1:i-1, ], pd_total_new,pd_top50[-(1:i-1), ])
    #     pd_top50 <- head(pd_top50,-1)
    #     break}
    #   if (i>1&N_within==pd_top50[i,3]&pd_total<=pd_top50[i,2]){
    #     pd_top50 <- rbind(pd_top50[1:i-1, ], pd_total_new,pd_top50[-(1:i-1), ])
    #     pd_top50 <- head(pd_top50,-1)
    #     break}
    # }}

for (i in 1:50){
  if (i==1&pd_total<pd_top50[i,2]){
    pd_top50 <- rbind(pd_total_new,pd_top50[(1:50), ])
    pd_top50 <- head(pd_top50,-1)
    break}
  if (i>1&pd_total<=pd_top50[i,2]){
    pd_top50 <- rbind(pd_top50[1:i-1, ], pd_total_new,pd_top50[-(1:i-1), ])
    pd_top50 <- head(pd_top50,-1)
    break}
}}


write.csv(pd_top50,"PD_total_0506.csv")


