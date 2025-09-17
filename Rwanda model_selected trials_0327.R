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
trials <- 20000
trial_select <- 50
# Specify the number of sub-groups we are going to run
sgroup <-25
# Command specifing the number of monthly time periods in the model (432 monthly time periods from 2012 - 2020)
# Note that the generation shall be divisible by 12. Otherwise, errors would occur in calibration
generations <- 324
deathmulti<-1
initialvf<-0.25
# Read the file that stores model parameter inputs
# source("~/Dropbox/VCU_PhD_Year 3&4/Paper 1/R code/v1.0/Import_data_v1.0.R")
# Read the file that stores model parameter inputs
# Set the pathway to read importing files
# Please change the pathway to shared drive when using the VCU computer
# If error mesage "incomplete final line" occurs, open csv file with text edit and add a new empty line
# setwd("/gpfs_fs/home/panz/US model/")
# setwd("~/Dropbox/VCU_PhD_Year 3&4/Paper 1/R code/Data/")
#setwd("A:/Common/Personnel work/Rasnick/Rwanda Model/Rwanda Calibration Model 1.0")
#setwd("C:/Users/rrasn/Desktop/Rwanda Calibration Model 1.0")
setwd("~/Dropbox/VCU_PhD_Year 3&4/GRA/Rwanda-Model/Full model/Data/")

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

pd_top50 <- read.csv(file='PD_total_0506.csv')
list_trial <- pd_top50$trial

# Set up empty variables to store the value for the parameter inputs 
# Set up empty variables to store the value for the parameter inputs 
# Note that the initial value of lambda shall set to 1 for separate rates otherwise error may occur (since we are multiplying the values)
lambda_inf <- array(1,dim=c(trial_select, sgroup, generations))
lambda_t_r <- array(0,dim=c(trial_select, sgroup, generations))

# Sets up variables to store the value for each compartment
# The initial value is set as null.
S <- array(NA,dim=c(trial_select, sgroup, generations))
S_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
S_total <- array(0,dim=c(trial_select, 1, generations))
S_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

S_F <- array(NA,dim=c(trial_select, sgroup, generations))
S_F_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
S_F_total <- array(0,dim=c(trial_select, 1, generations))
S_F_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

S_M <- array(NA,dim=c(trial_select, sgroup, generations))
S_M_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
S_M_total <- array(0,dim=c(trial_select, 1, generations))
S_M_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

S_F_young <- array(NA,dim=c(trial_select, sgroup, generations))
S_F_young_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
S_F_young_total <- array(0,dim=c(trial_select, 1, generations))
S_F_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

S_M_young <- array(NA,dim=c(trial_select, sgroup, generations))
S_M_young_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
S_M_young_total <- array(0,dim=c(trial_select, 1, generations))
S_M_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

D <- array(0,dim=c(trial_select, sgroup, generations))
X1 <- array(NA,dim=c(trial_select, sgroup, generations))
X2 <- array(NA,dim=c(trial_select, sgroup, generations))
X3 <- array(NA,dim=c(trial_select, sgroup, generations))
X4 <- array(NA,dim=c(trial_select, sgroup, generations))
X5 <- array(NA,dim=c(trial_select, sgroup, generations))
X6 <- array(NA,dim=c(trial_select, sgroup, generations))
X7 <- array(NA,dim=c(trial_select, sgroup, generations))
X8 <- array(NA,dim=c(trial_select, sgroup, generations))
X9 <- array(NA,dim=c(trial_select, sgroup, generations))
X10 <- array(NA,dim=c(trial_select, sgroup, generations))
X11 <- array(NA,dim=c(trial_select, sgroup, generations))
X12 <- array(NA,dim=c(trial_select, sgroup, generations))
X13 <- array(NA,dim=c(trial_select, sgroup, generations))
X14 <- array(NA,dim=c(trial_select, sgroup, generations))
X15 <- array(NA,dim=c(trial_select, sgroup, generations))
X16 <- array(NA,dim=c(trial_select, sgroup, generations))
X17 <- array(NA,dim=c(trial_select, sgroup, generations))
X18 <- array(NA,dim=c(trial_select, sgroup, generations))
X19 <- array(NA,dim=c(trial_select, sgroup, generations))
X20 <- array(NA,dim=c(trial_select, sgroup, generations))
X21 <- array(NA,dim=c(trial_select, sgroup, generations))
X22 <- array(NA,dim=c(trial_select, sgroup, generations))
X23 <- array(NA,dim=c(trial_select, sgroup, generations))
X24 <- array(NA,dim=c(trial_select, sgroup, generations))


N_lr_u_f <- array(NA,dim=c(trial_select, sgroup, generations))
N_hr_f <- array(NA,dim=c(trial_select, sgroup, generations))
N_lr_r_f <- array(NA,dim=c(trial_select, sgroup, generations))
N_u_m <- array(NA,dim=c(trial_select, sgroup, generations))
N_r_m <- array(NA,dim=c(trial_select, sgroup, generations))

# This command sets up variables to store the value of model parameters.
# Rate at which people enter the susceptible population compartment (population growth).
rho <- array(NA,dim=c(trial_select, sgroup, generations))

# Rate of natural history HIV disease pregression.
delta_1 <- array(NA,dim=c(trial_select, sgroup, generations))
delta_2 <- array(NA,dim=c(trial_select, sgroup, generations))
delta_3 <- array(NA,dim=c(trial_select, sgroup, generations))
delta_4 <- array(NA,dim=c(trial_select, sgroup, generations))
delta_5 <- array(NA,dim=c(trial_select, sgroup, generations))
delta_6 <- array(NA,dim=c(trial_select, sgroup, generations))
m_delta <- array(NA,dim=c(trial_select, sgroup, generations))

# Rate of getting diagnosed with HIV.
alpha_1 <- array(NA,dim=c(trial_select, sgroup, generations))
alpha_2 <- array(NA,dim=c(trial_select, sgroup, generations))
alpha_3 <- array(NA,dim=c(trial_select, sgroup, generations))
alpha_4 <- array(NA,dim=c(trial_select, sgroup, generations))

# Rate of being linked.
sigma_1 <- array(NA,dim=c(trial_select, sgroup, generations))
sigma_2 <- array(NA,dim=c(trial_select, sgroup, generations))
sigma_3 <- array(NA,dim=c(trial_select, sgroup, generations))
sigma_4 <- array(NA,dim=c(trial_select, sgroup, generations))
# sigma_1_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# sigma_2_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# sigma_3_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# sigma_4_ipv <- array(NA,dim=c(trial_select, sgroup, generations))

# Rate of getting lost from care.
gamma_1 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_2 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_3 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_4 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_5 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_6 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_7 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_8 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_9 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_10 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_11 <- array(NA,dim=c(trial_select, sgroup, generations))
gamma_12 <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_1_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_2_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_3_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_4_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_5_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_6_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_7_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_8_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_9_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_10_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_11_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# gamma_12_ipv <- array(NA,dim=c(trial_select, sgroup, generations))

# Rate of failure to maintain viral suppression.
psi_1 <- array(NA,dim=c(trial_select, sgroup, generations))
psi_2 <- array(NA,dim=c(trial_select, sgroup, generations))
psi_3 <- array(NA,dim=c(trial_select, sgroup, generations))
psi_4 <- array(NA,dim=c(trial_select, sgroup, generations))

# Rate of being on ART and suppressed.
theta_1 <- array(NA,dim=c(trial_select, sgroup, generations))
theta_2 <- array(NA,dim=c(trial_select, sgroup, generations))
theta_3 <- array(NA,dim=c(trial_select, sgroup, generations))
theta_4 <- array(NA,dim=c(trial_select, sgroup, generations))
# theta_1_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# theta_2_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# theta_3_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# theta_4_ipv <- array(NA,dim=c(trial_select, sgroup, generations))

# Rate of return to ART and being suppressed.
tau_1 <- array(NA,dim=c(trial_select, sgroup, generations))
tau_2 <- array(NA,dim=c(trial_select, sgroup, generations))
tau_3 <- array(NA,dim=c(trial_select, sgroup, generations))
tau_4 <- array(NA,dim=c(trial_select, sgroup, generations))

# Mortality rate.
mu_0 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_1 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_2 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_3 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_4 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_5 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_6 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_7 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_8 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_9 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_10 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_11 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_12 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_13 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_14 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_15 <- array(NA,dim=c(trial_select, sgroup, generations))
mu_16 <- array(NA,dim=c(trial_select, sgroup, generations))
m_mu <- array(NA,dim=c(trial_select, sgroup, generations))

# Percentage of individuas consistently using a condom in sub-group r.
c_r <- array(NA,dim=c(trial_select, sgroup, generations))
#c_r_ipv <- array(NA,dim=c(trial_select, sgroup, generations))
# Reduction in probability of HIV transmission when using a condom. 
epsilon <- array(NA,dim=c(trial_select, sgroup, generations))
# Weight appleid for condom use in sub-group r.
omega_r <- array(NA,dim=c(trial_select, sgroup, generations))
# Reduction in probability of HIV transmission when on ART and viral suppressed.
upsilon <- array(NA,dim=c(trial_select, sgroup, generations))
# Weight applied for viral suppreessed individuals
kappa <- array(NA,dim=c(trial_select, sgroup, generations))
# probability of HIV transmission per unprotected sex contact when not virally suppressed.
beta <- array(NA,dim=c(trial_select, sgroup, generations))
# Average number of sexual acts per time period in sub-group r
n_r <- array(NA,dim=c(trial_select, sgroup, generations))

# Total number of individuals in sub-group r
Total <- array(NA,dim=c(trial_select, sgroup, generations))
initial_diag <- array(NA,dim=c(trial_select, sgroup, generations))
initial_link <- array(NA,dim=c(trial_select, sgroup, generations))
initial_ltfu <- array(NA,dim=c(trial_select, sgroup, generations))
initial_artvs <- array(NA,dim=c(trial_select, sgroup, generations))
initial_artvf <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_unlink_1 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_unlink_2 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_unlink_3 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_unlink_4 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_link_1 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_link_2 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_link_3 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_link_4 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_artvs_1 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_artvs_2 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_artvs_3 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_artvs_4 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_artvf_1 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_artvf_2 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_artvf_3 <- array(NA,dim=c(trial_select, sgroup, generations))
initial_cd4_artvf_4 <- array(NA,dim=c(trial_select, sgroup, generations))

# Set up variables to store model outputs
# Set up variables to store the number of individuals infected.
N <- array(0,dim=c(trial_select, sgroup, generations))
N_y <- array(0,dim=c(trial_select, sgroup, generations/12))
N_total <- array(0,dim=c(trial_select, 1, generations))
N_total_y <- matrix(0,nrow=trial_select, ncol=generations/12)
N_F <- array(0,dim=c(trial_select, sgroup, generations))
N_F_y <- array(0,dim=c(trial_select, sgroup, generations/12))
N_F_total <- array(0,dim=c(trial_select, 1, generations))
N_F_total_y <- matrix(0,nrow=trial_select, ncol=generations/12)
N_M <- array(0,dim=c(trial_select, sgroup, generations))
N_M_y <- array(0,dim=c(trial_select, sgroup, generations/12))
N_M_total <- array(0,dim=c(trial_select, 1, generations))
N_M_total_y <- matrix(0,nrow=trial_select, ncol=generations/12)
N_young <- array(0,dim=c(trial_select, sgroup, generations))
N_young_y <- array(0,dim=c(trial_select, sgroup, generations/12))
N_young_total <- array(0,dim=c(trial_select, 1, generations))
N_young_total_y <- matrix(0,nrow=trial_select, ncol=generations/12)
N_F_young_total <- array(0,dim=c(trial_select, 1, generations))
N_F_young_total_y <- matrix(0,nrow=trial_select, ncol=generations/12)
N_M_young_total <- array(0,dim=c(trial_select, 1, generations))
N_M_young_total_y <- matrix(0,nrow=trial_select, ncol=generations/12)
N_adult <- array(0,dim=c(trial_select, sgroup, generations))
N_adult_y <- array(0,dim=c(trial_select, sgroup, generations/12))
N_adult_total <- array(0,dim=c(trial_select, 1, generations))
N_adult_total_y <- matrix(0,nrow=trial_select, ncol=generations/12)

S_F_total<- array(0,dim=c(trial_select, 1, generations))
S_M_total<- array(0,dim=c(trial_select, 1, generations))
S_young_total<- array(0,dim=c(trial_select, 1, generations))
S_F_young_total<- array(0,dim=c(trial_select, 1, generations))
S_M_young_total<- array(0,dim=c(trial_select, 1, generations))

#This command sets up variables to store the number of individuals virally not suppresed.
Z <- array(0,dim=c(trial_select, sgroup, generations))
Z_total <- array(0,dim=c(trial_select, 1, generations))
Z_total_y <- matrix(0,nrow=trial_select, ncol=generations/12)
#This command sets up variables to store the number of individuals virally suppressed.
V <- array(0,dim=c(trial_select, sgroup, generations))
V_total <- array(0,dim=c(trial_select, 1, generations))
V_total_y <- matrix(0,nrow=trial_select, ncol=generations/12)

# This command sets up variables to store the percent of individuals virally not suppressed.
p_r <- array(0,dim=c(trial_select, sgroup, generations))

# This command sets up variables to store the percent of individuals virally suppressed.
p_r_v <- array(0,dim=c(trial_select, sgroup, generations))

# This command sets up variables to store HIV treatment engagement status in the cohort. 
Undiagnosed <- array(NA,dim=c(trial_select, sgroup, generations))
Undiagnosed_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Undiagnosed_total <- array(0,dim=c(trial_select, 1, generations))
Undiagnosed_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

Diagnosed <- array(NA,dim=c(trial_select, sgroup, generations))
Diagnosed_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Diagnosed_total <- array(0,dim=c(trial_select, 1, generations))
Diagnosed_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Diagnosed_F_total <- array(0,dim=c(trial_select, 1, generations))
Diagnosed_F_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Diagnosed_M_total <- array(0,dim=c(trial_select, 1, generations))
Diagnosed_M_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Diagnosed_F_young_total <- array(0,dim=c(trial_select, 1, generations))
Diagnosed_F_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Diagnosed_M_young_total <- array(0,dim=c(trial_select, 1, generations))
Diagnosed_M_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Diagnosed_total_py <- matrix(NA,nrow=trial_select, ncol=generations/12)

Rate_diagnosed <- array(NA,dim=c(trial_select, sgroup, generations))
Rate_diagnosed_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Rate_diagnosed_total <- array(0,dim=c(trial_select, 1, generations))
Rate_diagnosed_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Rate_diagnosed_F_total <- array(0,dim=c(trial_select, 1, generations))
Rate_diagnosed_F_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Rate_diagnosed_M_total <- array(0,dim=c(trial_select, 1, generations))
Rate_diagnosed_M_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Rate_diagnosed_F_young_total <- array(0,dim=c(trial_select, 1, generations))
Rate_diagnosed_F_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Rate_diagnosed_M_young_total <- array(0,dim=c(trial_select, 1, generations))
Rate_diagnosed_M_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

New_Diagnosed <- array(NA,dim=c(trial_select, sgroup, generations))
New_Diagnosed_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
New_Diagnosed_total <- array(0,dim=c(trial_select, 1, generations))
New_Diagnosed_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

New_Diagnosed_MSM_total <- array(0,dim=c(trial_select, 1, generations))
New_Diagnosed_MSM_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
New_Diagnosed_MSMIDU_total <- array(0,dim=c(trial_select, 1, generations))
New_Diagnosed_MSMIDU_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
New_Diagnosed_IDUM_total <- array(0,dim=c(trial_select, 1, generations))
New_Diagnosed_IDUM_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
New_Diagnosed_IDUF_total <- array(0,dim=c(trial_select, 1, generations))
New_Diagnosed_IDUF_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
New_Diagnosed_HETM_total <- array(0,dim=c(trial_select, 1, generations))
New_Diagnosed_HETM_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
New_Diagnosed_HETF_total <- array(0,dim=c(trial_select, 1, generations))
New_Diagnosed_HETF_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

Rate_new_Diagnosed <- array(NA,dim=c(trial_select, sgroup, generations))
Rate_new_Diagnosed_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Rate_new_Diagnosed_total <- array(0,dim=c(trial_select, 1, generations))
Rate_new_Diagnosed_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

Linked <- array(NA,dim=c(trial_select, sgroup, generations))
Linked_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Linked_total <- array(0,dim=c(trial_select, 1, generations))
Linked_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

Lost <- array(NA,dim=c(trial_select, sgroup, generations))
Lost_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Lost_total <- array(0,dim=c(trial_select, 1, generations))
Lost_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

Suppressed <- array(NA,dim=c(trial_select, sgroup, generations))
Suppressed_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Suppressed_total <- array(0,dim=c(trial_select, 1, generations))
Suppressed_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Suppressed_total_F <- array(0,dim=c(trial_select, 1, generations))
Suppressed_total_F_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Suppressed_total_M <- array(0,dim=c(trial_select, 1, generations))
Suppressed_total_M_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Suppressed_total_F_young <- array(0,dim=c(trial_select, 1, generations))
Suppressed_total_F_young_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Suppressed_total_M_young <- array(0,dim=c(trial_select, 1, generations))
Suppressed_total_M_young_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_ART <- matrix(NA,nrow=trial_select, ncol=generations/12)
P_Target_ART <- matrix(NA,nrow=trial_select, ncol=generations/12)

P_suppressed <- array(NA,dim=c(trial_select, sgroup, generations))
P_suppressed_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
P_suppressed_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
P_suppressed_F_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_F_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
P_suppressed_M_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_M_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
P_suppressed_F_young_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_F_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
P_suppressed_M_young_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_M_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

P_suppressed_art <- array(NA,dim=c(trial_select, sgroup, generations))
P_suppressed_art_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
P_suppressed_art_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_art_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
P_suppressed_art_F_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_art_F_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
P_suppressed_art_M_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_art_M_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
P_suppressed_art_F_young_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_art_F_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
P_suppressed_art_M_young_total <- array(0,dim=c(trial_select, 1, generations))
P_suppressed_art_M_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)


Suppressed_onART_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Suppressed_onART_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VSonART <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VS_art_F<- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VS_art_M<- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VS_art_F_young<- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VS_art_M_young<- matrix(NA,nrow=trial_select, ncol=generations/12)

Target_VS <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VS_F<- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VS_M<- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VS_F_young<- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VS_M_young<- matrix(NA,nrow=trial_select, ncol=generations/12)

Not_Suppressed <- array(NA,dim=c(trial_select, sgroup, generations))
Not_Suppressed_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Not_Suppressed_total <- array(0,dim=c(trial_select, 1, generations))
Not_Suppressed_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

On_ART <- array(NA,dim=c(trial_select, sgroup, generations))
On_ART_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
On_ART_p <- array(0,dim=c(trial_select, 1, generations))
On_ART_total <- array(0,dim=c(trial_select, 1, generations))
On_ART_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
On_ART_total_F <- array(0,dim=c(trial_select, 1, generations))
On_ART_total_F_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
On_ART_total_M <- array(0,dim=c(trial_select, 1, generations))
On_ART_total_M_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
On_ART_total_F_young <- array(0,dim=c(trial_select, 1, generations))
On_ART_total_F_young_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
On_ART_total_M_young <- array(0,dim=c(trial_select, 1, generations))
On_ART_total_M_young_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
On_ART_grand_total<- array(0,dim=c(trial_select, 1, generations))
On_ART_grand_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

On_ART_p_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
On_ART_total_ny <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_onART_N <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_onART_P <- matrix(NA,nrow=trial_select, ncol=generations/12)

On_ART_p_F <- array(0,dim=c(trial_select, 1, generations))
On_ART_p_F_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
On_ART_p_M <- array(0,dim=c(trial_select, 1, generations))
On_ART_p_M_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

On_ART_p_F_young <- array(0,dim=c(trial_select, 1, generations))
On_ART_p_F_young_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
On_ART_p_M_young <- array(0,dim=c(trial_select, 1, generations))
On_ART_p_M_young_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

# This command sets up variables to store HIV incidence
Incidence <- array(NA,dim=c(trial_select, sgroup, generations))
Incidence_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Incidence_total <- array(0,dim=c(trial_select, 1, generations))
Incidence_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_incidence <- matrix(NA,nrow=trial_select, ncol=generations/12)
Incidence_F <- array(NA,dim=c(trial_select, sgroup, generations))
Incidence_F_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Incidence_F_total <- array(0,dim=c(trial_select, 1, generations))
Incidence_F_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_incidence_F <- matrix(NA,nrow=trial_select, ncol=generations/12)
Incidence_M <- array(NA,dim=c(trial_select, sgroup, generations))
Incidence_M_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Incidence_M_total <- array(0,dim=c(trial_select, 1, generations))
Incidence_M_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_incidence_M <- matrix(NA,nrow=trial_select, ncol=generations/12)
Incidence_F_young <- array(NA,dim=c(trial_select, sgroup, generations))
Incidence_F_young_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Incidence_F_young_total <- array(0,dim=c(trial_select, 1, generations))
Incidence_F_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_incidence_F_young  <- matrix(NA,nrow=trial_select, ncol=generations/12)
Incidence_M_young  <- array(NA,dim=c(trial_select, sgroup, generations))
Incidence_M_young_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Incidence_M_young_total <- array(0,dim=c(trial_select, 1, generations))
Incidence_M_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_incidence_M_young  <- matrix(NA,nrow=trial_select, ncol=generations/12)

# This command sets up variables to store HIV prevalence.
Prevalence <- array(NA,dim=c(trial_select, sgroup, generations))
Prevalence_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
Prevalence_total <- array(0,dim=c(trial_select, 1, generations))
Prevalence_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_F <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_M <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_young <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_F_young <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_M_young <- matrix(NA,nrow=trial_select, ncol=generations/12)

DHS_Prevalence <- array(NA,dim=c(trial_select, sgroup, generations))
DHS_Prevalence_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
DHS_Prevalence_total <- array(0,dim=c(trial_select, 1, generations))
DHS_Prevalence_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
DHS_Target_prev <- matrix(NA,nrow=trial_select, ncol=generations/12)

DHS_Prevalence_F <- array(NA,dim=c(trial_select, sgroup, generations))
DHS_Prevalence_y_F<- array(NA,dim=c(trial_select, sgroup, generations/12))
DHS_Prevalence_F_total <- array(0,dim=c(trial_select, 1, generations))
DHS_Prevalence_F_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
DHS_Target_prev_F <- matrix(NA,nrow=trial_select, ncol=generations/12)

DHS_Prevalence_M <- array(NA,dim=c(trial_select, sgroup, generations))
DHS_Prevalence_y_M <- array(NA,dim=c(trial_select, sgroup, generations/12))
DHS_Prevalence_M_total <- array(0,dim=c(trial_select, 1, generations))
DHS_Prevalence_M_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
DHS_Target_prev_M <- matrix(NA,nrow=trial_select, ncol=generations/12)

DHS_Prevalence_F_young <- array(NA,dim=c(trial_select, sgroup, generations))
DHS_Prevalence_y_F_young <- array(NA,dim=c(trial_select, sgroup, generations/12))
DHS_Prevalence_F_total_young <- array(0,dim=c(trial_select, 1, generations))
DHS_Prevalence_F_total_y_young <- matrix(NA,nrow=trial_select, ncol=generations/12)
DHS_Target_prev_F_young <- matrix(NA,nrow=trial_select, ncol=generations/12)

DHS_Prevalence_M_young <- array(NA,dim=c(trial_select, sgroup, generations))
DHS_Prevalence_y_M_young <- array(NA,dim=c(trial_select, sgroup, generations/12))
DHS_Prevalence_M_total_young <- array(0,dim=c(trial_select, 1, generations))
DHS_Prevalence_M_total_y_young <- matrix(NA,nrow=trial_select, ncol=generations/12)
DHS_Target_prev_M_young <- matrix(NA,nrow=trial_select, ncol=generations/12)






Prevalence_F_total <- array(0,dim=c(trial_select, 1, generations))
Prevalence_F_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_F <- matrix(NA,nrow=trial_select, ncol=generations/12)


Prevalence_M_total <- array(0,dim=c(trial_select, 1, generations))
Prevalence_M_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_M <- matrix(NA,nrow=trial_select, ncol=generations/12)

Prevalence_young_total <- array(0,dim=c(trial_select, 1, generations))
Prevalence_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_young <- matrix(NA,nrow=trial_select, ncol=generations/12)

Prevalence_F_young_total <- array(0,dim=c(trial_select, 1, generations))
Prevalence_F_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_F_young <- matrix(NA,nrow=trial_select, ncol=generations/12)


Prevalence_M_young_total <- array(0,dim=c(trial_select, 1, generations))
Prevalence_M_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_prev_M_young <- matrix(NA,nrow=trial_select, ncol=generations/12)

# This command sets up variables to store the number of new infection.
NI <- array(NA,dim=c(trial_select, sgroup, generations))
NI_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
NI_total <- array(0,dim=c(trial_select, 1, generations))
NI_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
NI_F <- array(NA,dim=c(trial_select, sgroup, generations))
NI_F_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
NI_F_total <- array(0,dim=c(trial_select, 1, generations))
NI_F_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
NI_M <- array(NA,dim=c(trial_select, sgroup, generations))
NI_M_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
NI_M_total <- array(0,dim=c(trial_select, 1, generations))
NI_M_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
NI_young <- array(NA,dim=c(trial_select, sgroup, generations))
NI_young_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
NI_young_total <- array(0,dim=c(trial_select, 1, generations))
NI_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
NI_M_young <- array(NA,dim=c(trial_select, sgroup, generations))
NI_M_young_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
NI_M_young_total <- array(0,dim=c(trial_select, 1, generations))
NI_M_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
NI_F_young <- array(NA,dim=c(trial_select, sgroup, generations))
NI_F_young_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
NI_F_young_total <- array(0,dim=c(trial_select, 1, generations))
NI_F_young_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)
NI_adult <- array(NA,dim=c(trial_select, sgroup, generations))
NI_adult_y <- array(NA,dim=c(trial_select, sgroup, generations/12))
NI_adult_total <- array(0,dim=c(trial_select, 1, generations))
NI_adult_total_y <- matrix(NA,nrow=trial_select, ncol=generations/12)

Target_N <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_N_F <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_N_M <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_N_young <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_N_adult <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_NI <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_NI_F <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_NI_M <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_NI_young <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_NI_adult <- matrix(NA,nrow=trial_select, ncol=generations/12)


# This command sets up variables to store total unprotected sex for heterosexual female.
N_het_f_u<- array(0,dim=c(trial_select, sgroup, generations))
N_het_f_r<- array(0,dim=c(trial_select, sgroup, generations))
# This command sets up variables to store total unprotected sex for heterosexual male.
N_het_m_u<- array(0,dim=c(trial_select, sgroup, generations))
N_het_m_r<- array(0,dim=c(trial_select, sgroup, generations))

# Model calibration goodness of fit parameters 
# number on ART
pd_n_art_total<- numeric(trial_select)

# percent diagnosed
pd_p_diagnosed<-numeric(trial_select)
pd_p_diagnosed_F<-numeric(trial_select)
pd_p_diagnosed_M<-numeric(trial_select)
pd_p_diagnosed_F_young<-numeric(trial_select)
pd_p_diagnosed_M_young<-numeric(trial_select)  

# Percent on Treatment
pd_art_p <- numeric(trial_select)
pd_art_p_F <- numeric(trial_select)
pd_art_p_M <- numeric(trial_select)
pd_art_p_F_young <- numeric(trial_select)
pd_art_p_M_young <- numeric(trial_select)

# Percent virally suppressed
pd_vs_total <- numeric(trial_select)
pd_vs_total_F <- numeric(trial_select)
pd_vs_total_M <- numeric(trial_select)
pd_vs_total_M_young <- numeric(trial_select)
pd_vs_total_F_young <- numeric(trial_select)


# Percent virally suppressed on ART
pd_vs_ART_total <- numeric(trial_select)
pd_vs_ART_total_F <- numeric(trial_select)
pd_vs_ART_total_M <- numeric(trial_select)
pd_vs_ART_total_M_young <- numeric(trial_select)
pd_vs_ART_total_F_young <- numeric(trial_select)


# HIV Prevalence
pd_prev_total <- numeric(trial_select)
pd_prev_total_F <- numeric(trial_select)
pd_prev_total_M <- numeric(trial_select)
pd_prev_total_young <- numeric(trial_select)
pd_prev_total_F_young <- numeric(trial_select)
pd_prev_total_M_young <- numeric(trial_select)
pd_incidence_total <- numeric(trial_select)
pd_incidence_F_total <- numeric(trial_select)
pd_incidence_M_total <- numeric(trial_select)
pd_incidence_F_young_total <- numeric(trial_select)
pd_incidence_M_young_total <- numeric(trial_select)

Target_N_diag <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_R_diag <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_N_HIV <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_R_HIV_F <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_R_HIV_M <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_R_HIV_F_young <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_R_HIV_M_young <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_R_HIV <- matrix(NA,nrow=trial_select, ncol=generations/12)
Target_VS_art <- matrix(NA,nrow=trial_select, ncol=generations/12)

pd_total=numeric(trial_select)

#Initialize the parameter inputs by random draw
# Initial population
# Susceptible population
i=1
for (t in 1:trials) {
  # Total population for each sub-population
  TotalUrbanWomen=sum(Pop[1:10,"Total"])
  TotalHighRisk = runif(1,min=11000, max=45000)
  TotalLowRiskUrban=TotalUrbanWomen-TotalHighRisk
  if (t %in% list_trial ){
    # High risk women
    Total[i,6,1]=TotalHighRisk*(157+470)/1338 #15-24
    Total[i,7,1]=TotalHighRisk*(342+288/2)/1338 #25-34
    Total[i,8,1]=TotalHighRisk*(288/2+81/5)/1338 #35-44
    Total[i,9,1]=TotalHighRisk*(81/5+81/5)/1338 #45-54
    Total[i,10,1]=TotalHighRisk*(81/5+81/5)/1338 #55-64
    # low risk urban women
    women=c(535205,460559,327619,256033,201099,209695,154203,124472,89338,67466)
    totalwomen=sum(women)
    Total[i,1,1]=TotalUrbanWomen*(sum(535205,460559)/totalwomen)-Total[i,6,1] #15-24
    Total[i,2,1]=TotalUrbanWomen*(sum(327619,256033)/totalwomen)-Total[i,7,1] #25-34
    Total[i,3,1]=TotalUrbanWomen*(sum(201099,209695)/totalwomen)-Total[i,8,1] #35-44
    Total[i,4,1]=TotalUrbanWomen*(sum(154203,124472)/totalwomen)-Total[i,9,1]#45-54
    Total[i,5,1]=TotalUrbanWomen*(sum(89338,67466)/totalwomen)-Total[i,10,1] #55-64
    
    for (s in 11:sgroup) {
      Total[i,s,1]=Pop[s,"Total"]
    }}
  # HIV prevalence
  for (l in 1:25) {
    a=Pop_HIVPrev[l,"Alpha"]
    b=Pop_HIVPrev[l,"Beta"]
    prev_ini=rbeta(1,a,b)
    if (t %in% list_trial ){Prevalence[i,l,1]=prev_ini}
  }
  # Susceptible population
  if (t %in% list_trial ){
    for (s in 1:sgroup) {S[i,s,1]=Total[i,s,1]*(1-Prevalence[i,s,1])}}
  
  # By infectious compartment
  # Diagnosis
  # By age
  multi_diag=rtri(1,0.15,0.8,0.3)
  
  a=Pop_comp_diag[1,"Alpha"]
  b=Pop_comp_diag[1,"Beta"]
  #mode=Pop_comp_diag[1,"Mode"]
  initial_diag_s=rbeta(1,a,b)*multi_diag
  
  # Apply to all sub-population
  if (t %in% list_trial ){
    for (s in 1:sgroup) {initial_diag[i,s,1]=initial_diag_s}}
  
  
  # Linkage
  # By age
  a=Pop_comp_link[1,"Alpha"]
  b=Pop_comp_link[1,"Beta"]
  initial_link_s=rbeta(1,a,b)
  
  # Apply to all sub-population
  if (t %in% list_trial ){
    for (s in 1:sgroup) {initial_link[i,s,1]=initial_link_s}}
  
  # LTFU (No one is LTFU)
  
  # On ART and suppressed
  # By age
  a=Pop_comp_artvs[1,"Alpha"]
  b=Pop_comp_artvs[1,"Beta"]
  initial_artvs_s=rbeta(1,a,b)
  
  # Apply to all sub-population
  if (t %in% list_trial ){
    for (s in 1:sgroup) {initial_artvs[i,s,1]=initial_artvs_s}}
  
  # On ART and not suppressed
  # By age
  a=Pop_comp_artvf[1,"Alpha"]
  b=Pop_comp_artvf[1,"Beta"]
  initialvf_s=rbeta(1,a,b)
  
  # Apply to all sub-population
  if (t %in% list_trial ){
    for (s in 1:sgroup) {initial_artvf[i,s,1]=initialvf_s}}
  
  # Distribution by CD4
  # Not linked, CD4>500
  #a=Pop_CD4[1,"Alpha"]
  #b=Pop_CD4[1,"Beta"]
  if (t %in% list_trial){
    initial_cd4_unlink_1[i,1,1]=Pop_CD4[1,"Mean"]
    # Not linked, CD4 350-500
    #a=Pop_CD4[2,"Alpha"]
    #b=Pop_CD4[2,"Beta"]
    initial_cd4_unlink_2[i,1,1]=Pop_CD4[2,"Mean"]
    # Not linked, CD4 200-350
    #a=Pop_CD4[3,"Alpha"]
    #b=Pop_CD4[3,"Beta"]
    initial_cd4_unlink_3[i,1,1]=Pop_CD4[3,"Mean"]
    # Not linked, CD4 <200
    initial_cd4_unlink_4[i,1,1]=1-initial_cd4_unlink_1[i,1,1]-initial_cd4_unlink_2[i,1,1]-initial_cd4_unlink_3[i,1,1]
    # Linked, CD4>500
    #a=Pop_CD4[4,"Alpha"]
    #b=Pop_CD4[4,"Beta"]
    initial_cd4_link_1[i,1,1]=Pop_CD4[4,"Mean"]
    # Linked, CD4 350-500
    #a=Pop_CD4[5,"Alpha"]
    #b=Pop_CD4[5,"Beta"]
    initial_cd4_link_2[i,1,1]=Pop_CD4[5,"Mean"]
    # Linked, CD4 200-350
    a=Pop_CD4[6,"Alpha"]
    b=Pop_CD4[6,"Beta"]
    initial_cd4_link_3[i,1,1]=Pop_CD4[6,"Mean"]
    # Linked, CD4 <200
    initial_cd4_link_4[i,1,1]=1-initial_cd4_link_1[i,1,1]-initial_cd4_link_2[i,1,1]-initial_cd4_link_3[i,1,1]
    # ARTVS, CD4>500
    #a=Pop_CD4[7,"Alpha"]
    #b=Pop_CD4[7,"Beta"]
    initial_cd4_artvs_1[i,1,1]=Pop_CD4[7,"Mean"]
    # ARTVS, CD4 350-500
    #a=Pop_CD4[8,"Alpha"]
    #b=Pop_CD4[8,"Beta"]
    initial_cd4_artvs_2[i,1,1]=Pop_CD4[8,"Mean"]
    # ARTVS, CD4 200-350
    #a=Pop_CD4[9,"Alpha"]
    #b=Pop_CD4[9,"Beta"]
    initial_cd4_artvs_3[i,1,1]=Pop_CD4[9,"Mean"]
    # ARTVS, CD4 <200
    initial_cd4_artvs_4[i,1,1]=1-initial_cd4_artvs_1[i,1,1]-initial_cd4_artvs_2[i,1,1]-initial_cd4_artvs_3[i,1,1]
    # ARTVF, CD4>500
    # a=Pop_CD4[7,"Alpha"]
    #  b=Pop_CD4[7,"Beta"]
    initial_cd4_artvf_1[i,1,1]=Pop_CD4[10,"Mean"]
    # ARTVF, CD4 350-500
    # a=Pop_CD4[8,"Alpha"]
    #  b=Pop_CD4[8,"Beta"]
    initial_cd4_artvf_2[i,1,1]=Pop_CD4[11,"Mean"]
    # ARTVF, CD4 200-350
    # a=Pop_CD4[9,"Alpha"]
    #  b=Pop_CD4[9,"Beta"]
    initial_cd4_artvf_3[i,1,1]=Pop_CD4[12,"Mean"]
    # ARTVF, CD4 <200
    initial_cd4_artvf_4[i,1,1]=1-initial_cd4_artvf_1[i,1,1]-initial_cd4_artvf_2[i,1,1]-initial_cd4_artvf_3[i,1,1]
    
    # Apply to all sub-population
    for (s in 1:sgroup){
      initial_cd4_unlink_1[i,s,1]=initial_cd4_unlink_1[i,1,1]
      initial_cd4_unlink_2[i,s,1]=initial_cd4_unlink_2[i,1,1]
      initial_cd4_unlink_3[i,s,1]=initial_cd4_unlink_3[i,1,1]
      initial_cd4_unlink_4[i,s,1]=initial_cd4_unlink_4[i,1,1]
      initial_cd4_link_1[i,s,1]=initial_cd4_link_1[i,1,1]
      initial_cd4_link_2[i,s,1]=initial_cd4_link_2[i,1,1]
      initial_cd4_link_3[i,s,1]=initial_cd4_link_3[i,1,1]
      initial_cd4_link_4[i,s,1]=initial_cd4_link_4[i,1,1]
      initial_cd4_artvs_1[i,s,1]=initial_cd4_artvs_1[i,1,1]
      initial_cd4_artvs_2[i,s,1]=initial_cd4_artvs_2[i,1,1]
      initial_cd4_artvs_3[i,s,1]=initial_cd4_artvs_3[i,1,1]
      initial_cd4_artvs_4[i,s,1]=initial_cd4_artvs_4[i,1,1]
      initial_cd4_artvf_1[i,s,1]=initial_cd4_artvf_1[i,1,1]
      initial_cd4_artvf_2[i,s,1]=initial_cd4_artvf_2[i,1,1]
      initial_cd4_artvf_3[i,s,1]=initial_cd4_artvf_3[i,1,1]
      initial_cd4_artvf_4[i,s,1]=initial_cd4_artvf_4[i,1,1]
    }
  }
  
  # Infectious compartment
  if (t %in% list_trial){
    for (s in 1:sgroup){
      X1[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*(1-initial_diag[i,s,1])*initial_cd4_unlink_1[i,s,1]
      X2[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*(1-initial_diag[i,s,1])*initial_cd4_unlink_2[i,s,1]
      X3[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*(1-initial_diag[i,s,1])*initial_cd4_unlink_3[i,s,1]
      X4[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*(1-initial_diag[i,s,1])*initial_cd4_unlink_4[i,s,1]
      X5[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*(1-initial_link[i,s,1])*initial_cd4_link_1[i,s,1]
      X6[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*(1-initial_link[i,s,1])*initial_cd4_link_2[i,s,1]
      X7[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*(1-initial_link[i,s,1])*initial_cd4_link_3[i,s,1]
      X8[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*(1-initial_link[i,s,1])*initial_cd4_link_4[i,s,1]
      X9[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*(1-initial_artvs[i,s,1])*initial_cd4_link_1[i,s,1]
      X10[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*(1-initial_artvs[i,s,1])*initial_cd4_link_2[i,s,1]
      X11[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*(1-initial_artvs[i,s,1])*initial_cd4_link_3[i,s,1]
      X12[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*(1-initial_artvs[i,s,1])*initial_cd4_link_4[i,s,1]
      X13[i,s,1]=0
      X14[i,s,1]=0
      X15[i,s,1]=0
      X16[i,s,1]=0
      X17[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*initial_artvs[i,s,1]*(1-initial_artvf[i,s,1])*initial_cd4_artvs_1[i,s,1]
      X18[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*initial_artvs[i,s,1]*(1-initial_artvf[i,s,1])*initial_cd4_artvs_2[i,s,1]
      X19[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*initial_artvs[i,s,1]*(1-initial_artvf[i,s,1])*initial_cd4_artvs_3[i,s,1]
      X20[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*initial_artvs[i,s,1]*(1-initial_artvf[i,s,1])*initial_cd4_artvs_4[i,s,1]
      X21[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*initial_artvs[i,s,1]*initial_artvf[i,s,1]*initial_cd4_artvf_1[i,s,1]
      X22[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*initial_artvs[i,s,1]*initial_artvf[i,s,1]*initial_cd4_artvf_2[i,s,1]
      X23[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*initial_artvs[i,s,1]*initial_artvf[i,s,1]*initial_cd4_artvf_3[i,s,1]
      X24[i,s,1]=Total[i,s,1]*Prevalence[i,s,1]*initial_diag[i,s,1]*initial_link[i,s,1]*initial_artvs[i,s,1]*initial_artvf[i,s,1]*initial_cd4_artvf_4[i,s,1]
    }}
  
  # Probability of HIV transmission (heterosexual, converted to monthly probability)
  # FSW
  min=p_transmission[1,"Min"]
  max=p_transmission[1,"Max"]
  mode=p_transmission[1,"Mode"]
  # random draw
  beta_fsw=rtri(1,min,max,mode)
  # Heterosexual female
  min=p_transmission[2,"Min"]
  max=p_transmission[2,"Max"]
  mode=p_transmission[2,"Mode"]
  # random draw
  beta_hetf=rtri(1,min,max,mode)
  # Heterosexual male
  min=p_transmission[3,"Min"]
  max=p_transmission[3,"Max"]
  mode=p_transmission[3,"Mode"]
  # random draw
  beta_hetm=rtri(1,min,max,mode)
  # Apply the random draw value to all sub-population
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      if (s>=1&s<=5||s>=11&s<=15) {beta[i,s,1]=beta_hetf}
      if (s>=6&s<=10){beta[i,s,1]=beta_fsw}
      if (s>=16&s<=25){beta[i,s,1]=beta_hetm}
      for (g in 2:generations) {beta[i,s,g]=beta[i,s,1]}}}

  # Condom effectiveness
  # HET
  min=p_condomEffect[1,"Min"]
  max=p_condomEffect[1,"Max"]
  mode=p_condomEffect[1,"Mode"]
  # random draw
  epsilon_s=rtri(1,min,max,mode)
  # Apply the random draw value to all sub-population
  if (t %in% list_trial){
    for (s in 1:sgroup) {epsilon[i,s,1]=epsilon_s
    for (g in 2:generations) {
      epsilon[i,s,g]=epsilon[i,s,1]}}}

  # Proportion of consistent condom use
  # Low-risk urban women
  for (l in 1:5){
    a=p_condomUse[l,"Alpha"]
    b=p_condomUse[l,"Beta"]
    c_r_uw1=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l,1]=c_r_uw1}
    
    a=p_condomUse[l+5,"Alpha"]
    b=p_condomUse[l+5,"Beta"]
    c_r_uw2=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l,73]=c_r_uw2}
    
    a=p_condomUse[l+10,"Alpha"]
    b=p_condomUse[l+10,"Beta"]
    c_r_uw3=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l,133]=c_r_uw3}
    
    # Low-risk rural women
    a=p_condomUse[l+15,"Alpha"]
    b=p_condomUse[l+15,"Beta"]
    c_r_rw1=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+10,1]=c_r_rw1}
    
    a=p_condomUse[l+20,"Alpha"]
    b=p_condomUse[l+20,"Beta"]
    c_r_rw2=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+10,73]=c_r_rw2}
    
    a=p_condomUse[l+25,"Alpha"]
    b=p_condomUse[l+25,"Beta"]
    c_r_rw3=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+10,133]=c_r_rw3}
    
    # Urban men
    a=p_condomUse[l+30,"Alpha"]
    b=p_condomUse[l+30,"Beta"]
    c_r_um1=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+15,1]=c_r_um1}
    
    a=p_condomUse[l+35,"Alpha"]
    b=p_condomUse[l+35,"Beta"]
    c_r_um2=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+15,73]=c_r_um2}
    
    a=p_condomUse[l+40,"Alpha"]
    b=p_condomUse[l+40,"Beta"]
    c_r_um3=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+15,133]=c_r_um3}
    
    # Rural men
    a=p_condomUse[l+45,"Alpha"]
    b=p_condomUse[l+45,"Beta"]
    c_r_rm1=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+20,1]=c_r_rm1}
    
    a=p_condomUse[l+50,"Alpha"]
    b=p_condomUse[l+50,"Beta"]
    c_r_rm2=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+20,73]=c_r_rm2}
    
    a=p_condomUse[l+55,"Alpha"]
    b=p_condomUse[l+55,"Beta"]
    c_r_rm3=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+20,133]=c_r_rm3}
    
    # High-risk women
    a=p_condomUse[l+60,"Alpha"]
    b=p_condomUse[l+60,"Beta"]
    c_r_fsw=rbeta(1,a,b)
    if (t %in% list_trial){c_r[i,l+5,1]=c_r_fsw
    c_r[i,l+5,73]=c_r[i,l+5,1]
    c_r[i,l+5,133]=c_r[i,l+5,1]}}
  
  #Apply to all generations
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      for (g in 2:generations) {
        if (g>=1&g<=72) {
          c_r[i,s,g]=c_r[i,s,1]
          #c_r_ipv[i,s,g]=c_r_ipv[i,s,1]
        }
        if (g>=73&g<=132) {
          c_r[i,s,g]=c_r[i,s,73]
          #c_r_ipv[i,s,g]=c_r_ipv[i,s,49]
        }
        if (g>=133&g<=generations) {
          c_r[i,s,g]=c_r[i,s,133]
          #c_r_ipv[i,s,g]=c_r_ipv[i,s,133]
        }}}
  }
  
  # Effectiveness of viral suppression
  # Sexual
  min=p_VSEffect[1,"Min"]
  max=p_VSEffect[1,"Max"]
  mode=p_VSEffect[1,"Mode"]
  # random draw
  upsilon_s=rtri(1,min,max,mode)
  if (t %in% list_trial){upsilon[i,1,1]=upsilon_s}
  
  # Apply the random draw value to all sub-population
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      upsilon[i,s,1]=upsilon[i,1,1]
      for (g in 2:generations) {
        upsilon[i,s,g]=upsilon[i,s,1]
      }}
  }
 
  # Number of sex acts
  # Low risk
  min=n_sexact[1,"Min"]
  max=n_sexact[1,"Max"]
  mode=n_sexact[1,"Mode"] 
  n_r_lr=rtri(1,min,max,mode)
  if (t %in% list_trial){n_r[i,1,1]=n_r_lr
  n_r[i,2,1]=n_r[i,1,1]}
  # High risk
  min=n_sexact[2,"Min"]
  max=n_sexact[2,"Max"]
  mode=n_sexact[2,"Mode"] 
  n_r_hr=rtri(1,min,max,mode)*runif(n=1,min=1,max=5)
  if (t %in% list_trial){n_r[i,6,1]=n_r_hr}
  
  #Youth multiplier
  min=0.2
  max=1.0
  mode=0.5
  youth_multi_f=rtri(1,min,max,mode)
  youth_multi_m=rtri(1,min,max,mode)
  
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      # Assign value to all other sub-populations 
      if (s==1||s==11){n_r[i,s,1]=n_r[i,2,1]*youth_multi_f}
      if (s==16||s==21){n_r[i,s,1]=n_r[i,2,1]*youth_multi_m}
      if (s>=2&s<=5||s>=12&s<=15||s>=17&s<=20||s>=22&s<=25) {n_r[i,s,1]=n_r[i,2,1]}
      if (s>=6&s<=10) {n_r[i,s,1]=n_r[i,6,1]}
      #Apply to all generations
      for (g in 2:generations) {n_r[i,s,g]=n_r[i,s,1]}
    }
  }
 
  # Parameter inputs for infectious stages
  # Population growth
  # for each trial, draw randomly from the input distribution (monthly)
  for(s in 1:sgroup){
    u=p_popgrowth[s,"Mean"]
    sd=p_popgrowth[s,"SD"]
    # Random draw (Note that rho can be either positive or negative)
    rho_s <- rnorm(1,mean=u,sd=sd)
    if (t %in% list_trial){rho[i,s,1]=rho_s}
  }
  
  # For sub-population with same characters, we assume that the parameters are the same
  # For example, for all individuals that are white and aged 14-25 would share the same population growth rate
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      # We assume that the population growth does not change overtime
      for (g in 2:generations){rho[i,s,g]=rho[i,s,1]}
    }
  }
 
  # Natural history
  # for each trial, draw randomly from the input distribution (monthly)
  # Assume that the natural history does not differ by engagement in care
  if (t %in% list_trial){m_delta[i,1,1]=1}
  
  # Disease progression from CD4>500 to CD4>350-500
  a=p_DisProgress[1,"Alpha"]
  b=p_DisProgress[1,"Beta"]
  # Random draw
  delta_1_f <- rbeta(1,a,b)
  if (t %in% list_trial){delta_1[i,1,1] <- 1-exp(log(1-delta_1_f)/12)}
  
  a=p_DisProgress[4,"Alpha"]
  b=p_DisProgress[4,"Beta"]
  # Random draw
  delta_1_m <- rbeta(1,a,b)
  if (t %in% list_trial){delta_1[i,16,1] <- 1-exp(log(1-delta_1_m)/12)}
  
  # Step length is set to 1 since sub-populations would have same natural history with the same CD4 level
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      if (s>=1 & s<=15) {
        delta_1[i,s,1] <- delta_1[i,1,1]
        delta_4[i,s,1] <- delta_1[i,1,1]*m_delta[i,1,1]}
      if (s>=16 & s<=25) {
        delta_1[i,s,1] <- delta_1[i,16,1]
        delta_4[i,s,1] <- delta_1[i,16,1]*m_delta[i,1,1]
      }
      # We assume that the natural history does not change overtime
      for (g in 2:generations){
        delta_1[i,s,g]=delta_1[i,s,1]
        delta_4[i,s,g]=delta_4[i,s,1]}
    }}
  
  # Disease progression from CD4>350-500 to CD4>200-350
  a=p_DisProgress[2,"Alpha"]
  b=p_DisProgress[2,"Beta"]
  # Random draw
  delta_2_f <- rbeta(1,a,b)
  if (t %in% list_trial){ delta_2[i,1,1] <- 1-exp(log(1-delta_2_f)/12)}
  
  a=p_DisProgress[5,"Alpha"]
  b=p_DisProgress[5,"Beta"]
  # Random draw
  delta_2_m <- rbeta(1,a,b)
  if (t %in% list_trial){delta_2[i,16,1] <- 1-exp(log(1-delta_2_m)/12)}
  
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      # Step length is set to 1 since sub-populations would have same natural history with the same CD4 level
      if (s>=1 & s<=15) {
        delta_2[i,s,1] <- delta_2[i,1,1]
        delta_5[i,s,1] <- delta_2[i,1,1]*m_delta[i,1,1]}
      if (s>=16 & s<=25) {
        delta_2[i,s,1] <- delta_2[i,16,1]
        delta_5[i,s,1] <- delta_2[i,16,1]*m_delta[i,1,1]}
      # We assume that the natural history does not change overtime
      for (g in 2:generations){
        delta_2[i,s,g]=delta_2[i,s,1]
        delta_5[i,s,g]=delta_5[i,s,1]}
    }}
  
  # Disease progression from CD4>200-350 to CD4<200
  a=p_DisProgress[3,"Alpha"]
  b=p_DisProgress[3,"Beta"]
  # Random draw
  delta_3_f <- rbeta(1,a,b)
  if (t %in% list_trial){delta_3[i,1,1] <- 1-exp(log(1-delta_3_f)/12)}
  
  a=p_DisProgress[6,"Alpha"]
  b=p_DisProgress[6,"Beta"]
  # Random draw
  delta_3_m <- rbeta(1,a,b)
  if (t %in% list_trial){delta_3[i,16,1] <- 1-exp(log(1-delta_3_m)/12)}
  
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      # Step length is set to 1 since sub-populations would have same natural history with the same CD4 level
      if (s>=1 & s<=15) {
        delta_3[i,s,1] <- delta_3[i,1,1]
        delta_6[i,s,1] <- delta_3[i,1,1]*m_delta[i,1,1]}
      if (s>=16 & s<=25) {
        delta_3[i,s,1] <- delta_3[i,16,1]
        delta_6[i,s,1] <- delta_3[i,16,1]*m_delta[i,1,1]}
      # We assume that the natural history does not change overtime
      for (g in 2:generations){
        delta_3[i,s,g]=delta_3[i,s,1]
        delta_6[i,s,g]=delta_6[i,s,1]}
    }
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
    alpha_1_uw=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l,1]=-log(1-alpha_1_uw)
      alpha_2[i,l,1]=alpha_1[i,l,1]
      alpha_3[i,l,1]=alpha_1[i,l,1]
      alpha_4[i,l,1]=alpha_1[i,l,1]*cd4_multi}
    
    a=p_diagnosis[l+5,"Alpha"]
    b=p_diagnosis[l+5,"Beta"]
    alpha_1_uw2=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l,73]=-log(1-alpha_1_uw2)
      alpha_2[i,l,73]=alpha_1[i,l,73]
      alpha_3[i,l,73]=alpha_1[i,l,73]
      alpha_4[i,l,73]=alpha_1[i,l,73]*cd4_multi}
    
    a=p_diagnosis[l+10,"Alpha"]
    b=p_diagnosis[l+10,"Beta"]
    alpha_1_uw3=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l,133]=-log(1-alpha_1_uw3)
      alpha_2[i,l,133]=alpha_1[i,l,133]
      alpha_3[i,l,133]=alpha_1[i,l,133]
      alpha_4[i,l,133]=alpha_1[i,l,133]*cd4_multi}
    
    # Low-risk rural women
    a=p_diagnosis[l+15,"Alpha"]
    b=p_diagnosis[l+15,"Beta"]
    alpha_1_rw=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l+10,1]=-log(1-alpha_1_rw)
      alpha_2[i,l+10,1]=alpha_1[i,l+10,1]
      alpha_3[i,l+10,1]=alpha_1[i,l+10,1]
      alpha_4[i,l+10,1]=alpha_1[i,l+10,1]*cd4_multi}
    
    a=p_diagnosis[l+20,"Alpha"]
    b=p_diagnosis[l+20,"Beta"]
    alpha_1_rw2=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l+10,73]=-log(1-alpha_1_rw2)
      alpha_2[i,l+10,73]=alpha_1[i,l+10,73]
      alpha_3[i,l+10,73]=alpha_1[i,l+10,73]
      alpha_4[i,l+10,73]=alpha_1[i,l+10,73]*cd4_multi}
    
    a=p_diagnosis[l+25,"Alpha"]
    b=p_diagnosis[l+25,"Beta"]
    alpha_1_rw3=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l+10,133]=-log(1-alpha_1_rw3)
      alpha_2[i,l+10,133]=alpha_1[i,l+10,133]
      alpha_3[i,l+10,133]=alpha_1[i,l+10,133]
      alpha_4[i,l+10,133]=alpha_1[i,l+10,133]*cd4_multi}
    
    # Urban men
    a=p_diagnosis[l+30,"Alpha"]
    b=p_diagnosis[l+30,"Beta"]
    alpha_1_um=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l+15,1]=-log(1-alpha_1_um)
      alpha_2[i,l+15,1]=alpha_1[i,l+15,1]
      alpha_3[i,l+15,1]=alpha_1[i,l+15,1]
      alpha_4[i,l+15,1]=alpha_1[i,l+15,1]*cd4_multi}
    
    a=p_diagnosis[l+35,"Alpha"]
    b=p_diagnosis[l+35,"Beta"]
    alpha_1_um2=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l+15,73]=-log(1-alpha_1_um2)
      alpha_2[i,l+15,73]=alpha_1[i,l+15,73]
      alpha_3[i,l+15,73]=alpha_1[i,l+15,73]
      alpha_4[i,l+15,73]=alpha_1[i,l+15,73]*cd4_multi}
    
    a=p_diagnosis[l+40,"Alpha"]
    b=p_diagnosis[l+40,"Beta"]
    alpha_1_um3=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l+15,133]=-log(1-alpha_1_um3)
      alpha_2[i,l+15,133]=alpha_1[i,l+15,133]
      alpha_3[i,l+15,133]=alpha_1[i,l+15,133]
      alpha_4[i,l+15,133]=alpha_1[i,l+15,133]*cd4_multi}
    
    # Rural men
    a=p_diagnosis[l+45,"Alpha"]
    b=p_diagnosis[l+45,"Beta"]
    alpha_1_rm=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l+20,1]=-log(1-alpha_1_rm)
      alpha_2[i,l+20,1]=alpha_1[i,l+20,1]
      alpha_3[i,l+20,1]=alpha_1[i,l+20,1]
      alpha_4[i,l+20,1]=alpha_1[i,l+20,1]*cd4_multi}
    
    a=p_diagnosis[l+50,"Alpha"]
    b=p_diagnosis[l+50,"Beta"]
    alpha_1_rm2=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l+20,73]=-log(1-alpha_1_rm2)
      alpha_2[i,l+20,73]=alpha_1[i,l+20,73]
      alpha_3[i,l+20,73]=alpha_1[i,l+20,73]
      alpha_4[i,l+20,73]=alpha_1[i,l+20,73]*cd4_multi}
    
    a=p_diagnosis[l+55,"Alpha"]
    b=p_diagnosis[l+55,"Beta"]
    alpha_1_rm3=rbeta(1,a,b)
    if (t %in% list_trial){
      alpha_1[i,l+20,133]=-log(1-alpha_1_rm3)
      alpha_2[i,l+20,133]=alpha_1[i,l+20,133]
      alpha_3[i,l+20,133]=alpha_1[i,l+20,133]
      alpha_4[i,l+20,133]=alpha_1[i,l+20,133]*cd4_multi}
  
  # High-risk women
  a=p_diagnosis[61,"Alpha"]
  b=p_diagnosis[61,"Beta"]
  alpha_1_fsw=rbeta(1,a,b)
  if (t %in% list_trial){
    alpha_1[i,6,1]=-log(1-alpha_1_fsw)
    alpha_2[i,6,1]=alpha_1[i,6,1]
    alpha_3[i,6,1]=alpha_1[i,6,1]
    alpha_4[i,6,1]=alpha_1[i,6,1]*cd4_multi}
  }
  
  if (t %in% list_trial){
    for (s in 1:sgroup){
      if (s>=6&s<=10){
        alpha_1[i,s,1]=alpha_1[i,6,1]
        alpha_2[i,s,1]=alpha_2[i,6,1]
        alpha_3[i,s,1]=alpha_3[i,6,1]
        alpha_4[i,s,1]=alpha_4[i,6,1]
        alpha_1[i,s,73]=alpha_1[i,6,1]
        alpha_2[i,s,73]=alpha_2[i,6,1]
        alpha_3[i,s,73]=alpha_3[i,6,1]
        alpha_4[i,s,73]=alpha_4[i,6,1]
        alpha_1[i,s,133]=alpha_1[i,6,1]
        alpha_2[i,s,133]=alpha_2[i,6,1]
        alpha_3[i,s,133]=alpha_3[i,6,1]
        alpha_4[i,s,133]=alpha_4[i,6,1]}}}
  
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      if (s<=15){
        alpha_1[i,s,1]=alpha_1[i,s,1]*correct_f
        alpha_2[i,s,1]=alpha_2[i,s,1]*correct_f
        alpha_3[i,s,1]=alpha_3[i,s,1]*correct_f
        alpha_4[i,s,1]=alpha_4[i,s,1]*correct_f
        alpha_1[i,s,73]=alpha_1[i,s,73]*correct_f
        alpha_2[i,s,73]=alpha_2[i,s,73]*correct_f
        alpha_3[i,s,73]=alpha_3[i,s,73]*correct_f
        alpha_4[i,s,73]=alpha_4[i,s,73]*correct_f
        alpha_1[i,s,133]=alpha_1[i,s,133]*correct_f
        alpha_2[i,s,133]=alpha_2[i,s,133]*correct_f
        alpha_3[i,s,133]=alpha_3[i,s,133]*correct_f
        alpha_4[i,s,133]=alpha_4[i,s,133]*correct_f
      }
      if (s>15){
        alpha_1[i,s,1]=alpha_1[i,s,1]*correct_m
        alpha_2[i,s,1]=alpha_2[i,s,1]*correct_m
        alpha_3[i,s,1]=alpha_3[i,s,1]*correct_m
        alpha_4[i,s,1]=alpha_4[i,s,1]*correct_m
        alpha_1[i,s,73]=alpha_1[i,s,73]*correct_m
        alpha_2[i,s,73]=alpha_2[i,s,73]*correct_m
        alpha_3[i,s,73]=alpha_3[i,s,73]*correct_m
        alpha_4[i,s,73]=alpha_4[i,s,73]*correct_m
        alpha_1[i,s,133]=alpha_1[i,s,133]*correct_m
        alpha_2[i,s,133]=alpha_2[i,s,133]*correct_m
        alpha_3[i,s,133]=alpha_3[i,s,133]*correct_m
        alpha_4[i,s,133]=alpha_4[i,s,133]*correct_m
      }}
    
    for (s in 1:sgroup){
      # Convert to probability
      alpha_1[i,s,1]=1-exp(-alpha_1[i,s,1])
      alpha_2[i,s,1]=1-exp(-alpha_2[i,s,1])
      alpha_3[i,s,1]=1-exp(-alpha_3[i,s,1])
      alpha_4[i,s,1]=1-exp(-alpha_4[i,s,1])
      alpha_1[i,s,73]=1-exp(-alpha_1[i,s,73])
      alpha_2[i,s,73]=1-exp(-alpha_2[i,s,73])
      alpha_3[i,s,73]=1-exp(-alpha_3[i,s,73])
      alpha_4[i,s,73]=1-exp(-alpha_4[i,s,73])
      alpha_1[i,s,133]=1-exp(-alpha_1[i,s,133])
      alpha_2[i,s,133]=1-exp(-alpha_2[i,s,133])
      alpha_3[i,s,133]=1-exp(-alpha_3[i,s,133])
      alpha_4[i,s,133]=1-exp(-alpha_4[i,s,133])}
    
    #Apply to all generations
    for (s in 1:sgroup) {
      for (g in 2:generations) {
        if (g>=1&g<=72) {
          alpha_1[i,s,g]=alpha_1[i,s,1]
          alpha_2[i,s,g]=alpha_2[i,s,1]
          alpha_3[i,s,g]=alpha_3[i,s,1]
          alpha_4[i,s,g]=alpha_4[i,s,1]}
        if (g>=73&g<=132) {
          alpha_1[i,s,g]=alpha_1[i,s,73]
          alpha_2[i,s,g]=alpha_2[i,s,73]
          alpha_3[i,s,g]=alpha_3[i,s,73]
          alpha_4[i,s,g]=alpha_4[i,s,73]}
        if (g>=133&g<=generations) {
          alpha_1[i,s,g]=alpha_1[i,s,133]
          alpha_2[i,s,g]=alpha_2[i,s,133]
          alpha_3[i,s,g]=alpha_3[i,s,133]
          alpha_4[i,s,g]=alpha_4[i,s,133]}}}
  }

   
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
  sigma_1_fsw=rtri(1,min,max,mode)
  if (t %in% list_trial){
    sigma_1[i,6,1]=-log(1-sigma_1_fsw)
    sigma_2[i,6,1]=sigma_1[i,6,1]
    sigma_3[i,6,1]=sigma_1[i,6,1]
    sigma_4[i,6,1]=sigma_1[i,6,1]*cd4_multi}
  # Low-risk
  min=p_link[2,"Min"]
  max=p_link[2,"Max"]
  mode=p_link[2,"Mode"]
  sigma_1_lr=rtri(1,min,max,mode)
  if (t %in% list_trial){
    sigma_1[i,1,1]=-log(1-sigma_1_lr)
    sigma_2[i,1,1]=sigma_1[i,1,1]
    sigma_3[i,1,1]=sigma_1[i,1,1]
    sigma_4[i,1,1]=sigma_1[i,1,1]*cd4_multi}
  
  if (t %in% list_trial){
    # Convert to probabilities
    sigma_1[i,1,1]=1-exp(-sigma_1[i,1,1])
    sigma_2[i,1,1]=1-exp(-sigma_2[i,1,1])
    sigma_3[i,1,1]=1-exp(-sigma_3[i,1,1])
    sigma_4[i,1,1]=1-exp(-sigma_4[i,1,1])
    sigma_1[i,6,1]=1-exp(-sigma_1[i,6,1])
    sigma_2[i,6,1]=1-exp(-sigma_2[i,6,1])
    sigma_3[i,6,1]=1-exp(-sigma_3[i,6,1])
    sigma_4[i,6,1]=1-exp(-sigma_4[i,6,1])
    
    for (s in 1:sgroup){
      if (s>=6&s<=10){
        sigma_1[i,s,1]=sigma_1[i,6,1]
        sigma_2[i,s,1]=sigma_2[i,6,1]
        sigma_3[i,s,1]=sigma_3[i,6,1]
        sigma_4[i,s,1]=sigma_4[i,6,1]
      }else{sigma_1[i,s,1]=sigma_1[i,1,1]
      sigma_2[i,s,1]=sigma_2[i,1,1]
      sigma_3[i,s,1]=sigma_3[i,1,1]
      sigma_4[i,s,1]=sigma_4[i,1,1]}
      # We assume that the probability of diagnosis will change every two years (24 months)
      for (g in 2:generations){
        sigma_1[i,s,g]=sigma_1[i,s,1]
        sigma_2[i,s,g]=sigma_2[i,s,1]
        sigma_3[i,s,g]=sigma_3[i,s,1]
        sigma_4[i,s,g]=sigma_4[i,s,1]}}
  }
  
  # Probability of LTFU
  # For each trial, draw randomly from the distribution 
  a=p_LTFU[1,"Alpha"]
  b=p_LTFU[1,"Beta"]
  #mode=p_LTFU[1,"Mode"]
  gamma_1_s=rbeta(1,a,b)
  if (t %in% list_trial){gamma_1[i,1,1]=1-exp(log(1-gamma_1_s)/12)}
  
  a=p_LTFU[2,"Alpha"]
  b=p_LTFU[2,"Beta"]
  #mode=p_LTFU[2,"Mode"]
  gamma_2_s=rbeta(1,a,b)
  if (t %in% list_trial){gamma_2[i,1,1]=1-exp(log(1-gamma_2_s)/12)}
  
  a=p_LTFU[3,"Alpha"]
  b=p_LTFU[3,"Beta"]
  #mode=p_LTFU[3,"Mode"]
  gamma_3_s=rbeta(1,a,b)
  if (t %in% list_trial){gamma_3[i,1,1]=1-exp(log(1-gamma_3_s)/12)}
  
  a=p_LTFU[4,"Alpha"]
  b=p_LTFU[4,"Beta"]
  #mode=p_LTFU[4,"Mode"]
  gamma_4_s=rbeta(1,a,b)
  if (t %in% list_trial){gamma_4[i,1,1]=1-exp(log(1-gamma_4_s)/12)}
  
  a=p_LTFU[5,"Alpha"]
  b=p_LTFU[5,"Beta"]
  #mode=p_LTFU[5,"Mode"]
  gamma_5_s=rbeta(1,a,b)
  if (t %in% list_trial){gamma_5[i,1,1]=1-exp(log(1-gamma_5_s)/12)
  gamma_9[i,1,1]=gamma_5[i,1,1]}
  
  a=p_LTFU[6,"Alpha"]
  b=p_LTFU[6,"Beta"]
  #mode=p_LTFU[6,"Mode"]
  gamma_6_s=rbeta(1,a,b)
  if (t %in% list_trial){
    gamma_6[i,1,1]=1-exp(log(1-gamma_6_s)/12)
    gamma_10[i,1,1]=gamma_6[i,1,1]}
  
  a=p_LTFU[7,"Alpha"]
  b=p_LTFU[7,"Beta"]
  #mode=p_LTFU[7,"Mode"]
  gamma_7_s=rbeta(1,a,b)
  if (t %in% list_trial){
    gamma_7[i,1,1]=1-exp(log(1-gamma_7_s)/12)
    gamma_11[i,1,1]=gamma_7[i,1,1]}
  
  a=p_LTFU[8,"Alpha"]
  b=p_LTFU[8,"Beta"]
  #mode=p_LTFU[8,"Mode"]
  gamma_8_s=rbeta(1,a,b)
  if (t %in% list_trial){
    gamma_8[i,1,1]=1-exp(log(1-gamma_8_s)/12)
    gamma_12[i,1,1]=gamma_8[i,1,1]}
  
  if (t %in% list_trial){
    for (s in 1:sgroup){
      gamma_1[i,s,1]=gamma_1[i,1,1]
      gamma_2[i,s,1]=gamma_2[i,1,1]
      gamma_3[i,s,1]=gamma_3[i,1,1]
      gamma_4[i,s,1]=gamma_4[i,1,1]
      gamma_5[i,s,1]=gamma_5[i,1,1]
      gamma_6[i,s,1]=gamma_6[i,1,1]
      gamma_7[i,s,1]=gamma_7[i,1,1]
      gamma_8[i,s,1]=gamma_8[i,1,1]
      gamma_9[i,s,1]=gamma_9[i,1,1]
      gamma_10[i,s,1]=gamma_10[i,1,1]
      gamma_11[i,s,1]=gamma_11[i,1,1]
      gamma_12[i,s,1]=gamma_12[i,1,1]
      # We assume that the probability of diagnosis will be constant
      for (g in 2:generations){
        gamma_1[i,s,g]=gamma_1[i,s,1]
        gamma_2[i,s,g]=gamma_2[i,s,1]
        gamma_3[i,s,g]=gamma_3[i,s,1]
        gamma_4[i,s,g]=gamma_4[i,s,1]
        gamma_5[i,s,g]=gamma_5[i,s,1]
        gamma_6[i,s,g]=gamma_6[i,s,1]
        gamma_7[i,s,g]=gamma_7[i,s,1]
        gamma_8[i,s,g]=gamma_8[i,s,1]
        gamma_9[i,s,g]=gamma_9[i,s,1]
        gamma_10[i,s,g]=gamma_10[i,s,1]
        gamma_11[i,s,g]=gamma_11[i,s,1]
        gamma_12[i,s,g]=gamma_12[i,s,1]
      }
    }   
  }
 
  # Probability of On ART and Suppressed
  # For each trial, draw randomly from the distribution 
  a=p_ARTVS[1,"Alpha"]
  b=p_ARTVS[1,"Beta"]
  #mode=p_ARTVS[1,"Mode"]
  theta_1_s=rbeta(1,a,b)
  if (t %in% list_trial){theta_1[i,1,1]=1-exp(log(1-theta_1_s)/12)}
  
  a=p_ARTVS[2,"Alpha"]
  b=p_ARTVS[2,"Beta"]
  #mode=p_ARTVS[2,"Mode"]
  theta_2_s=rbeta(1,a,b)
  if (t %in% list_trial){theta_2[i,1,1]=1-exp(log(1-theta_2_s)/12)}
  
  a=p_ARTVS[3,"Alpha"]
  b=p_ARTVS[3,"Beta"]
  #mode=p_ARTVS[3,"Mode"]
  theta_3_s=rbeta(1,a,b)
  if (t %in% list_trial){theta_3[i,1,1]=1-exp(log(1-theta_3_s)/12)}
  
  a=p_ARTVS[4,"Alpha"]
  b=p_ARTVS[4,"Beta"]
  #mode=p_ARTVS[4,"Mode"]
  theta_4_s=rbeta(1,a,b)
  if (t %in% list_trial){theta_4[i,1,1]=1-exp(log(1-theta_4_s)/12)}
  
  a=p_ARTVS[5,"Alpha"]
  b=p_ARTVS[5,"Beta"]
  #mode=p_ARTVS[5,"Mode"]
  theta_1_s2=rbeta(1,a,b)
  if (t %in% list_trial){theta_1[i,1,37]=1-exp(log(1-theta_1_s2)/12)}
  
  a=p_ARTVS[6,"Alpha"]
  b=p_ARTVS[6,"Beta"]
  #mode=p_ARTVS[6,"Mode"]
  theta_2_s2=rbeta(1,a,b)
  if (t %in% list_trial){theta_2[i,1,37]=1-exp(log(1-theta_2_s2)/12)}
  
  
  a=p_ARTVS[7,"Alpha"]
  b=p_ARTVS[7,"Beta"]
  #mode=p_ARTVS[7,"Mode"]
  theta_3_s2=rbeta(1,a,b)
  if (t %in% list_trial){theta_3[i,1,37]=1-exp(log(1-theta_3_s2)/12)}
  
  a=p_ARTVS[8,"Alpha"]
  b=p_ARTVS[8,"Beta"]
  #mode=p_ARTVS[8,"Mode"]
  theta_4_s2=rbeta(1,a,b)
  if (t %in% list_trial){theta_4[i,1,37]=1-exp(log(1-theta_4_s2)/12)}
  
  a=p_ARTVS[9,"Alpha"]
  b=p_ARTVS[9,"Beta"]
  #mode=p_ARTVS[9,"Mode"]
  theta_1_s3=rbeta(1,a,b)
  if (t %in% list_trial){theta_1[i,1,97]=1-exp(log(1-theta_1_s3)/12)}
  
  a=p_ARTVS[10,"Alpha"]
  b=p_ARTVS[10,"Beta"]
  #mode=p_ARTVS[10,"Mode"]
  theta_2_s3=rbeta(1,a,b)
  if (t %in% list_trial){theta_2[i,1,97]=1-exp(log(1-theta_2_s3)/12)}
  
  a=p_ARTVS[11,"Alpha"]
  b=p_ARTVS[11,"Beta"]
  #mode=p_ARTVS[11,"Mode"]
  theta_3_s3=rbeta(1,a,b)
  if (t %in% list_trial){theta_3[i,1,97]=1-exp(log(1-theta_3_s3)/12)}
  
  a=p_ARTVS[12,"Alpha"]
  b=p_ARTVS[12,"Beta"]
  #mode=p_ARTVS[12,"Mode"]
  theta_4_s3=rbeta(1,a,b)
  if (t %in% list_trial){theta_4[i,1,97]=1-exp(log(1-theta_4_s3)/12)}
  
  a=p_ARTVS[13,"Alpha"]
  b=p_ARTVS[13,"Beta"]
  #mode=p_ARTVS[13,"Mode"]
  theta_1_s4=rbeta(1,a,b)
  if (t %in% list_trial){theta_1[i,1,145]=1-exp(log(1-theta_1_s4)/12)}
  
  a=p_ARTVS[14,"Alpha"]
  b=p_ARTVS[14,"Beta"]
  #mode=p_ARTVS[14,"Mode"]
  theta_2_s4=rbeta(1,a,b)
  if (t %in% list_trial){theta_2[i,1,145]=1-exp(log(1-theta_2_s4)/12)}
  
  a=p_ARTVS[15,"Alpha"]
  b=p_ARTVS[15,"Beta"]
  #mode=p_ARTVS[15,"Mode"]
  theta_3_s4=rbeta(1,a,b)
  if (t %in% list_trial){theta_3[i,1,145]=1-exp(log(1-theta_3_s4)/12)}
  
  a=p_ARTVS[16,"Alpha"]
  b=p_ARTVS[16,"Beta"]
  #mode=p_ARTVS[16,"Mode"]
  theta_4_s4=rbeta(1,a,b)
  if (t %in% list_trial){theta_4[i,1,145]=1-exp(log(1-theta_4_s4)/12)}
  
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      theta_1[i,s,1]=theta_1[i,1,1]
      theta_2[i,s,1]=theta_2[i,1,1]
      theta_3[i,s,1]=theta_3[i,1,1]
      theta_4[i,s,1]=theta_4[i,1,1]
      theta_1[i,s,37]=theta_1[i,1,37]
      theta_2[i,s,37]=theta_2[i,1,37]
      theta_3[i,s,37]=theta_3[i,1,37]
      theta_4[i,s,37]=theta_4[i,1,37]
      theta_1[i,s,97]=theta_1[i,1,97]
      theta_2[i,s,97]=theta_2[i,1,97]
      theta_3[i,s,97]=theta_3[i,1,97]
      theta_4[i,s,97]=theta_4[i,1,97]
      theta_1[i,s,145]=theta_1[i,1,145]
      theta_2[i,s,145]=theta_2[i,1,145]
      theta_3[i,s,145]=theta_3[i,1,145]
      theta_4[i,s,145]=theta_4[i,1,145]
      if (s>=6&s<11){
        theta_1[i,s,1]=theta_1[i,1,1]#*runif(1,min=0.2,max=1)
        theta_2[i,s,1]=theta_2[i,1,1]#*runif(1,min=0.2,max=1)
        theta_3[i,s,1]=theta_3[i,1,1]#*runif(1,min=0.2,max=1)
        theta_4[i,s,1]=theta_4[i,1,1]#*runif(1,min=0.2,max=1)
        theta_1[i,s,37]=theta_1[i,1,37]#*runif(1,min=0.2,max=1)
        theta_2[i,s,37]=theta_2[i,1,37]#*runif(1,min=0.2,max=1)
        theta_3[i,s,37]=theta_3[i,1,37]#*runif(1,min=0.2,max=1)
        theta_4[i,s,37]=theta_4[i,1,37]#*runif(1,min=0.2,max=1)
        theta_1[i,s,97]=theta_1[i,1,97]#*runif(1,min=0.2,max=1)
        theta_2[i,s,97]=theta_2[i,1,97]#*runif(1,min=0.2,max=1)
        theta_3[i,s,97]=theta_3[i,1,97]#*runif(1,min=0.2,max=1)
        theta_4[i,s,97]=theta_4[i,1,97]#*runif(1,min=0.2,max=1)
        theta_1[i,s,145]=theta_1[i,1,145]#*runif(1,min=0.2,max=1)
        theta_2[i,s,145]=theta_2[i,1,145]#*runif(1,min=0.2,max=1)
        theta_3[i,s,145]=theta_3[i,1,145]#*runif(1,min=0.2,max=1)
        theta_4[i,s,145]=theta_4[i,1,145]#*runif(1,min=0.2,max=1)
      }
      # # We assume that the probability of diagnosis will change every years (12 months)
      for (g in 2:generations){
        if (g<=36) {
          theta_1[i,s,g]=theta_1[i,s,1]
          theta_2[i,s,g]=theta_2[i,s,1]
          theta_3[i,s,g]=theta_3[i,s,1]
          theta_4[i,s,g]=theta_4[i,s,1]}
        if (g>=37&g<=96) {
          theta_1[i,s,g]=theta_1[i,s,37]
          theta_2[i,s,g]=theta_2[i,s,37]
          theta_3[i,s,g]=theta_3[i,s,37]
          theta_4[i,s,g]=theta_4[i,s,37]}
        if (g>=97&g<=144) {
          theta_1[i,s,g]=theta_1[i,s,97]
          theta_2[i,s,g]=theta_2[i,s,97]
          theta_3[i,s,g]=theta_3[i,s,97]
          theta_4[i,s,g]=theta_4[i,s,97]}
        if (g>=145&g<=generations) {
          theta_1[i,s,g]=theta_1[i,s,145]
          theta_2[i,s,g]=theta_2[i,s,145]
          theta_3[i,s,g]=theta_3[i,s,145]
          theta_4[i,s,g]=theta_4[i,s,145]}}
    }   
  }

  if (t %in% list_trial){
    for (s in 1:sgroup){
      for (g in 1:generations){
        theta_1[i,s,g]=-log(1-theta_1[i,s,g])
        theta_2[i,s,g]=-log(1-theta_2[i,s,g])
        theta_3[i,s,g]=-log(1-theta_3[i,s,g])
        theta_4[i,s,g]=-log(1-theta_4[i,s,g])}}}

  correct_female=rtri(1,1,2,1.5)
  correct_old=rtri(1,1,2,1.5)
  if (t %in% list_trial){
    for (s in 1:sgroup){
      for (g in 1:generations){
        if (s>=2&s<=5||s>=7&s<=10||s>=12&s<=15){
          theta_1[i,s,g]=theta_1[i,s,g]*correct_female*correct_old
          theta_2[i,s,g]=theta_2[i,s,g]*correct_female*correct_old
          theta_3[i,s,g]=theta_3[i,s,g]*correct_female*correct_old
          theta_4[i,s,g]=theta_4[i,s,g]*correct_female*correct_old}
        if (s>=17&s<=20||s>=22&s<=25){
          theta_1[i,s,g]=theta_1[i,s,g]*correct_old
          theta_2[i,s,g]=theta_2[i,s,g]*correct_old
          theta_3[i,s,g]=theta_3[i,s,g]*correct_old
          theta_4[i,s,g]=theta_4[i,s,g]*correct_old}
      }
    }
  }

  if (t %in% list_trial){
    for (s in 1:sgroup){
      for (g in 1:generations){
        theta_1[i,s,g]=1-exp(-theta_1[i,s,g])
        theta_2[i,s,g]=1-exp(-theta_2[i,s,g])
        theta_3[i,s,g]=1-exp(-theta_3[i,s,g])
        theta_4[i,s,g]=1-exp(-theta_4[i,s,g])}}}
 
  # Probability of On ART and Not suppressed
  # For each trial, draw randomly from the distribution 
  # Extract the baseline value
  a=p_ARTnotVS[1,"Alpha"]
  b=p_ARTnotVS[1,"Beta"]
  #mode=p_ARTnotVS[1,"Mode"]
  # Random draw
  psi_1_s=rbeta(1,a,b)
  if (t %in% list_trial){psi_1[i,1,1]=psi_1_s/12}
  
  a=p_ARTnotVS[2,"Alpha"]
  b=p_ARTnotVS[2,"Beta"]
  #mode=p_ARTnotVS[2,"Mode"]
  # Random draw
  psi_2_s=rbeta(1,a,b)
  if (t %in% list_trial){psi_2[i,1,1]=psi_2_s/12}
  
  a=p_ARTnotVS[3,"Alpha"]
  b=p_ARTnotVS[3,"Beta"]
  #mode=p_ARTnotVS[3,"Mode"]
  # Random draw
  psi_3_s=rbeta(1,a,b)
  if (t %in% list_trial){psi_3[i,1,1]=psi_3_s/12}
  
  a=p_ARTnotVS[4,"Alpha"]
  b=p_ARTnotVS[4,"Beta"]
  #mode=p_ARTnotVS[4,"Mode"]
  # Random draw
  psi_4_s=rbeta(1,a,b)
  if (t %in% list_trial){psi_4[i,1,1]=psi_4_s/12}
  
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
  
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      if (s>=1&s<=15){
        # Apply baseline value to all sub-population
        psi_1[i,s,1]=psi_1[i,1,1]*multi_female
        psi_2[i,s,1]=psi_2[i,1,1]*multi_female
        psi_3[i,s,1]=psi_3[i,1,1]*multi_female
        psi_4[i,s,1]=psi_4[i,1,1]*multi_female}
      if (s>=16&s<=25){
        # Apply baseline value to all sub-population
        psi_1[i,s,1]=psi_1[i,1,1]*multi_male
        psi_2[i,s,1]=psi_2[i,1,1]*multi_male
        psi_3[i,s,1]=psi_3[i,1,1]*multi_male
        psi_4[i,s,1]=psi_4[i,1,1]*multi_male}}
    
    for (s in 1:sgroup) {
      # Convert rate to probability
      psi_1[i,s,1]=1-exp(-psi_1[i,s,1])
      psi_2[i,s,1]=1-exp(-psi_2[i,s,1])
      psi_3[i,s,1]=1-exp(-psi_3[i,s,1])
      psi_4[i,s,1]=1-exp(-psi_4[i,s,1])
      # We assume that the probability is constant overtime
      for (g in 2:generations){
        psi_1[i,s,g]=psi_1[i,s,1]
        psi_2[i,s,g]=psi_2[i,s,1]
        psi_3[i,s,g]=psi_3[i,s,1]
        psi_4[i,s,g]=psi_4[i,s,1]}
    }
  }
 
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
    mu_9_f=rtri(1,min,max,mode)
    if (t %in% list_trial){
      mu_9[i,l,1]=mu_9_f
      mu_1[i,l,1]=mu_9_f*pre_cd4_1
      mu_5[i,l,1]=mu_1[i,l,1]
      mu_13[i,l,1]=mu_1[i,l,1]*deathmulti
      mu_1[i,l+5,1]=mu_1[i,l,1]
      mu_5[i,l+5,1]=mu_5[i,l,1]
      mu_9[i,l+5,1]=mu_9[i,l,1]
      mu_13[i,l+5,1]=mu_13[i,l,1]
      mu_1[i,l+10,1]=mu_1[i,l,1]
      mu_5[i,l+10,1]=mu_5[i,l,1]
      mu_9[i,l+10,1]=mu_9[i,l,1]
      mu_13[i,l+10,1]=mu_13[i,l,1]}
    # 
    min=p_death[l+5,"Min"]
    max=p_death[l+5,"Max"]
    mode=p_death[l+5,"Mode"]
    mu_10_f=rtri(1,min,max,mode)
    if (t %in% list_trial){
      mu_10[i,l,1]=mu_10_f
      mu_2[i,l,1]=mu_10_f*pre_cd4_2
      mu_6[i,l,1]=mu_2[i,l,1]
      mu_14[i,l,1]=mu_2[i,l,1]*deathmulti
      mu_2[i,l+5,1]=mu_2[i,l,1]
      mu_6[i,l+5,1]=mu_6[i,l,1]
      mu_10[i,l+5,1]=mu_10[i,l,1]
      mu_14[i,l+5,1]=mu_14[i,l,1]
      mu_2[i,l+10,1]=mu_2[i,l,1]
      mu_6[i,l+10,1]=mu_6[i,l,1]
      mu_10[i,l+10,1]=mu_10[i,l,1]
      mu_14[i,l+10,1]=mu_14[i,l,1]}
    # 
    min=p_death[l+10,"Min"]
    max=p_death[l+10,"Max"]
    mode=p_death[l+10,"Mode"]
    mu_11_f=rtri(1,min,max,mode)
    if (t %in% list_trial){
      mu_11[i,l,1]=mu_11_f
      mu_3[i,l,1]=mu_11_f*pre_cd4_3
      mu_7[i,l,1]=mu_3[i,l,1]
      mu_15[i,l,1]=mu_3[i,l,1]*deathmulti
      mu_3[i,l+5,1]=mu_3[i,l,1]
      mu_7[i,l+5,1]=mu_7[i,l,1]
      mu_11[i,l+5,1]=mu_11[i,l,1]
      mu_15[i,l+5,1]=mu_15[i,l,1]
      mu_3[i,l+10,1]=mu_3[i,l,1]
      mu_7[i,l+10,1]=mu_7[i,l,1]
      mu_11[i,l+10,1]=mu_11[i,l,1]
      mu_15[i,l+10,1]=mu_15[i,l,1]}
    # mu_3[i,l+15,1]=mu_3[i,l,1]
    # mu_7[i,l+15,1]=mu_7[i,l,1]
    # mu_11[i,l+15,1]=mu_11[i,l,1]
    # mu_15[i,l+15,1]=mu_15[i,l,1]
    # 
    min=p_death[l+15,"Min"]
    max=p_death[l+15,"Max"]
    mode=p_death[l+15,"Mode"]
    mu_12_f=rtri(1,min,max,mode)
    if (t %in% list_trial){
      mu_12[i,l,1]=mu_12_f
      mu_4[i,l,1]=mu_12_f*pre_cd4_4
      mu_8[i,l,1]=mu_4[i,l,1]
      mu_16[i,l,1]=mu_4[i,l,1]*deathmulti
      mu_4[i,l+5,1]=mu_4[i,l,1]
      mu_8[i,l+5,1]=mu_8[i,l,1]
      mu_12[i,l+5,1]=mu_12[i,l,1]
      mu_16[i,l+5,1]=mu_16[i,l,1]
      mu_4[i,l+10,1]=mu_4[i,l,1]
      mu_8[i,l+10,1]=mu_8[i,l,1]
      mu_12[i,l+10,1]=mu_12[i,l,1]
      mu_16[i,l+10,1]=mu_16[i,l,1]}
    # mu_4[i,l+15,1]=mu_4[i,l,1]
    # mu_8[i,l+15,1]=mu_8[i,l,1]
    # mu_12[i,l+15,1]=mu_12[i,l,1]
    # mu_16[i,l+15,1]=mu_16[i,l,1]
    
    #Men
    min=p_death[l+20,"Min"]
    max=p_death[l+20,"Max"]
    mode=p_death[l+20,"Mode"]
    mu_9_m=rtri(1,min,max,mode)
    if (t %in% list_trial){
      mu_9[i,l+15,1]=mu_9_m
      mu_1[i,l+15,1]=mu_9_m*pre_cd4_1
      mu_5[i,l+15,1]=mu_1[i,l+15,1]
      mu_13[i,l+15,1]=mu_1[i,l+15,1]*deathmulti
      mu_1[i,l+20,1]=mu_1[i,l+15,1]
      mu_5[i,l+20,1]=mu_5[i,l+15,1]
      mu_9[i,l+20,1]=mu_9[i,l+15,1]
      mu_13[i,l+20,1]=mu_13[i,l+15,1]}
    
    min=p_death[l+25,"Min"]
    max=p_death[l+25,"Max"]
    mode=p_death[l+25,"Mode"]
    mu_10_m=rtri(1,min,max,mode)
    if (t %in% list_trial){
      mu_10[i,l+15,1]=mu_10_m
      mu_2[i,l+15,1]=mu_10_m*pre_cd4_2
      mu_6[i,l+15,1]=mu_2[i,l+15,1]
      mu_14[i,l+15,1]=mu_2[i,l+15,1]*deathmulti
      mu_2[i,l+20,1]=mu_2[i,l+15,1]
      mu_6[i,l+20,1]=mu_6[i,l+15,1]
      mu_10[i,l+20,1]=mu_10[i,l+15,1]
      mu_14[i,l+20,1]=mu_14[i,l+15,1]}
    
    min=p_death[l+30,"Min"]
    max=p_death[l+30,"Max"]
    mode=p_death[l+30,"Mode"]
    mu_11_m=rtri(1,min,max,mode)
    if (t %in% list_trial){
      mu_11[i,l+15,1]=mu_11_m
      mu_3[i,l+15,1]=mu_11_m*pre_cd4_3
      mu_7[i,l+15,1]=mu_3[i,l+15,1]
      mu_15[i,l+15,1]=mu_3[i,l+15,1]*deathmulti
      mu_3[i,l+20,1]=mu_3[i,l+15,1]
      mu_7[i,l+20,1]=mu_7[i,l+15,1]
      mu_11[i,l+20,1]=mu_11[i,l+15,1]
      mu_15[i,l+20,1]=mu_15[i,l+15,1]}
    
    min=p_death[l+35,"Min"]
    max=p_death[l+35,"Max"]
    mode=p_death[l+35,"Mode"]
    mu_12_m=rtri(1,min,max,mode)
    if (t %in% list_trial){
      mu_12[i,l+15,1]=mu_12_m
      mu_4[i,l+15,1]=mu_12_m*pre_cd4_4
      mu_8[i,l+15,1]=mu_4[i,l+15,1]
      mu_16[i,l+15,1]=mu_4[i,l+15,1]*deathmulti
      mu_4[i,l+20,1]=mu_4[i,l+15,1]
      mu_8[i,l+20,1]=mu_8[i,l+15,1]
      mu_12[i,l+20,1]=mu_12[i,l+15,1]
      mu_16[i,l+20,1]=mu_16[i,l+15,1]}
  }
  
  if (t %in% list_trial){
    for (s in 1:sgroup) {
      # Convert to probabilities
      mu_1[i,s,1]=1-exp(-mu_1[i,s,1]/12)
      mu_2[i,s,1]=1-exp(-mu_2[i,s,1]/12)
      mu_3[i,s,1]=1-exp(-mu_3[i,s,1]/12)
      mu_4[i,s,1]=1-exp(-mu_4[i,s,1]/12)
      mu_5[i,s,1]=1-exp(-mu_5[i,s,1]/12)
      mu_6[i,s,1]=1-exp(-mu_6[i,s,1]/12)
      mu_7[i,s,1]=1-exp(-mu_7[i,s,1]/12)
      mu_8[i,s,1]=1-exp(-mu_8[i,s,1]/12)
      mu_9[i,s,1]=1-exp(-mu_9[i,s,1]/12)
      mu_10[i,s,1]=1-exp(-mu_10[i,s,1]/12)
      mu_11[i,s,1]=1-exp(-mu_11[i,s,1]/12)
      mu_12[i,s,1]=1-exp(-mu_12[i,s,1]/12)
      mu_13[i,s,1]=1-exp(-mu_13[i,s,1]/12)
      mu_14[i,s,1]=1-exp(-mu_14[i,s,1]/12)
      mu_15[i,s,1]=1-exp(-mu_15[i,s,1]/12)
      mu_16[i,s,1]=1-exp(-mu_16[i,s,1]/12)}
    
    # We assume that death does not change overtime
    
    for (g in 1:generations) {
      for (s in 1:sgroup){
        #undiagnosed, diagnosed but not linked to care, LTFU
        mu_1[i,s,g]=mu_1[i,s,1]
        mu_2[i,s,g]=mu_2[i,s,1]
        mu_3[i,s,g]=mu_3[i,s,1]
        mu_4[i,s,g]=mu_4[i,s,1]
        # linked to care
        mu_5[i,s,g]=mu_5[i,s,1]
        mu_6[i,s,g]=mu_6[i,s,1]
        mu_7[i,s,g]=mu_7[i,s,1]
        mu_8[i,s,g]=mu_8[i,s,1]
        # On ART & Suppressed
        mu_9[i,s,g]=mu_9[i,s,1]
        mu_10[i,s,g]=mu_10[i,s,1]
        mu_11[i,s,g]=mu_11[i,s,1]
        mu_12[i,s,g]=mu_12[i,s,1]
        # On ART & Not suppressed
        mu_13[i,s,g]=mu_13[i,s,1]
        mu_14[i,s,g]=mu_14[i,s,1]
        mu_15[i,s,g]=mu_15[i,s,1]
        mu_16[i,s,g]=mu_16[i,s,1]}
    }
  }
 
  # Probability of re-engagement in care
  # For each trial, draw randomly from the distribution 
  min=p_reengage[1,"Min"]
  max=p_reengage[1,"Max"]
  mode=p_reengage[1,"Mode"]
  tau_1_s=1-exp(-rtri(1,min,max,mode)/12)
  tau_2_s=1-exp(-rtri(1,min,max,mode)/12)
  tau_3_s=1-exp(-rtri(1,min,max,mode)/12)
  tau_4_s=1-exp(-rtri(1,min,max,mode)/12)
  # Random draw
  if (t %in% list_trial){
    tau_1[i,1,1]=tau_1_s
    tau_2[i,1,1]=tau_2_s
    tau_3[i,1,1]=tau_3_s
    tau_4[i,1,1]=tau_4_s
# p_reengage_trial[i,1]=tau_1[i,1,1]
# We assume that it is the same for all sub-population and overtime
for (s in 1:sgroup) {
  tau_1[i,s,1]=tau_1[i,1,1]
  tau_2[i,s,1]=tau_2[i,1,1]
  tau_3[i,s,1]=tau_3[i,1,1]
  tau_4[i,s,1]=tau_4[i,1,1]
  for (g in 2:generations) {
    tau_1[i,s,g]=tau_1[i,s,1]
    tau_2[i,s,g]=tau_2[i,s,1]
    tau_3[i,s,g]=tau_3[i,s,1]
    tau_4[i,s,g]=tau_4[i,s,1]
  }
}
  }
  
  if (t %in% list_trial ){i=i+1}
}


# Initalize baseline population characteristics
for (t in 1:(trial_select)) {
  #Initalize parameters for each sub-group
  for (s in 1:(sgroup)) {
    # Calculate the initial HIV prevalence 
    
    # Total HIV infected individuals
    N[t,s,1] <- X1[t,s,1]+X2[t,s,1]+X3[t,s,1]+X4[t,s,1]+X5[t,s,1]+X6[t,s,1]+X7[t,s,1]+X8[t,s,1]+X9[t,s,1]+X10[t,s,1]+X11[t,s,1]+X12[t,s,1]+X13[t,s,1]+X14[t,s,1]+X15[t,s,1]+X16[t,s,1]+X17[t,s,1]+X18[t,s,1]+X19[t,s,1]+X20[t,s,1]+X21[t,s,1]+X22[t,s,1]+X23[t,s,1]+X24[t,s,1]
    
    # Total population of individuals who are not viral suppressed
    Z[t,s,1] <- X1[t,s,1]+X2[t,s,1]+X3[t,s,1]+X4[t,s,1]+X5[t,s,1]+X6[t,s,1]+X7[t,s,1]+X8[t,s,1]+X9[t,s,1]+X10[t,s,1]+X11[t,s,1]+X12[t,s,1]+X13[t,s,1]+X14[t,s,1]+X15[t,s,1]+X16[t,s,1]+X21[t,s,1]+X22[t,s,1]+X23[t,s,1]+X24[t,s,1]
    
    # Total population of individuals who are viral suppressed
    V[t,s,1] <- X17[t,s,1]+X18[t,s,1]+X19[t,s,1]+X20[t,s,1]
    
    # Total population of individuals infected but undiagnosed 
    # Undiagnosed[t,s,1] <- X1[t,s,1]+X2[t,s,1]+X3[t,s,1]+X4[t,s,1]
    
    # Total population of individuals diagnosed with HIV
    Diagnosed[t,s,1] <- X5[t,s,1]+X6[t,s,1]+X7[t,s,1]+X8[t,s,1]+X9[t,s,1]+X10[t,s,1]+X11[t,s,1]+X12[t,s,1]+X13[t,s,1]+X14[t,s,1]+X15[t,s,1]+X16[t,s,1]+X17[t,s,1]+X18[t,s,1]+X19[t,s,1]+X20[t,s,1]+X21[t,s,1]+X22[t,s,1]+X23[t,s,1]+X24[t,s,1]
    
    # Total population of individuals newly diagnosed
    New_Diagnosed[t,s,1] <- X1[t,s,1]*alpha_1[t,s,g]+X2[t,s,1]*alpha_2[t,s,g]+X3[t,s,1]*alpha_3[t,s,g]+X4[t,s,1]*alpha_4[t,s,g]
    
    # Total population of inidivudlas linked to care
    # Linked[t,s,1] <- X9[t,s,1]+X10[t,s,1]+X11[t,s,1]+X12[t,s,1]
    
    # Total population of individuals lost
    # Lost[t,s,1] <-  X13[t,s,1]+X14[t,s,1]+X15[t,s,1]+X16[t,s,1]
    
    # Total population of individuals on ART and suppressed
    Suppressed[t,s,1] <- X17[t,s,1]+X18[t,s,1]+X19[t,s,1]+X20[t,s,1]
    
    # Total population of individuals on ART and not suppressed
    # Not_Suppressed[t,s,1] <- X21[t,s,1]+X22[t,s,1]+X23[t,s,1]+X24[t,s,1]
    
    # Total population of individuals on ART
    On_ART[t,s,1] <- X17[t,s,1]+X18[t,s,1]+X19[t,s,1]+X20[t,s,1]+X21[t,s,1]+X22[t,s,1]+X23[t,s,1]+X24[t,s,1]
    
    if (N[t,s,1]==0) {
      
      p_r[t,s,1] <- 1
      p_r_v[t,s,1] <- 0
      
    } else {
      # Prevalence of individuals who are viral supprressed among people living with HIV
      p_r_v[t,s,1] <- V[t,s,1]/N[t,s,1]  
      
      # percentage of individuals who are not viral supprressed among people living with HIV
      p_r[t,s,1] <- 1-p_r_v[t,s,1]
    }
    
    # initial value for force of infection and new infection (set to null)
    lambda_t_r[t,s,1] <- NA
    NI[t,s,1] <- NA
    N_total[t,1,1]<- NA
    #lambda_inf[t,s,1] <- 0
    # lambda_t_r[t,s,1] <- 0
    # NI[t,s,1] <- 0
    # N_total[t,1,1]<- 0
  }
}

# Start of the loop for each trial
for (t in 1:(trial_select)){
  # Start the loop by generations
  # The loop of generations goes before loops of sub-groups since operations in time period t+1 depends on the outcomes from all sub-groups in time period t.
  # The generations starts from time 2 so that the equations in R are in line with the difference equation. 
  for (g in 2:(generations))
  {
    for (s in 1:sgroup) {
      #for (j in 1:sgroup) {
      
      # total population of low risk urban women from last generation, who will mix with urban men 
      N_lr_u_f[t,s,g-1] <- S[t,1,g-1]+S[t,2,g-1]+S[t,3,g-1]+S[t,4,g-1]+S[t,5,g-1]+N[t,1,g-1]+N[t,2,g-1]+N[t,3,g-1]+N[t,4,g-1]+N[t,5,g-1]
      # total population of high risk urban women from last generation, who will mix with urban men 
      N_hr_f[t,s,g-1] <- S[t,6,g-1]+S[t,7,g-1]+S[t,8,g-1]+S[t,9,g-1]+S[t,10,g-1]+N[t,6,g-1]+N[t,7,g-1]+N[t,8,g-1]+N[t,9,g-1]+N[t,10,g-1]
      # total population of low risk rural women from last generation, who will mix with rural men 
      N_lr_r_f[t,s,g-1] <- S[t,11,g-1]+S[t,12,g-1]+S[t,13,g-1]+S[t,14,g-1]+S[t,15,g-1]+N[t,11,g-1]+N[t,12,g-1]+N[t,13,g-1]+N[t,14,g-1]+N[t,15,g-1]
      # total population of urban men from last generation, who will mix with high risk women and low risk urban women
      N_u_m[t,s,g-1] <- S[t,16,g-1]+S[t,17,g-1]+S[t,18,g-1]+S[t,19,g-1]+S[t,20,g-1]+N[t,16,g-1]+N[t,17,g-1]+N[t,18,g-1]+N[t,19,g-1]+N[t,20,g-1]
      # total population of rural men from last generation, who will mix with high risk women and low risk rural women
      N_r_m[t,s,g-1] <- S[t,21,g-1]+S[t,22,g-1]+S[t,23,g-1]+S[t,24,g-1]+S[t,25,g-1]+N[t,21,g-1]+N[t,22,g-1]+N[t,23,g-1]+N[t,24,g-1]+N[t,25,g-1]
      #}
      # Growth rate in this generation
      # rho_2[t,s,g-1] <- rho[t,s,g]*(1+N[t,s,g]/S[t,s,g])
      
      # low risk urban women mixed with urban men
      if (s>=1&s<=5) 
      {
        for (j in 16:20) 
        {
          # force of infection low risk urban women
          #lambda_t_r[t,s,g] <-  lambda_t_r[t,s,g]+((N[t,j,g-1]+S[t,j,g-1])/N_u_m[t,s,g-1])*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*(1-(1-(omega_r[t,s,g]*beta[t,s,g]*(kappa[t,s,g]*p_r_v[t,j,g-1]+p_r[t,j,g-1])))^n_r[t,s,g])
          lambda_inf[t,s,g] <-  lambda_inf[t,s,g]*(1-((N[t,j,g-1]+S[t,j,g-1])/N_u_m[t,s,g-1])*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*beta[t,s,g]*((1-upsilon[t,s,g])*p_r_v[t,j,g-1]+p_r[t,j,g-1]))^((1-epsilon[t,s,g]*c_r[t,s,g])*n_r[t,s,g])
          
        } 
      }
      
      # high risk women mixed with urban men
      if (s>=6&s<=10) 
      {
        for (j in 16:20) 
        {
          # force of infection (high risk urban women)
          #lambda_t_r[t,s,g] <-  lambda_t_r[t,s,g]+((N[t,j,g-1]+S[t,j,g-1])/N_u_m[t,s,g-1])*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*(1-(1-(omega_r[t,s,g]*beta[t,s,g]*(kappa[t,s,g]*p_r_v[t,j,g-1]+p_r[t,j,g-1])))^n_r[t,s,g])
          lambda_inf[t,s,g] <-  lambda_inf[t,s,g]*(1-((N[t,j,g-1]+S[t,j,g-1])/N_u_m[t,s,g-1])*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*beta[t,s,g]*((1-upsilon[t,s,g])*p_r_v[t,j,g-1]+p_r[t,j,g-1]))^((1-epsilon[t,s,g]*c_r[t,s,g])*n_r[t,s,g])
          
        }
      }
      
      # low risk rural women mixed with rural men
      if (s>=11&s<=15) 
      {
        for (j in 21:25) {
          # force of infection (low risk rural women)
          #lambda_t_r[t,s,g] <-  lambda_t_r[t,s,g]+((N[t,j,g-1]+S[t,j,g-1])/N_u_m[t,s,g-1])*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*(1-(1-(omega_r[t,s,g]*beta[t,s,g]*(kappa[t,s,g]*p_r_v[t,j,g-1]+p_r[t,j,g-1])))^n_r[t,s,g])
          lambda_inf[t,s,g] <-  lambda_inf[t,s,g]*(1-((N[t,j,g-1]+S[t,j,g-1])/N_r_m[t,s,g-1])*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*beta[t,s,g]*((1-upsilon[t,s,g])*p_r_v[t,j,g-1]+p_r[t,j,g-1]))^((1-epsilon[t,s,g]*c_r[t,s,g])*n_r[t,s,g])
          
        }
      }
      
      # urban men mixed with low risk urban women and high risk women
      if (s>=16&s<=20) 
      {
        for (j in 1:10) 
        {
          # force of infection (urban men)
          #lambda_t_r[t,s,g] <-  lambda_t_r[t,s,g]+((N[t,j,g-1]+S[t,j,g-1])/N_u_m[t,s,g-1])*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*(1-(1-(omega_r[t,s,g]*beta[t,s,g]*(kappa[t,s,g]*p_r_v[t,j,g-1]+p_r[t,j,g-1])))^n_r[t,s,g])
          lambda_inf[t,s,g] <-  lambda_inf[t,s,g]*(1-((N[t,j,g-1]+S[t,j,g-1])/(N_lr_u_f[t,s,g-1]+N_hr_f[t,s,g-1]))*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*beta[t,s,g]*((1-upsilon[t,s,g])*p_r_v[t,j,g-1]+p_r[t,j,g-1]))^((1-epsilon[t,s,g]*c_r[t,s,g])*n_r[t,s,g])
          
        }
      }
      
      # rural men mixed with low risk rural women
      if (s>=21&s<=25) 
      {
        for (j in 11:15) 
        {
          # force of infection (rural men)
          #lambda_t_r[t,s,g] <-  lambda_t_r[t,s,g]+((N[t,j,g-1]+S[t,j,g-1])/N_u_m[t,s,g-1])*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*(1-(1-(omega_r[t,s,g]*beta[t,s,g]*(kappa[t,s,g]*p_r_v[t,j,g-1]+p_r[t,j,g-1])))^n_r[t,s,g])
          lambda_inf[t,s,g] <-  lambda_inf[t,s,g]*(1-((N[t,j,g-1]+S[t,j,g-1])/N_lr_r_f[t,s,g-1])*(N[t,j,g-1]/(N[t,j,g-1]+S[t,j,g-1]))*beta[t,s,g]*((1-upsilon[t,s,g])*p_r_v[t,j,g-1]+p_r[t,j,g-1]))^((1-epsilon[t,s,g]*c_r[t,s,g])*n_r[t,s,g])
          
        }
      }
      
      # Calculate new infections
      #NI[t,s,g] <- lambda_t_r[t,s,g]*S[t,s,g-1]
      lambda_t_r[t,s,g]<-(1-lambda_inf[t,s,g])
      NI[t,s,g] <- (1-lambda_inf[t,s,g])*S[t,s,g-1]
      
      # Difference equation for each compartments
      # Difference equation for susceptable population
      S[t,s,g] <- S[t,s,g-1]+rho[t,s,g]*S[t,s,g-1]-lambda_t_r[t,s,g]*S[t,s,g-1]
      
      # Difference equations for compartments X1 - X4 (undiagnosed individuals)
      
      X1[t,s,g] <- X1[t,s,g-1]+lambda_t_r[t,s,g]*S[t,s,g-1]-delta_1[t,s,g]*X1[t,s,g-1]-alpha_1[t,s,g]*X1[t,s,g-1]-mu_1[t,s,g]*X1[t,s,g-1]
      
      X2[t,s,g] <- X2[t,s,g-1]+delta_1[t,s,g]*X1[t,s,g-1]-delta_2[t,s,g]*X2[t,s,g-1]-alpha_2[t,s,g]*X2[t,s,g-1]-mu_2[t,s,g]*X2[t,s,g-1]
      
      X3[t,s,g] <- X3[t,s,g-1]+delta_2[t,s,g]*X2[t,s,g-1]-delta_3[t,s,g]*X3[t,s,g-1]-alpha_3[t,s,g]*X3[t,s,g-1]-mu_3[t,s,g]*X3[t,s,g-1]
      
      X4[t,s,g] <- X4[t,s,g-1]+delta_3[t,s,g]*X3[t,s,g-1]-alpha_4[t,s,g]*X4[t,s,g-1]-mu_4[t,s,g]*X4[t,s,g-1]
      
      
      
      # Difference equations for compartments X5 - X8 (diagnosed individuals)
      
      X5[t,s,g] <- X5[t,s,g-1]+alpha_1[t,s,g]*X1[t,s,g-1]-delta_1[t,s,g]*X5[t,s,g-1]-sigma_1[t,s,g]*X5[t,s,g-1]-mu_1[t,s,g]*X5[t,s,g-1]
      
      X6[t,s,g] <- X6[t,s,g-1]+alpha_2[t,s,g]*X2[t,s,g-1]+delta_1[t,s,g]*X5[t,s,g-1]-delta_2[t,s,g]*X6[t,s,g-1]-sigma_2[t,s,g]*X6[t,s,g-1]-mu_2[t,s,g]*X6[t,s,g-1]
      
      X7[t,s,g] <- X7[t,s,g-1]+alpha_3[t,s,g]*X3[t,s,g-1]+delta_2[t,s,g]*X6[t,s,g-1]-delta_3[t,s,g]*X7[t,s,g-1]-sigma_3[t,s,g]*X7[t,s,g-1]-mu_3[t,s,g]*X7[t,s,g-1]
      
      X8[t,s,g] <- X8[t,s,g-1]+alpha_4[t,s,g]*X4[t,s,g-1]+delta_3[t,s,g]*X7[t,s,g-1]-sigma_4[t,s,g]*X8[t,s,g-1]-mu_4[t,s,g]*X8[t,s,g-1]
      
      
      
      # Difference equations for compartments X9 - X12 (individuals linked to care)
      
      X9[t,s,g] <- X9[t,s,g-1]+sigma_1[t,s,g]*X5[t,s,g-1]-delta_4[t,s,g]*X9[t,s,g-1]-theta_1[t,s,g]*X9[t,s,g-1]-gamma_1[t,s,g]*X9[t,s,g-1]-mu_5[t,s,g]*X9[t,s,g-1]
      
      X10[t,s,g] <- X10[t,s,g-1]+sigma_2[t,s,g]*X6[t,s,g-1]+delta_4[t,s,g]*X9[t,s,g-1]-delta_5[t,s,g]*X10[t,s,g-1]-theta_2[t,s,g]*X10[t,s,g-1]-gamma_2[t,s,g]*X10[t,s,g-1]-mu_6[t,s,g]*X10[t,s,g-1]
      
      X11[t,s,g] <- X11[t,s,g-1]+sigma_3[t,s,g]*X7[t,s,g-1]+delta_5[t,s,g]*X10[t,s,g-1]-delta_6[t,s,g]*X11[t,s,g-1]-theta_3[t,s,g]*X11[t,s,g-1]-gamma_3[t,s,g]*X11[t,s,g-1]-mu_7[t,s,g]*X11[t,s,g-1]
      
      X12[t,s,g] <- X12[t,s,g-1]+sigma_4[t,s,g]*X8[t,s,g-1]+delta_6[t,s,g]*X11[t,s,g-1]-theta_4[t,s,g]*X12[t,s,g-1]-gamma_4[t,s,g]*X12[t,s,g-1]-mu_8[t,s,g]*X12[t,s,g-1]
      
      
      
      # Difference equations for compartments X13 - X16 (individuals lost from care)
      
      X13[t,s,g] <- X13[t,s,g-1]+gamma_1[t,s,g]*X9[t,s,g-1]+gamma_5[t,s,g]*X17[t,s,g-1]+gamma_9[t,s,g]*X21[t,s,g-1]-delta_1[t,s,g]*X13[t,s,g-1]-tau_1[t,s,g]*X13[t,s,g-1]-mu_1[t,s,g]*X13[t,s,g-1]
      
      X14[t,s,g] <- X14[t,s,g-1]+gamma_2[t,s,g]*X10[t,s,g-1]+gamma_6[t,s,g]*X18[t,s,g-1]+gamma_10[t,s,g]*X22[t,s,g-1]+delta_1[t,s,g]*X13[t,s,g-1]-tau_2[t,s,g]*X14[t,s,g-1]-delta_2[t,s,g]*X14[t,s,g-1]-mu_2[t,s,g]*X14[t,s,g-1]
      
      X15[t,s,g] <- X15[t,s,g-1]+gamma_3[t,s,g]*X11[t,s,g-1]+gamma_7[t,s,g]*X19[t,s,g-1]+gamma_11[t,s,g]*X23[t,s,g-1]+delta_2[t,s,g]*X14[t,s,g-1]-tau_3[t,s,g]*X15[t,s,g-1]-delta_3[t,s,g]*X15[t,s,g-1]-mu_3[t,s,g]*X15[t,s,g-1]
      
      X16[t,s,g] <- X16[t,s,g-1]+gamma_4[t,s,g]*X12[t,s,g-1]+gamma_8[t,s,g]*X20[t,s,g-1]+gamma_12[t,s,g]*X24[t,s,g-1]+delta_3[t,s,g]*X15[t,s,g-1]-tau_4[t,s,g]*X16[t,s,g-1]-mu_4[t,s,g]*X16[t,s,g-1]
      
      
      
      # Difference equations for compartments X17 - X20 (individuals who are viral suppressed)
      
      X17[t,s,g] <- X17[t,s,g-1]+theta_1[t,s,g]*X9[t,s,g-1]+tau_1[t,s,g]*X13[t,s,g-1]-gamma_5[t,s,g]*X17[t,s,g-1]-psi_1[t,s,g]*X17[t,s,g-1]-mu_9[t,s,g]*X17[t,s,g-1]
      
      X18[t,s,g] <- X18[t,s,g-1]+theta_2[t,s,g]*X10[t,s,g-1]+tau_2[t,s,g]*X14[t,s,g-1]-gamma_6[t,s,g]*X18[t,s,g-1]-psi_2[t,s,g]*X18[t,s,g-1]-mu_10[t,s,g]*X18[t,s,g-1]
      
      X19[t,s,g] <- X19[t,s,g-1]+theta_3[t,s,g]*X11[t,s,g-1]+tau_3[t,s,g]*X15[t,s,g-1]-gamma_7[t,s,g]*X19[t,s,g-1]-psi_3[t,s,g]*X19[t,s,g-1]-mu_11[t,s,g]*X19[t,s,g-1]
      
      X20[t,s,g] <- X20[t,s,g-1]+theta_4[t,s,g]*X12[t,s,g-1]+tau_4[t,s,g]*X16[t,s,g-1]-gamma_8[t,s,g]*X20[t,s,g-1]-psi_4[t,s,g]*X20[t,s,g-1]-mu_12[t,s,g]*X20[t,s,g-1]
      
      
      # Difference equations compartments X21 - X24 (individuals who are not viral suppressed)
      
      X21[t,s,g] <- X21[t,s,g-1]+psi_1[t,s,g]*X17[t,s,g-1]-gamma_9[t,s,g]*X21[t,s,g-1]-mu_13[t,s,g]*X21[t,s,g-1]
      
      X22[t,s,g] <- X22[t,s,g-1]+psi_2[t,s,g]*X18[t,s,g-1]-gamma_10[t,s,g]*X22[t,s,g-1]-mu_14[t,s,g]*X22[t,s,g-1]
      
      X23[t,s,g] <- X23[t,s,g-1]+psi_3[t,s,g]*X19[t,s,g-1]-gamma_11[t,s,g]*X23[t,s,g-1]-mu_15[t,s,g]*X23[t,s,g-1]
      
      X24[t,s,g] <- X24[t,s,g-1]+psi_4[t,s,g]*X20[t,s,g-1]-gamma_12[t,s,g]*X24[t,s,g-1]-mu_16[t,s,g]*X24[t,s,g-1]
      
      # Difference equations for compartment D (individuals who have died)
      
      D[t,s,g] <- D[t,s,g-1]+mu_1[t,s,g]*(X1[t,s,g-1]+X5[t,s,g-1]+X13[t,s,g-1])+mu_2[t,s,g]*(X2[t,s,g-1]+X6[t,s,g-1]+X14[t,s,g-1])+mu_3[t,s,g]*(X3[t,s,g-1]+X7[t,s,g-1]+X15[t,s,g-1])+ mu_4[t,s,g]*(X4[t,s,g-1]+X8[t,s,g-1]+X16[t,s,g-1])+mu_5[t,s,g]*X9[t,s,g-1]+mu_6[t,s,g]*X10[t,s,g-1]+mu_7[t,s,g]*X11[t,s,g-1]+mu_8[t,s,g]*X12[t,s,g-1]+mu_9[t,s,g]*X17[t,s,g-1]+mu_10[t,s,g]*X18[t,s,g-1]+mu_11[t,s,g]*X19[t,s,g-1]+mu_12[t,s,g]*X20[t,s,g-1]+mu_13[t,s,g]*X21[t,s,g-1]+mu_14[t,s,g]*X22[t,s,g-1]+mu_15[t,s,g]*X23[t,s,g-1]+mu_16[t,s,g]*X24[t,s,g-1]
      
      
      # Total population of HIV infected individuals
      N[t,s,g] <- X1[t,s,g]+X2[t,s,g]+X3[t,s,g]+X4[t,s,g]+X5[t,s,g]+X6[t,s,g]+X7[t,s,g]+X8[t,s,g]+X9[t,s,g]+X10[t,s,g]+X11[t,s,g]+X12[t,s,g]+X13[t,s,g]+X14[t,s,g]+X15[t,s,g]+X16[t,s,g]+X17[t,s,g]+X18[t,s,g]+X19[t,s,g]+X20[t,s,g]+X21[t,s,g]+X22[t,s,g]+X23[t,s,g]+X24[t,s,g]
      
      # Total population of individuals diagnosed with HIV
      Diagnosed[t,s,g] <- X5[t,s,g]+X6[t,s,g]+X7[t,s,g]+X8[t,s,g]+X9[t,s,g]+X10[t,s,g]+X11[t,s,g]+X12[t,s,g]+X13[t,s,g]+X14[t,s,g]+X15[t,s,g]+X16[t,s,g]+X17[t,s,g]+X18[t,s,g]+X19[t,s,g]+X20[t,s,g]+X21[t,s,g]+X22[t,s,g]+X23[t,s,g]+X24[t,s,g]#+
      #X1[t,s,g]*alpha_1[t,s,g]+X2[t,s,g]*alpha_2[t,s,g]+X3[t,s,g]*alpha_3[t,s,g]+X4[t,s,g]*alpha_4[t,s,g]
      
      # Total population of individuals newly diagnosed
      New_Diagnosed[t,s,g] <- X1[t,s,g]*alpha_1[t,s,g]+X2[t,s,g]*alpha_2[t,s,g]+X3[t,s,g]*alpha_3[t,s,g]+X4[t,s,g]*alpha_4[t,s,g]
      
      # Total population of individuals on ART
      On_ART[t,s,g] <- X17[t,s,g]+X18[t,s,g]+X19[t,s,g]+X20[t,s,g]+X21[t,s,g]+X22[t,s,g]+X23[t,s,g]+X24[t,s,g]
      
      # Total population of individuals who are not viral suppressed
      Z[t,s,g] <- X1[t,s,g]+X2[t,s,g]+X3[t,s,g]+X4[t,s,g]+X5[t,s,g]+X6[t,s,g]+X7[t,s,g]+X8[t,s,g]+X9[t,s,g]+X10[t,s,g]+X11[t,s,g]+X12[t,s,g]+X13[t,s,g]+X14[t,s,g]+X15[t,s,g]+X16[t,s,g]+X21[t,s,g]+X22[t,s,g]+X23[t,s,g]+X24[t,s,g]
      
      # Total population of individuals who are viral suppressed
      V[t,s,g] <- X17[t,s,g]+X18[t,s,g]+X19[t,s,g]+X20[t,s,g]
      
      if (N[t,s,g]==0) {
        
        p_r[t,s,g] <- 1
        p_r_v[t,s,g] <- 0
        
      } else {
        
        # Prevalence of HIV in the population based on individuals who are not viral supprressed
        p_r[t,s,g] <- Z[t,s,g]/N[t,s,g]
        # Prevalence of HIV in the population based on individuals who are viral supprressed
        p_r_v[t,s,g] <- V[t,s,g]/N[t,s,g]  
      }
    }
  }
}
#}
#}
# Model calibration_targets&GOF
# Calibration targets #
# Import calibration targets #
#files <- list.files(path="~/Dropbox/VCU_PhD_Year 3&4/Paper 1/R code/Calibration/", pattern = "csv")
#setwd("~/Dropbox/VCU_PhD_Year 3&4/Paper 1/R code/Calibration")
Target_all <- read.csv("Targets_00_15-54_new.csv")
Target_new <- read.csv("Targets_01.csv")

# Generate predicted outcomes
for (t in (1:trial_select)){  
  for (g in (1:generations)) {
    N_F_total[t,1,g] <- 0
    N_M_total[t,1,g] <- 0
    N_young_total[t,1,g] <- N[t,1,g]+N[t,6,g]+N[t,11,g]+N[t,16,g]+N[t,21,g]
    N_F_young_total[t,1,g] <- N[t,1,g]+N[t,6,g]+N[t,11,g]
    N_M_young_total[t,1,g] <- N[t,16,g]+N[t,21,g]
    N_adult_total[t,1,g] <- sum(N[t,1:3,g]+N[t,6:8,g]+N[t,11:13,g]+N[t,16:18,g]+N[t,21:23,g])
    NI_total[t,1,g] <- 0
    NI_F_total[t,1,g] <- 0
    NI_M_total[t,1,g] <- 0
    NI_young_total[t,1,g] <-  NI[t,1,g]+NI[t,6,g]+NI[t,11,g]+NI[t,16,g]+NI[t,21,g]
    NI_F_young_total[t,1,g] <- NI[t,1,g]+NI[t,6,g]+NI[t,11,g]
    NI_M_young_total[t,1,g] <- NI[t,16,g]+NI[t,21,g]
    NI_adult_total[t,1,g] <- sum(NI[t,1:3,g]+NI[t,6:8,g]+NI[t,11:13,g]+NI[t,16:18,g]+NI[t,21:23,g])
    N_total[t,1,g] <- 0
    Z_total[t,1,g] <- 0
    V_total[t,1,g] <- 0
    S_total[t,1,g] <- 0
    S_F_total[t,1,g] <- 0
    S_M_total[t,1,g] <- 0
    S_F_young_total[t,1,g] <- S[t,1,g]+S[t,6,g]+S[t,11,g]
    S_M_young_total[t,1,g] <- S[t,16,g]+S[t,21,g]
    
    Diagnosed_total[t,1,g] <- 0
    Diagnosed_F_total[t,1,g] <- 0
    Diagnosed_M_total[t,1,g] <- 0
    Diagnosed_F_young_total[t,1,g] <- Diagnosed[t,1,g]+Diagnosed[t,6,g]+Diagnosed[t,11,g]
    Diagnosed_M_young_total[t,1,g] <- Diagnosed[t,16,g]+Diagnosed[t,21,g]
    
    
    New_Diagnosed_total[t,1,g] <- 0
    
    Suppressed_total[t,1,g] <- 0
    Suppressed_total_F[t,1,g] <- 0
    Suppressed_total_M[t,1,g] <- 0
    Suppressed_total_F_young[t,1,g] <- V[t,1,g]+V[t,6,g]+V[t,11,g]
    Suppressed_total_M_young[t,1,g] <- V[t,16,g]+V[t,21,g]
    
    On_ART_total[t,1,g] <- 0
    On_ART_total_F[t,1,g] <- 0
    On_ART_total_M[t,1,g] <- 0
    On_ART_total_F_young[t,1,g] <- On_ART[t,1,g]+On_ART[t,6,g]+On_ART[t,11,g]
    On_ART_total_M_young[t,1,g] <- On_ART[t,16,g]+On_ART[t,21,g]
    
    
    Linked_total[t,1,g] <- 0
    for(s in c(1:3, 6:8, 11:13)){ # for female subpopulations
      N_F_total[t,1,g] <- N_F_total[t,1,g]+N[t,s,g]
      NI_F_total[t,1,g] <-  NI_F_total[t,1,g]+NI[t,s,g]
      S_F_total[t,1,g] = S_F_total[t,1,g]+S[t,s,g]
      On_ART_total_F[t,1,g]=On_ART_total_F[t,1,g]+On_ART[t,s,g]
      Suppressed_total_F[t,1,g]=Suppressed_total_F[t,1,g]+V[t,s,g]
      Diagnosed_F_total[t,1,g] = Diagnosed_F_total[t,1,g]+Diagnosed[t,s,g]
    }
    for(s in c(16:18, 21:23)){ # for male subpopulations
      N_M_total[t,1,g] <- N_M_total[t,1,g]+N[t,s,g]
      NI_M_total[t,1,g] <- NI_M_total[t,1,g]+NI[t,s,g]
      S_M_total[t,1,g]=S_M_total[t,1,g]+S[t,s,g]
      On_ART_total_M[t,1,g]=On_ART_total_M[t,1,g]+On_ART[t,s,g]
      Suppressed_total_M [t,1,g]=Suppressed_total_M[t,1,g]+V[t,s,g] 
      Diagnosed_M_total[t,1,g] = Diagnosed_M_total[t,1,g]+Diagnosed[t,s,g]
    }
    
    for (s in c(1:3, 6:8, 11:13, 16:18, 21:23)){
      # Total population    
      
      if(is.na(NI[t,s,g])){
        NI[t,s,g]=0
      }
      NI_total[t,1,g]=NI_total[t,1,g]+NI[t,s,g]
      N_total[t,1,g]=N_total[t,1,g]+N[t,s,g]
      S_total[t,1,g]=S_total[t,1,g]+S[t,s,g]
      Diagnosed_total[t,1,g]=Diagnosed_total[t,1,g]+Diagnosed[t,s,g]
      New_Diagnosed_total[t,1,g]=New_Diagnosed_total[t,1,g]+New_Diagnosed[t,s,g]
      Suppressed_total[t,1,g]=Suppressed_total[t,1,g]+V[t,s,g]
      Linked_total[t,1,g]=Linked_total[t,1,g]+Linked[t,s,g]
      On_ART_total[t,1,g]=On_ART_total[t,1,g]+On_ART[t,s,g]
    }
    for (s in c(1:3, 6:8, 11:13, 16:18, 21:23)){
      On_ART_grand_total[t,1,g]=On_ART_grand_total[t,1,g]+On_ART[t,s,g]
    }
  }
}  

N_F_total_all <- array(0, dim=c(trial_select,1,generations))
N_M_total_all<- array(0, dim=c(trial_select,1,generations))
N_adult_total_all<- array(0, dim=c(trial_select,1,generations))
S_total_all<- array(0, dim=c(trial_select,1,generations))
S_F_total_all<- array(0, dim=c(trial_select,1,generations))
S_M_total_all<- array(0, dim=c(trial_select,1,generations))
Suppressed_total_all<- array(0, dim=c(trial_select,1,generations))
Suppressed_total_F_all<- array(0, dim=c(trial_select,1,generations))
Suppressed_total_M_all<- array(0, dim=c(trial_select,1,generations))
N_total_all<- array(0, dim=c(trial_select,1,generations))

for (t in (1:trial_select)){  
  for (g in (1:generations)) {
    N_F_total_all[t,1,g] <- 0
    N_M_total_all[t,1,g] <- 0
    N_total_all[t,1,g] <- 0
    S_total_all[t,1,g] <- 0
    S_F_total_all[t,1,g] <- 0
    S_M_total_all[t,1,g] <- 0
    
    Suppressed_total_all[t,1,g] <- 0
    Suppressed_total_F_all[t,1,g] <- 0
    Suppressed_total_M_all[t,1,g] <- 0
    
    for(s in c(1:15)){ # for female subpopulations
      N_F_total_all[t,1,g] <- N_F_total_all[t,1,g]+N[t,s,g]
      S_F_total_all[t,1,g] = S_F_total_all[t,1,g]+S[t,s,g]
      Suppressed_total_F_all[t,1,g]=Suppressed_total_F_all[t,1,g]+V[t,s,g]
    }
    for(s in c(16:20, 21:25)){ # for male subpopulations
      N_M_total_all[t,1,g] <- N_M_total_all[t,1,g]+N[t,s,g]
      S_M_total_all[t,1,g]=S_M_total_all[t,1,g]+S[t,s,g]
      Suppressed_total_M_all [t,1,g]=Suppressed_total_M_all[t,1,g]+V[t,s,g] 
    }
    
    for (s in c(1:25)){
      # Total population    
      N_total_all[t,1,g]=N_total_all[t,1,g]+N[t,s,g]
      S_total_all[t,1,g]=S_total_all[t,1,g]+S[t,s,g]
      Suppressed_total_all[t,1,g]=Suppressed_total_all[t,1,g]+V[t,s,g]
    }
  }
}  

# Save the generated parameters
# save.image("Data_trial500.RData")

pd_prev_total_all <- numeric(trial_select)
pd_prev_total_F_all <- numeric(trial_select)
pd_prev_total_M_all <- numeric(trial_select)
pd_vs_total_all <- numeric(trial_select)
pd_vs_total_F_all <- numeric(trial_select)
pd_vs_total_M_all <- numeric(trial_select)

for (t in (1:trial_select)) {
  # number on ART
  pd_n_art_total[t]<- 0
  
  # percent diagnosed
  pd_p_diagnosed[t]<-0
  pd_p_diagnosed_F[t]<-0
  pd_p_diagnosed_M[t]<-0
  pd_p_diagnosed_F_young[t]<-0
  pd_p_diagnosed_M_young[t]<-0  
  
  # Percent on Treatment
  pd_art_p[t] <- 0
  pd_art_p_F[t] <- 0
  pd_art_p_M[t] <- 0
  pd_art_p_F_young[t] <- 0
  pd_art_p_M_young[t] <- 0
  
  # Percent virally suppressed
  pd_vs_total[t] <- 0
  pd_vs_total_F[t] <- 0
  pd_vs_total_M[t] <- 0
  pd_vs_total_M_young[t] <- 0
  pd_vs_total_F_young[t] <- 0
  pd_vs_total_all[t] <- 0
  pd_vs_total_F_all[t] <- 0
  pd_vs_total_M_all[t] <- 0
  
  
  # Percent virally suppressed on ART
  pd_vs_ART_total[t] <- 0
  pd_vs_ART_total_F[t] <- 0
  pd_vs_ART_total_M[t] <- 0
  pd_vs_ART_total_M_young[t] <- 0
  pd_vs_ART_total_F_young[t] <- 0
  
  
  # HIV Prevalence
  pd_prev_total[t] <- 0
  pd_prev_total_F[t] <- 0
  pd_prev_total_M[t] <- 0
  pd_prev_total_young[t] <- 0
  pd_prev_total_F_young[t] <- 0
  pd_prev_total_M_young[t] <- 0
  pd_prev_total_all[t] <- 0
  pd_prev_total_F_all[t] <- 0
  pd_prev_total_M_all[t] <- 0
  
  # HIV Incidence
  pd_incidence_total[t] <- 0
  pd_incidence_F_total[t] <- 0
  pd_incidence_M_total[t] <- 0
  pd_incidence_F_young_total[t] <- 0
  pd_incidence_M_young_total[t] <- 0
  
  
  

    Prevalence_total[t,1,] <- (N_total[t,1,])/(N_total[t,1,]+S_total[t,1,])
    Prevalence_total_y[t,] <- Prevalence_total[t,1,][seq(1,length(Prevalence_total[t,1,]),12)]
    for (y in (1:length(Prevalence_total_y[t,]))) {
    # Prevalence
    Target_prev[t,y]<- Target_all[y,"Prevalence"]
    if (is.na(Target_prev[t,y])||is.na(Prevalence_total_y[t,y]))
    {
      pd_prev_total[t] <- pd_prev_total[t]
    } else {
      pd_prev_total[t] <- pd_prev_total[t]+abs(Prevalence_total_y[t,y]-Target_prev[t,y])/(Target_prev[t,y])
    }
  }
  
    
    
    # HIV prevalence female
      Prevalence_F_total[t,1,] <- (N_F_total[t,1,])/(N_F_total[t,1,]+S_F_total[t,1,])
      Prevalence_F_total_y[t,] <- Prevalence_F_total[t,1,][seq(1,length(Prevalence_F_total[t,1,]),12)]
  for (y in (1:length(Prevalence_total_y[t,]))) {
    # Prevalence Female 
    Target_prev_F[t,y]<- Target_all[y,"Prevalence_F"]
    if (is.na(Target_prev_F[t,y])||is.na(Prevalence_F_total_y[t,y]))
    {
      pd_prev_total_F[t] <- pd_prev_total_F[t]
    } else {
      pd_prev_total_F[t] <- pd_prev_total_F[t]+abs(Prevalence_F_total_y[t,y]-Target_prev_F[t,y])/(Target_prev_F[t,y])
    }
  }
      
  # HIV prevalence male
  Prevalence_M_total[t,1,] <- (N_M_total[t,1,])/(N_M_total[t,1,]+S_M_total[t,1,])
  Prevalence_M_total_y[t,] <- Prevalence_M_total[t,1,][seq(1,length(Prevalence_M_total[t,1,]),12)]
  for (y in (1:length(Prevalence_total_y[t,]))) {
    # Prevalence male
    Target_prev_M[t,y]<- Target_all[y,"Prevalence_M"]
    if (is.na(Target_prev_M[t,y])||is.na(Prevalence_M_total_y[t,y]))
    {
      pd_prev_total_M[t] <- pd_prev_total_M[t]
    } else {
      pd_prev_total_M[t] <- pd_prev_total_M[t]+abs(Prevalence_M_total_y[t,y]-Target_prev_M[t,y])/(Target_prev_M[t,y])
    }
  }
  
  # HIV prevalence young female
  Prevalence_F_young_total[t,1,] <- (N_F_young_total[t,1,])/(N_F_young_total[t,1,]+S_F_young_total[t,1,])
  Prevalence_F_young_total_y[t,] <- Prevalence_F_young_total[t,1,][seq(1,length(Prevalence_F_young_total[t,1,]),12)]
  for (y in (1:length(Prevalence_total_y[t,]))) {
    # Prevalence F_youngemale 
    Target_prev_F_young[t,y]<- Target_all[y,"Prevalence_young_F"]
    if (is.na(Target_prev_F_young[t,y])||is.na(Prevalence_F_young_total_y[t,y]))
    {
      pd_prev_total_F_young[t] <- pd_prev_total_F_young[t]
    } else {
      pd_prev_total_F_young[t] <- pd_prev_total_F_young[t]+abs(Prevalence_F_young_total_y[t,y]-Target_prev_F_young[t,y])/(Target_prev_F_young[t,y])
    }
  }
  
  # HIV prevalence young male
  Prevalence_M_young_total[t,1,] <- (N_M_young_total[t,1,])/(N_M_young_total[t,1,]+S_M_young_total[t,1,])
  Prevalence_M_young_total_y[t,] <- Prevalence_M_young_total[t,1,][seq(1,length(Prevalence_M_young_total[t,1,]),12)]
  for (y in (1:length(Prevalence_total_y[t,]))) {
    # Prevalence male
    Target_prev_M_young[t,y]<- Target_all[y,"Prevalence_young_M"]
    if (is.na(Target_prev_M_young[t,y])||is.na(Prevalence_M_young_total_y[t,y]))
    {
      pd_prev_total_M_young[t] <- pd_prev_total_M_young[t]
    } else {
      pd_prev_total_M_young[t] <- pd_prev_total_M_young[t]+abs(Prevalence_M_young_total_y[t,y]-Target_prev_M_young[t,y])/(Target_prev_M_young[t,y])
    }
  }
  
  Prevalence_total[t,1,] <- (N_total_all[t,1,])/(N_total_all[t,1,]+S_total_all[t,1,])
  Prevalence_total_y[t,] <- Prevalence_total[t,1,][seq(1,length(Prevalence_total[t,1,]),12)]
  for (y in (1:length(Prevalence_total_y[t,]))) {
    # Prevalence
    Target_prev[t,y]<- Target_new[y,"Prevalence"]
    if (is.na(Target_prev[t,y])||is.na(Prevalence_total_y[t,y]))
    {
      pd_prev_total_all[t] <- pd_prev_total_all[t]
    } else {
      pd_prev_total_all[t] <- pd_prev_total_all[t]+abs(Prevalence_total_y[t,y]-Target_prev[t,y])/(Target_prev[t,y])
    }
  }
  
    Prevalence_F_total[t,1,] <- (N_F_total_all[t,1,])/(N_F_total_all[t,1,]+S_F_total_all[t,1,])
    Prevalence_F_total_y[t,] <- Prevalence_F_total[t,1,][seq(1,length(Prevalence_F_total[t,1,]),12)]
  for (y in (1:length(Prevalence_total_y[t,]))) {
    # Prevalence Female 
    Target_prev_F[t,y]<- Target_new[y,"Prevalence_F"]
    if (is.na(Target_prev_F[t,y])||is.na(Prevalence_F_total_y[t,y]))
    {
      pd_prev_total_F_all[t] <- pd_prev_total_F_all[t]
    } else {
      pd_prev_total_F_all[t] <- pd_prev_total_F_all[t]+abs(Prevalence_F_total_y[t,y]-Target_prev_F[t,y])/(Target_prev_F[t,y])
    }
  }
  
  
    Prevalence_M_total[t,1,] <- (N_M_total_all[t,1,])/(N_M_total_all[t,1,]+S_M_total_all[t,1,])
    Prevalence_M_total_y[t,] <- Prevalence_M_total[t,1,][seq(1,length(Prevalence_M_total[t,1,]),12)]
  for (y in (1:length(Prevalence_total_y[t,]))) {
    # Prevalence male
    Target_prev_M[t,y]<- Target_new[y,"Prevalence_M"]
    if (is.na(Target_prev_M[t,y])||is.na(Prevalence_M_total_y[t,y]))
    {
      pd_prev_total_M_all[t] <- pd_prev_total_M_all[t]
    } else {
      pd_prev_total_M_all[t] <- pd_prev_total_M_all[t]+abs(Prevalence_M_total_y[t,y]-Target_prev_M[t,y])/(Target_prev_M[t,y])
    }
  }
  
    # HIV incidence
      NI_total_y[t,] <- rollapply(NI_total[t,1,],12,sum,by=12)
      S_total_y[t,] <- S_total[t,1,][seq(1, length(S_total[t,1,]), 12)]
      N_total_y[t,] <- N_total[t,1,][seq(1, length(N_total[t,1,]), 12)]
      Incidence_total_y[t,] <- NI_total_y[t,]/S_total_y[t,]
  for (y in (1:length(Incidence_total_y[t,]))) {
    #Incidence
    Target_incidence[t,y]<- Target_all[y,"Incidence"]
    if (is.na(Target_incidence[t,y])||is.na(Incidence_total_y[t,y]))
    {
      pd_incidence_total[t] <- pd_incidence_total[t]
    } else {
      pd_incidence_total[t] <- pd_incidence_total[t]+abs(Incidence_total_y[t,y]-Target_incidence[t,y])/ (Target_incidence[t,y])
    }
  }
  
      
      # HIV incidence female
        NI_F_total_y[t,] <- rollapply(NI_F_total[t,1,],12,sum,by=12)
        S_F_total_y[t,] <- S_F_total[t,1,][seq(1, length(S_F_total[t,1,]), 12)]
        N_F_total_y[t,] <- N_F_total[t,1,][seq(1, length(N_F_total[t,1,]), 12)]
        Incidence_F_total_y[t,] <- NI_F_total_y[t,]/S_F_total_y[t,]
  for (y in (1:length(Incidence_F_total_y[t,]))) {
    #Incidence female
    Target_incidence_F[t,y]<- Target_all[y,"Incidence_F"]
    if (is.na(Target_incidence_F[t,y])||is.na(Incidence_F_total_y[t,y]))
    {
      pd_incidence_F_total[t] <- pd_incidence_F_total[t]
    } else {
      pd_incidence_F_total[t] <- pd_incidence_F_total[t]+abs(Incidence_F_total_y[t,y]-Target_incidence_F[t,y])/ (Target_incidence_F[t,y])
    }
  }
        
        
   # HIV incidence male
   NI_M_total_y[t,] <- rollapply(NI_M_total[t,1,],12,sum,by=12)
   S_M_total_y[t,] <- S_M_total[t,1,][seq(1, length(S_M_total[t,1,]), 12)]
   N_M_total_y[t,] <- N_M_total[t,1,][seq(1, length(N_M_total[t,1,]), 12)]
   Incidence_M_total_y[t,] <- NI_M_total_y[t,]/S_M_total_y[t,]
  for (y in (1:length(Incidence_M_total_y[t,]))) {
    #Incidence male
    Target_incidence_M[t,y]<- Target_all[y,"Incidence_M"]
    if (is.na(Target_incidence_M[t,y])||is.na(Incidence_M_total_y[t,y]))
    {
      pd_incidence_M_total[t] <- pd_incidence_M_total[t]
    } else {
      pd_incidence_M_total[t] <- pd_incidence_M_total[t]+abs(Incidence_M_total_y[t,y]-Target_incidence_M[t,y])/ (Target_incidence_M[t,y])
    }
  }
   
   
   # Number on ART 
     On_ART_grand_total_y[t,] <- On_ART_grand_total[t,1,][seq(1,length(On_ART_grand_total[t,1,]),12)]
  for (y in (1:length(P_suppressed_total_y[t,]))) {
    # number on ART
    Target_ART[t,y]<- Target_all[y,"N_ART"]
    if (is.na(Target_ART[t,y])||is.na(On_ART_grand_total_y[t,y]))
    {
      pd_n_art_total[t] <- pd_n_art_total[t]
    } else {
      pd_n_art_total[t] <- pd_n_art_total[t]+abs(On_ART_grand_total_y[t,y]-Target_ART[t,y])/ (Target_ART[t,y])
    }
  }
  
     
  # Percent on ART 
  Diagnosed_total_y[t,] <- Diagnosed_total[t,1,][seq(1, length(Diagnosed_total[t,1,]), 12)]
  On_ART_p_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]/Diagnosed_total_y[t,]*100
  for (y in (1:length(P_suppressed_total_y[t,]))) {
    # percent on ART
    P_Target_ART[t,y]<- Target_new[y,"P_ART_All"]
    if (is.na(P_Target_ART[t,y])||is.na(On_ART_p_y[t,y]))
    {
      pd_art_p[t] <- pd_art_p[t]
    } else {
      pd_art_p[t] <- pd_art_p[t]+abs(On_ART_p_y[t,y]-P_Target_ART[t,y])/ (P_Target_ART[t,y])
    }
  }
  
       
  # Percent on ART female
  Diagnosed_F_total_y[t,] <- Diagnosed_F_total[t,1,][seq(1, length(Diagnosed_F_total[t,1,]), 12)]
  On_ART_p_F_y[t,] <- On_ART_total_F[t,1,][seq(1,length(On_ART_total_F[t,1,]),12)]/Diagnosed_F_total_y[t,]*100
  for (y in (1:length(P_suppressed_total_y[t,]))) {
    # percent on ART female
    P_Target_ART[t,y]<- Target_new[y,"P_ART_All_F"]
    if (is.na(P_Target_ART[t,y])||is.na(On_ART_p_F_y[t,y]))
    {
      pd_art_p_F[t] <- pd_art_p_F[t]
    } else {
      pd_art_p_F[t] <- pd_art_p_F[t]+abs(On_ART_p_F_y[t,y]-P_Target_ART[t,y])/ (P_Target_ART[t,y])
    }
  }
  
         
  # Percent on ART male
  Diagnosed_M_total_y[t,] <- Diagnosed_M_total[t,1,][seq(1, length(Diagnosed_M_total[t,1,]), 12)]
  On_ART_p_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]/Diagnosed_M_total_y[t,]*100
  for (y in (1:length(P_suppressed_total_y[t,]))) {
    # percent on ART male
    P_Target_ART[t,y]<- Target_new[y,"P_ART_All_M"]
    if (is.na(P_Target_ART[t,y])||is.na(On_ART_p_M_y[t,y]))
    {
      pd_art_p_M[t] <- pd_art_p_M[t]
    } else {
      pd_art_p_M[t] <- pd_art_p_M[t]+abs(On_ART_p_M_y[t,y]-P_Target_ART[t,y])/ (P_Target_ART[t,y])
    }
  }
  
  # Percent viral suppression (among ART)
    On_ART_total_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
    Suppressed_total_y[t,] <- Suppressed_total[t,1,][seq(1,length(Suppressed_total[t,1,]),12)]
    P_suppressed_art_total_y[t,] <- Suppressed_total_y[t,]/On_ART_total_y[t,]*100
  for (y in (1:length(P_suppressed_art_total_y[t,]))) {
    # % VS among linked
    Target_VS_art[t,y]<- Target_all[y,"Percent_VS_on_ART_total"]
    if (is.na(Target_VS_art[t,y])||is.na(P_suppressed_art_total_y[t,y]))
    {
      pd_vs_ART_total[t] <- pd_vs_ART_total[t]
    } else {
      pd_vs_ART_total[t] <- pd_vs_ART_total[t]+abs(P_suppressed_art_total_y[t,y]-Target_VS_art[t,y])/ (Target_VS_art[t,y])
    }
  }
  
    
    # Percent viral suppression female (among ART)
      On_ART_total_F_y[t,] <- On_ART_total_F[t,1,][seq(1,length(On_ART_total_F[t,1,]),12)]
      Suppressed_total_F_y[t,] <- Suppressed_total_F[t,1,][seq(1,length(Suppressed_total_F[t,1,]),12)]
      P_suppressed_art_F_total_y[t,] <- Suppressed_total_F_y[t,]/On_ART_total_F_y[t,]*100
  for (y in (1:length(P_suppressed_art_total_y[t,]))) {
    # % VS among linked females
    Target_VS_art_F[t,y]<- Target_all[y,"Percent_VS_on_ART_Female"]
    if (is.na(Target_VS_art_F[t,y])||is.na(P_suppressed_art_F_total_y[t,y]))
    {
      pd_vs_ART_total_F[t] <- pd_vs_ART_total_F[t]
    } else {
      pd_vs_ART_total_F[t] <- pd_vs_ART_total_F[t]+abs(P_suppressed_art_F_total_y[t,y]-Target_VS_art_F[t,y])/ (Target_VS_art_F[t,y])
    }
  }
  
      # Percent viral suppression male (among ART)
        On_ART_total_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]
        Suppressed_total_M_y[t,] <- Suppressed_total_M[t,1,][seq(1,length(Suppressed_total_M[t,1,]),12)]
        P_suppressed_art_M_total_y[t,] <- Suppressed_total_M_y[t,]/On_ART_total_M_y[t,]*100
  for (y in (1:length(P_suppressed_art_total_y[t,]))) {
    # % VS among linked males
    Target_VS_art_M[t,y]<- Target_all[y,"Percent_VS_on_ART_Male"]
    if (is.na(Target_VS_art_M[t,y])||is.na(P_suppressed_art_M_total_y[t,y]))
    {
      pd_vs_ART_total_M[t] <- pd_vs_ART_total_M[t]
    } else {
      pd_vs_ART_total_M[t] <- pd_vs_ART_total_M[t]+abs(P_suppressed_art_M_total_y[t,y]-Target_VS_art_M[t,y])/ (Target_VS_art_M[t,y])
    }
  }
  
    # Percent viral suppression (among ALL)
      Suppressed_total_y[t,] <- Suppressed_total[t,1,][seq(1,length(Suppressed_total[t,1,]),12)]
      P_suppressed_total_y[t,] <- Suppressed_total_y[t,]/N_total_y[t,]*100
    
    for (y in (1:length(P_suppressed_art_total_y[t,]))) {
    # % VS linked
    Target_VS[t,y]<- Target_all[y,"Percent_VS_total"]
    if (is.na(Target_VS[t,y])||is.na(P_suppressed_total_y[t,y]))
    {
      pd_vs_total[t] <- pd_vs_total[t]
    } else {
      pd_vs_total[t] <- pd_vs_total[t]+abs(P_suppressed_total_y[t,y]-Target_VS[t,y])/ (Target_VS[t,y])
    }
  }
  
      # Percent viral suppression female
        Suppressed_total_F_y[t,] <- Suppressed_total_F[t,1,][seq(1,length(Suppressed_total_F[t,1,]),12)]
        P_suppressed_F_total_y[t,] <- Suppressed_total_F_y[t,]/N_F_total_y[t,]*100
  for (y in (1:length(P_suppressed_art_total_y[t,]))) {
    # % VS among females
    Target_VS_F[t,y]<- Target_all[y,"Percent_VS_Female"]
    if (is.na(Target_VS_F[t,y])||is.na(P_suppressed_F_total_y[t,y]))
    {
      pd_vs_total_F[t] <- pd_vs_total_F[t]
    } else {
      pd_vs_total_F[t] <- pd_vs_total_F[t]+abs(P_suppressed_F_total_y[t,y]-Target_VS_F[t,y])/ (Target_VS_F[t,y])
    }
  }
  
  # Percent viral suppression male
  Suppressed_total_M_y[t,] <- Suppressed_total_M[t,1,][seq(1,length(Suppressed_total_M[t,1,]),12)]
  P_suppressed_M_total_y[t,] <- Suppressed_total_M_y[t,]/N_M_total_y[t,]*100
  for (y in (1:length(P_suppressed_art_total_y[t,]))) {
    # % VS among males
    Target_VS_M[t,y]<- Target_all[y,"Percent_VS_Male"]
    if (is.na(Target_VS_M[t,y])||is.na(P_suppressed_M_total_y[t,y]))
    {
      pd_vs_total_M[t] <- pd_vs_total_M[t]
    } else {
      pd_vs_total_M[t] <- pd_vs_total_M[t]+abs(P_suppressed_M_total_y[t,y]-Target_VS_M[t,y])/ (Target_VS_M[t,y])
    }
  }
  

  N_total_y[t,] <- N_total_all[t,1,][seq(1, length(N_total_all[t,1,]), 12)]
  Suppressed_total_y[t,] <- Suppressed_total_all[t,1,][seq(1,length(Suppressed_total_all[t,1,]),12)]
  P_suppressed_total_y[t,] <- Suppressed_total_y[t,]/N_total_y[t,]*100
  for (y in (1:length(P_suppressed_art_total_y[t,]))) {
    # % VS linked
    Target_VS[t,y]<- Target_new[y,"Percent_VS_total"]
    if (is.na(Target_VS[t,y])||is.na(P_suppressed_total_y[t,y]))
    {
      pd_vs_total_all[t] <- pd_vs_total_all[t]
    } else {
      pd_vs_total_all[t] <- pd_vs_total_all[t]+abs(P_suppressed_total_y[t,y]-Target_VS[t,y])/ (Target_VS[t,y])
    }
  }
  

  N_F_total_y[t,] <- N_F_total_all[t,1,][seq(1, length(N_F_total_all[t,1,]), 12)]
  Suppressed_total_F_y[t,] <- Suppressed_total_F_all[t,1,][seq(1,length(Suppressed_total_F_all[t,1,]),12)]
  P_suppressed_F_total_y[t,] <- Suppressed_total_F_y[t,]/N_F_total_y[t,]*100
  for (y in (1:length(P_suppressed_art_total_y[t,]))) {
    # % VS among females
    Target_VS_F[t,y]<- Target_new[y,"Percent_VS_Female"]
    if (is.na(Target_VS_F[t,y])||is.na(P_suppressed_F_total_y[t,y]))
    {
      pd_vs_total_F_all[t] <- pd_vs_total_F_all[t]
    } else {
      pd_vs_total_F_all[t] <- pd_vs_total_F_all[t]+abs(P_suppressed_F_total_y[t,y]-Target_VS_F[t,y])/ (Target_VS_F[t,y])
    }
  }
  

  N_M_total_y[t,] <- N_M_total_all[t,1,][seq(1, length(N_M_total_all[t,1,]), 12)]
  Suppressed_total_M_y[t,] <- Suppressed_total_M_all[t,1,][seq(1,length(Suppressed_total_M_all[t,1,]),12)]
  P_suppressed_M_total_y[t,] <- Suppressed_total_M_y[t,]/N_M_total_y[t,]*100
  for (y in (1:length(P_suppressed_art_total_y[t,]))) {
    # % VS among males
    Target_VS_M[t,y]<- Target_new[y,"Percent_VS_Male"]
    if (is.na(Target_VS_M[t,y])||is.na(P_suppressed_M_total_y[t,y]))
    {
      pd_vs_total_M_all[t] <- pd_vs_total_M_all[t]
    } else {
      pd_vs_total_M_all[t] <- pd_vs_total_M_all[t]+abs(P_suppressed_M_total_y[t,y]-Target_VS_M[t,y])/ (Target_VS_M[t,y])
    }
  }
  
  pd_total[t]= pd_n_art_total[t]+
    pd_art_p[t] + pd_art_p_F[t] + pd_art_p_M[t]+
    pd_vs_total[t] + pd_vs_total_F[t] + pd_vs_total_M[t] + 
    pd_vs_ART_total[t] + pd_vs_ART_total_F[t] + pd_vs_ART_total_M[t] + 
    pd_prev_total[t] + pd_prev_total_F[t] + pd_prev_total_M[t] + 
    pd_incidence_total[t]+pd_incidence_F_total[t]+pd_incidence_M_total[t]+
    pd_prev_total_all[t] + pd_prev_total_F_all[t] + pd_prev_total_M_all[t]+
    pd_vs_total_all[t] + pd_vs_total_F_all[t] + pd_vs_total_M_all[t]+
    pd_prev_total_F_young[t] + pd_prev_total_M_young[t]
  
}

pdf(file="Calibration_top50_15-45_0506.pdf")

# Assign a vector to store the year variables
year <- seq(from=2004, to=2030, by=1)



for (t in 1:trials){
  Prevalence_total[t, ] <- N_total/(N_total+S_total)
  Prevalence_total_y[t,] <- Prevalence_total[t, seq(1,length(Prevalence_total[t,1,]),12)]}
df_prev <- as.data.frame(Prevalence_total_y)
colnames(df_prev) <- c(year)
df_prev_long <- df_prev %>%
  pivot_longer(
    cols=1:27,
    names_to="Year",
    values_to="Prev"
  )
ggplot()+
  geom_point(data=df_prev_long, mapping=aes(x=Year,y=Prev), color = "darkgreen", size = 1, alpha = 1)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Prevalence,ymin=Prevalence_LB,ymax=Prevalence_UB),color="red", size = 0.5)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=DHS_Prevalence,ymin=DHS_Prevalence_LB,ymax=DHS_Prevalence_UB),color="purple", size = 0.5)+
  labs(title="HIV prevalence RPHIA (red) DHS (purple)")

# HIV prevalence female
for (t in 1:trial_select){
  Prevalence_F_total[t,1,] <- (N_F_total[t,1,])/(N_F_total[t,1,]+S_F_total[t,1,])
  Prevalence_F_total_y[t,] <- Prevalence_F_total[t,1,][seq(1,length(Prevalence_F_total[t,1,]),12)]}
df_prev_F <- as.data.frame(Prevalence_F_total_y)
colnames(df_prev_F) <- c(year)
df_prev_F_long <- df_prev_F %>%
  pivot_longer(
    cols=1:27,
    names_to="Year",
    values_to="Prev"
  )
ggplot()+
  geom_point(data=df_prev_F_long, mapping=aes(x=Year,y=Prev), color = "darkgreen", size = 1, alpha = 1)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Prevalence_F,ymin=Prevalence_LB_F,ymax=Prevalence_UB_F),color="red", size = 0.5)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=DHS_Prevalence_F,ymin=DHS_Prevalence_LB_F,ymax=DHS_Prevalence_UB_F),color="purple", size = 0.5)+
  labs(title="HIV prevalence females RPHIA (red) DHS (purple)")

# HIV prevalence male
for (t in 1:trial_select){
  Prevalence_M_total[t,1,] <- (N_M_total[t,1,])/(N_M_total[t,1,]+S_F_total[t,1,])
  Prevalence_M_total_y[t,] <- Prevalence_M_total[t,1,][seq(1,length(Prevalence_M_total[t,1,]),12)]}
df_prev_M <- as.data.frame(Prevalence_M_total_y)
colnames(df_prev_M) <- c(year)
df_prev_M_long <- df_prev_M %>%
  pivot_longer(
    cols=1:27,
    names_to="Year",
    values_to="Prev"
  )
ggplot()+
  geom_point(data=df_prev_M_long, mapping=aes(x=Year,y=Prev), color = "darkgreen", size = 1, alpha = 1)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Prevalence_M,ymin=Prevalence_LB_M,ymax=Prevalence_UB_M),color="red", size = 0.5)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=DHS_Prevalence_M,ymin=DHS_Prevalence_LB_M,ymax=DHS_Prevalence_UB_M),color="purple", size = 0.5)+
  labs(title="HIV prevalence males RPHIA (red) DHS (purple)")

for (t in 1:trial_select){
  Prevalence_total[t,1,] <- (N_total_all[t,1,])/(N_total_all[t,1,]+S_total_all[t,1,])
  Prevalence_total_y[t,] <- Prevalence_total[t,1,][seq(1,length(Prevalence_total[t,1,]),12)]}
df_prev <- as.data.frame(Prevalence_total_y)
colnames(df_prev) <- c(year)
df_prev_long <- df_prev %>%
  pivot_longer(
    cols=1:27,
    names_to="Year",
    values_to="Prev"
  )
ggplot()+
  geom_point(data=df_prev_long, mapping=aes(x=Year,y=Prev), color = "darkgreen", size = 1, alpha = 1)+
  geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=Prevalence,ymin=Prevalence_LB,ymax=Prevalence_UB),color="red", size = 0.5)+
  labs(title="HIV prevalence RPHIA (15-64)")

# HIV prevalence female
for (t in 1:trial_select){
  Prevalence_F_total[t,1,] <- (N_F_total_all[t,1,])/(N_F_total_all[t,1,]+S_F_total_all[t,1,])
  Prevalence_F_total_y[t,] <- Prevalence_F_total[t,1,][seq(1,length(Prevalence_F_total[t,1,]),12)]}
df_prev_F <- as.data.frame(Prevalence_F_total_y)
colnames(df_prev_F) <- c(year)
df_prev_F_long <- df_prev_F %>%
  pivot_longer(
    cols=1:27,
    names_to="Year",
    values_to="Prev"
  )
ggplot()+
  geom_point(data=df_prev_F_long, mapping=aes(x=Year,y=Prev), color = "darkgreen", size = 1, alpha = 1)+
  geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=Prevalence_F,ymin=Prevalence_LB_F,ymax=Prevalence_UB_F),color="red", size = 0.5)+
  labs(title="HIV prevalence females RPHIA (15-64)")

# HIV prevalence male
for (t in 1:trial_select){
  Prevalence_M_total[t,1,] <- (N_M_total_all[t,1,])/(N_M_total_all[t,1,]+S_F_total_all[t,1,])
  Prevalence_M_total_y[t,] <- Prevalence_M_total[t,1,][seq(1,length(Prevalence_M_total[t,1,]),12)]}
df_prev_M <- as.data.frame(Prevalence_M_total_y)
colnames(df_prev_M) <- c(year)
df_prev_M_long <- df_prev_M %>%
  pivot_longer(
    cols=1:27,
    names_to="Year",
    values_to="Prev"
  )
ggplot()+
  geom_point(data=df_prev_M_long, mapping=aes(x=Year,y=Prev), color = "darkgreen", size = 1, alpha = 1)+
  geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=Prevalence_M,ymin=Prevalence_LB_M,ymax=Prevalence_UB_M),color="red", size = 0.5)+
  labs(title="HIV prevalence males RPHIA (15-64)")

# HIV prevalence 15-24
# for (t in 1:trial_select){
#   Prevalence_young_total[t,1,] <- (N_F_young_total[t,1,]+N_M_young_total[t,1,])/(N_F_young_total[t,1,]+N_M_young_total[t,1,]+S_F_young_total[t,1,]+S_M_young_total[t,1,])
#   Prevalence_young_total_y[t,] <- Prevalence_young_total[t,1,][seq(1,length(Prevalence_young_total[t,1,]),12)]}
# df_prev_young <- as.data.frame(Prevalence_young_total_y)
# colnames(df_prev_young) <- c(year)
# df_prev_long <- df_prev %>%
#   pivot_longer(
#     cols=1:27,
#     names_to="Year",
#     values_to="Prev"
#   )
# ggplot()+
#   geom_point(data=df_prev_long, mapping=aes(x=Year,y=Prev), color = "darkgreen", size = 1, alpha = 1)+
#   geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Prevalence_young,ymin=Prevalence_LB_young,ymax=Prevalence_UB_young),color="red", size = 0.5)+
#   labs(title="HIV prevalence 15-24 RPHIA")

# HIV prevalence female 15-24
for (t in 1:trial_select){
  Prevalence_F_young_total[t,1,] <- (N_F_young_total[t,1,])/(N_F_young_total[t,1,]+S_F_young_total[t,1,])
  Prevalence_F_young_total_y[t,] <- Prevalence_F_young_total[t,1,][seq(1,length(Prevalence_F_young_total[t,1,]),12)]}
df_prev_F_young <- as.data.frame(Prevalence_F_young_total_y)
colnames(df_prev_F_young) <- c(year)
df_prev_F_young_long <- df_prev_F_young %>%
  pivot_longer(
    cols=1:27,
    names_to="Year",
    values_to="Prev"
  )
ggplot()+
  geom_point(data=df_prev_F_young_long, mapping=aes(x=Year,y=Prev), color = "darkgreen", size = 1, alpha = 1)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Prevalence_young_F,ymin=Prevalence_LB_young_F,ymax=Prevalence_UB_young_F),color="red", size = 0.5)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=DHS_Prevalence_young_F,ymin=DHS_Prevalence_LB_young_F,ymax=DHS_Prevalence_UB_young_F),color="purple", size = 0.5)+
  labs(title="HIV prevalence female 15-24 RPHIA (red) DHS (purple)")

# HIV prevalence male 15-24
for (t in 1:trial_select){
  Prevalence_M_young_total[t,1,] <- (N_M_young_total[t,1,])/(N_M_young_total[t,1,]+S_F_young_total[t,1,])
  Prevalence_M_young_total_y[t,] <- Prevalence_M_young_total[t,1,][seq(1,length(Prevalence_M_young_total[t,1,]),12)]}
df_prev_M <- as.data.frame(Prevalence_M_total_y)
colnames(df_prev_M) <- c(year)
df_prev_M_long <- df_prev_M %>%
  pivot_longer(
    cols=1:27,
    names_to="Year",
    values_to="Prev"
  )
ggplot()+
  geom_point(data=df_prev_M_long, mapping=aes(x=Year,y=Prev), color = "darkgreen", size = 1, alpha = 1)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Prevalence_young_M,ymin=Prevalence_LB_young_M,ymax=Prevalence_UB_young_M),color="red", size = 0.5)+
  geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=DHS_Prevalence_young_M,ymin=DHS_Prevalence_LB_young_M,ymax=DHS_Prevalence_UB_young_M),color="purple", size = 0.5)+
  labs(title="HIV prevalence Male 15-24 RPHIA (red) DHS (purple)")




# HIV incidence
for (t in 1:trial_select) {
  NI_total_y[t,] <- rollapply(NI_total[t,1,],12,sum,by=12)
  S_total_y[t,] <- S_total[t,1,][seq(1, length(S_total[t,1,]), 12)]
  N_total_y[t,] <- N_total[t,1,][seq(1, length(N_total[t,1,]), 12)]
  Incidence_total_y[t,] <- NI_total_y[t,]/S_total_y[t,]}
  df_inci <- as.data.frame(Incidence_total_y)
  colnames(df_inci) <- c(year)
  df_inci_long <- df_inci %>%
    pivot_longer(
      cols=1:27,
      names_to="Year",
      values_to="Incidence"
    )
  ggplot()+
    geom_point(data=df_inci_long, mapping=aes(x=Year,y=Incidence), color = "darkgreen", size = 1, alpha = 1)+
    geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Incidence,ymin=Incidence_LB,ymax=Incidence_UB),color="red", size = 0.5)+
    labs(title="HIV incidence")
  
  # HIV incidence female
  for (t in 1:trial_select) {
    NI_F_total_y[t,] <- rollapply(NI_F_total[t,1,],12,sum,by=12)
    S_F_total_y[t,] <- S_F_total[t,1,][seq(1, length(S_F_total[t,1,]), 12)]
    N_F_total_y[t,] <- N_F_total[t,1,][seq(1, length(N_F_total[t,1,]), 12)]
    Incidence_F_total_y[t,] <- NI_F_total_y[t,]/S_F_total_y[t,]}
    df_inci <- as.data.frame(Incidence_F_total_y)
    colnames(df_inci) <- c(year)
    df_inci_long <- df_inci %>%
      pivot_longer(
        cols=1:27,
        names_to="Year",
        values_to="Incidence"
      )
    ggplot()+
      geom_point(data=df_inci_long, mapping=aes(x=Year,y=Incidence), color = "darkgreen", size = 1, alpha = 1)+
      geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Incidence_F,ymin=Incidence_F_LB,ymax=Incidence_F_UB),color="red", size = 0.5)+
      labs(title="HIV incidence female")
    
    # HIV incidence male
    for (t in 1:trial_select) {
      NI_M_total_y[t,] <- rollapply(NI_M_total[t,1,],12,sum,by=12)
      S_M_total_y[t,] <- S_M_total[t,1,][seq(1, length(S_M_total[t,1,]), 12)]
      N_M_total_y[t,] <- N_M_total[t,1,][seq(1, length(N_M_total[t,1,]), 12)]
      Incidence_M_total_y[t,] <- NI_M_total_y[t,]/S_M_total_y[t,]}
      df_inci <- as.data.frame(Incidence_M_total_y)
      colnames(df_inci) <- c(year)
      df_inci_long <- df_inci %>%
        pivot_longer(
          cols=1:27,
          names_to="Year",
          values_to="Incidence"
        )
      ggplot()+
        geom_point(data=df_inci_long, mapping=aes(x=Year,y=Incidence), color = "darkgreen", size = 1, alpha = 1)+
        geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Incidence_M,ymin=Incidence_M_LB,ymax=Incidence_M_UB),color="red", size = 0.5)+
        labs(title="HIV incidence male")
      
      
      
      # HIV incidence female 15-24
      for (t in 1:trial_select) {
        NI_F_young_total_y[t,] <- rollapply(NI_F_young_total[t,1,],12,sum,by=12)
        S_F_young_total_y[t,] <- S_F_young_total[t,1,][seq(1, length(S_F_young_total[t,1,]), 12)]
        N_F_young_total_y[t,] <- N_F_young_total[t,1,][seq(1, length(N_F_young_total[t,1,]), 12)]
        Incidence_F_young_total_y[t,] <- NI_F_young_total_y[t,]/S_F_young_total_y[t,]}
        df_inci <- as.data.frame(Incidence_F_young_total_y)
        colnames(df_inci) <- c(year)
        df_inci_long <- df_inci %>%
          pivot_longer(
            cols=1:27,
            names_to="Year",
            values_to="Incidence"
          )
        ggplot()+
          geom_point(data=df_inci_long, mapping=aes(x=Year,y=Incidence), color = "darkgreen", size = 1, alpha = 1)+
          geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Incidence_F_young,ymin=Incidence_F_young_LB,ymax=Incidence_F_young_UB),color="red", size = 0.5)+
          labs(title="HIV incidence female 15-24")
        
        # HIV incidence male 15-24 
        for (t in 1:trial_select) {
          NI_M_young_total_y[t,] <- rollapply(NI_M_young_total[t,1,],12,sum,by=12)
          S_M_young_total_y[t,] <- S_M_young_total[t,1,][seq(1, length(S_M_young_total[t,1,]), 12)]
          N_M_young_total_y[t,] <- N_M_young_total[t,1,][seq(1, length(N_M_young_total[t,1,]), 12)]
          Incidence_M_young_total_y[t,] <- NI_M_young_total_y[t,]/S_M_young_total_y[t,]}
          df_inci <- as.data.frame(Incidence_M_young_total_y)
          colnames(df_inci) <- c(year)
          df_inci_long <- df_inci %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="Incidence"
            )
          ggplot()+
            geom_point(data=df_inci_long, mapping=aes(x=Year,y=Incidence), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Incidence_M_young,ymin=Incidence_M_young_LB,ymax=Incidence_M_young_UB),color="red", size = 0.5)+
            labs(title="HIV incidence male 15-24")
          
          
          # Rate of people living with diagnosed HIV
          for (t in 1:trial_select) {
            Diagnosed_total_y[t,] <- Diagnosed_total[t,1,][seq(1, length(Diagnosed_total[t,1,]), 12)]
            Rate_diagnosed_total_y[t,] <- Diagnosed_total_y[t,]/(N_total_y[t,])*100}
          df_rate_hiv <- as.data.frame(Rate_diagnosed_total_y)
          colnames(df_rate_hiv) <- c(year)
          df_rate_hiv_long <- df_rate_hiv %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="HIV_rate"
            )
          ggplot()+
            geom_point(data=df_rate_hiv_long, mapping=aes(x=Year,y=HIV_rate), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_diagnosed_total,ymin=Percent_diagnosed_total_LB,ymax=Percent_diagnosed_total_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of people with known HIV")
          
          
          # Rate of people living with diagnosed HIV female
          for (t in 1:trial_select) {
            Diagnosed_F_total_y[t,] <- Diagnosed_F_total[t,1,][seq(1, length(Diagnosed_F_total[t,1,]), 12)]
            Rate_diagnosed_F_total_y[t,] <- Diagnosed_F_total_y[t,]/(N_F_total_y[t,])*100}
          df_rate_hiv_F <- as.data.frame(Rate_diagnosed_F_total_y)
          colnames(df_rate_hiv_F) <- c(year)
          df_rate_hiv_F_long <- df_rate_hiv_F %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="HIV_rate"
            )
          ggplot()+
            geom_point(data=df_rate_hiv_F_long, mapping=aes(x=Year,y=HIV_rate), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_diagnosed_Female,ymin=Percent_diagnosed_Female_LB,ymax=Percent_diagnosed_Female_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of females with known HIV")
          
          
          # Rate of people living with diagnosed HIV male
          for (t in 1:trial_select) {
            Diagnosed_M_total_y[t,] <- Diagnosed_M_total[t,1,][seq(1, length(Diagnosed_M_total[t,1,]), 12)]
            Rate_diagnosed_M_total_y[t,] <- Diagnosed_M_total_y[t,]/(N_M_total_y[t,])*100}
          df_rate_hiv_M <- as.data.frame(Rate_diagnosed_M_total_y)
          colnames(df_rate_hiv_M) <- c(year)
          df_rate_hiv_M_long <- df_rate_hiv_M %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="HIV_rate"
            )
          ggplot()+
            geom_point(data=df_rate_hiv_M_long, mapping=aes(x=Year,y=HIV_rate), color = "darkgreen", size = 1, alpha = 1)+
            geom_point(data=Target_all, mapping=aes(x=1:27, y=Percent_diagnosed_Male,ymin=Percent_diagnosed_Male_LB,ymax=Percent_diagnosed_Male_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of males with known HIV")
          
          # Rate of people living with diagnosed HIV female 15-24
          for (t in 1:trial_select) {
            Diagnosed_F_young_total_y[t,] <- Diagnosed_F_young_total[t,1,][seq(1, length(Diagnosed_F_young_total[t,1,]), 12)]
            Rate_diagnosed_F_young_total_y[t,] <- Diagnosed_F_young_total_y[t,]/(N_F_young_total_y[t,])*100}
          df_rate_hiv_F_young <- as.data.frame(Rate_diagnosed_F_young_total_y)
          colnames(df_rate_hiv_F_young) <- c(year)
          df_rate_hiv_F_young_long <- df_rate_hiv_F_young %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="HIV_rate"
            )
          ggplot()+
            geom_point(data=df_rate_hiv_F_young_long, mapping=aes(x=Year,y=HIV_rate), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_diagnosed_Female_YOUNG,ymin=Percent_diagnosed_Female_YOUNG_LB,ymax=Percent_diagnosed_Female_YOUNG_UB),color="red", size = 1)+ylim(0,100)+
            labs(title="% of females 15-24 with known HIV")
          
          
          # Rate of people living with diagnosed HIV male 15-24
          for (t in 1:trial_select) {
            Diagnosed_M_young_total_y[t,] <- Diagnosed_M_young_total[t,1,][seq(1, length(Diagnosed_M_young_total[t,1,]), 12)]
            Rate_diagnosed_M_young_total_y[t,] <- Diagnosed_M_young_total_y[t,]/(N_M_young_total_y[t,])*100}
          df_rate_hiv_M_young <- as.data.frame(Rate_diagnosed_M_young_total_y)
          colnames(df_rate_hiv_M_young) <- c(year)
          df_rate_hiv_M_young_long <- df_rate_hiv_M_young %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="HIV_rate"
            )
          ggplot()+
            geom_point(data=df_rate_hiv_M_young_long, mapping=aes(x=Year,y=HIV_rate), color = "darkgreen", size = 1, alpha = 1)+
            # geom_point(data=Target_all, mapping=aes(x=1:27, y=Percent_diagnosed_Male_YOUNG),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of males 15-24 with known HIV")
          
          
          # Number on ART 
          for (t in 1:trial_select){
            #On_ART_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            On_ART_grand_total_y[t,] <- On_ART_grand_total[t,1,][seq(1,length(On_ART_grand_total[t,1,]),12)]}
          df_p_vs <- as.data.frame(On_ART_grand_total_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="ART"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=ART), color = "darkgreen", size = 1, alpha = 1)+
            geom_point(data=Target_all, mapping=aes(x=1:27, y=N_ART),color="red", size = 1)+
            labs(title="Number on ART 15-49")
          
          # Percent on ART
          for (t in 1:trial_select){
            #On_ART_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            On_ART_p_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]/ N_total_y[t,]*100}
          df_p_vs <- as.data.frame(On_ART_p_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="ART"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=ART), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=P_ART_All,ymin=P_ART_All_LB,ymax=P_ART_All_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% on ART Among All living with HIV")
          
          
          # Percent on ART female
          for (t in 1:trial_select){
            #On_ART_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            On_ART_p_F_y[t,] <- On_ART_total_F[t,1,][seq(1,length(On_ART_total_F[t,1,]),12)]/ N_F_total_y[t,]*100}
          df_p_vs <- as.data.frame(On_ART_p_F_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="ART"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=ART), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=P_ART_All_F,ymin=P_ART_All_F_LB,ymax=P_ART_All_F_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% on ART Among females living with HIV")
          
          
          # Percent on ART male
          for (t in 1:trial_select){
            #On_ART_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            On_ART_p_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]/ N_M_total_y[t,]*100}
          df_p_vs <- as.data.frame(On_ART_p_M_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="ART"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=ART), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=P_ART_All_M,ymin=P_ART_All_M_LB,ymax=P_ART_All_M_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% on ART Among males living with HIV")
          
          
          
          # Percent on ART female 15-24
          for (t in 1:trial_select){
            #On_ART_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            On_ART_p_F_young_y[t,] <- On_ART_total_F_young[t,1,][seq(1,length(On_ART_total_F_young[t,1,]),12)]/ N_F_young_total_y[t,]*100}
          df_p_vs <- as.data.frame(On_ART_p_F_young_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="ART"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=ART), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=P_ART_All_F_young,ymin=P_ART_All_F_young_LB,ymax=P_ART_All_F_young_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% on ART Among females 15-24 living with HIV")
          
          # Percent on ART among diagnosed
          for (t in 1:trial_select){
            #On_ART_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            On_ART_p_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]/Diagnosed_total_y[t,]*100}
          df_p_vs <- as.data.frame(On_ART_p_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="ART"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=ART), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=P_ART_All,ymin=P_ART_All_LB,ymax=P_ART_All_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% on ART Among Diagnosed")
          
          
          # Percent on ART female
          for (t in 1:trial_select){
            #On_ART_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            On_ART_p_F_y[t,] <- On_ART_total_F[t,1,][seq(1,length(On_ART_total_F[t,1,]),12)]/Diagnosed_F_total_y[t,]*100}
          df_p_vs <- as.data.frame(On_ART_p_F_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="ART"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=ART), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=P_ART_All_F,ymin=P_ART_All_F_LB,ymax=P_ART_All_F_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% on ART Among females diagnosed")
          
          
          # Percent on ART male
          for (t in 1:trial_select){
            #On_ART_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            On_ART_p_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]/Diagnosed_M_total_y[t,]*100}
          df_p_vs <- as.data.frame(On_ART_p_M_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="ART"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=ART), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=P_ART_All_M,ymin=P_ART_All_M_LB,ymax=P_ART_All_M_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% on ART Among males diagnosed")
          
          
          
          # Percent on ART female 15-24
          for (t in 1:trial_select){
            #On_ART_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            On_ART_p_F_young_y[t,] <- On_ART_total_F_young[t,1,][seq(1,length(On_ART_total_F_young[t,1,]),12)]/Diagnosed_F_young_total_y[t,]*100}
          df_p_vs <- as.data.frame(On_ART_p_F_young_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="ART"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=ART), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=P_ART_All_F_young,ymin=P_ART_All_F_young_LB,ymax=P_ART_All_F_young_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% on ART Among females 15-24 diagnosed")
          
          # Percent viral suppression (among ART)
          for (t in 1:trial_select){
            On_ART_total_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            Suppressed_total_y[t,] <- Suppressed_total[t,1,][seq(1,length(Suppressed_total[t,1,]),12)]
            P_suppressed_art_total_y[t,] <- Suppressed_total_y[t,]/On_ART_total_y[t,]*100}
          df_p_vs_art <- as.data.frame(P_suppressed_art_total_y)
          colnames(df_p_vs_art) <- c(year)
          df_p_vs_art_long <- df_p_vs_art %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS_art"
            )
          ggplot()+
            geom_point(data=df_p_vs_art_long, mapping=aes(x=Year,y=P_VS_art), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_on_ART_total,ymin=Percent_VS_on_ART_total_LB,ymax=Percent_VS_on_ART_total_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression among art")
          
          # Percent viral suppression female (among ART)
          for (t in 1:trial_select){
            On_ART_total_F_y[t,] <- On_ART_total_F[t,1,][seq(1,length(On_ART_total_F[t,1,]),12)]
            Suppressed_total_F_y[t,] <- Suppressed_total_F[t,1,][seq(1,length(Suppressed_total_F[t,1,]),12)]
            P_suppressed_art_F_total_y[t,] <- Suppressed_total_F_y[t,]/On_ART_total_F_y[t,]*100}
          df_p_vs_art_F <- as.data.frame(P_suppressed_art_F_total_y)
          colnames(df_p_vs_art_F) <- c(year)
          df_p_vs_art_F_long <- df_p_vs_art_F %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS_art"
            )
          ggplot()+
            geom_point(data=df_p_vs_art_F_long, mapping=aes(x=Year,y=P_VS_art), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_on_ART_Female,ymin=Percent_VS_on_ART_Female_LB,ymax=Percent_VS_on_ART_Female_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression female among art")
          
          # Percent viral suppression male (among ART)
          for (t in 1:trial_select){
            On_ART_total_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]
            Suppressed_total_M_y[t,] <- Suppressed_total_M[t,1,][seq(1,length(Suppressed_total_M[t,1,]),12)]
            P_suppressed_art_M_total_y[t,] <- Suppressed_total_M_y[t,]/On_ART_total_M_y[t,]*100}
          df_p_vs_art_M <- as.data.frame(P_suppressed_art_M_total_y)
          colnames(df_p_vs_art_M) <- c(year)
          df_p_vs_art_M_long <- df_p_vs_art_M %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS_art"
            )
          ggplot()+
            geom_point(data=df_p_vs_art_M_long, mapping=aes(x=Year,y=P_VS_art), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_on_ART_Male,ymin=Percent_VS_on_ART_Male_LB,ymax=Percent_VS_on_ART_Male_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression male among art")
          
          # Percent viral suppression female 15-25 (among ART)
          for (t in 1:trial_select){
            On_ART_total_F_young_y[t,] <- On_ART_total_F_young[t,1,][seq(1,length(On_ART_total_F_young[t,1,]),12)]
            Suppressed_total_F_young_y[t,] <- Suppressed_total_F_young[t,1,][seq(1,length(Suppressed_total_F_young[t,1,]),12)]
            P_suppressed_art_F_young_total_y[t,] <- Suppressed_total_F_young_y[t,]/On_ART_total_F_young_y[t,]*100}
          df_p_vs_art_F_young <- as.data.frame(P_suppressed_art_F_young_total_y)
          colnames(df_p_vs_art_F_young) <- c(year)
          df_p_vs_art_F_young_long <- df_p_vs_art_F_young %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS_art"
            )
          ggplot()+
            geom_point(data=df_p_vs_art_F_young_long, mapping=aes(x=Year,y=P_VS_art), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_on_ART_Female_YOUNG,ymin=Percent_VS_on_ART_Female_YOUNG_LB,ymax=Percent_VS_on_ART_Female_YOUNG_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression female 15-24 among art")
          
          # Percent viral suppression male 15-24 (among ART)
          for (t in 1:trial_select){
            On_ART_total_M_young_y[t,] <- On_ART_total_M_young[t,1,][seq(1,length(On_ART_total_M_young[t,1,]),12)]
            Suppressed_total_M_young_y[t,] <- Suppressed_total_M_young[t,1,][seq(1,length(Suppressed_total_M_young[t,1,]),12)]
            P_suppressed_art_M_young_total_y[t,] <- Suppressed_total_M_young_y[t,]/On_ART_total_M_young_y[t,]*100}
          df_p_vs_art_M_young <- as.data.frame(P_suppressed_art_M_young_total_y)
          colnames(df_p_vs_art_M_young) <- c(year)
          df_p_vs_art_M_young_long <- df_p_vs_art_M_young %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS_art"
            )
          ggplot()+
            geom_point(data=df_p_vs_art_M_young_long, mapping=aes(x=Year,y=P_VS_art), color = "darkgreen", size = 1, alpha = 1)+
            #geom_point(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_on_ART_Male_YOUNG),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression male 15-24 among art")
          
          
          
          
          
          # Percent viral suppression (among ALL)
          for (t in 1:trial_select){
            #On_ART_total_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            Suppressed_total_y[t,] <- Suppressed_total[t,1,][seq(1,length(Suppressed_total[t,1,]),12)]
            P_suppressed_total_y[t,] <- Suppressed_total_y[t,]/N_total_y[t,]*100}
          df_p_vs <- as.data.frame(P_suppressed_total_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=P_VS), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_total, ymin=Percent_VS_total_LB, ymax=Percent_VS_total_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression")
          
          # Percent viral suppression female
          for (t in 1:trial_select){
            #On_ART_total_F_y[t,] <- On_ART_total_F[t,1,][seq(1,length(On_ART_total_F[t,1,]),12)]
            Suppressed_total_F_y[t,] <- Suppressed_total_F[t,1,][seq(1,length(Suppressed_total_F[t,1,]),12)]
            P_suppressed_F_total_y[t,] <- Suppressed_total_F_y[t,]/N_F_total_y[t,]*100}
          df_p_vs_F <- as.data.frame(P_suppressed_F_total_y)
          colnames(df_p_vs_F) <- c(year)
          df_p_vs_F_long <- df_p_vs_F %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS"
            )
          ggplot()+
            geom_point(data=df_p_vs_F_long, mapping=aes(x=Year,y=P_VS), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_Female, ymin=Percent_VS_Female_LB, ymax=Percent_VS_Female_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression female")
          
          # Percent viral suppression male
          for (t in 1:trial_select){
            #On_ART_total_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]
            Suppressed_total_M_y[t,] <- Suppressed_total_M[t,1,][seq(1,length(Suppressed_total_M[t,1,]),12)]
            P_suppressed_M_total_y[t,] <- Suppressed_total_M_y[t,]/N_M_total_y[t,]*100}
          df_p_vs_M <- as.data.frame(P_suppressed_M_total_y)
          colnames(df_p_vs_M) <- c(year)
          df_p_vs_M_long <- df_p_vs_M %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS"
            )
          ggplot()+
            geom_point(data=df_p_vs_M_long, mapping=aes(x=Year,y=P_VS), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_Male, ymin=Percent_VS_Male_LB, ymax=Percent_VS_Male_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression male")
          
          # Percent viral suppression (among ALL)
          for (t in 1:trial_select){
            #On_ART_total_y[t,] <- On_ART_total[t,1,][seq(1,length(On_ART_total[t,1,]),12)]
            N_total_y[t,] <- N_total_all[t,1,][seq(1, length(N_total_all[t,1,]), 12)]
            Suppressed_total_y[t,] <- Suppressed_total_all[t,1,][seq(1,length(Suppressed_total_all[t,1,]),12)]
            P_suppressed_total_y[t,] <- Suppressed_total_y[t,]/N_total_y[t,]*100}
          df_p_vs <- as.data.frame(P_suppressed_total_y)
          colnames(df_p_vs) <- c(year)
          df_p_vs_long <- df_p_vs %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS"
            )
          ggplot()+
            geom_point(data=df_p_vs_long, mapping=aes(x=Year,y=P_VS), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=Percent_VS_total, ymin=Percent_VS_total_LB, ymax=Percent_VS_total_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression (15-64)")
          
          # Percent viral suppression female
          for (t in 1:trial_select){
            N_F_total_y[t,] <- N_F_total_all[t,1,][seq(1, length(N_F_total_all[t,1,]), 12)]
            Suppressed_total_F_y[t,] <- Suppressed_total_F[t,1,][seq(1,length(Suppressed_total_F[t,1,]),12)]
            P_suppressed_F_total_y[t,] <- Suppressed_total_F_y[t,]/N_F_total_y[t,]*100}
          df_p_vs_F <- as.data.frame(P_suppressed_F_total_y)
          colnames(df_p_vs_F) <- c(year)
          df_p_vs_F_long <- df_p_vs_F %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS"
            )
          ggplot()+
            geom_point(data=df_p_vs_F_long, mapping=aes(x=Year,y=P_VS), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=Percent_VS_Female, ymin=Percent_VS_Female_LB, ymax=Percent_VS_Female_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression female (15-64)")
          
          # Percent viral suppression male
          for (t in 1:trial_select){
            #On_ART_total_M_y[t,] <- On_ART_total_M[t,1,][seq(1,length(On_ART_total_M[t,1,]),12)]
            N_M_total_y[t,] <- N_M_total_all[t,1,][seq(1, length(N_M_total_all[t,1,]), 12)]
            Suppressed_total_M_y[t,] <- Suppressed_total_M[t,1,][seq(1,length(Suppressed_total_M[t,1,]),12)]
            P_suppressed_M_total_y[t,] <- Suppressed_total_M_y[t,]/N_M_total_y[t,]*100}
          df_p_vs_M <- as.data.frame(P_suppressed_M_total_y)
          colnames(df_p_vs_M) <- c(year)
          df_p_vs_M_long <- df_p_vs_M %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS"
            )
          ggplot()+
            geom_point(data=df_p_vs_M_long, mapping=aes(x=Year,y=P_VS), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_new, mapping=aes(x=1:27, y=Percent_VS_Male, ymin=Percent_VS_Male_LB, ymax=Percent_VS_Male_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression male (15-64)")
          
          # Percent viral suppression female 15-25
          for (t in 1:trial_select){
            #On_ART_total_F_young_y[t,] <- On_ART_total_F_young[t,1,][seq(1,length(On_ART_total_F_young[t,1,]),12)]
            Suppressed_total_F_young_y[t,] <- Suppressed_total_F_young[t,1,][seq(1,length(Suppressed_total_F_young[t,1,]),12)]
            P_suppressed_F_young_total_y[t,] <- Suppressed_total_F_young_y[t,]/N_F_young_total_y[t,]*100}
          df_p_vs_F_young <- as.data.frame(P_suppressed_F_young_total_y)
          colnames(df_p_vs_F_young) <- c(year)
          df_p_vs_F_young_long <- df_p_vs_F_young %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS"
            )
          ggplot()+
            geom_point(data=df_p_vs_F_young_long, mapping=aes(x=Year,y=P_VS), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_Female_YOUNG, ymin=Percent_VS_Female_YOUNG_LB, ymax=Percent_VS_Female_YOUNG_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression female 15-24")
          
          # Percent viral suppression male 15-24
          for (t in 1:trial_select){
            #On_ART_total_M_young_y[t,] <- On_ART_total_M_young[t,1,][seq(1,length(On_ART_total_M_young[t,1,]),12)]
            Suppressed_total_M_young_y[t,] <- Suppressed_total_M_young[t,1,][seq(1,length(Suppressed_total_M_young[t,1,]),12)]
            P_suppressed_M_young_total_y[t,] <- Suppressed_total_M_young_y[t,]/N_M_young_total_y[t,]*100}
          df_p_vs_M_young <- as.data.frame(P_suppressed_M_young_total_y)
          colnames(df_p_vs_M_young) <- c(year)
          df_p_vs_M_young_long <- df_p_vs_M_young %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="P_VS"
            )
          ggplot()+
            geom_point(data=df_p_vs_M_young_long, mapping=aes(x=Year,y=P_VS), color = "darkgreen", size = 1, alpha = 1)+
            # geom_point(data=Target_all, mapping=aes(x=1:27, y=Percent_VS_Male_YOUNG, ymin=Percent_VS_Male_YOUNG_LB, ymax=Percent_VS_Male_YOUNG_UB),color="red", size = 1)+
            ylim(0,100)+
            labs(title="% of viral suppression male 15-24")
          
          
          # new infections
          for (t in 1:trial_select) {
            NI_total_y[t,] <- rollapply(NI_total[t,1,],12,sum,by=12)}
          df_NI_HIV_ALL <- as.data.frame(NI_total_y)
          colnames(df_NI_HIV_ALL) <- c(year)
          df_NI_HIV_ALL_long <- df_NI_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="NI"
            )
          ggplot()+
            geom_point(data=df_NI_HIV_ALL_long, mapping=aes(x=Year,y=NI), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_NEW_ALL,ymin=N_NEW_ALL_LB,ymax=N_NEW_ALL_UB),color="red", size = 0.5)+
            labs(title="New Infections (corroborration)")
          
          # new infections females
          for (t in 1:trial_select) {
            NI_F_total_y[t,] <- rollapply(NI_F_total[t,1,],12,sum,by=12)}
          df_NI_F_HIV_ALL <- as.data.frame(NI_F_total_y)
          colnames(df_NI_F_HIV_ALL) <- c(year)
          df_NI_F_HIV_ALL_long <- df_NI_F_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="NI"
            )
          ggplot()+
            geom_point(data=df_NI_F_HIV_ALL_long, mapping=aes(x=Year,y=NI), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_NEW_Female,ymin=N_NEW_Female_LB,ymax=N_NEW_Female_UB),color="red", size = 0.5)+
            labs(title="New Infections females with HIV (corroborration)")
          
          
          # new infections males
          for (t in 1:trial_select) {
            NI_M_total_y[t,] <- rollapply(NI_M_total[t,1,],12,sum,by=12)}
          df_NI_M_HIV_ALL <- as.data.frame(NI_M_total_y)
          colnames(df_NI_M_HIV_ALL) <- c(year)
          df_NI_M_HIV_ALL_long <- df_NI_M_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="NI"
            )
          ggplot()+
            geom_point(data=df_NI_M_HIV_ALL_long, mapping=aes(x=Year,y=NI), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_NEW_Male,ymin=N_NEW_Male_LB,ymax=N_NEW_Male_UB),color="red", size = 0.5)+
            labs(title="New Infections males with HIV (corroborration)")
          
          # new infections young (15-24)
          for (t in 1:trial_select) {
            NI_young_total_y[t,] <- rollapply(NI_young_total[t,1,],12,sum,by=12)}
          df_NI_young_HIV_ALL <- as.data.frame(NI_young_total_y)
          colnames(df_NI_young_HIV_ALL) <- c(year)
          df_NI_young_HIV_ALL_long <- df_NI_young_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="NI"
            )
          ggplot()+
            geom_point(data=df_NI_young_HIV_ALL_long, mapping=aes(x=Year,y=NI), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_NEW_Young,ymin=N_NEW_Young_LB,ymax=N_NEW_Young_UB),color="red", size = 0.5)+
            labs(title="New Infections age 15-24 with HIV (corroborration)")
          
          
          # new infections adult (15-44)
          for (t in 1:trial_select) {
            NI_adult_total_y[t,] <- rollapply(NI_adult_total[t,1,],12,sum,by=12)}
          df_NI_adult_HIV_ALL <- as.data.frame(NI_adult_total_y)
          colnames(df_NI_adult_HIV_ALL) <- c(year)
          df_NI_adult_HIV_ALL_long <- df_NI_adult_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="NI"
            )
          ggplot()+
            geom_point(data=df_NI_adult_HIV_ALL_long, mapping=aes(x=Year,y=NI), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_NEW_Adults,ymin=N_NEW_Adults_LB,ymax=N_NEW_Adults_UB),color="red", size = 0.5)+
            labs(title="New Infections age 15-44 with HIV (corroborration)")
          
          # Number with HIV
          for (t in 1:trial_select) {
            N_total_y[t,] <- N_total[t,1,][seq(1,length(N_total[t,1,]),12)]}
          df_N_HIV_ALL <- as.data.frame(N_total_y)
          colnames(df_N_HIV_ALL) <- c(year)
          df_N_HIV_ALL_long <- df_N_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="N"
            )
          ggplot()+
            geom_point(data=df_N_HIV_ALL_long, mapping=aes(x=Year,y=N), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_HIV_ALL,ymin=N_HIV_ALL_LB,ymax=N_HIV_ALL_UB),color="red", size = 0.5)+
            labs(title="Total with HIV (corroborration)")
          
          
          # Number with HIV females
          for (t in 1:trial_select) {
            N_F_total_y[t,] <- N_F_total[t,1,][seq(1,length(N_F_total[t,1,]),12)]}
          df_N_F_HIV_ALL <- as.data.frame(N_F_total_y)
          colnames(df_N_F_HIV_ALL) <- c(year)
          df_N_F_HIV_ALL_long <- df_N_F_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="N"
            )
          ggplot()+
            geom_point(data=df_N_F_HIV_ALL_long, mapping=aes(x=Year,y=N), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_HIV_Female,ymin=N_HIV_Female_LB,ymax=N_HIV_Female_UB),color="red", size = 0.5)+
            labs(title="Total females with HIV (corroborration)")
          
          
          # Number with HIV males
          for (t in 1:trial_select) {
            N_M_total_y[t,] <- N_M_total[t,1,][seq(1,length(N_M_total[t,1,]),12)]}
          df_N_M_HIV_ALL <- as.data.frame(N_M_total_y)
          colnames(df_N_M_HIV_ALL) <- c(year)
          df_N_M_HIV_ALL_long <- df_N_M_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="N"
            )
          ggplot()+
            geom_point(data=df_N_M_HIV_ALL_long, mapping=aes(x=Year,y=N), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_HIV_Male,ymin=N_HIV_Male_LB,ymax=N_HIV_Male_UB),color="red", size = 0.5)+
            labs(title="Total males with HIV (corroborration)")
          
          # Number with HIV young (15-24)
          for (t in 1:trial_select) {
            N_young_total_y[t,] <- N_young_total[t,1,][seq(1,length(N_young_total[t,1,]),12)]}
          df_N_young_HIV_ALL <- as.data.frame(N_young_total_y)
          colnames(df_N_young_HIV_ALL) <- c(year)
          df_N_young_HIV_ALL_long <- df_N_young_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="N"
            )
          ggplot()+
            geom_point(data=df_N_young_HIV_ALL_long, mapping=aes(x=Year,y=N), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_HIV_Young,ymin=N_HIV_Young_LB,ymax=N_HIV_Young_UB),color="red", size = 0.5)+
            labs(title="Total age 15-24 with HIV (corroborration)")
          
          
          # Number with HIV adult (15-44)
          for (t in 1:trial_select) {
            N_adult_total_y[t,] <- N_adult_total[t,1,][seq(1,length(N_adult_total[t,1,]),12)]}
          df_N_adult_HIV_ALL <- as.data.frame(N_adult_total_y)
          colnames(df_N_adult_HIV_ALL) <- c(year)
          df_N_adult_HIV_ALL_long <- df_N_adult_HIV_ALL %>%
            pivot_longer(
              cols=1:27,
              names_to="Year",
              values_to="N"
            )
          ggplot()+
            geom_point(data=df_N_adult_HIV_ALL_long, mapping=aes(x=Year,y=N), color = "darkgreen", size = 1, alpha = 1)+
            geom_pointrange(data=Target_all, mapping=aes(x=1:27, y=N_HIV_Adults,ymin=N_HIV_Adults_LB,ymax=N_HIV_Adults_UB),color="red", size = 0.5)+
            labs(title="Total age 15-44 with HIV (corroborration)")
          
dev.off()

          