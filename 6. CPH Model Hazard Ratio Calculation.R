# For cox hazard ratio for IDH
rm( list=ls() )
library( ggplot2 )
library( readxl )
library( ggpubr )
library(naniar)

library(emmeans)
library(tidyverse)

library(tidyr)

library(survminer)
library(ggfortify)

library("maxstat")
library("survival")
library("forestmodel")


setwd("D:/IDH paper 11-8-23/R_analysis")

# sink("Cutoff_point_for_survival_141_23-9-23.txt")
df <- read_excel("df_1_ratio_fd.xlsx")

colnames(df)

# for removing NAN values from fractal dimension
df1 <- df %>% drop_na("ed_meanfd","et_meanfd","Survival_months","ncr_net_meanfd","IDH_status")

df1_WT<-df1[df1$IDH_status == "WT",]
df1_Mut<-df1[df1$IDH_status == "Mutant",] 


# For DF1 IDH Fractal dimension
# for survival plots

df1 <- df1 %>% mutate("FD_Enhancing" = ifelse(et_meanfd < 0.69, "< 0.69",">= 0.69" ))
df1$"FD_Enhancing" <- factor(df1$"FD_Enhancing")


df1 <- df1 %>% mutate("FD_Nonenhancing" = ifelse(ncr_net_meanfd < 1.2, "< 1.2", ">= 1.2" ))
df1$"FD_Nonenhancing" <- factor(df1$"FD_Nonenhancing")

df1 <- df1 %>% mutate("FD_Edema" = ifelse(ed_meanfd < 1.8, "< 1.8",  ">= 1.8"))
df1$"FD_Edema" <- factor(df1$"FD_Edema")

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = df1$Survival_months, event = df1$Vital_status_1_dead)
surv_object
fit1<- survfit( surv_object~  df1$"FD_Enhancing", data=df1 )
ggsurvplot(fit1, data = df1,xlim=c(0,175), ylim=c(0,1) ,pval = TRUE)

# Fit a Cox proportional hazards model
fit.coxph <- coxph(Surv(Survival_months,Vital_status_1_dead) ~ FD_Enhancing + FD_Nonenhancing + FD_Edema, data = df1)
tiff("Coxph_multi_fd.tiff",units="in", width=8.5, height=3.5, res=600)
ggforest(fit.coxph, data = df1,fontsize=1,cpositions = c(0.005, 0.26, 0.4))
dev.off()

# Fit a Cox proportional hazards model
fit.coxph <- coxph(Surv(Survival_months,Vital_status_1_dead) ~FD_Enhancing, data = df1)
tiff("Coxph_univ_ET_fd.tiff",units="in", width=8.5, height=2, res=600)
ggforest(fit.coxph, data = df1,fontsize=1,cpositions = c(0.005, 0.26, 0.4))
dev.off()

# for median overall survival
FD_ETSurv<- survfit(Surv(Survival_months,Vital_status_1_dead) ~FD_Enhancing, data = df1)
print(FD_ETSurv)

# Fit a Cox proportional hazards model
fit.coxph <- coxph(Surv(Survival_months,Vital_status_1_dead) ~FD_Nonenhancing, data = df1)
tiff("Coxph_univ_NET_fd.tiff",units="in", width=8.5, height=2, res=600)
ggforest(fit.coxph, data = df1,fontsize=1,cpositions = c(0.005, 0.26, 0.4))
dev.off()

# for median overall survival
FD_NETSurv<- survfit(Surv(Survival_months,Vital_status_1_dead) ~FD_Nonenhancing, data = df1)
print(FD_NETSurv)

# Fit a Cox proportional hazards model
fit.coxph <- coxph(Surv(Survival_months,Vital_status_1_dead) ~FD_Edema, data = df1)
tiff("Coxph_univ_ED_fd.tiff",units="in", width=8.5, height=2, res=600)
ggforest(fit.coxph, data = df1,fontsize=1,cpositions = c(0.005, 0.26, 0.4))
dev.off()

# for median overall survival
FD_EDSurv<- survfit(Surv(Survival_months,Vital_status_1_dead) ~FD_Edema, data = df1)
print(FD_EDSurv)


# For DF2 IDH lacunarity
# for cox proportion hazard ratio plots

#For removing the NAN values from Lacunarity 
df2 <- df %>% drop_na("ncr_net_meanlac", "et_meanlac","Survival_months", "ed_meanlac", "IDH_status")

df2_WT<-df2[df2$IDH_status == "WT",]
df2_Mut<-df2[df2$IDH_status == "Mutant",] 

df2 <- df2 %>% mutate("Lacunarity_Enhancing" = ifelse(et_meanlac < 3.52, "< 3.52",">= 3.52" ))
df2$"Lacunarity_Enhancing" <- factor(df2$"Lacunarity_Enhancing")

df2 <- df2 %>% mutate("Lacunarity_Nonenhancing" = ifelse(ncr_net_meanlac < 1.48, "< 1.48", ">= 1.48" ))
df2$"Lacunarity_Nonenhancing" <- factor(df2$"Lacunarity_Nonenhancing")

df2 <- df2 %>% mutate("Lacunarity_Edema" = ifelse(ed_meanlac < 0.97, "< 0.97",  ">= 0.97"))
df2$"Lacunarity_Edema" <- factor(df2$"Lacunarity_Edema")

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = df2$Survival_months, event = df2$Vital_status_1_dead)
surv_object
# fit2<- survfit( surv_object~  `Lacunarity Enhancing`, data=df2 )
# ggsurvplot(fit2, data = df2,xlim=c(0,175), ylim=c(0,1) ,pval = TRUE)


# Fit a Cox proportional hazards model
fit.coxph <- coxph(Surv(Survival_months,Vital_status_1_dead) ~Lacunarity_Enhancing + Lacunarity_Nonenhancing + Lacunarity_Edema, data = df2)
tiff("Coxph_lac_multi.tiff",units="in", width=8.5, height=3.5, res=600)
ggforest(fit.coxph, data = df2,fontsize=1,cpositions = c(0.005, 0.26, 0.4))
dev.off()

# Fit a Cox proportional hazards model Enhancing 
fit.coxph <- coxph(Surv(Survival_months,Vital_status_1_dead) ~Lacunarity_Enhancing , data = df2)
tiff("Coxph_lac_univ_ET.tiff",units="in", width=8.5, height=2, res=600)
ggforest(fit.coxph, data = df2,fontsize=1,cpositions = c(0.005, 0.26, 0.4))
dev.off()

# for median overall survival
Lac_ETSurv<- survfit(Surv(Survival_months,Vital_status_1_dead) ~ Lacunarity_Enhancing, data = df2)
print(Lac_ETSurv)

# Fit a Cox proportional hazards model
fit.coxph <- coxph(Surv(Survival_months,Vital_status_1_dead) ~ Lacunarity_Nonenhancing , data = df2)
tiff("Coxph_lac_univ_NET.tiff",units="in", width=8.5, height=2, res=600)
ggforest(fit.coxph, data = df2,fontsize=1,cpositions = c(0.005, 0.26, 0.4))
dev.off()

# for median overall survival
Lac_NETSurv<- survfit(Surv(Survival_months,Vital_status_1_dead) ~ Lacunarity_Nonenhancing, data = df2)
print(Lac_NETSurv)

# Fit a Cox proportional hazards model
fit.coxph <- coxph(Surv(Survival_months,Vital_status_1_dead) ~ Lacunarity_Edema, data = df2)
tiff("Coxph_lac_univ_ED.tiff",units="in", width=8.5, height=2, res=600)
ggforest(fit.coxph, data = df2,fontsize=1,cpositions = c(0.005, 0.26, 0.4))
dev.off()

# for median overall survival
Lac_EDSurv<- survfit(Surv(Survival_months,Vital_status_1_dead) ~ Lacunarity_Edema, data = df2)
print(Lac_EDSurv)

#Sink()
