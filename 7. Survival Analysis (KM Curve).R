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


setwd("D:/IDH paper 11-8-23/R_analysis")

# sink("Cutoff_point_for_survival_IDH_whole_23-9-23.txt")
df <- read_excel("df_1_ratio_fd.xlsx")

colnames(df)

# for removing NAN values from fractal dimension
df1 <- df %>% drop_na("ed_meanfd","et_meanfd","Survival_months","ncr_net_meanfd","IDH_status")

#For removing the NAN values from Lacunarity 
df2 <- df %>% drop_na("ncr_net_meanlac", "et_meanlac","Survival_months", "ed_meanlac", "IDH_status")


print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
print(" For df1 Fractal Dimension")
print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for FD Enhancing Tumor")
  mtHL_et <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~et_meanfd,
                          data=df1, smethod="LogRank", pmethod="HL")
  print(mtHL_et)
  
  tiff("mtHL_et_meanfd_IDH.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_et,xlim=c(0,2),ylim=c(0,5),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_et <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ et_meanfd + ncr_net_meanfd +ed_meanfd ,
                             data=df1, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -FD Enhancing Tumor")
  print( pAdj_HL_et)
  plot(pAdj_HL_et,  xlab="Fractal Dimension",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for FD Non-Enhancing Tumor")
  mtHL_net <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ncr_net_meanfd,
                           data=df1, smethod="LogRank", pmethod="HL")
  print(mtHL_net)
  
  tiff("mtHL_ncr_net_meanfd_IDH.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_net,xlim=c(0,2),ylim=c(0,5),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_net <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ ncr_net_meanfd +et_meanfd + ed_meanfd ,
                              data=df1, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -FD Non-Enhancing Tumor ")
  print( pAdj_HL_net)
  plot(pAdj_HL_net,  xlab="Fractal Dimension",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for FD Edema")
  mtHL_ed <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ ed_meanfd,
                          data=df1, smethod="LogRank", pmethod="HL")
  print(mtHL_ed)
  
  tiff("mtHL_ed_meanfd_IDH.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_ed,xlim=c(1.8,2),ylim=c(0,5),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_ed <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~  ed_meanfd +ncr_net_meanfd +et_meanfd ,
                             data=df1, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -FD Edema ")
  print( pAdj_HL_ed)
  plot(pAdj_HL_ed,  xlab="Fractal Dimension",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}




print("")
print("")
print("==================================================================================================")
print("For df2 Lacunarity")
print("==================================================================================================")

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for Lacunarity Enhancing Tumor")
  mtHL_etL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~et_meanlac,
                           data=df2, smethod="LogRank", pmethod="HL")
  print(mtHL_etL)
  
  tiff("mtHL_et_meanlac_df2_2.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_etL,xlim=c(1,8),ylim=c(0,6),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_etL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ et_meanlac + ncr_net_meanlac +ed_meanlac ,
                              data=df2, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value Lacunarity Enhancing Tumor ")
  print( pAdj_HL_etL)
  plot(pAdj_HL_etL,  xlab="Lacunarity",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for Lacunarity Non-Enhancing Tumor")
  mtHL_netL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ncr_net_meanlac,
                            data=df2, smethod="LogRank", pmethod="HL")
  print(mtHL_netL)
  
  tiff("mtHL_ncr_net_meanlac_df2_2.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_netL,xlim=c(1,6),ylim=c(0,6),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_netL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ ncr_net_meanlac +et_meanlac + ed_meanlac ,
                               data=df2, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -Lacunarity Non-Enhancing Tumo ")
  print( pAdj_HL_netL)
  plot(pAdj_HL_netL,  xlab="Lacunarity",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for Lacunarity Edema")
  mtHL_edL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ ed_meanlac,
                           data=df2, smethod="LogRank", pmethod="HL")
  print(mtHL_edL)
  
  tiff("mtHL_ed_meanlac_df2_2.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_edL,xlim=c(0.6,2),ylim=c(0,6),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_edL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~  ed_meanlac +ncr_net_meanlac +et_meanlac ,
                              data=df2, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -Lacunarity Edema ")
  print( pAdj_HL_edL)
  plot(pAdj_HL_edL,  xlab="Lacunarity",ylab="Log-rank statistic" )
  print("-------------------------------------------------------------------------------------------------")
}



print("")
print("")
print("==================================================================================================")
print("For df3  total subjects")
print("==================================================================================================")


#For removing the NAN values from Fractal dimension total 
df3_FD <- df %>% drop_na("ed_meanfd","et_meanfd","Survival_months","ncr_net_meanfd")

#For removing the NAN values from Lacunarity total 
df3_lac <- df %>% drop_na("ncr_net_meanlac", "et_meanlac","Survival_months", "ed_meanlac")



print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
print(" For df3_FD Fractal Dimension")
print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for FD Enhancing Tumor")
  mtHL_et <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~et_meanfd,
                          data=df3_FD, smethod="LogRank", pmethod="HL")
  print(mtHL_et)
  
  tiff("mtHL_et_meanfd_total.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_et,xlim=c(0,2),ylim=c(0,5),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_et <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ et_meanfd + ncr_net_meanfd +ed_meanfd ,
                             data=df3_FD, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -FD Enhancing Tumor")
  print( pAdj_HL_et)
  plot(pAdj_HL_et,  xlab="Fractal Dimension",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for FD Non-Enhancing Tumor")
  mtHL_net <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ncr_net_meanfd,
                           data=df3_FD, smethod="LogRank", pmethod="HL")
  print(mtHL_net)
  
  tiff("mtHL_ncr_net_meanfd_total.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_net,xlim=c(0,2),ylim=c(0,5),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_net <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ ncr_net_meanfd +et_meanfd + ed_meanfd ,
                              data=df3_FD, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -FD Non-Enhancing Tumor ")
  print( pAdj_HL_net)
  plot(pAdj_HL_net,  xlab="Fractal Dimension",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for FD Edema")
  mtHL_ed <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ ed_meanfd,
                          data=df3_FD, smethod="LogRank", pmethod="HL")
  print(mtHL_ed)
  
  tiff("mtHL_ed_meanfd_total.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_ed,xlim=c(1.8,2),ylim=c(0,5),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_ed <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~  ed_meanfd +ncr_net_meanfd +et_meanfd ,
                             data=df3_FD, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -FD Edema ")
  print( pAdj_HL_ed)
  plot(pAdj_HL_ed,  xlab="Fractal Dimension",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}


print("")
print("==================================================================================================")
print("For df3_lac Lacunarity")
print("==================================================================================================")

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for Lacunarity Enhancing Tumor")
  mtHL_etL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~et_meanlac,
                           data=df3_lac, smethod="LogRank", pmethod="HL")
  print(mtHL_etL)
  
  tiff("mtHL_et_meanlac_df3_lac_total.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_etL,xlim=c(1,8),ylim=c(0,6),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_etL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ et_meanlac + ncr_net_meanlac +ed_meanlac ,
                              data=df3_lac, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value Lacunarity Enhancing Tumor ")
  print( pAdj_HL_etL)
  plot(pAdj_HL_etL,  xlab="Lacunarity",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for Lacunarity Non-Enhancing Tumor")
  mtHL_netL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ncr_net_meanlac,
                            data=df3_lac, smethod="LogRank", pmethod="HL")
  print(mtHL_netL)
  
  tiff("mtHL_ncr_net_meanlac_df3_lac_total.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_netL,xlim=c(1,6),ylim=c(0,6),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_netL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ ncr_net_meanlac +et_meanlac + ed_meanlac ,
                               data=df3_lac, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -Lacunarity Non-Enhancing Tumo ")
  print( pAdj_HL_netL)
  plot(pAdj_HL_netL,  xlab="Lacunarity",ylab="Log-rank statistic" )
  print("-----------------------------------------------------------------------------------------------")
}

{
  print("----------------------------------------------------------------------------------------------")
  print( "Cut off for Lacunarity Edema")
  mtHL_edL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~ ed_meanlac,
                           data=df3_lac, smethod="LogRank", pmethod="HL")
  print(mtHL_edL)
  
  tiff("mtHL_ed_meanlac_df3_lac_total.tiff",units="in", width=6, height=6, res=600)
  plot(mtHL_edL,xlim=c(0.6,2),ylim=c(0,6),pch=19,col="black",cex=1,axes=FALSE,lwd=4,
       xlab="",ylab="" , bty="n")
  axis(side = 1, lwd = 3,cex.axis=2)
  axis(side = 2, lwd = 3, cex.axis=2)
  dev.off()
  
  pAdj_HL_edL <- maxstat.test(Surv(Survival_months,Vital_status_1_dead) ~  ed_meanlac +ncr_net_meanlac +et_meanlac ,
                              data=df3_lac, smethod="LogRank", pmethod="HL", abseps=0.01)
  print("for Adjusted p value -Lacunarity Edema ")
  print( pAdj_HL_edL)
  plot(pAdj_HL_edL,  xlab="Lacunarity",ylab="Log-rank statistic" )
  print("-------------------------------------------------------------------------------------------------")
}

