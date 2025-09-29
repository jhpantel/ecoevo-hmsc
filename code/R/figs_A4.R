## Code for Figures of Appendix4
# Author. J.H. Pantel &  R.J. Hermann
# Date 29.07.2025
## libraries
library(Hmsc)
library(tidyverse)
library(patchwork)
library(bayesplot)
bayesplot_theme_set(new = theme_default())
library(reshape2)
library(cowplot)
library(MASS)
library(gridGraphics)
library(ggpubr)


#### Step 1. Read in the data
load("./data/emp/hmsc_emp.RData")
################################################
#Figures for the raw data
data_abun_l <- pivot_longer(data_abun,cols=5:8,names_to = "Species",values_to = "Count")
data_abun_lm <- data_abun_l %>%
  group_by(Week, C.Light,C.Nutrients,Species) %>%
  summarise(Count = mean(Count))
data_abun_lm$C.Light <- as.factor(data_abun_lm$C.Light)
data_abun_lm$C.Nutrients <- as.factor(data_abun_lm$C.Nutrients)

f1<-ggplot(data_abun_lm,aes(x=Week,y=log(Count),color=C.Light,linetype = C.Nutrients))+
  geom_line()+
  ylab("Ln(Count)")+
  facet_wrap(~Species)+
  scale_x_continuous(breaks = c(3,5,7,9,11),labels=c(3,5,7,9,11))+
  theme_pubr()

data_pheno_dl <- pivot_longer(data_pheno_d,cols=5:8,names_to = "Species",values_to = "Deltatrait")
data_pheno_dlm <- data_pheno_dl %>%
  group_by(Week, C.Light,C.Nutrients,Species) %>%
  summarise(Deltatrait = mean(Deltatrait))
data_pheno_dlm$C.Light <- as.factor(data_pheno_dlm$C.Light)
data_pheno_dlm$C.Nutrients <- as.factor(data_pheno_dlm$C.Nutrients)

f2<-ggplot(data_pheno_dlm,aes(x=Week,y=log(Deltatrait),color=C.Light,linetype = C.Nutrients))+
  geom_line()+
  ylab("Ln(Change in frond size)")+
  facet_wrap(~Species)+
  scale_x_continuous(breaks = c(3,5,7,9,11),labels=c(3,5,7,9,11))+
  theme_pubr()

f3<-ggplot(data_pheno_dlm,aes(x=Week,y=Deltatrait,color=C.Light,linetype = C.Nutrients))+
  geom_line()+
  ylab("Change in frond size")+
  facet_wrap(~Species)+
  scale_x_continuous(breaks = c(3,5,7,9,11),labels=c(3,5,7,9,11))+
  theme_pubr()

data_pheno_l <- pivot_longer(data_pheno,cols=5:8,names_to = "Species",values_to = "Trait")
data_pheno_m <- data_pheno_l %>%
  group_by(Week, C.Light,C.Nutrients,Species) %>%
  summarise(Trait = mean(Trait))

data_pheno_m$C.Nutrients <- as.factor(data_pheno_m$C.Nutrients)
data_pheno_m$C.Light <- as.factor(data_pheno_m$C.Light)

f4<-ggplot(data_pheno_m,aes(x=Week,y=Trait,color=C.Light,linetype = C.Nutrients))+
  geom_line()+
  facet_wrap(~Species)+
  ylab("Frond size (cm^2)")+
  scale_x_continuous(breaks = c(1,3,5,7,9,11),labels=c(1,3,5,7,9,11))+
  theme_pubr()

d2<-ggarrange(f1, f2, f3, f4, ncol=2, nrow=2, common.legend = TRUE, legend="top",labels="auto")

ggsave("./output/fig_S4_1.pdf",d2,width=10,height=7,limitsize=F)
###################################################################################
###################################################################################

data_abun_lm$Deltatrait <- data_pheno_dlm$Deltatrait
library(patchwork)


data <- pivot_wider(data_abun_lm,names_from = Species, values_from = c(Count,Deltatrait))


p1<-ggplot(data,aes(x=Deltatrait_Rn,y=log(Count_Rn),color=C.Light,shape = C.Nutrients))+
  geom_point()+
  ylab("Ln(count Rn)")+
  xlab("dRn")+
  theme_pubr()

p2<-ggplot(data,aes(x=Deltatrait_Wc,y=log(Count_Lm),color=C.Light,shape = C.Nutrients))+
  geom_point()+
  theme_pubr()+
  ylab("Ln(count Lm)")+
  xlab("dWc")

p3<-ggplot(data,aes(x=Deltatrait_Lm,y=log(Count_Sp),color=C.Light,shape = C.Nutrients))+
  geom_point()+
  theme_pubr()+
  ylab("Ln(count Sp)")+
  xlab("dLm")

p4<-ggplot(data,aes(x=Deltatrait_Sp,y=log(Count_Sp),color=C.Light,shape = C.Nutrients))+
  geom_point()+
  theme_pubr()+
  ylab("Ln(count Sp)")+
  xlab("dSp")

d3<-ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="top",labels="auto")

ggsave("./output/fig_S4_2.pdf",d3,width=10,height=7,limitsize=F)
##################################################################################
#
m.post = Hmsc::convertToCodaObject(m.1_evo)
m.df <- as.data.frame(rbind(m.post$Beta[[1]],m.post$Beta[[2]],m.post$Beta[[3]],m.post$Beta[[4]]))

#Checking the beta effect of dLm on Sp
#With pre-checking we could find that both species only co-occur in patches patch 4,7,41, 46
par(mai=c(0,0,0,0),oma=c(2,2,0.5,0.5),mar=c(0,0,0,0))
#For beta of |deltax| species 15 on species 12
#dRn on Rn
p5<-bayesplot::mcmc_areas(as.data.frame(m.df$`B[dRn (C7), Rn (S2)]` )) + labs(x="Posterior distribution",y="Density")+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  ggtitle("Estimated beta effect of dRn on density of Rn")

#dWc on Lm
p6<-bayesplot::mcmc_areas(as.data.frame(m.df$`B[dWc (C9), Lm (S1)]`)) + labs(x="Posterior distribution",y="Density")+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  ggtitle("Estimated beta effect of dWc on density of Lm")

#dLm on Sp
p7<-bayesplot::mcmc_areas(as.data.frame(m.df$`B[dLm (C6), Sp (S3)]`)) + labs(x="Posterior distribution",y="Density")+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  ggtitle("Estimated beta effect of dLm on density of Sp")

#dSp on Sp
p8<-bayesplot::mcmc_areas(as.data.frame(m.df$`B[dSp (C8), Sp (S3)]`)) + labs(x="Posterior distribution",y="Density")+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  ggtitle("Estimated beta effect of dSp on density of Sp")


d4<-ggarrange(p5, p6, p7, p8, ncol=2, nrow=2, common.legend = TRUE, legend="top",labels="auto")

ggsave("./output/fig_S4_3.pdf",d4,width=10,height=7,limitsize=F)