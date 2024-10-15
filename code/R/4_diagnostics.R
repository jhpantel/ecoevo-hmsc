# Title. diagnostics.R
# Author. R. Hermann
# Date 25-08-2023
# Description. This script contains analyses for the HMSC_ecoevo project and associated manuscript Hermann & Pantel 202x.

#setwd("./output")
#load("h_01_d_minus3_hmsc_UPDATE.RData")

#### Libraries ---------------------------------------------------------------
library(Hmsc)
library(tidyverse)
library(bayesplot)

#### Checking the MCMC chains --------------
#### Step 1 Create the necessary data.frame --------------

#To determine the convergence of the chains we need to check two things, first off how we will check how each chain converges on itself (hence for example no autocorrelation) and second how well they converge both chains converge with each other. This can be done with the effective sample size and potential scale reduction factors.
m.post <- convertToCodaObject(m.1)
m.df <- as.data.frame(rbind(m.post$Beta[[1]],m.post$Beta[[2]]))
ess.beta<- effectiveSize(m.post$Beta)
psrf.beta <- gelman.diag(m.post$Beta, multivariate=FALSE)$psrf

s <- ncol(m.1$Y)
npred <- ncol(m.1$XData)+1
sp<-colnames(m.1$Y)
sp_n<-as.numeric(sub("y","",colnames(m.1$Y)))
VP <- computeVariancePartitioning(m.1, group = c(1,2,rep(3,s)), groupnames = c("Intercept","Env","deltaX"))
pdf("./output/vp2.pdf",width=8,height=4)
plotVariancePartitioning(m.1, VP = VP, args.legend=list(cex=0.75,bg="transparent"))
dev.off()

#### Step 2 Creatin simpple histograms --------------
#Creating a histogram of the effective sample size (ess) of all Beta parameters
#ess shows how many data points of the MCMC chain were independently and randomly drawn and as such can be effectively used
hist(ess.beta,xlab="Effective sample site",main="")

#Plotting the potential scale reduction factors (psrf) of all beta parameters
#psrf from the Gelman-Rubin diagnostic is based a comparison of within-chain and between-chain variances, and is similar to a classical analysis of variance
hist(psrf.beta,xlab="Potential scale reduction factos",main="")

#### Step 3 Beta separation by fixed effects --------------

#Per Species we have 47 fixed effects, which are split up accordingly:
#1 = intercept
#2-16 = impact of $N_{j,t-1}$
#17 = impact of environment
#18-32 = impact of $x_{j,t-1}$
#33-47 = impact of $x_{j,t-1} * N_{j,t-1}$

#For ess
ess_effects <- data.frame(ess=0,effect="0")
for (i in 1:s) {
  x <- i*npred
  y <- x-(npred-1)
  ess_effects[y,1] <- unname(ess.beta)[y]
  ess_effects[y,2] <- "intercept"
  ess_effects[(y+1),1] <- unname(ess.beta)[(y+1)]
  ess_effects[(y+1),2] <- "environment"
  ess_effects[(y+2):(y+s+1),1] <- unname(ess.beta)[(y+2):(y+s+1)]
  ess_effects[(y+2):(y+s+1),2] <- rep("x_t-1",s)
}


ggplot( ess_effects,aes(x=effect,y=ess))+
  geom_boxplot(aes(fill=effect))+
  ylab("Effective sanple size")+xlab("Fixed effects")

#For psrf
psrf_effects <- data.frame(psrf=0,effect="0")
for (i in 1:s) {
  x <- i*npred
  y <- x-(npred-1)
  psrf_effects[y,1] <- unname(psrf.beta)[y,1]
  psrf_effects[y,2] <- "intercept"
  psrf_effects[(y+1),1] <- unname(psrf.beta)[(y+1),1]
  psrf_effects[(y+1),2] <- "environment"
  psrf_effects[(y+2):(y+s+1),1] <- unname(psrf.beta)[(y+2):(y+s+1),1]
  psrf_effects[(y+2):(y+s+1),2] <- rep("x_t-1",s)
}

ggplot(psrf_effects,aes(x=effect,y=psrf))+
  geom_boxplot(aes(fill=effect))+
  ylab("Potenital scale reduction factors")+xlab("Fixed effects")

#### Step 4 Beta separation by species --------------

#We have 15 species with each 47 different effects, so in order to plot the Beta diagnostics we use the following code.

#For ess
ess_species <- data.frame(ess=0,species="0")
for (i in 1:s) {
  x <- i*npred
  y <- x-(npred-1)
  ess_species[y:x,1]<- unname(ess.beta)[y:x]
  if (sp_n[i]<10)
    ess_species[y:x,2]<- rep(paste("y",0,sp_n[i],sep=""),length(m.post$Beta[y:x]))
  if (sp_n[i]>9)
    ess_species[y:x,2]<- rep(paste("y",sp_n[i],sep=""),length(m.post$Beta[y:x]))
  
}

ggplot(ess_species,aes(x=species,y=ess))+
  geom_boxplot(aes(fill=species))+
  ylab("Effective sample size")+xlab("Species")

#For psrf
psrf_species <- data.frame(psrf=0,species="0")
for (i in 1:s) {
  x <- i*npred
  y <- x-(npred-1)
  psrf_species[y:x,1]<- unname(psrf.beta)[y:x,1]
  if (i<10)
    psrf_species[y:x,2]<- rep(paste("y",0,sp_n[i],sep=""),length(m.post$Beta[y:x]))
  if (i>9)
    psrf_species[y:x,2]<- rep(paste("y",sp_n[i],sep=""),length(m.post$Beta[y:x]))
}

ggplot(psrf_species,aes(x=species,y=psrf))+
  geom_boxplot(aes(fill=species))+
  ylab("Potenital scale reduction factors")+xlab("Species")

#### Step 5 Mean and SD of the Beta --------------
B.mean = NULL
B.sd = NULL
m.df <- as.data.frame(rbind(m.post$Beta[[1]],m.post$Beta[[2]]))
B.N <- m.df[,2:(1+s)]
for(i in 1:(s-1)){
  B.N <- cbind(B.N,m.df[,2+(i*npred):(s+(i*npred-1))])
}
for (i in 1:s) {
  x <- B.N[,(seq(i,dim(B.N)[2],by=s))]
  m <- apply(x,2,mean)
  st <- apply(x,2,sd)
  n<-names(m)
  for (j in 1:s){
    names(m)[j] <- str_split_fixed(n[j]," ",n=4)[4]
    names(st)[j] <- str_split_fixed(n[j]," ",n=4)[4]
  }
  me<-data.frame(Sp_X=names(m),Beta=unname(m))
  sta<-data.frame(Sp_X=names(st),Beta=unname(st))
  me$Sp_y  <- paste("y",sp_n[i],sep="")
  sta$Sp_y <- paste("y",sp_n[i],sep="")
  B.mean <- rbind(B.mean,me)
  B.sd <- rbind(B.sd,sta)
}


ggplot(B.mean,aes(x=Sp_y,y=Beta))+
  geom_boxplot()+
  scale_x_discrete(limits=sp)+
  ggtitle("Mean of Beta")


ggplot(B.sd,aes(x=Sp_y,y=Beta))+
  geom_boxplot()+
  scale_x_discrete(limits=sp)+
  ggtitle("SD of Beta")

#### Step 6 Further detailed checks --------------
#To check autocorrelation on specific predictors
#Environmental influence on all species
mcmc_acf(m.post$Beta[,seq(2,(s*npred),by=npred)])
#Intercept
mcmc_acf(m.post$Beta[,seq(1,(s*npred),by=npred)])
#Plotting xt-1 from species 15 on all other species
mcmc_acf(m.post$Beta[,seq((s+2),(s*npred),by=npred)])

#To plot specific MCMC iterations
#Environmental influence on all species
plot(m.post$Beta[,seq(2,(s*npred),by=npred)])
#Intercept
plot(m.post$Beta[,seq(1,(s*npred),by=npred)])
#Plotting xt-1 from species 15 on all other species
plot(m.post$Beta[,seq((s+2),(s*npred),by=npred)])


#To plot the Nt-1 Betas for a specific species
#B.sp1 <- B.N[,(seq(1,dim(B.N)[2],by=s))]
#bayesplot::mcmc_areas(B.sp1)

