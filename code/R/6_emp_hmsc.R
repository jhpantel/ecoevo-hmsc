# Title. 6_emp_hmsc.R
# Author. J.H. Pantel &  R.J. Hermann
# Date 29.07.2025
# Description. This script contains analyses for empircial data of the HMSC_ecoevo project and associated manuscript Pantel & Hermann 2025.

#### Libraries ---------------------------------------------------------------
library(tidyverse)
library(Hmsc)
library(tictoc)
library(statgenSTA)
library(statgenGxE)

#### Step 0. Seeting a random seed to make the result reproducable
set.seed(77)

#### Step 1. Loading the empirical data from the publication Jewell & Bell 2023 --------------
dat_ab <- read.csv("./data/emp/Phase1_abundance.csv")
dat_ph <- read.csv("./data/emp/Phase1_phenotypes.csv")

#Deleting abundance data of Week 14 and Week 0, because they did not measure the FrondArea on that Weeks
dat_ab <- subset(dat_ab,Week<14)
dat_ab <- subset(dat_ab,Week>0)

#Changing the H(for high), M (for moderate) and L (for low) for the actual values of the experiment
#C.Light is the amount of layers of 50% shade cloth (described as percentage in the publication)
dat_ph$C.Light[dat_ph$C.Light=="L"]<- 3
dat_ph$C.Light[dat_ph$C.Light=="M"]<- 12
dat_ph$C.Light[dat_ph$C.Light=="H"]<- 50

dat_ab$C.Light[dat_ab$C.Light=="L"]<- 3
dat_ab$C.Light[dat_ab$C.Light=="M"]<- 12
dat_ab$C.Light[dat_ab$C.Light=="H"]<- 50

#C.Nutrients is the concentrations of added nutrients, we use her the µg/L concentrations for dissolved phosphor.
dat_ph$C.Nutrients[dat_ph$C.Nutrients=="L"]<- 10
dat_ph$C.Nutrients[dat_ph$C.Nutrients=="M"]<- 40
dat_ph$C.Nutrients[dat_ph$C.Nutrients=="H"]<- 160

dat_ab$C.Nutrients[dat_ab$C.Nutrients=="L"]<- 10
dat_ab$C.Nutrients[dat_ab$C.Nutrients=="M"]<- 40
dat_ab$C.Nutrients[dat_ab$C.Nutrients=="H"]<- 160

#Setting them to numeric, as they were characters before
dat_ab$C.Light <- as.numeric(dat_ab$C.Light)
dat_ab$C.Nutrients <- as.numeric(dat_ab$C.Nutrients)

dat_ph$C.Light <- as.numeric(dat_ph$C.Light)
dat_ph$C.Nutrients <- as.numeric(dat_ph$C.Nutrients)

#Calculating mean phenotype (FrondArea) for each species for each sample
dat_ph_m <-dat_ph %>%
  group_by(Week,C.Light,C.Nutrients,C.Rep,Species) %>%
  summarise(Frondmean = mean(FrondArea))

#Sorting the data matrices to all have the same order
ord_dat_ab <- dat_ab[order(dat_ab$Week,dat_ab$C.Light,dat_ab$C.Nutrients,dat_ab$Species,dat_ab$C.Rep),]
ord_dat_ph_m <- dat_ph_m[order(dat_ph_m$Week,dat_ph_m$C.Light,dat_ph_m$C.Nutrients,dat_ph_m$Species,dat_ph_m$C.Rep),]

data_abun<-pivot_wider(ord_dat_ab ,names_from = Species,values_from=Individuals)
data_pheno<-pivot_wider(ord_dat_ph_m,names_from = Species,values_from=Frondmean)

#### Step 2. Broad-sense heritability estimates --------------
## Load in reciprocal transplant data
dat_rt <- read.csv("./data/emp/ReciprocalTransplant.csv")
## Community and Environment as factors
dat_rt$C.Env <- as.factor(dat_rt$C.Env)
dat_rt$E.Env <- as.factor(dat_rt$E.Env)
## Create environmental value
dat_rt$Env.Val <- NA
Env.Val <- by(dat_rt,factor(dat_rt$E.Env),function(x){mean(x$FrondArea)})
for(i in 1:length(levels(factor(dat_rt$E.Env)))){
  dat_rt$Env.Val[dat_rt$E.Env==levels(factor(dat_rt$E.Env))[i]] <- as.numeric(Env.Val[i])
}
## Calculate phenotypic variance components in a linear mixed model
dat_rt$Env.Val <- as.factor(dat_rt$Env.Val)
dat_rt$Env.Val <- as.ordered(dat_rt$Env.Val)
Wc <- subset(dat_rt,grepl("Wc", Species))
Sp <- subset(dat_rt,grepl("Sp", Species))
Lm <- subset(dat_rt,grepl("Lm", Species))
Rn <- subset(dat_rt,grepl("Rn", Species))
## Broad-sense heritability estimates by species
#Wc
Wc.TD <- statgenSTA::createTD(data = Wc, genotype = "C.Env", trial = "Env.Val")
drops.Wc <- gxeVarComp(TD = Wc.TD, trait = "FrondArea")
var.Wc <- vc(drops.Wc)
h2.Wc <- var.Wc$Component[1] / (var.Wc$Component[1] + (var.Wc$Component[2] / 9) + (var.Wc$Component[3]/(9*2)))
#Sp
Sp.TD <- statgenSTA::createTD(data = Sp, genotype = "C.Env", trial = "Env.Val")
drops.Sp <- gxeVarComp(TD = Sp.TD, trait = "FrondArea")
var.Sp <- vc(drops.Sp)
h2.Sp <- var.Sp$Component[1] / (var.Sp$Component[1] + (var.Sp$Component[2] / 9) + (var.Sp$Component[3]/(9*2)))
#Lm
Lm.TD <- statgenSTA::createTD(data = Lm, genotype = "C.Env", trial = "Env.Val")
drops.Lm <- gxeVarComp(TD = Lm.TD, trait = "FrondArea")
var.Lm <- vc(drops.Lm)
h2.Lm <- var.Lm$Component[1] / (var.Lm$Component[1] + (var.Lm$Component[2] / 9) + (var.Lm$Component[3]/(9*2)))
#Rn
Rn.TD <- statgenSTA::createTD(data = Rn, genotype = "C.Env", trial = "Env.Val")
drops.Rn <- gxeVarComp(TD = Rn.TD, trait = "FrondArea")
var.Rn <- vc(drops.Rn)
h2.Rn <- var.Rn$Component[1] / (var.Rn$Component[1] + (var.Rn$Component[2] / 9) + (var.Rn$Component[3]/(9*2)))


#Calculating the Deltatrait value, as each week has exactly 18 datapoints (rows), I substract the trait value with its corresponding value 18 rows above.
data_pheno_d <- data_pheno
for (i in 2:6){
  y=i*18-17
  z=i*18
  data_pheno_d$Lm[y:z] = abs(data_pheno$Lm[y:z] - data_pheno$Lm[(y-18):(z-18)])
  data_pheno_d$Rn[y:z] = abs(data_pheno$Rn[y:z] - data_pheno$Rn[(y-18):(z-18)])
  data_pheno_d$Sp[y:z] = abs(data_pheno$Sp[y:z] - data_pheno$Sp[(y-18):(z-18)])
  data_pheno_d$Wc[y:z] = abs(data_pheno$Wc[y:z] - data_pheno$Wc[(y-18):(z-18)])
}

#Now multiplying the trait-change values with the heritability value to estimate the impact of evolution
data_pheno_d$Lm <- data_pheno_d$Lm*h2.Lm
data_pheno_d$Rn <- data_pheno_d$Rn*h2.Rn
data_pheno_d$Sp <- data_pheno_d$Sp*h2.Sp
data_pheno_d$Wc <- data_pheno_d$Wc*h2.Wc

#Squaring the environmental terms
data_pheno_d$Lightsq <- data_abun$C.Light^2
data_pheno_d$Nutrientsq <- data_abun$C.Nutrients^2

#Setting the first 18 values to NA as Deltatrait cannot be calculated for them
data_pheno_d$Lm[1:18] = NA
data_pheno_d$Rn[1:18] = NA
data_pheno_d$Sp[1:18] = NA
data_pheno_d$Wc[1:18] = NA

colnames(data_pheno_d)[5:8] <- c("dLm","dRn","dSp","dWc")

#Taking out the first week, as we do not have deltatrait for those
data_pheno_d<-subset(data_pheno_d, Week>1)
data_abun <- subset(data_abun,Week>1)

#Checking if both dataframes are ordered in the same way
nrow(data_abun)==nrow(data_pheno_d)
all(data_abun$C.Rep == data_pheno_d$C.Rep)
all(data_abun$C.Light == data_pheno_d$C.Light)
all(data_abun$C.Nutrients == data_pheno_d$C.Nutrients)
all(data_abun$Week == data_pheno_d$Week)

#Setting the Y data, studydesign and XData
Y <- data_abun[,5:8]
Y_log <- log(Y)
studyDesign <- data_abun[,1:4]
studyDesign$Week <- as.factor(data_abun$Week)
studyDesign$C.Light <- as.factor(data_abun$C.Light)
studyDesign$C.Nutrients <- as.factor(data_abun$C.Nutrients)
studyDesign$C.Rep <- as.factor(data_abun$C.Rep)
XData <- data_pheno_d[,c(2,3,5:10)]

studyDesign <- as.data.frame(studyDesign)
Y_log <- as.data.frame(Y_log)
XData <- as.data.frame(XData)

#Creating the rand effects and the 
rL.time <- Hmsc::HmscRandomLevel(units = unique(data_abun$Week))
rL.rep <- Hmsc::HmscRandomLevel(units = unique(data_abun$C.Rep))

#Determining the length of the MCMC
thin = 5000
samples = 250
transient = 125*thin
verbose = 5000
chains = 4
nParallel = 4

#Setting and running the evolution model
XFormula= ~C.Light + Lightsq + C.Nutrients + Nutrientsq + dLm + dRn + dSp + dWc
m <- Hmsc::Hmsc(Y=Y_log,XData=XData,studyDesign=studyDesign,XFormula=XFormula,ranLevels=list("Week"=rL.time,"C.Rep"=rL.rep),distr="normal",XScale=TRUE)

tic()
m.1_evo <- Hmsc::sampleMcmc(m,samples=samples,transient=transient,nChains=chains,verbose=verbose,thin=thin,updater = list(GammaEta = FALSE),nParallel=nParallel)
a <- toc()

save.image(file="./data/emp/hmsc_emp.RData")

#Setting and running the no-evolution model
XData_noevo<- XData[,c(1,2,7,8)]
XFormula_noevo= ~C.Light + Lightsq + C.Nutrients + Nutrientsq
m_noevo <- Hmsc::Hmsc(Y=Y_log,XData=XData_noevo,studyDesign=studyDesign,XFormula=XFormula_noevo,ranLevels=list("Week"=rL.time,"C.Rep"=rL.rep),distr="normal",XScale=TRUE)

m.1_noevo <- Hmsc::sampleMcmc(m_noevo,samples=samples,transient=transient,nChains=chains,verbose=verbose,thin=thin,updater = list(GammaEta = FALSE),nParallel=nParallel)
save.image(file="./data/emp/hmsc_emp.RData")

#Setting and running the null model
XData_zero<- XData_noevo[,-c(1:4)]
XFormula_zero= ~1
m_zero <- Hmsc::Hmsc(Y=Y_log,XData=XData_zero,studyDesign=studyDesign,XFormula=XFormula_zero,ranLevels=list("Week"=rL.time,"C.Rep"=rL.rep),distr="normal",XScale=TRUE)

m.1_zero <- Hmsc::sampleMcmc(m_zero,samples=samples,transient=transient,nChains=chains,verbose=verbose,thin=thin,updater = list(GammaEta = FALSE),nParallel=nParallel)
save.image(file="./data/emp/hmsc_emp.RData")
