# Title. emp_hmsc_v01.R
# Author. R.J. Hermann
# Date 24-08-2023
# Description. This script contains analyses for enmpircial data of the HMSC_ecoevo project and associated manuscript Hermann & Pantel 202x.

#### Libraries ---------------------------------------------------------------
library(tidyverse)
library(Hmsc)
library(tictoc)

#### Step 1. Loading the empirical data from the publication Jewell & Bell 2023 --------------
dat_ab <- read.csv("empirical data/Phase1_abundance.csv")
dat_ph <- read.csv("empirical data/Phase1_phenotypes.csv")

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

#Now I need to get the mean phenotype (FrondArea) for each species for each sample
dat_ph_m <-dat_ph %>%
  group_by(Week,C.Light,C.Nutrients,C.Rep,Species) %>%
  summarise(Frondmean = mean(FrondArea))

#Saving the data files
write.csv(dat_ab,"empirical data/abundance_ordered_numeric.csv",row.names = F)
write.csv(dat_ph,"empirical data/phenotype_ordered_numeric.csv",row.names = F)
write.csv(dat_ph_m,"empirical data/phenotype_ordered_m_numeric.csv",row.names = F)

#Here I ordered all of the files the same way in excel, as this was the fastest way for me to do it.
#They were sorted with the following sorting priorities: Week, C.Light, C.Nutrients, Species and C.Rep

#Now we read in the ordered files
dat_ab<-read.csv("empirical data/abundance_ordered_numeric.csv")
dat_ph<-read.csv("empirical data/phenotype_ordered_numeric.csv")
dat_ph_m<-read.csv("empirical data/phenotype_ordered_m_numeric.csv")

data_abun<-pivot_wider(dat_ab,names_from = Species,values_from=Individuals)
data_pheno<-pivot_wider(dat_ph_m,names_from = Species,values_from=Frondmean)

#Here I calculate the Deltatrait value, as each week has exactly 18 datapoints (rows), I substract the trait value with its corresponding value 18 rows above.
data_pheno_d <- data_pheno
for (i in 2:6){
  y=i*18-17
  z=i*18
  data_pheno_d$Lm[y:z] = abs(data_pheno$Lm[y:z] - data_pheno$Lm[(y-18):(z-18)])
  data_pheno_d$Rn[y:z] = abs(data_pheno$Rn[y:z] - data_pheno$Rn[(y-18):(z-18)])
  data_pheno_d$Sp[y:z] = abs(data_pheno$Sp[y:z] - data_pheno$Sp[(y-18):(z-18)])
  data_pheno_d$Wc[y:z] = abs(data_pheno$Wc[y:z] - data_pheno$Wc[(y-18):(z-18)])
}
#Setting the first 18 values to NA as I cannot calculate Deltatrait for them
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
XData <- data_pheno_d[,2:8]
XData <- select(XData, -C.Rep)

studyDesign <- as.data.frame(studyDesign)
YData <- as.data.frame(YData)
XData <- as.data.frame(XData)

#Creating the rand effects and the 
rL.time <- Hmsc::HmscRandomLevel(units = unique(data_abun$Week))
rL.rep <- Hmsc::HmscRandomLevel(units = unique(data_abun$C.Rep))

#Determing the length of the MCMC
thin = 5000
samples = 250
transient = 125*thin
verbose = 5000
chains = 4
nParallel = 4

#Setting and running the evolution model
XFormula= ~C.Light + C.Nutrients + dLm + dRn + dSp + dWc
m <- Hmsc::Hmsc(Y=Y_log,XData=XData,studyDesign=studyDesign,XFormula=XFormula,ranLevels=list("Week"=rL.time,"C.Rep"=rL.rep),distr="normal",XScale=TRUE)

tic()
m.1_evo <- Hmsc::sampleMcmc(m,samples=samples,transient=transient,nChains=chains,verbose=verbose,thin=thin,updater = list(GammaEta = FALSE),nParallel=nParallel)
a <- toc()

save.image(file="data/hmsc_emp_coll.RData")

#Setting and running the no-evolution model
XData_noevo<- XData[,1:2]
XFormula_noevo= ~C.Light + C.Nutrients
m <- Hmsc::Hmsc(Y=Y_log,XData=XData_noevo,studyDesign=studyDesign,XFormula=XFormula_noevo,ranLevels=list("Week"=rL.time,"C.Rep"=rL.rep),distr="normal",XScale=TRUE)

m.1_noevo <- Hmsc::sampleMcmc(m,samples=samples,transient=transient,nChains=chains,verbose=verbose,thin=thin,updater = list(GammaEta = FALSE),nParallel=nParallel)
save.image(file="data/hmsc_emp_coll.RData")
