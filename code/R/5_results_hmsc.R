# Title. 5_results_hmsc.R
# Author. J.H. Pantel &  R.J. Hermann
# Date 29.07.2025
# Description. This script contains analyses for the HMSC_ecoevo project and associated manuscript Pantel & Hermann 2025.

#### Libraries ---------------------------------------------------------------
library(bayesplot)
library(tidyverse)
library(gridExtra)
library(Hmsc)

#### Step 0 Load the data --------------
d_lev <- c("d_zero","d_minus9","d_minus8","d_minus7","d_minus6","d_minus5","d_minus4","d_minus3","d_minus2","d_01","d_02","d_03","d_04","d_05","d_06","d_07","d_08","d_09","d_10")
h_lev <- c("h_0_","h_01_","h_02_","h_03_","h_04_","h_05_","h_06_","h_07_","h_08_","h_09_","h_10_")

# HMSC for single results condition
z <- 4
f <- 8

print(paste(h_lev[z],d_lev[f],sep=""))
result <- paste(h_lev[z],d_lev[f],sep="")
load(paste("./data/mc/",result,"_hmsc.RData",sep=""))

#### Step 1. Evaluate posteriors --------------
m.post = Hmsc::convertToCodaObject(m.1)
m.df <- do.call(rbind, m.post$Beta)

#### Create Variance Partitioning plot
s <- ncol(m.1$Y)
VP <- computeVariancePartitioning(m.1, group = c(1,2,rep(3,s)), groupnames = c("Intercept","Env","deltaX"))
plotVariancePartitioning(m.1, VP = VP, args.legend=list(cex=0.75,bg="transparent"))

# All beta coefficients where HPDI doesn't include 0 ####
a <- apply(m.df,2,function(x) quantile(x,probs=c(0.025,0.975)))
b <- apply(a,2,function(x) prod(x[1],x[2]) < 0)
B.not0 <- m.df[,b==FALSE]
not0 <- colnames(B.not0)
bayesplot::mcmc_areas(B.not0)
not0 <- colnames(B.not0)

# Separate coefficient estimate posteriors into groups
s <- ncol(m.1$Y)
npred <- 1 + length(colnames(m.1$XData))

# Intercept values ##########
B.int <- m.df[,seq(1,length(colnames(m.df)),npred)]
# All values, colored by overlap of HPDI with 0
xc<-NULL
for(i in 1:ncol(B.int)) xc[i]<-as.numeric(sum(colnames(B.int)[i]==not0))+1
theme_set(theme_grey())
#Renaming it to E and the according species number
for(i in 1:(s-1)){
  #B.int <- cbind(B.int,m.df[,seq(1,length(colnames(m.df)),npred)])
  colnames(B.int)[i+1] <- paste0("Int", (i+1), sep = "")
}
colnames(B.int)[1] <- "Int1"
bayesplot::mcmc_areas(B.int,rhat=xc) + legend_none() +  theme(text = element_text(size = 20))

# N ~ Environment ##########
B.E <- m.df[,seq(2,length(colnames(m.df)),npred)]
# All values, colored by overlap of HPDI with 0
xc<-NULL
for(i in 1:ncol(B.E)) xc[i]<-as.numeric(sum(colnames(B.E)[i]==not0))+1
theme_set(theme_grey())
#Renaming it to E and the according species number
for(i in 1:(s-1)){
  #B.E <- cbind(B.E,m.df[,2+(i*npred)])
  colnames(B.E)[i+1] <- paste0("E", (i+1), sep = "")
}
colnames(B.E)[1] <- "E1"
bayesplot::mcmc_areas(B.E,rhat=xc) + legend_none() +  theme(text = element_text(size = 20))


# N ~ |deltaX| values ####
B.dX <- m.df[,3:(2+s)]
for(i in 1:(s-1)){
  B.dX <- cbind(B.dX,m.df[,((3+(i*npred)):(3+(i*npred)+s-1))])
}
bayesplot::mcmc_areas(B.dX)
# All values, colored by overlap of HPDI with 0
xd=NULL
for(k in 1:s) xd[k]<-as.numeric(sum(colnames(B.dX)[k]==not0))+1
for(i in 1:(s-1)){
  test <- m.df[,((3+(i*npred)):(3+(i*npred)+s-1))]
  for(k in 1:s) xc[k]<-as.numeric(sum(colnames(test)[k]==not0))+1
  xd <- c(xd,xc)
}
bayesplot::mcmc_areas(B.dX,rhat=xd)+legend_none()

# by species
id <- as.numeric(sub("y","",colnames(m.1$Y)))
plist <- list()
r=NULL
for(i in 1:s){
  blah <- B.dX[,(seq(i,dim(B.dX)[2],by=s))]
  assign(paste0("B.", id[i], sep = ""), blah)
  for(k in 1:s) r[k]<-as.numeric(sum(colnames(blah)[k]==not0))+1
  plist[[i]] <- bayesplot::mcmc_areas(eval(parse(text=paste0("B.", id[i], sep = ""))),rhat=r) +legend_none()
}
gridExtra::grid.arrange(grobs=plist,nrow=round(s/3))

# All large coefficients where HPDI doesn't include 0 ####
B.giant <- B.not0[,colMeans(B.not0) > 0.5 | colMeans(B.not0) < -0.5]
bayesplot::mcmc_areas(B.giant)
