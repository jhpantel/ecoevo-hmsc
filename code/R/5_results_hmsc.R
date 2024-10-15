# Title. results.R
# Author. J.H. Pantel
# Date 28-08-2023
# Description. This script contains analyses for the HMSC_ecoevo project and associated manuscript Hermann & Pantel 202x.

#### Libraries ---------------------------------------------------------------
library(bayesplot)
library(tidyverse)
library(gridExtra)

#### Step 1. Evaluate posteriors --------------
m.post = Hmsc::convertToCodaObject(m.1)
m.df <- as.data.frame(rbind(m.post$Beta[[1]],m.post$Beta[[2]]))

#### Create Variance Paritioning plot
s <- ncol(m.1$Y)
VP <- computeVariancePartitioning(m.1, group = c(1,2,rep(3,s)), groupnames = c("Intercept","Env","deltaX"))
#pdf("./output/vp2.pdf",width=8,height=4)
plotVariancePartitioning(m.1, VP = VP, args.legend=list(cex=0.75,bg="transparent"))
#dev.off()

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
# All values
bayesplot::mcmc_areas(B.int)
# All values, colored by overlap of HPDI with 0
xc<-NULL
for(i in 1:ncol(B.int)) xc[i]<-as.numeric(sum(colnames(B.int)[i]==not0))+1
theme_set(theme_grey())
bayesplot::mcmc_areas(B.int,rhat=xc) + legend_none() +  theme(text = element_text(size = 20))

# N ~ Nt-1 values ##########
# All values
B.N <- m.df[,2:(1+s)]
for(i in 1:(s-1)){
  B.N <- cbind(B.N,m.df[,2+(i*npred):(s+(i*npred-1))])
}
bayesplot::mcmc_areas(B.N)
# All values, colored by overlap of HPDI with 0
xd=NULL
for (k in 1:s) xd[k] <- as.numeric(sum(colnames(B.N)[k]==not0))+1
for(i in 1:(s-1)){
  test <- m.df[,2+(i*npred):(s+(i*npred-1))]
  for(k in 1:s) xc[k]<-as.numeric(sum(colnames(test)[k]==not0))+1
  xd <- c(xd,xc)
}
bayesplot::mcmc_areas(B.N,rhat=xd) + legend_none()

# by species
id <- as.numeric(sub("y","",colnames(m.1$Y)))
plist <- list()
for(i in 1:s){
  blah <- B.N[,(seq(i,dim(B.N)[2],by=s))]
  assign(paste0("B.", id[i], sep = ""), blah)
  plist[[i]] <- bayesplot::mcmc_areas(eval(parse(text=paste0("B.", id[i], sep = "")))) +legend_none()
}
gridExtra::grid.arrange(grobs=plist,nrow=round(s/3))

# N ~ Environment ##########
B.E <- m.df[,13]
for(i in 1:(s-1)){
  B.E <- cbind(B.E,m.df[,13+(i*npred)])
  colnames(B.E)[i+1] <- paste0("E", (i+1), sep = "")
}
colnames(B.E)[1] <- "E1"
bayesplot::mcmc_areas(B.E)

# N ~ |deltaX| values ####
B.dX <- m.df[,14:(14+s-1)]
for(i in 1:(s-1)){
  B.dX <- cbind(B.dX,m.df[,((14+(i*npred)):(14+(i*npred)+s-1))])
}
bayesplot::mcmc_areas(B.dX)
# All values, colored by overlap of HPDI with 0
xd=NULL
for(k in 1:s) xd[k]<-as.numeric(sum(colnames(B.dX)[k]==not0))+1
for(i in 1:(s-1)){
  test <- m.df[,((14+(i*npred)):(14+(i*npred)+s-1))]
  for(k in 1:s) xc[k]<-as.numeric(sum(colnames(test)[k]==not0))+1
  xd <- c(xd,xc)
}
bayesplot::mcmc_areas(B.dX,rhat=xd)+legend_none()

# by species
plist.dX <- list()
for(i in 1:s){
  blah <- B.dX[,(seq(i,dim(B.dX)[2],by=s))]
  assign(paste0("B.dX.", id[i], sep = ""), blah)
  plist.dX[[i]] <- bayesplot::mcmc_areas(eval(parse(text=paste0("B.dX.", id[i], sep = ""))))
}
gridExtra::grid.arrange(grobs=plist.dX,nrow=round(s/3))
# All values, colored by overlap of HPDI with 0
id <- as.numeric(sub("y","",colnames(m.1$Y)))
plist <- list()
r=NULL
for(i in 1:s){
  blah <- B.N[,(seq(i,dim(B.N)[2],by=s))]
  assign(paste0("B.", id[i], sep = ""), blah)
  for(k in 1:s) r[k]<-as.numeric(sum(colnames(blah)[k]==not0))+1
  plist[[i]] <- bayesplot::mcmc_areas(eval(parse(text=paste0("B.", id[i], sep = ""))),rhat=r) +legend_none()
}
gridExtra::grid.arrange(grobs=plist,nrow=round(s/3))

# N ~ Nt-1 x |deltaX| values ####
B.NdX <- m.df[,25:(25+s-1)]
for(i in 1:(s-1)){
  B.NdX <- cbind(B.NdX,m.df[,((25+(i*npred)):(25+(i*npred)+s-1))])
}
bayesplot::mcmc_areas(B.NdX)
# by species
plist.NdX <- list()
for(i in 1:s){
  blah <- B.NdX[,(seq(i,dim(B.NdX)[2],by=s))]
  assign(paste0("B.NdX.", id[i], sep = ""), blah)
  plist.NdX[[i]] <- bayesplot::mcmc_areas(eval(parse(text=paste0("B.NdX.", id[i], sep = ""))))
}
gridExtra::grid.arrange(grobs=plist.NdX,nrow=round(s/3))
# All values, colored by overlap of HPDI with 0
xd=NULL
for(k in 1:s) xd[k]<-as.numeric(sum(colnames(B.NdX)[k]==not0))+1
for(i in 1:(s-1)){
  test <- m.df[,((14+(i*npred)):(14+(i*npred)+s-1))]
  for(k in 1:s) xc[k]<-as.numeric(sum(colnames(test)[k]==not0))+1
  xd <- c(xd,xc)
}
bayesplot::mcmc_areas(B.NdX,rhat=xd) + legend_none()
# by species
## RH, I cannot get this to work. Error with R. ##
plist.NdX <- list()
for(i in 1:s){
  blah <- B.NdX[,(seq(i,dim(B.NdX)[2],by=s))]
  assign(paste0("B.NdX.", id[i], sep = ""), blah)
  for(k in 1:(s)) r[k]<-as.numeric(sum(colnames(blah)[k]==not0))+1
  plist.NdX[[i]] <- bayesplot::mcmc_areas(eval(parse(text=paste0("B.NdX.", id[i], sep = ""))),rhat=r) + legend_none()
}
gridExtra::grid.arrange(grobs=plist.NdX,nrow=round(s/3))
## RH, I cannot get this to work. Error with R. ##

# All large coefficients where HPDI doesn't include 0 ####
B.giant <- B.not0[,colMeans(B.not0) > 0.5 | colMeans(B.not0) < -0.5]
bayesplot::mcmc_areas(B.giant)

#### Step 2. Evaluate posteriors --------------
m.post = Hmsc::convertToCodaObject(m.1)