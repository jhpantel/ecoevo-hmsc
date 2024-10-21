## Code for Figure 1
## J.H. Pantel
## 13-10-2024

## libraries
library(ecoevor)
library(Hmsc)
library(ggplot2)
library(patchwork)
library(bayesplot)
library(reshape2)
library(cowplot)

## set random seed
set.seed(8)

# 1. Sim: Two species, logistic growth + competition, environmental covariate ---------
### Sim1. Two species, Leslie-Gower competition, environmental covariate
# Initial conditions
N0 <- c(10,10)
r <- c(1.7,1.7)
alpha.11 <- 0.01
alpha.22 <- 0.01
alpha.12 <- 0.005
alpha.21 <- 0.01
alpha <- matrix(c(alpha.11,alpha.21,alpha.12,alpha.22),nrow=2,byrow=FALSE) # careful to make sure you have correct interaction matrix
E.0 <- 0.8
x <- c(0.6,0.8)
# Simulation of model for t time steps
t <- 40
N <- array(NA,dim=c(t,length(r)))
N <- as.data.frame(N)
colnames(N) <- paste0("N",1:length(r))
N[1,] <- N0
E <- rep(NA, t)
E[1] <- E.0
for (i in 2:t) {
  N[i,] <- disc_LV_E(r=r,N0=N[i-1,],alpha=alpha,E=E[i-1],x=x)
  E[i] <- E[i-1] + rnorm(1,0,0.1)
}

N$time <- 1:t
dat <- melt(N, id.vars="time")
ggplot2::ggplot(dat, aes(time, value, col=variable)) + geom_point() + geom_hline(yintercept = ((r[1] - 1)/alpha.11),linetype = "dashed", color = "gray")

### Mod1: HMSC model fit 1 ###
# prepare data in HMSC format
dat <- as.data.frame(log(N))
dat$time <- 1:t
df <- data.frame(dat[(2:t),-3])
Y <- as.matrix(df)
XData <- data.frame(dat[1:(t-1),-3],E[1:(t-1)],E[1:(t-1)]^2)
colnames(XData) <- c(paste0("n",1:length(r)),"E","Esq")
studyDesign = data.frame(sample = as.factor(1:(t-1)))
rL = HmscRandomLevel(units = studyDesign$sample)
m.1.hmsc = Hmsc(Y = Y, XData = XData, XFormula = ~E+Esq, studyDesign = studyDesign, ranLevels = list(sample = rL))
# Bayesian model parameters
nChains <- 2
thin <- 5
samples <- 1000
transient <- 500*thin
verbose <- 500*thin
# sample MCMC
m.1.sample <- sampleMcmc(m.1.hmsc,thin=thin,sample=samples,transient=transient,nChains=nChains,verbose=verbose)
# prepare for plot
pred <- predict(m.1.sample)
hold1 <- list()
hold2 <- list()
for(i in 1:length(r)){
  hold1[[i]] <- as.numeric(Y[,i])
  data <- lapply(pred, function(x) x[, colnames(df)[i], drop = FALSE])
  hold2[[i]] <- matrix(unlist(data), nrow=length(pred), byrow=TRUE)
}
y <- unlist(hold1) # vector of each species' observed N across time points
yrep <- as.matrix(data.frame(hold2))
### Mod1 plot ###
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(2,2,1,1))
color_scheme_set("brightblue")
#ppc_dens_overlay(y, yrep[1:50, ])
a <- ppc_intervals(y,yrep,x=rep(dat$time[1:(t-1)],times=2),prob = 0.95) +
  labs(x = "time",y = "ln(N)",) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
### Gradient plot ###
Gradient <- constructGradient(m.1.sample,focalVariable="E",non.focalVariables=list(Esq=list(2)),ngrid=39)
Gradient$XDataNew$Esq <- Gradient$XDataNew$E^2
predY <- predict(m.1.sample,XData=Gradient$XDataNew,expected=TRUE)
plotGradient(m.1.sample,Gradient,pred=predY,showData=T,measure="Y",index=1,main="",xlab="E_t",ylab="predicted N1_t+1",showPosteriorSupport=FALSE,cex.axis=0.75)
b <- recordPlot()
plotGradient(m.1.sample,Gradient,pred=predY,showData=T,measure="Y",index=2,main="",xlab="E_t",ylab="predicted N2_t+1",showPosteriorSupport=FALSE,cex.axis=0.75)
c <- recordPlot()
### Beta posterior plot ###
m.post = Hmsc::convertToCodaObject(m.1.sample)
m.beta <- as.data.frame(rbind(m.post$Beta[[1]],m.post$Beta[[2]]))
colnames(m.beta) <- c("Int1","E1","Esq1","Int2","E2","Esq2")
d <- bayesplot::mcmc_areas(m.beta) + scale_x_continuous(limits=c(-20,30))
#postBeta = getPostEstimate(m.1.sample, parName = "Beta")
#plotBeta(m.1.sample, post = postBeta, param = "Mean", supportLevel = 0.7)
#d <- recordPlot()
### variance partitioning ###
VP <- computeVariancePartitioning(m.1.sample, group = c(1, 1, 1), groupnames = "Env")
plotVariancePartitioning(m.1.sample, VP, args.legend = list(cex = 0.6, bg = "transparent"),cex.axis=0.75)
e <- recordPlot()
### ### ### ### ### ### ###
### all plots together ###
### ### ### ### ### ### ###
plot_grid(a,plot_grid(b,c,nrow=2),d,e,nrow=1,rel_widths=c(2,1,1,1))

# 2. Sim: Two species, logistic growth + competition, environmental covariate, evolution ---------
### Sim2. Two species, Leslie-Gower competition, environmental covariate, evolution
# Initial conditions
N0 <- c(10,10)
alpha.11 <- 0.01
alpha.22 <- 0.01
alpha.12 <- 0.005
alpha.21 <- 0.01
alpha <- matrix(c(alpha.11,alpha.21,alpha.12,alpha.22),nrow=2,byrow=FALSE) # careful to make sure you have correct interaction matrix
E.0 <- 0.8
x <- c(0.1,0.8)
P <- 1
w <- 2
Wmax <- 2
h2 <- 1
k <- (w + (1 - h2) * P)/(P + w)





