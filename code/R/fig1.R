## Code for Figure 1
## J.H. Pantel
## 13-10-2024

## libraries
library(ecoevor)
library(Hmsc)
library(ggplot2)
library(patchwork)
library(bayesplot)

## set random seed
set.seed(42)

# 1a. Sim: One species, logistic growth + environmental covariate ---------
### Sim1. One species, logistic growth ###
# Simulate initial species population growth
N1.0 <- 10
r1.0 <- 1.67
alpha.11 <- 0.00125
# Simulation of model for t time steps
t <- 40
N <- rep(NA, t)
N[1] <- N1.0
for (i in 2:t) {
  N[i] <- disc_log(r=r1.0, N0=N[i-1], alpha=alpha.11)
}
### Mod1: HMSC model fit 1 ###
# prepare data in HMSC format
# Plot simulation: ggplot
dat <- as.data.frame(N)
dat$time <- 1:t
Y <- as.matrix(log(dat$N[2:t]))
XData <- data.frame(x=log(dat$N[1:(t-1)]))
m.1.hmsc <- Hmsc(Y=Y,XData=XData,XFormula=~x)
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
y <- as.numeric(XData[,1])
yrep <- matrix(unlist(pred), nrow=length(pred), byrow=TRUE)
## Sim2. One species, logistic growth, fluctuating environmental covariate
# Simulate initial species population growth with environment fluctuations
N1.0 <- 10
r1.0 <- 1.67
alpha.11 <- 0.00125
E.0 <- 0.8
x1.0 <- 0.8
# Simulation of model for t time steps
t <- 40
N2 <- rep(NA, t)
N2[1] <- N1.0
E <- rep(NA, t)
E[1] <- E.0
for (i in 2:t) {
  N2[i] <- disc_log_E(r=r1.0, N0=N2[i-1], alpha=alpha.11, E=E[i-1], x=x1.0)
  E[i] <- E[i-1] + rnorm(1,0,.1)
}
### Mod2: HMSC model fit 1 ###
# prepare data in HMSC format
dat <- as.data.frame(cbind(N2,E))
dat$time <- 1:t
df <- data.frame(cbind(log(dat$N2[2:t]),log(dat$N2[1:(t-1)]),E[1:(t-1)],E[1:(t-1)]^2))
colnames(df) <- c("Nt1","Nt","E","Esq")
Y <- as.matrix(log(dat$N2[2:t]))
XData <- df
m.2.hmsc <- Hmsc(Y=Y,XData=XData,XFormula=~Nt+E+Esq)
# Bayesian model parameters
nChains <- 2
thin <- 5
samples <- 1000
transient <- 500*thin
verbose <- 500*thin
# sample MCMC
m.2.sample <- sampleMcmc(m.2.hmsc,thin=thin,sample=samples,transient=transient,nChains=nChains,verbose=verbose)
# prepare for plot
pred <- predict(m.2.sample)
y <- c(y,as.numeric(XData[,1]))
yrep2 <- matrix(unlist(pred), nrow=length(pred), byrow=TRUE)
yrep <- cbind(yrep,yrep2)
### Mod1 and 2 plot ###
color_scheme_set("brightblue")
#ppc_dens_overlay(y, yrep[1:50, ])
ppc_intervals(y,yrep,x=rep(dat$time[1:(t-1)],times=2),prob = 0.95) +
  labs(x = "time",y = "N",)





