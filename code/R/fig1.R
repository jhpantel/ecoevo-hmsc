## Code for Figure 1
## J.H. Pantel
## 13-10-2024

## libraries
library(ecoevor)
library(Hmsc)
library(ggplot2)
library(patchwork)

## set random seed
set.seed(42)

# 1a. Sim: One species, logistic growth + environmental covariate ---------
## Sim1. One species, logistic growth
# Simulate initial species population growth
N1.0 <- 10
r1.0 <- 1.67
alpha.11 <- 0.00125
# Simulation of model for t time steps
t <- 30
N <- rep(NA, t)
N[1] <- N1.0
for (i in 2:t) {
  N[i] <- disc_log(r=r1.0, N0=N[i-1], alpha=alpha.11)
}
## Sim2. One species, logistic growth, fluctuating environmental covariate
# Simulate initial species population growth with environment fluctuations
N1.0 <- 10
r1.0 <- 1.67
alpha.11 <- 0.00125
E.0 <- 0.8
x1.0 <- 0.8
# Simulation of model for t time steps
t <- 40
N <- rep(NA, t)
N[1] <- N1.0
E <- rep(NA, t)
E[1] <- E.0
for (i in 2:t) {
  N[i] <- disc_log_E(r=r1.0, N0=N[i-1], alpha=alpha.11, E=E[i-1], x=x1.0)
  E[i] <- E[i-1] + rnorm(1,0,.1)
}



