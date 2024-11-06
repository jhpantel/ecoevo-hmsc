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
library(MASS)

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

### Plot of time series ###
gdat <- N
gdat$time <- 1:t
dat <- melt(gdat, id.vars="time")
pa1 <- ggplot2::ggplot(dat, aes(time, value, col=variable)) + geom_point() + geom_hline(yintercept = ((r[1] - 1)/alpha.11),linetype = "dashed", color = "gray")
### ###

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
par(mgp=c(2,0.45,0), tcl=-0.4, mar=c(1.3,1.2,0,0))
color_scheme_set("brightblue")
#ppc_dens_overlay(y, yrep[1:50, ])
p1 <- ppc_intervals(y,yrep,x=rep(dat$time[1:(t-1)],times=2),prob = 0.95) +
  labs(x = "time",y = "ln(N)",) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
### Gradient plot ###
Gradient <- constructGradient(m.1.sample,focalVariable="E",non.focalVariables=list(Esq=list(2)),ngrid=39)
Gradient$XDataNew$Esq <- Gradient$XDataNew$E^2
predY <- predict(m.1.sample,XData=Gradient$XDataNew,expected=TRUE)
plotGradient(m.1.sample,Gradient,pred=predY,showData=T,measure="Y",index=1,main="",xlab="E_t",ylab="predicted N1_t+1",showPosteriorSupport=FALSE,cex.axis=0.75,bty='n')
p2 <- recordPlot()
plotGradient(m.1.sample,Gradient,pred=predY,showData=T,measure="Y",index=2,main="",xlab="E_t",ylab="predicted N2_t+1",showPosteriorSupport=FALSE,cex.axis=0.75,bty='n')
p3 <- recordPlot()
### Beta posterior plot ###
m.post = Hmsc::convertToCodaObject(m.1.sample)
m.beta <- as.data.frame(rbind(m.post$Beta[[1]],m.post$Beta[[2]]))
colnames(m.beta) <- c("Int1","E1","Esq1","Int2","E2","Esq2")
p4 <- bayesplot::mcmc_areas(m.beta) + scale_x_continuous(limits=c(-20,30)) +
  theme(axis.text=element_text(size=8))
### variance partitioning ###
par(mgp=c(2,0.45,0), tcl=-0.4, mar=c(1.3,1.2,0.5,0.5))
VP <- computeVariancePartitioning(m.1.sample, group = c(1, 1, 1), groupnames = "Env")
plotVariancePartitioning(m.1.sample, VP, args.legend = list(cex = 0.6, bg = "transparent"),cex.axis=0.6,main="",cols=c("white","orange"))
p5 <- recordPlot()
### ### ### ### ### ### ###
### all plots together ###
### ### ### ### ### ### ###
a <- plot_grid(p1,plot_grid(p2,p3,nrow=2),p4,p5,nrow=1,rel_widths=c(2,1,1,1))

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
x.0 <- c(0.1,0.8)
P <- 1
w <- 2
Wmax <- 2
h2 <- 1
k <- (w + (1 - h2) * P)/(P + w)
# Simulation of model for t time steps
# Simulation of model for t time steps
t <- 40
N <- array(NA,dim=c(t,length(r)))
N <- as.data.frame(N)
colnames(N) <- paste0("N",1:length(r))
N[1,] <- N0
E <- rep(NA, t)
E[1] <- E.0
x <- array(NA,dim=c(t,length(r)))
x <- as.data.frame(x)
colnames(x) <- paste0("x",1:length(r))
x[1,] <- x.0
r <- array(NA,dim=c(t,length(r)))
r <- as.data.frame(r)
colnames(r) <- paste0("r",1:length(r))
What <- Wmax*sqrt(w/(P+w))
r[1,] <- What*exp((-(((w+(1-h2)*P)/(P+w))*(E[1]-x[1,]))^2)/(2*(P+w)))
for (i in 2:t) {
  res <- disc_LV_evol(N0=N[i-1,],alpha=alpha,E=E[i-1],x=x[i-1,],P=P,w=w,Wmax=Wmax,h2=h2)
  N[i,] <- res$Nt1
  r[i,] <- res$r
  # trait change
  d <- E[i-1] - x[i-1,]
  d1 <- k*d
  x[i,] <- E[i-1] - d1
  # environmental change
  E[i] <- E[i-1] + abs(rnorm(1, 0, 0.05))
}

### Plot of time series ###
gdat <- N
gdat$time <- 1:t
dat <- melt(gdat, id.vars="time")
pa2 <- ggplot2::ggplot(dat, aes(time, value, col=variable)) + geom_point()
### ###

### Mod2: HMSC model fit 2 ###
# prepare data in HMSC format
dat <- as.data.frame(cbind(log(N),x))
dat$time <- 1:t
df <- data.frame(dat[(2:t),-5])
colnames(df) <- c("Nt1","Nt2","xt1","xt2")
df$dx1 <- abs(dat$x1[2:t] - dat$x1[1:(t-1)])
df$dx2 <- abs(dat$x2[2:t] - dat$x2[1:(t-1)])
Y <- as.matrix(cbind(df$Nt1, df$Nt2))
XData <- data.frame(cbind(E[1:(t-1)],E[1:(t-1)]^2,abs(dat$x1[2:t] - dat$x1[1:(t-1)]),abs(dat$x2[2:t] - dat$x2[1:(t-1)])))
colnames(XData) <- c("E","Esq","dx1","dx2")
studyDesign = data.frame(sample = as.factor(1:(t-1)))
rL = HmscRandomLevel(units = studyDesign$sample)
m.2.hmsc = Hmsc(Y = Y, XData = XData, XFormula = ~E+Esq+dx1+dx2, studyDesign = studyDesign, ranLevels = list(sample = rL))
# Bayesian model parameters
nChains <- 2
thin <- 5
samples <- 2000
transient <- 1000*thin
verbose <- 500*thin
# sample MCMC
m.2.sample <- sampleMcmc(m.2.hmsc,thin=thin,sample=samples,transient=transient,nChains=nChains,verbose=verbose)
# prepare for plot
pred2 <- predict(m.2.sample)
hold1 <- list()
hold2 <- list()
for(i in 1:length(r)){
  hold1[[i]] <- as.numeric(Y[,i])
  data <- lapply(pred2, function(x) x[, colnames(x)[i], drop = FALSE])
  hold2[[i]] <- matrix(unlist(data), nrow=length(pred2), byrow=TRUE)
}
y <- unlist(hold1) # vector of each species' observed N across time points
yrep <- as.matrix(data.frame(hold2))
### Mod2 plot ###
par(mgp=c(2,0.45,0), tcl=-0.4, mar=c(1.3,1.2,0,0))
color_scheme_set("brightblue")
#ppc_dens_overlay(y, yrep[1:50, ])
p6 <- ppc_intervals(y,yrep,x=rep(dat$time[1:(t-1)],times=2),prob = 0.95) +
  labs(x = "time",y = "ln(N)",) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none")
### Gradient plot ###
Gradient <- constructGradient(m.2.sample,focalVariable="E",non.focalVariables=list(Esq=list(2)),ngrid=39)
Gradient$XDataNew$Esq <- Gradient$XDataNew$E^2
predY <- predict(m.2.sample,XData=Gradient$XDataNew,expected=TRUE)
plotGradient(m.2.sample,Gradient,pred=predY,showData=T,measure="Y",index=1,main="",xlab="E_t",ylab="predicted N1_t+1",showPosteriorSupport=FALSE,cex.axis=0.75,mar=c(1,1,0,0),bty='n')
p7 <- recordPlot()
plotGradient(m.2.sample,Gradient,pred=predY,showData=T,measure="Y",index=2,main="",xlab="E_t",ylab="predicted N2_t+1",showPosteriorSupport=FALSE,cex.axis=0.75,mar=c(1,1,0,0),bty='n')
p8 <- recordPlot()
### Beta posterior plot ###
m.post = Hmsc::convertToCodaObject(m.2.sample)
m.beta <- as.data.frame(rbind(m.post$Beta[[1]],m.post$Beta[[2]]))
colnames(m.beta) <- c("Int1","E1","Esq1","dx1_1","dx1_2","Int2","E2","Esq2","dx2_1","dx2_2")
p9 <- bayesplot::mcmc_areas(m.beta) + scale_x_continuous(limits=c(-20,30)) + theme(axis.text=element_text(size=8))
### variance partitioning ###
par(mgp=c(2,0.45,0), tcl=-0.4, mar=c(1.3,1.2,0.5,0.5))
VP <- computeVariancePartitioning(m.2.sample,group=c(1,1,1,2,3),groupnames=c("Env","Sp1","Sp2"))
plotVariancePartitioning(m.2.sample,VP,cols=c("white","skyblue","darkgrey","orange"),args.legend=list(cex=0.6,bg="transparent"),cex.axis=0.6,main="")
p10 <- recordPlot()
### ### ### ### ### ### ###
### all plots together ###
### ### ### ### ### ### ###
b <- plot_grid(p6,plot_grid(p7,p8,nrow=2),p9,p10,nrow=1,rel_widths=c(2,1,1,1))

# 3. Sim: Two species, logistic growth + competition, environmental covariate in spatially structured environment ---------
### Sim3. Two species, Leslie-Gower competition, environmental covariate, evolution, multiple sites
set.seed(42)
# Initial conditions
# Simulation of model for t time steps, i sites
j <- 10
# random site locations
xycoords = matrix(runif(2*j), ncol = 2)
rownames(xycoords) <- 1:j
t <- 40
N <- array(NA,dim=c(j*t,length(r)))
N <- as.data.frame(N)
colnames(N) <- paste0("N",1:length(r))
x <- array(NA,dim=c(j*t,length(r)))
x <- as.data.frame(x)
colnames(x) <- paste0("x",1:length(r))
r <- array(NA,dim=c(j*t,length(r)))
r <- as.data.frame(r)
colnames(r) <- paste0("r",1:length(r))
N[1:j,] <- rep(N0,each=j)
E <- rep(NA, t)
E[1] <- E.0
x[1:j,] <- rep(x.0,each=j)
What <- Wmax*sqrt(w/(P+w))
r0 <- as.numeric(What*exp((-(((w+(1-h2)*P)/(P+w))*(E[1]-x[1,]))^2)/(2*(P+w))))
r[1:j,] <- rep(r0,each=j)
# Keep track of site and time
study <- array(NA,dim=c(j*t,3))
study <- as.data.frame(study)
colnames(study) <- c("obs","site","time")
study$obs <- 1:dim(study)[1]
study$site <- rep(1:j,times=t)
study$time <- rep(1:t,each=j)

for (i in 2:t) {
  for(z in 1:j){
    res <- disc_LV_evol(N0=N[study$site==z & study$time==i-1,],alpha=alpha,E=E[i-1],x=x[study$site==z & study$time==i-1,],P=P,w=w,Wmax=Wmax,h2=h2)
    N[study$site==z & study$time==i,] <- res$Nt1
    r[study$site==z & study$time==i,] <- res$r
    # trait change
    d <- E[i-1] - x[study$site==z & study$time==i-1,]
    d1 <- k*d
    x[study$site==z & study$time==i,] <- E[i-1] - d1
  }
  # environmental change
  E[i] <- E[i-1] + abs(rnorm(1, 0, 0.05))
}
# site-level random effect
sigma <- 0
sigma.spatial <- 4
alpha.spatial <- 0.5
Sigma = sigma.spatial^2*exp(-as.matrix(dist(xycoords))/alpha.spatial)
# draw from covariance matrix
a = mvrnorm(mu=rep(0,j), Sigma = Sigma)
for(i in 1:t){
  N$N1[study$time==i] <- N$N1[study$time==i] + a
}
a = mvrnorm(mu=rep(0,j), Sigma = Sigma)
for(i in 1:t){
  N$N2[study$time==i] <- N$N2[study$time==i] + a
}

### Plot of time series ###
gdat <- N
gdat$site <- study$site
gdat$time <- study$time
dat <- melt(gdat, id.vars=c("site","time"))
pa3 <- ggplot2::ggplot(dat, aes(time, value, col=variable)) + geom_point() + facet_wrap(~site,nrow=2)
### ###

### Mod3: HMSC model fit 3 ###
# prepare data in HMSC format
#N[N <0] <- 0.01
# prepare data in HMSC format
dat <- as.data.frame(cbind(log(N),x))
dat$time <- study$time
df <- data.frame(cbind(dat$N1[dat$time !=1], dat$N2[dat$time !=1], dat$x1[dat$time !=1], dat$x2[dat$time !=1]))
colnames(df) <- c("Nt1", "Nt2", "xt1", "xt2")
Y <- as.matrix(cbind(df$Nt1, df$Nt2))
XData <- data.frame(E=rep(E[1:(t-1)],each=j))
XData$Esq <- XData$E^2
XData$dx1 <- abs(dat$x1[dat$time != 1] - dat$x1[dat$time != t])
XData$dx2 <- abs(dat$x2[dat$time != 1] - dat$x2[dat$time != t])
# Update study design to include spatial random effects
samp1 <- 1:(j*t)
samp1 <- as.data.frame(samp1)
samp1 <- subset(samp1,samp1>10)
rownames(XData) <- row.names(samp1)
studyDesign = data.frame(sample=1:(j*(t-1)),site=study$site[study$time > 1])
studyDesign$sample <- as.factor(studyDesign$sample)
studyDesign$site <- as.factor(studyDesign$site)
rL.spatial = HmscRandomLevel(sData = xycoords)  
rL.sample = HmscRandomLevel(units = studyDesign$sample)
# Fit HMSC model
m.3.hmsc = Hmsc(Y = Y, XData = XData, XFormula = ~E+Esq+dx1+dx2, studyDesign = studyDesign, ranLevels = list("sample"=rL.sample,"site"=rL.spatial))
# Bayesian model parameters
nChains <- 2
thin <- 5
samples <- 2000
transient <- 1000*thin
verbose <- 500*thin
# sample MCMC
m.3.sample <- sampleMcmc(m.3.hmsc,thin=thin,sample=samples,transient=transient,nChains=nChains,verbose=verbose)
### Gradient plot ###
par(mgp=c(2,0.45,0), tcl=-0.4, mar=c(1.3,1.2,0,0))
Gradient <- constructGradient(m.3.sample,focalVariable="dx2",non.focalVariables=list(E=list(2),Esq=list(2),dx1=list(2)),ngrid=39)
Gradient$XDataNew$Esq <- Gradient$XDataNew$E^2
predY <- predict(m.3.sample,XData=Gradient$XDataNew,expected=TRUE)
plotGradient(m.3.sample,Gradient,pred=predY,showData=T,measure="Y",index=1,main="",xlab="E_t",ylab="predicted N1_t+1",showPosteriorSupport=FALSE,cex.axis=0.75,mar=c(0,0,0,0),bty='n',xaxs = "r", yaxs = "r", tcl=-0.4)
p11 <- recordPlot()
plotGradient(m.3.sample,Gradient,pred=predY,showData=T,measure="Y",index=2,main="",xlab="E_t",ylab="predicted N2_t+1",showPosteriorSupport=FALSE,cex.axis=0.75,mar=c(0,0,0,0),bty='n',xaxs = "r", yaxs = "r", tcl=-0.4)
p12 <- recordPlot()
### Beta posterior plot ###
m3.post.hmsc <- convertToCodaObject(m.3.sample)
m.beta <- as.data.frame(rbind(m3.post.hmsc$Beta[[1]],m3.post.hmsc$Beta[[2]]))
colnames(m.beta) <- c("Int1","E1","Esq1","dx1_1","dx1_2","Int2","E2","Esq2","dx2_1","dx2_2")
p13 <- bayesplot::mcmc_trace(m3.post.hmsc$Beta)
p14 <- bayesplot::mcmc_areas(m.beta) + scale_x_continuous(limits=c(-20,30)) + theme(axis.text=element_text(size=8))
par(mgp=c(2,0.45,0), tcl=-0.4, mar=c(1.3,1.2,0.5,0.5))
VP <- computeVariancePartitioning(m.3.sample,group=c(1,1,1,2,3),groupnames=c("Env","Sp1","Sp2"))
plotVariancePartitioning(m.3.sample,VP,cols=c("white","skyblue","darkgrey","orange","pink"),args.legend=list(cex=0.6,bg="transparent"),main="",cex.axis=0.6)
p15 <- recordPlot()

# prepare for plot
pred3 <- predict(m.3.sample)
hold1 <- list()
hold2 <- list()
for(i in 1:length(r)){
  hold1[[i]] <- as.numeric(Y[,i])
  data <- lapply(pred3, function(x) x[, colnames(x)[i], drop = FALSE])
  hold2[[i]] <- matrix(unlist(data), nrow=length(pred3), byrow=TRUE)
}
y <- unlist(hold1) # vector of each species' observed N across time points
yrep <- as.matrix(data.frame(hold2))
### Mod3 plot ###
gg_dat <- data.frame(y=c(Y[,1],Y[,2]))
gg_dat$species <- rep(c(1,2),each=390)
gg_dat$site <- rep(study$site[study$time !=1],times=2)
gg_dat$time <- rep(study$time[study$time !=1],times=2)
gg_dat$ci_025 <- apply(yrep,2,quantile,probs=.025)
gg_dat$ci_975 <- apply(yrep,2,quantile,probs=.975)
p16 <- ggplot2::ggplot(gg_dat, aes(time, y, col=species)) + geom_point(size = 0.5) + facet_wrap(~site,nrow=2) + geom_errorbar(aes(ymin=ci_025, ymax=ci_975), alpha=0.5) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  theme(legend.position = "none",strip.text = element_text(size=6))

### ### ### ### ### ### ###
### all plots together ###
### ### ### ### ### ### ###
c <- plot_grid(p16,plot_grid(p11,p12,nrow=2),p14,p15,nrow=1,rel_widths=c(3,1,1,1),rel_heights=c(2,1,1,1))


### ### ### ### ### ### ###
### all plots together ###
### ### ### ### ### ### ###
d <- plot_grid(p1,plot_grid(p2,p3,nrow=2),p4,p5,p6,plot_grid(p7,p8,nrow=2),p9,p10,p16,plot_grid(p11,p12,nrow=2),p14,p15,nrow=3,ncol=4,rel_widths=c(3,1,1,1,3,1,1,1,3,1,1,1))

pdf("./output/fig1.pdf",width=8,height=4)
d
dev.off()


