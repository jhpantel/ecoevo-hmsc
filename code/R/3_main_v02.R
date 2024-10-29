# Title. main_v02.R
# Author. J.H. Pantel
# Date 28-10-2024
# Description. This script contains analyses for the HMSC_ecoevo project and associated manuscript Hermann & Pantel 202x.

#### Libraries ---------------------------------------------------------------
library(Hmsc)
library(tictoc)

#### Step 1. Run MATLAB simulations --------------

### Scenario 1. ###
# Weak selection, fixed environment, all dispersal and heritability levels
# P = 1, w = 10, v[E] = 0
# /code/R/metacom_evol_sim_same_init_v01.m

### Scenario 2. ###
# Strong selection, fixed environment, d = 10^-3, h2 = 0.1
# P = 0.1, w = 0.2, v[E] = 0
# HMSC_ecoevo/code/simulation/metacom_evol_sim_same_init_v02.m

### Scenario 3. ###
# Weak selection, variable environment, d = 10^-3, h2 = 0.1
# P = 1, w = 10, v[E] = 0.01
# HMSC_ecoevo/code/simulation/metacom_evol_sim_same_init_v03.m

### Scenario 4. ###
# Strong selection, variable environment, d = 10^-3, h2 = 0.1
# P = 0.1, w = 0.2, v[E] = 0.01
# HMSC_ecoevo/code/simulation/metacom_evol_sim_same_init_v04.m

#### Step 2. Prepare for HMSC analyses --------------
d_lev <- c("d_zero","d_minus9","d_minus8","d_minus7","d_minus6","d_minus5","d_minus4","d_minus3","d_minus2","d_01","d_02","d_03","d_04","d_05","d_06","d_07","d_08","d_09","d_10")
h_lev <- c("h_0_","h_01_","h_02_","h_03_","h_04_","h_05_","h_06_","h_07_","h_08_","h_09_","h_10_")
#setwd("./data")

# HMSC for single results condition
z <- 2
f <- 4

print(paste(h_lev[z],d_lev[f],sep=""))
result <- paste(h_lev[z],d_lev[f],sep="")
load(paste("/Users/jhpantel/Nextcloud/Pantel/data/",result,"_res_v02.RData",sep="")) # Scenario 1
# res  <- R.matlab::readMat(paste(result,"_res_v02.mat",sep="")) # Scenario 2
# res  <- R.matlab::readMat(paste(result,"_res_v03.mat",sep="")) # Scenario 3
# res  <- R.matlab::readMat(paste(result,"_res_v04.mat",sep="")) # Scenario 4

# Address extinct species
a <- apply(N,3,function(x) colSums(x) > 0)
b <- apply(a,1,function(x) min(which(x==FALSE)))
c <- which(is.finite(b)) # identity of species that go extinct
d <- b[c] # generation of last occurrence for species that go extinct

# Arrange all data for HMSC
## site-by-species abundance data as a time series
Y <- array(NA,dim=c(50000,15),dimnames=list(NULL,c("y1","y2","y3","y4","y5","y6","y7","y8","y9","y10","y11","y12","y13","y14","y15")))
count <- seq(1,50000,by=50)
for(i in 1:1000){
  Y[(count[i]:(count[i]+49)),] <- N[,,i]
}

## site-by-species trait data as a time series
trait_x <- xt
trait_init <- trait_x[,,1]  # Eliminate my initial trait value of 9999 for absent species
trait_init[trait_init == 9999] <- NaN
trait_x[,,1] <- trait_init
## site-by-species CHANGE in trait value as a time series
## absolute value of change in trait value
delta_x <- array(NA,c(50,15,999),dimnames=list(NULL,c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15"),NULL))
for(i in 2:time){
  delta_x[,,i-1] <- abs(trait_x[,,i] - trait_x[,,i-1])
}
delta_x_time <- array(NA,dim=c(49950,15),dimnames=list(NULL,c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15")))
count <- seq(1,49950,by=50)
for(i in 1:999){
  delta_x_time[(count[i]:(count[i]+49)),] <- delta_x[,,i]
}
delta_x_time[is.finite(delta_x_time) == FALSE] <- 0

## Random effects of site and year
Random <- array(NA,dim=c(50000,2),dimnames=list(NULL,c("site","year")))
Random[,1] <- rep(1:50,1000)
Random[,2] <- rep(1:1000,each=50)
Random <- as.data.frame(Random)

## Fixed effects: environment, site abundance the year before, species trait value the year before
X <- Y[Random$year != 1000,]  # Use abundance in years 1-999 as fixed effects - the abundance the year before each response variable (which will be for Years 2-1000).
Y <- Y[Random$year != 1,]  # Trim Y to only include years 2-1000

## Scaling will be done in the HMSC command, the argument YScale=TRUE

Random <- Random[Random$year != 1,]  # Trim Random to only include years 2-1000
Random$site <- as.factor(Random$site)
Random$year <- as.factor(Random$year)
Random$sample <- as.factor(1:nrow(Random))

## Random effects
rownames(xy) <- 1:patch
rL.spatial <- Hmsc::HmscRandomLevel(sData = xy)
rL.spatial = Hmsc::setPriors(rL.spatial,nfMin=1,nfMax=5)

rL.time <- Hmsc::HmscRandomLevel(units = unique(Random$year))
rL.time = Hmsc::setPriors(rL.time,nfMin=1,nfMax=5)

## Fixed effects: environment
env_hold <- array(NA,dim=c(dim(Y)[1],dim(E)[2]))
for(i in 1:dim(E)[2]){ # for each environmental variable
  env_hold[,i] <- rep(E[,i],999)
  colnames(env_hold)[i] <- paste("env",i,sep="")
}

## Interaction effects: site abundance the year before X change in species trait value
X_delta_x <- X * delta_x_time
colnames(X_delta_x) <- paste("x",colnames(X_delta_x),sep="")

if(z==1){
  X <- env_hold
} else{
  X <- cbind(env_hold,delta_x_time)
  X <- as.data.frame(X)
}
X <- as.data.frame(X)

# Remove rare / extinct species
if (length(c) > 0) { # if there are any extinct species
  e <- c[d <= 50] # for species that go extinct by generation 50
  # Remove rare species
  Y <- Y[,-e]
  X <- X[]
  # colnames(X)[c(e+1)]
  X <- X[,-c(e+1)]
}

#### Step 4. Run HMSC analyses --------------
#XFormula = ~y1+y2+y3+y4+y5+y6+y7+y8+y9+y10+y11+y12+y13+y14+y15+env1+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+xy1+xy2+xy3+xy4+xy5+xy6+xy7+xy8+xy9+xy10+xy11+xy12+xy13+xy14+xy14
sp <- s[-e]
XFormula = ~env1+x1+x3+x4+x6+x7+x8+x9+x10+x12+x13+x14+x15
m <- Hmsc::Hmsc(Y=Y,XData=X,XFormula = XFormula,studyDesign=Random,ranLevels=list("site"=rL.spatial,"year"=rL.time),distr="normal",XScale=TRUE,YScale=TRUE)
tic()

nChains <- 4
thin <- 100
samples <- 250
transient <- 250*thin
verbose <- 50*thin

m.1 <- Hmsc::sampleMcmc(m,samples=samples,transient=transient,nChains=nChains,nParallel=nChains,verbose=verbose,thin=thin,updater = list(GammaEta = FALSE))
a <- toc()
save.image(file=file.path(paste("/Users/jhpantel/Nextcloud/Pantel/data/",h_lev[z],d_lev[f],"_hmsc_v02_short.RData",sep=""))) # Scenario 1
# save.image(file=file.path(paste(h_lev[z],d_lev[f],"_hmsc_v02.RData",sep=""))) # Scenario 2
# save.image(file=file.path(paste(h_lev[z],d_lev[f],"_hmsc_v03.RData",sep=""))) # Scenario 3
# save.image(file=file.path(paste(h_lev[z],d_lev[f],"_hmsc_v04.RData",sep=""))) # Scenario 4
tic()

#### Step 5. Evaluate HMSC diagnostics --------------
# load("../output/h_01_d_minus3_hmsc_v01.RData")
# /code/R/diagnostics.R

#### Step 6. Evaluate HMSC results --------------
# /code/R/results.R

