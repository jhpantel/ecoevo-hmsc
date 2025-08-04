# Title. 3_main.R
# Author. J.H. Pantel &  R.J. Hermann
# Date 29.07.2025
# Description. This script contains analyses for the HMSC_ecoevo project and associated manuscript Pantel & Hermann 2025.

#### Libraries ---------------------------------------------------------------
library(Hmsc)
library(tictoc)

#### Step 1. Run metacommunity simulations --------------

### Scenario ###
# Strong selection, fixed environment
# P = 0.1, w = 0.2, v[E] = 0
# ecoevo-hmsc/code/R/1_metacom_sim.R

#### Step 2. Prepare for HMSC analyses --------------
d_lev <- c("d_zero","d_minus9","d_minus8","d_minus7","d_minus6","d_minus5","d_minus4","d_minus3","d_minus2","d_01","d_02","d_03","d_04","d_05","d_06","d_07","d_08","d_09","d_10")
h_lev <- c("h_0_","h_01_","h_02_","h_03_","h_04_","h_05_","h_06_","h_07_","h_08_","h_09_","h_10_")

# HMSC for single results condition
z <- 10
f <- 1

print(paste(h_lev[z],d_lev[f],sep=""))
result <- paste(h_lev[z],d_lev[f],sep="")
load(paste("./data/mc/",result,"_res.RData",sep="")) # Scenario 1

# Address extinct species
a <- apply(N,3,function(x) colSums(x) > 0)
b <- apply(a,1,function(x) min(which(x==FALSE)))
c <- which(is.finite(b)) # identity of species that go extinct
d <- b[c] # generation of last occurrence for species that go extinct

# Arrange all data for HMSC
## site-by-species abundance data as a time series
spec_name = NULL
trait_name = NULL
for (i in 1:spec){
  spec_name[i] = paste("y",i,"",sep="")
  trait_name[i] = paste("x",i,"",sep="")
}

Y <- array(NA,dim=c((patch*cut),spec),dimnames=list(NULL,spec_name))
count <- seq(1,(patch*cut),by=patch)
for(i in 1:cut){
  Y[(count[i]:(count[i]+(patch-1))),] <- N[,,i]
}
Y[Y != 0] <- log(Y[Y != 0]) ## log(N)

## site-by-species trait data as a time series
trait_x <- xt
trait_init <- trait_x[,,1]  # Eliminate my initial trait value of 9999 for absent species
trait_init[trait_init == 9999] <- NaN
trait_x[,,1] <- trait_init
## site-by-species CHANGE in trait value as a time series
## absolute value of change in trait value
delta_x <- array(NA,c(patch,spec,(cut-1)),dimnames=list(NULL,trait_name,NULL))
for(i in 2:cut){
  delta_x[,,i-1] <- abs(trait_x[,,i] - trait_x[,,i-1])
}
delta_x_time <- array(NA,dim=c((patch*cut-patch),spec),dimnames=list(NULL,trait_name))
count <- seq(1,(patch*cut-patch),by=patch)
for(i in 1:(cut-1)){
  delta_x_time[(count[i]:(count[i]+(patch-1))),] <- delta_x[,,i]
}
delta_x_time[is.finite(delta_x_time) == FALSE] <- 0

## Random effects of site and year
Random <- array(NA,dim=c((patch*cut),2),dimnames=list(NULL,c("site","year")))
Random[,1] <- rep(1:patch,cut)
Random[,2] <- rep(1:cut,each=patch)
Random <- as.data.frame(Random)

## Fixed effects: environment, site abundance the year before, species trait value the year before
X <- Y[Random$year != cut,]  # Use abundance in years 1-999 as fixed effects - the abundance the year before each response variable (which will be for Years 2-1000).
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
env_names = NULL
for (i in 1:dim(E)[2]){
  env_names[i] = paste("env",i,"",sep="")
}
env_hold <- array(NA,dim=c(dim(Y)[1],dim(E)[2]),dimnames = list(NULL,env_names))
count <- seq(1,(patch*cut-patch),by=patch)
for(i in 1:(cut-1)){ # for each environmental variable
  env_hold[(count[i]:(count[i]+(patch-1))),] <- Et[,i+1]
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
#XFormula = ~env1+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15
sp <- s[-e]
sp
XFormula = ~env1+x1+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15
if (f == 1){
  m <- Hmsc::Hmsc(Y=Y,XData=X,XFormula = XFormula,studyDesign=Random,ranLevels=list("year"=rL.time),distr="normal",XScale=TRUE,YScale=TRUE)
  
} else {
  m <- Hmsc::Hmsc(Y=Y,XData=X,XFormula = XFormula,studyDesign=Random,ranLevels=list("site"=rL.spatial,"year"=rL.time),distr="normal",XScale=TRUE,YScale=TRUE)
}
tic()

nChains <- 4
thin <- 100
samples <- 250
transient <- 250*thin
verbose <- 50*thin

m.1 <- Hmsc::sampleMcmc(m,samples=samples,transient=transient,nChains=nChains,nParallel=nChains,verbose=verbose,thin=thin,updater = list(GammaEta = FALSE))
a <- toc()
save.image(file=file.path(paste("./data/mc/",h_lev[z],d_lev[f],"_hmsc.RData",sep=""))) # 
tic()
