## Evolving metacommunity simulation
## J.H. Pantel & Ruben Hermann
## June 26 2024
## Scenario 1. Weak selection, fixed environment, all dispersal and heritability levels
## # P = 1, w = 10, v[E] = 0
  
#set.seed(1)

###################################
## Step 0. Small functions    ##
###################################
# I needed to write some small function, which lets me run the dispersal lines smoother
# Here I check if the drawn runif number is smaller then the likelihood of dispersal
# if so then one individual with be migrating
# x = rlist (list of randomly generated numbers (runif) for each individual in a patch of species i)
# i = species number (is part of the loop for dispersal)
check <- function(x){
    x <=d[i]
}

# This function is the randome sample to chose amone sites where the species can immigrate,
# weighted by connectivitx.
# It takes simply a vector of all patches in, so:
# x = site/1:50 and i = species number (is part of the loop for dispersal)
disp_patch <- function(x){
    sample(x=site[-x],size=mig_count[x],replace=T,prob=connec[i,-x,x])
}
#Function to update the trait value, where the environmental value is substracted
#by the distance of the trait (which is changed by evo)
#x = The number of the species
#xt_up = the updated trait value, which is saved into xt
up_trait <- function(x){
  xt_up<-E-dt[,x,t]
}
###################################
## Step 1. Initialize species    ##
###################################
patch = 30   # Number of patches in the landscape
spec = 10    # Number of species in the community
time = 1000  # Number of time steps to run the simulation
N = array(numeric(),c(patch,spec,time))   # Record of abundances of all species over time patch x spec x time
P = array(1,c(1,spec))    # Phenotypic variances for all species
w = array(10,c(1,spec))   # Width of fitness curve
Wmax = 2    # Maximum population growth rate for all species
Wt = array(numeric(),c(patch,spec,time))   # Record of fitness vlaues over time, patch x spec x time
dt = array(numeric(),c(patch,spec,time))   # Record of distance from optimum value over time, patch x spec x time
xt = array(9999,c(patch,spec,time))     # Record of phenotype value over time, patch x spec x time

## Heritability
# Scenario 1: No evolution
#h2 = array(0,c(1,spec))

# Scenario 2: Evolution, same heritatibility, high
#h2 = array(0.9,c(1,spec))

extprob = 0.001    # Probability that species goes extinct when under critical density
repnum = 1000     # Totan number of replicates of simulations

What = Wmax*sqrt(w/(P+w))[1]
K = array(1000,c(patch,spec)) # Carrying capacity of each species in each patch

###################################
## Step 2. Initialize patches    ##
###################################

# Spatial connectivity and dispersal
##Random spatial location
xy = array(runif(patch*2),c(patch,2))
## Unifrom spatial locations
#x.1 = seq(0.1,1,0.1)
#x = array(x.1,c(length(x.1)*5,1))
#y.1 = seq(0.1,1,0.2)
#y = rep(y.1,each=10)
#xy = array(c(x,y),dim=c(50,2))
distance = as.matrix(dist(xy))   #Spatial distance
dkern = array(0.01,c(1,spec)) 

## Build connec, a 15x50x50 array with each species' connectivity matrix
# Not translating this part as the vector 'd' is missing - RH

## Build connec, a 15x50x50 array with each species' connectivity matrix
connec = array(0,c(spec,patch,patch))  # Species x patch x patch matrix of connectivtiy values

for (i in 1:spec) {
    connec[i,,] = exp((-1/dkern[i])*distance)
}

## Environmental properties
E = array(runif(patch),c(patch,1))    # Each patch has an environmental optimum value E. d, distance between species phenotype and patch optimum, will equal some phenotype value of the species x minus this E quantity. So d = x - E

## Initial occupancy
# Random initial occupancy
occ = array(rpois(50,0.75),c(1,50))
n = list()
for (i in 1:patch) {n[[i]] = rpois(occ[i],10)}
init = list()
for (i in 1:patch) {init[[i]] = sample(1:spec,occ[i])}

N[,,1]=0 
start = array(numeric())
for (i in 1:patch){
    if(length(n[[i]]!=0)){  # If there are species in this location
        N[i,init[[i]],1] = n[[i]] # Seed initial species occupancy at half the carrying capacity
        start.1 = matrix(c(rep(i,length(init[[i]])),init[[i]]),c(length(init[[i]]),2)) #this is from line 92 of the matlab script, saving initial species and the patch in which they start
        start = rbind(start,start.1)
    }
}
start = start[order(start[,2]),]
a = start[,1]
b = start[,2]
# Uniform initial occupancy, low N0
# N[,,1] = 10

## Degree of maladaption
s = 1:spec
s0 = s[(apply(N[,,1],2,sum))>0]   #  Return index of all species present in initial patches
B0 = array(numeric(),c(1,spec))
B0[s0] = array(rgamma(length(B0[s0]),shape=0.75,scale=1))  # Drawn from gamma distribution
d0.1 = sqrt(B0*(w+P))
x0.1 = min(E)-d0.1     # Depends on E value for site with lowest E #why do we create the distance matrix first and then overwrite it?

x0 = matrix(x0.1,patch,spec,byrow=T)   # Create a patch x spec matrix for starting distance to optimal trait
x0[N[,,1]==0]<-NA    # Set all values to NA where no species is present in initial patches
d0 = matrix(rep(d0.1,each=patch),c(patch,spec))
d0[N[,,1]==0]=NA
d0[N[,,1]>0] = E[a] - x0[N[,,1]>0]

## Initial distance from optimum phenotype and fitness
dt[,,1] = d0
xt[,,1] = x0

calc_wt <- function(x) (x^2/(2*(P+w)))
Wt[,,1] = What*exp(-t(apply(dt[,,1],1,calc_wt)))

#k = ((w + ((1-h2)*P))/(w+P));  #The quantity k determines the evolutionary inertia of the character, no evolution if k = 1, if k = w/(W + P), evolutionary equlibrium is approached at fastest possible speed

# Introduce intra- and inter-specific competition coefficients for
# Leslie-Gower growth model. Diagonal is alpha_ii, set to 0. Competition
# matrix is not symmetric, and has values drawn from a random uniform
# distribution from 0 - 0.15.
alpha =diag(0.00125,spec,spec)
alpha[upper.tri(alpha,diag=F)]=0 + (.0015-0)*runif(length(alpha[upper.tri(alpha,diag=F)]))
alpha[lower.tri(alpha,diag=F)]=0 + (.0015-0)*runif(length(alpha[lower.tri(alpha,diag=F)]))

d_vals = c(0, 10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
h2_vals = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
d_names = c("0", "minus9", "minus8", "minus7", "minus6", "minus5", "minus4", "minus3", "minus2", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10")
h2_names = c("0", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10")


h2 = array(h2_vals[2],c(1,spec)) #which h2_value do you want to use? - just using 2 to test
k = ((w+((1-h2)*P))/(w+P))
d = array(d_vals[8],c(1,spec))  #which d_value do you want to use? - just using 3 to test


save.image(file="C:/Users/Ruben Hermann/Nextcloud/Shared_with_me/ruben.hermann/Projects/project #1 HMSC/code/R/init_comm.Rdata")


load("C:/Users/Ruben Hermann/Nextcloud/Shared_with_me/ruben.hermann/Projects/project #1 HMSC/code/R/init_comm.Rdata")
###################################################
## Step 3. Evolving metacommunity simulations    ##
###################################################

for (t in 2:200) {
    t
    pop = N[,,t] # Write this generation's population size into a temporary record 
    pop[N[,,t-1]==0] = 0 # Carry over all extinct species

    ## Stochastic extinction due to critical population size
    # Any species at risk of stochastic extinction due to low population size are evaluated for potential extinction
    pop[N[,,t-1]>0 & N[,,t-1] < 100]  = runif(length(pop[N[,,t-1]>0 & N[,,t-1] < 100]))
    extinct<-which(pop>0&pop<extprob,arr.ind=T) # Record position of species that went stochastically extinct
    if (nrow(extinct)>0) {
        for (z in 1:nrow(extinct)){   #For each species that is stochastically extinct
            xt[extinct[z,1],extinct[z,2],t] =0
            dt[extinct[z,1],extinct[z,2],t] =0
        }
    }
    pop[pop>0&pop<extprob]=0  # Stochastically extinct populations set to zero for this generation    

    last=N[,,t-1]
    
    ## Carry over all populations that are above the critical population density
    pop[last>=100] = last[last>=100] #this had to be moved up, as R does not like the NAs, when accesing position information
    pop[pop>extprob & pop < 100] = last[N[,,t-1]>0 & N[,,t-1]<100 & pop>extprob] # Surviving populations

  

    ## Growth of all remaining populations
    dt[,,t] = sweep(dt[,,t-1],2,k,"*") # Update distance to optimum phenotype
    # Revise the value for the population that went stochastically extinct this generation
    if (nrow(extinct)>0) {
        for (X in 1:nrow(extinct)) {
            dt[extinct[X,1],extinct[X,2],t] =NA  #but this was already done ine line 154? - RH
        }
    }
    xt[,,t] =  sapply(1:spec,up_trait)  # update population average trait value
    Wt[,,t] = What*exp(-(sweep((dt[,,t]^2),2,(2*(P+w)),"/")))  # update fitness
    idx = pop==0

    num = Wt[,,t] * pop  # Growth formula numerator
    idx_XY = which(idx==0,arr.ind=T) #
    comp = rowSums(pop[idx_XY[,1],]*alpha[idx_XY[,2],])      # Get the pop change due to intra snd inter competition   
    denom = pop
    denom[idx==0] = 1 + comp
  

    pop[idx==0] = round(num[idx==0]/denom[idx==0])

    id_na = which(is.na(pop),arr.ind=T)
    if(nrow(id_na)>0){break}

    ## Competition - elimination of species not at carrying capacity
    # COMMENT ALL BELOW FOR LESLIE GOWER
    
    # I will skip this for now as it anyways is not run - RH

    ## COMMENT ALL ABOVE FOR LESLIE GOWER

    ## Dispersal
    # Compare each ind to a per-capita disp rate, generate a list of dispersers from each site, each species
    for (i in 1:spec) {   # for each species
    ## Generate number of dispersing individuals
        disp_count = pop[,i] #  Copy this species N record into a new variable
        rmat =array(0,dim=c(50,max(disp_count))) #Why is this line still there, rmat is not used? - RH
        rlist = sapply(disp_count[is.finite(disp_count)],runif) # Generate a random number 0-1 for each individual in the population at each site
        mig_list = lapply(rlist, check) # Checking in the list if the runif nunmber is smaller than the dispersal likelihood (check top of script for function)
        mig_list = lapply(mig_list,sum) # Summing the values in each element of the list
        mig_count = unlist(mig_list) # Turning the list into a vector of length 50 (patches) --> Number of migrants per patch

        # Subtract dispersers from source patch
        disp_count[is.finite(disp_count)] = disp_count[is.finite(disp_count)] - mig_count

        # Add dispersers to new patches
        id = which(is.finite(disp_count),arr.ind=T) # List of migrant source patches
        site = 1:patch # List of sites in region to sample from
        if (sum(connec[i,,]) == 0) {disp_site=NULL} # If there was no dispersal
        if (sum(connec[i,,]) == 50) {disp_site=NULL} # What is this line of code for? - RH
        if (sum(connec[i,,]) != 0 & sum(connec[i,,]) != 50){
            disp_site = sapply(site,disp_patch)
        }
        ## Count number of individuals going into new patches, add that to disp_count
        count = sapply(disp_site,table)    # Return number of individuals migrating into each new patch -

        ### Need to update xt and dt to reflect potential of same
        ### species colonizing from multiple other patches, species
        ### already present here (e.g. it needs to be weighted averages
        ### taking N into account).
        for (j in 1:length(count)){
            local_N = disp_count[unique(disp_site[[j]])]
            local_X = xt[unique(disp_site[[j]]),i,t]
            incom_N = unname(count[[j]])
            incom_X = xt[id[j],i,t]

            # Update record of N in each patch with migrants
            disp_count[unique(disp_site[[j]])] = local_N + incom_N

            # Update record of population mean trait x in each patch with migrants
            local_X[is.na(local_X)] = 1  #what is this suposed to do? - RH
            xt[unique(disp_site[[j]]),i,t] = ((local_X*local_N) + (incom_N*incom_X)) / (local_N+incom_N)

            # Update record of distance to optimum phenotype in new patch
            dt[unique(disp_site[[j]]),i,t] = E[unique(disp_site[[j]])] - xt[unique(disp_site[[j]]),i,t]
        
            # Update record of fitness in new patch
             Wt[unique(disp_site[[j]]),i,t] = What*exp(-(dt[unique(disp_site[[j]]),i,t]^2)/2*(P[i]+w[i]))
        }
        pop[,i] = disp_count # Write population of source reduced by emigrants and destination patch increased by immigrants
        # Check to see if this process gives NaN value
        Y = sum(is.na(pop))
        if(Y>0) break
    }
    N[,,t]= pop # Update population size this generation
    ## Check to see if population trait information remains for extinct populations
    A = which(pop==0 & N[,,t-1]>0,arr.ind=T)
    if(nrow(A)>0){
        for (z in 1:nrow(A)) { # For each species / site that went extinct for other reasons
            xt[A[z,1],A[z,2],t] = NA
            dt[A[z,1],A[z,2],t] = NA
            Wt[A[z,1],A[z,2],t] = NA
        }
    }
}

## use randsample to randomly sample among the available sites (the other 49) with weights according to connectivity
## subtract those individuals from source population, subject to extinction mortality probability, add remaining to new site
## Will have to figure out how to update dt and Wt at the new site (given the species has adapted to the old site in some way)
## update N matrix, end round

#### Step 4. Plot results
save.image("C:/Users/Ruben Hermann/Nextcloud/Shared_with_me/ruben.hermann/Projects/project #1 HMSC/code/R/h_02_d_08_res.RData")



