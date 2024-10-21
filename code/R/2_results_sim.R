# Title. sim_results.R
# Author. J.H. Pantel &  R.J. Hermann
# Date 02-07-2024
# Description. This script contains different figures to check the simulation results for any interesting outcomes.

#### Libraries ---------------------------------------------------------------
#library(R.matlab)
library(vegan)
library(tidyverse)
library(patchwork)

##
####  Step 1. Read in the data and get it into shape --------------
load("././data/h_01_d_minus3_res.RData") # Scenario 1
##
################### Extinction checks from JP ###################
# Address extinct species
a <- apply(N,3,function(x) colSums(x) > 0)
b <- apply(a,1,function(x) min(which(x==FALSE)))
c <- which(is.finite(b)) # identity of species that go extinct
d <- b[c] # generation of last occurrence for species that go extinct

##
################### Data manipulation for patches plots ###################
## absolute value of change in trait value
delta_x <- array(NA,c(dim(xt)[1],dim(xt)[2],dim(xt)[3]),dimnames=list(NULL,c(paste("x",1:dim(xt)[2],"",sep="")),NULL))
delta_x[,,1]=0
for(i in 3:dim(xt)[3]){
  delta_x[,,i-1] <- abs(xt[,,i] - xt[,,i-1])
}

## fitness of species: -abs(abs(E)-x) --> so that negative values represent low fitness
fit_x <-  array(NA,c(dim(xt)[1],dim(xt)[2],dim(xt)[3]),dimnames=list(NULL,c(paste("y",1:dim(xt)[2],"",sep="")),NULL))
for(i in 1:dim(xt)[3]){
  fit_x[,,i] <- -abs(sweep(abs(xt[,,i]),1,E,"-"))
}

## average fitness of patch: -abs(abs(E)-x) --> so that negative values represent low fitness
fit_x_m_p <-  array(NA,c(dim(xt)[1],1,dim(xt)[3]),dimnames=list(c(paste("p",1:dim(xt)[1],"",sep="")),NULL,NULL))
for(i in 1:dim(xt)[3]){
  fit_x_m_p[,,i] <- rowMeans(fit_x[,,i],na.rm=T)
}

## species richness per patch
sp_rich <- array(NA,c(dim(N)[1],1,dim(N)[3]),dimnames=list(c(paste("p",1:dim(xt)[1],"",sep="")),NULL,NULL))
for(i in 1:dim(xt)[3]){
  sp_rich[,,i] <- rowSums(N[,,i]>0)
}

## total population size per patch
patch_dens <- array(NA,c(dim(N)[1],1,dim(N)[3]),dimnames=list(c(paste("p",1:dim(xt)[1],"",sep="")),NULL,NULL))
for(i in 1:dim(xt)[3]){
  patch_dens[,,i] <- rowSums(N[,,i])
}

## alpha diversity size per patch
#the vegan package can calculate three different alpha diversities,
#depending on which one you want copy paste it to index:
#"shannon", "simpson" and "invsimpson"
alpha_div <- array(NA,c(dim(N)[1],1,dim(N)[3]),dimnames=list(c(paste("p",1:dim(xt)[1],"",sep="")),NULL,NULL))
for(i in 1:dim(xt)[3]){
  alpha_div[,,i] <- diversity(N[,,i],index="invsimpson")
}

################### Data manipulation for species plots ###################
## amount of patches occupied by species
occ_pa <- array(NA,c(1,dim(N)[2],dim(N)[3]),dimnames=list(NULL,c(paste("y",1:dim(N)[2],"",sep="")),NULL))
for(i in 1:dim(N)[3]){
  occ_pa[,,i] <- colSums(N[,,i]>0)
}

## average fitness of a species over the whole metacommunity
fit_x_m <- array(NA,c(1,dim(N)[2],dim(N)[3]),dimnames=list(NULL,c(paste("y",1:dim(N)[2],"",sep="")),NULL))
for(i in 1:dim(N)[3]){
  fit_x_m[,,i] <- colMeans(fit_x[,,i],na.rm=T)
}

##
################### Data manipulation for global plots ###################
## Total species richness
total_rich <- array(NA,dim(N)[3],dimnames=list(NULL))
for(i in 1:dim(N)[3]){
  total_rich[i] <- sum(colSums(N[,,i]>0)>0)
}

## Total species density
total_den <- array(NA,dim(N)[3],dimnames=list(NULL))
for(i in 1:dim(N)[3]){
  total_den[i] <- sum(colSums(N[,,i]))
}

## Beta diversity
beta_div <- array(NA,dim(N)[3],dimnames=list(NULL))
for(i in 1:dim(N)[3]){
  beta_div[i] <- mean(vegdist(N[,,i],binary=T),na.rm=T)
}

## Gamma diversity
gamma_div <- array(NA,dim(N)[3],dimnames=list(NULL))
for(i in 1:dim(N)[3]){
  gamma_div[i] <- diversity(colSums(N[,,i]),"shannon")
}

##
#### Step 2. Plots of the 50 patches over time --------------
mycols <- c("maroon4","coral3","peachpuff3","gray47","hotpink2","darkorange","cyan3","khaki1","firebrick1","cadetblue4","mediumslateblue","mediumspringgreen","darkorchid4","midnightblue","rosybrown1")
mycols2<-mycols[1:ncol(N)]

#The number of plots per page needs to be adjusted manually, the rest should work automatically independed of the number of patches, species (up to 15) or running time
## Plot of metacommunity population dynamics over time
par(mfrow=(c(10,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
#par(mai=c(0.1,0.1,0.1,0.1))
for(j in 1:nrow(N)){
  plot(N[j,1,],col=mycols2[1],ylim=c(0,max(N,na.rm=T)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11,16,21,26,31,36,41,46))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% c(46,47,48,49,50))
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:ncol(N)){
    points(N[j,i,],col=mycols2[i])
    text(x=(dim(N)[3]),y=(max(N,na.rm=T)/2),labels=j,pos=2,col="black")
  }
}
mtext("Population dynamics in each plot over time", side = 3, line = 0, outer = TRUE)

## Plot of metacommunity trait dynamics over time
par(mfrow=(c(10,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
for(j in 1:nrow(xt)){
  plot(xt[j,1,],col=mycols2[1],ylim=c(min(xt,na.rm=T),1),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11,16,21,26,31,36,41,46))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% c(46,47,48,49,50))
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:ncol(xt)){
    points(xt[j,i,],col=mycols2[i])
    text(x=(dim(xt)[3]),y=min(xt,na.rm=T)/2,labels=j,pos=2,col="black")
  }
}
mtext("Trait dynamics in each plot over time", side = 3, line = 0, outer = TRUE)

## Plot of metacommunity delta trait dynamics over time
par(mfrow=(c(10,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
for(j in 1:nrow(delta_x)){
  plot(delta_x[j,1,],col=mycols2[1],ylim=c(0,max(delta_x,na.rm=T)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")

  if (j %in% c(1,6,11,16,21,26,31,36,41,46))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% c(46,47,48,49,50))
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:ncol(delta_x)){
    points(delta_x[j,i,],col=mycols2[i])
    text(x=(dim(delta_x)[3]),y=(max(delta_x,na.rm=T)/2),labels=j,pos=2,col="black")
  }
}
mtext("Delta trait in each plot over time", side = 3, line = 0, outer = TRUE)


# Plot of species fitness for given patch over time
par(mfrow=(c(10,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
for(j in 1:nrow(fit_x)){
  plot(fit_x[j,1,],col=mycols2[1],ylim=c(-5.6,0),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11,16,21,26,31,36,41,46))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% c(46,47,48,49,50))
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:ncol(fit_x)){
    points(fit_x[j,i,],col=mycols2[i])
    text(x=(dim(fit_x)[3]),y=(min(fit_x,na.rm=T)/2),labels=j,pos=2,col="black")
  }
}
mtext("Trait distance to E value of each patch over time", side = 3, line = 0, outer = TRUE)

## Plot of species richness for given patch over time
par(mfrow=(c(10,5)))
par(mar=c(1,1,1,1))
for(j in 1:nrow(sp_rich)){
  plot(sp_rich[j,1,],ylim=c(0,max(sp_rich,na.rm=T)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11,16,21,26,31,36,41,46))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% c(46,47,48,49,50))
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:ncol(sp_rich)){
    points(sp_rich[j,i,])
    text(x=(dim(sp_rich)[3]),y=(max(sp_rich,na.rm=T)/2),labels=j,pos=2)
  }
}
mtext("Species richness in each patch over time", side = 3, line = 0, outer = TRUE)

## Plot of alpha diversity for given patch over time
par(mfrow=(c(10,5)))
par(mar=c(1,1,1,1))
for(j in 1:nrow(alpha_div)){
  plot(alpha_div[j,1,],ylim=c(0,2),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11,16,21,26,31,36,41,46))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% c(46,47,48,49,50))
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:ncol(alpha_div)){
    points(alpha_div[j,i,])
    text(x=(dim(alpha_div)[3]),y=(max(alpha_div,na.rm=T)/2),labels=j,pos=2)
  }
}
mtext("Shannon H Index in each patch over time", side = 3, line = 0, outer = TRUE)

#### Step 3. Plots of the 15 species over time --------------
#These are sadly without color, as it is probably impossible to create 50 different colors that can be differentiated efficiently
## Plot of metacommunity population dynamics over time
par(mfrow=(c(3,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
for(j in 1:ncol(N)){
  plot(N[1,j,],ylim=c(0,max(N,na.rm=T)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% 11:15)
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:nrow(N)){
    points(N[i,j,])
    text(x=(dim(N)[3]),y=(max(N,na.rm=T)/2),labels=j,pos=2)
  }
}
mtext("Population dynamics of each species in its given patch", side = 3, line = 0, outer = TRUE)

## Plot of metacommunity trait dynamics over time
par(mfrow=(c(3,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
for(j in 1:ncol(xt)){
  plot(xt[1,j,],ylim=c(min(xt,na.rm=T),1),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% 11:15)
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:nrow(xt)){
    points(xt[i,j,])
    text(x=(dim(xt)[3]),y=(min(xt,na.rm=T)/2),labels=j,pos=2)
  }
}
mtext("Trait dynamics of each species in its given patch", side = 3, line = 0, outer = TRUE)


## Plot of metacommunity delta trait dynamics over time
par(mfrow=(c(3,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
for(j in 1:ncol(delta_x)){
  plot(delta_x[1,j,],ylim=c(0,max(delta_x,na.rm=T)),xlim=c(0,dim(delta_x)[3]),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% 11:15)
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:nrow(delta_x)){
    points(delta_x[i,j,])
    text(x=(dim(delta_x)[3]),y=(max(delta_x,na.rm=T)/2),labels=j,pos=2)
  }
}
mtext("Delta trait of each species in its given patch", side = 3, line = 0, outer = TRUE)

## Plot of species fitness for given patch over time
par(mfrow=(c(3,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
for(j in 1:ncol(fit_x)){
  plot(fit_x[1,j,],ylim=c(min(fit_x,na.rm=T),0),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% 11:15)
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:nrow(fit_x)){
    points(fit_x[i,j,])
    text(x=(dim(fit_x)[3]),y=(min(fit_x,na.rm=T)/2),labels=j,pos=2)
  }
}
mtext("Species fitness of each species in its patch", side = 3, line = -0, outer = TRUE)

## Average species fitness for the whole metacommunity over time
par(mfrow=(c(3,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
for(j in 1:ncol(fit_x_m)){
  plot(fit_x_m[1,j,],ylim=c(min(fit_x_m,na.rm=T),0),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% 11:15)
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:nrow(fit_x_m)){
    points(fit_x_m[i,j,])
    text(x=(dim(fit_x_m)[3]),y=(min(fit_x_m,na.rm=T)/2),labels=j,pos=2)
  }
}
mtext("Average species fitness over the whole community", side = 3, line = 0, outer = TRUE)

## Plot of occupied patches patch over time
par(mfrow=(c(3,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
for(j in 1:ncol(occ_pa)){
  plot(occ_pa[1,j,],ylim=c(0,nrow(N)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11))
    axis(2,col="grey40",col.axis="grey20",at = c(0,15,30,45),labels=T)
  if (j %in% 11:15)
    axis(1,col="grey40",col.axis="grey20",at = c(0,250,500,750,1000),labels=T)
  for(i in 1:nrow(occ_pa)){
    points(occ_pa[i,j,])
    text(x=(dim(occ_pa)[3]),y=(max(occ_pa,na.rm=T)/2),labels=j,pos=2)
  }
}
mtext("Number of occupied patches of each species", side = 3, line = - 2, outer = TRUE)


###
#### Step 4. Different multi-plots NOT over time ####
mycols <- c("maroon4","coral3","peachpuff3","gray47","hotpink2","darkorange","cyan3","khaki1","firebrick1","cadetblue4","mediumslateblue","mediumspringgreen","darkorchid4","midnightblue","rosybrown1")
mycols2<-mycols[1:ncol(N)]

## Plot of metacommunity population dynamics over delta X
par(mfrow=(c(10,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
#par(mai=c(0.1,0.1,0.1,0.1))
for(j in 1:nrow(N)){
  plot(x=delta_x[j,1,],y=N[j,1,],col=mycols2[1],ylim=c(0,max(N,na.rm=T)),xlim=c(0,max(delta_x,na.rm=T)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11,16,21,26,31,36,41,46))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% c(46,47,48,49,50))
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:ncol(N)){
    points(x=delta_x[j,i,],y=N[j,i,],col=mycols2[i])
    text(x=max(delta_x,na.rm=T),y=(max(N,na.rm=T)/2),labels=j,pos=2,col="black")
  }
}
mtext("Population dynamics in each plot over delta X", side = 3, line = 0, outer = TRUE)


## Plot of metacommunity population dynamics over fitness 
par(mfrow=(c(10,5)))
par(mai=c(0.1,0.1,0.1,0.1),oma=c(2,2,2,1))
#par(mai=c(0.1,0.1,0.1,0.1))
for(j in 1:nrow(N)){
  plot(x=fit_x[j,1,],y=N[j,1,],col=mycols2[1],ylim=c(0,max(N,na.rm=T)),xlim=c(min(fit_x,na.rm=T),0),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(1,6,11,16,21,26,31,36,41,46))
    axis(2,col="grey40",col.axis="grey20",labels=T)
  if (j %in% c(46,47,48,49,50))
    axis(1,col="grey40",col.axis="grey20",labels=T)
  for(i in 1:ncol(N)){
    points(x=fit_x[j,i,],y=N[j,i,],col=mycols2[i])
    text(x=min(fit_x,na.rm=T),y=(max(N,na.rm=T)/2),labels=j,pos=4,col="black")
  }
}
mtext("Population dynamics in each plot over its fitness", side = 3, line = 0, outer = TRUE)



#### Step 5. Global Plots over time --------------
par(mfrow=(c(1,1)))
par(mar=c(5.1, 4.1, 4.1, 2.1),oma = c(0, 0, 0, 0))
#Total species richness of the metacommunity - to see if extinctions happened
plot(1:dim(N)[3],total_rich,main="Species richness of the metacommunity",ylab="Number of species",xlab="Timesteps")

#Total species density of the metacommunity
plot(1:dim(N)[3],total_den,main="Species density of the metacommunity",ylab="Number of species",xlab="Timesteps")

#Beta diversity of the whole metacommunity
plot(1:dim(N)[3],beta_div,main="Sorensen Beta diversity",ylab="Index of dissimilarity",xlab="Timesteps",ylim=c(0.4,1))

#Gamma diversity of the whole metacommunity
plot(1:dim(N)[3],gamma_div,main="Gamma diversity according to Shannon H index",ylab="Diversity",xlab="Timesteps")


#### Step 6. End time point plots ####
par(mfrow=(c(2,2)))
par(mar=c(5.1, 4.1, 4.1, 2.1),oma = c(0, 0, 0, 0))

#First page of plots in histograms
#Final species richness per plot
hist(sp_rich[,,dim(sp_rich)[3]],breaks=10,xlab="Species richness per plot",main="Species richness per plot at the end of the simulation")
#Final species density per plot
hist(patch_dens[,,dim(patch_dens)[3]],breaks=15,xlab="Population density per plot",main="Population density per plot at the end of the simulation")
#Final Shannon H  per plot
hist(alpha_div[,,dim(alpha_div)[3]],breaks=15,xlab="Shannon H per plot",main="Shannon H per plot at the end of the simulation")
#Final average fitness  per plot
hist(rowMeans(fit_x[,,dim(fit_x)[3]],na.rm=T),breaks=15,xlab="Average fitness per plot",main="Average fitness plot at the end of the simulation")

#Creating data frame for ggplot2
plot_dat <-data.frame(x=xy[,1],y=xy[,2],sp_rich=sp_rich[,,dim(sp_rich)[3]],patch_dens=patch_dens[,,dim(patch_dens)[3]],alpha_div=alpha_div[,,dim(alpha_div)[3]],fit_x=rowMeans(fit_x[,,dim(fit_x)[3]],na.rm=T))


#Second page of plots in xy-coordinates of plot

#Final species richness per plot
p1<-ggplot(plot_dat,aes(x=x,y=y,color=sp_rich))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_gradientn(colors=rainbow(5))+
  ylab("Y coordinate")+xlab("X Coordinate")+
  ggtitle("Final species richness per plot")
#Final species density per plot
p2<-ggplot(plot_dat,aes(x=x,y=y,color=alpha_div))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_gradientn(colors=rainbow(5))+
  ylab("Y coordinate")+xlab("X Coordinate")+
  ggtitle("Final Shannon H per plot")
#Final Shannon H  per plot
p3<-ggplot(plot_dat,aes(x=x,y=y,color=fit_x))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_gradientn(colors=rainbow(5))+
  ylab("Y coordinate")+xlab("X Coordinate")+
  ggtitle("Average final fitness per plot")
#Final average fitness  per plot
p4<-ggplot(plot_dat,aes(x=x,y=y,color=patch_dens))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_gradientn(colors=rainbow(5))+
  ylab("Y coordinate")+xlab("X Coordinate")+
  ggtitle("Final patch density per plot")
plot <- p1+p2+p3+p4
plot

hist(fit_x_m[,,dim(fit_x_m)[3]],breaks=15,xlab="Average fitness of each species",main="Average fitness of each species at the end of the simulation")

hist(fit_x_m[,,1],breaks=15,xlab="Average fitness of each species",main="Average fitness of each species at the start of the simulation")

hist((abs(fit_x_m[,,1])-abs(fit_x_m[,,dim(fit_x_m)[3]])),breaks=15,xlab="Average fitness increase of each species",main="Average fitness increase of each species during the simulation")

par(mfrow=c(1,1))
#### Step 7. Specific plots
#Here you can plot specific plots that you want to check







