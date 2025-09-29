# Code for Box1 Figure
# Author. J.H. Pantel &  R.J. Hermann
# Date 29.07.2025

#### Libraries ---------------------------------------------------------------
library(bayesplot)
library(tidyverse)
library(gridExtra)
library(cowplot)

#### Load the data
load("data/mc/h_01_d_minus3_hmsc.RData")
m.post = Hmsc::convertToCodaObject(m.1)
m.df <- as.data.frame(do.call(rbind, m.post$Beta))

mycols <- c("maroon4","coral3","peachpuff3","gray47","hotpink2","darkorange","cyan3","azure2","firebrick1","azure2","mediumslateblue","mediumspringgreen","darkorchid4","azure2","rosybrown1")
############################################################
#Checking the beta effect of species 15 on 12
#With pre-checking we could find that both species only co-occur in patches patch 4,7,41, 46
par(mai=c(0,0,0,0),oma=c(2,2,0.5,0.5),mar=c(0,0,0,0))
#For beta of |deltax| species 15 on species 12
p1<-bayesplot::mcmc_areas(as.data.frame(m.df$`B[x15 (C14), y12 (S9)]`)) + labs(x="Posterior distribution",y="Density")+theme_classic()+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())+ggtitle(expression("|"~Delta~"x12,15|"))

for(j in c(4)){
  plot(N[j,1,],col=mycols[1],ylim=c(0,max(N,na.rm=T)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(4))
    axis(2,col="grey40",col.axis="grey20",labels=c(0,100,200,300,400,500),at=c(0,100,200,300,400,500),las = 1)
  if (j %in% c(4))
    axis(1,col="grey40",col.axis="grey20",labels=c(0,100,200),at=c(0,100,200),las = 1)
  for(i in 1:ncol(N)){
    points(N[j,i,],col=mycols[i])
    text(x=(dim(N)[3]),y=(max(N,na.rm=T)/5),labels=j,pos=2,col="black")
  }
  for(i in c(12,15)){
    points(N[j,i,],col=mycols[i])
  }
}
p2 <- recordPlot()
for(j in c(7)){
  plot(N[j,1,],col=mycols[1],ylim=c(0,max(N,na.rm=T)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(7))
    axis(2,col="grey40",col.axis="grey20",labels=c(0,100,200,300,400,500),at=c(0,100,200,300,400,500),las = 1)
  if (j %in% c(7))
    axis(1,col="grey40",col.axis="grey20",labels=c(0,100,200),at=c(0,100,200),las = 1)
  for(i in 1:ncol(N)){
    points(N[j,i,],col=mycols[i])
    text(x=(dim(N)[3]),y=(max(N,na.rm=T)/5),labels=j,pos=2,col="black")
  }
  for(i in c(12,15)){
    points(N[j,i,],col=mycols[i])
  }
}
p3 <- recordPlot()


for(j in c(4)){
  plot(delta_x[j,1,],col=mycols[1],ylim=c(0,0.11),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(4))
    axis(2,col="grey40",col.axis="grey20",labels=T,las = 1)
  if (j %in% c(4))
    axis(1,col="grey40",col.axis="grey20",labels=c(0,100,200),at=c(0,100,200),las = 1)
  for(i in 15){
    points(delta_x[j,i,],col=mycols[i])
    text(x=(dim(N)[3]),y=(max(delta_x,na.rm=T)/5),labels=j,pos=2,col="black")
  }
}
p4 <- recordPlot()

for(j in c(7)){
  plot(delta_x[j,1,],col=mycols[1],ylim=c(0,0.11),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(7))
    axis(2,col="grey40",col.axis="grey20",labels=T,las = 1)
  if (j %in% c(7))
    axis(1,col="grey40",col.axis="grey20",labels=c(0,100,200),at=c(0,100,200),las = 1)
  for(i in 15){
    points(delta_x[j,i,],col=mycols[i])
    text(x=(dim(N)[3]),y=(max(delta_x,na.rm=T)/5),labels=j,pos=2,col="black")
  }
}
p5 <- recordPlot()


######################################################################
mycols <- c("maroon4","coral3","peachpuff3","gray47","hotpink2","darkorange","cyan3","azure2","azure2","azure2","mediumslateblue","mediumspringgreen","darkorchid4","midnightblue","azure2")
#Checking the beta effect of species 3 on 4
#With pre-chcking we could find that both species co-occur in patches 2, 26, 31, 32, 33,30, 43, 48 from which patches 2 & 30 are shown
#For beta of |deltax| species 3 on species 4
load("data/mc/h_09_d_minus3_hmsc.RData")
m.post = Hmsc::convertToCodaObject(m.1)
m.df <- as.data.frame(do.call(rbind, m.post$Beta))


p6<-bayesplot::mcmc_areas(as.data.frame(m.df$`B[x3 (C4), y4 (S3)]`)) + labs(x="Posterior distribution",y="Density")+theme_classic()+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())+ggtitle(expression("|"~Delta~"x4,3|"))



for(j in c(2)){
  plot(N[j,1,],col=mycols[1],ylim=c(0,max(N,na.rm=T)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(2))
    axis(2,col="grey40",col.axis="grey20",labels=c(0,100,200,300,400,500),at=c(0,100,200,300,400,500),las = 1)
  if (j %in% c(2))
    axis(1,col="grey40",col.axis="grey20",labels=c(0,100,200),at=c(0,100,200),las = 1)
  for(i in 1:ncol(N)){
    points(N[j,i,],col=mycols[i])
    text(x=(dim(N)[3]),y=(max(N,na.rm=T)/2),labels=j,pos=2,col="black")
  }
  for(i in 3:4){
    points(N[j,i,],col=mycols[i])
  }
}
p7 <- recordPlot()

for(j in c(26)){
  plot(N[j,1,],col=mycols[1],ylim=c(0,max(N,na.rm=T)),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(26))
    axis(2,col="grey40",col.axis="grey20",labels=c(0,100,200,300,400,500),at=c(0,100,200,300,400,500),las = 1)
  if (j %in% c(26))
    axis(1,col="grey40",col.axis="grey20",labels=c(0,100,200),at=c(0,100,200),las = 1)
  for(i in 1:ncol(N)){
    points(N[j,i,],col=mycols[i])
    text(x=(dim(N)[3]),y=(max(N,na.rm=T)/2),labels=j,pos=2,col="black")
  }
  for(i in 3:4){
    points(N[j,i,],col=mycols[i])
  }
}
p8 <- recordPlot()


for(j in c(2)){
  plot(delta_x[j,1,],col=mycols[1],ylim=c(0,0.11),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(2))
    axis(2,col="grey40",col.axis="grey20",labels=T,las = 1)
  if (j %in% c(2))
    axis(1,col="grey40",col.axis="grey20",labels=c(0,100,200),at=c(0,100,200),las = 1)
  for(i in 3){
    points(delta_x[j,i,],col=mycols[i])
    text(x=(dim(delta_x)[3]),y=(max(delta_x,na.rm=T)/20),labels=j,pos=2,col="black")
  }
}
p9 <- recordPlot()

for(j in c(26)){
  plot(delta_x[j,1,],col=mycols[1],ylim=c(0,0.11),axes=FALSE,type="n",ylab=NA,xlab=NA,xaxs="i",yaxs="i")
  if (j %in% c(26))
    axis(2,col="grey40",col.axis="grey20",labels=T,las = 1)
  if (j %in% c(26))
    axis(1,col="grey40",col.axis="grey20",labels=c(0,100,200),at=c(0,100,200),las = 1)
  for(i in 3){
    points(delta_x[j,i,],col=mycols[i])
    text(x=(dim(delta_x)[3]),y=(max(delta_x,na.rm=T)/20),labels=j,pos=2,col="black")
  }
}
p10 <- recordPlot()


p1.g <- ggdraw(p1)
p2.g <- ggdraw(p2)
p3.g <- ggdraw(p3)
p4.g <- ggdraw(p4)
p5.g <- ggdraw(p5)
p6.g <- ggdraw(p6)
p7.g <- ggdraw(p7)
p8.g <- ggdraw(p8)
p9.g <- ggdraw(p9)
p10.g <- ggdraw(p10)


#Putting plot together
d<-plot_grid(p1.g,plot_grid(p2.g,p3.g,nrow=2,align="v"),plot_grid(p4.g,p5.g,nrow=2,align="v"),p6.g,plot_grid(p7.g,p8.g,nrow=2,align="v"),plot_grid(p9.g,p10.g,nrow=2,align="v"),nrow=2,ncol=3,align="h",axis="b",labels = c('a','b','c','d','e','f'),label_size = 12,hjust=c(-1,1,-0.5),rel_widths = c(1.5,1,1))

ggsave(file="./output/box1_fig.pdf",plot=d,width=297,height=210,units="mm")
