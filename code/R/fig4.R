## Code for Figure 4
# Author. J.H. Pantel &  R.J. Hermann
# Date 29.07.2025

## libraries
library(tidyverse)
library(Hmsc)
library(cowplot)
library(gridExtra)
library(corrplot)
library(ggpubr)

#Loading the empirical data
load("./data/emp/hmsc_emp.RData")

#Calculating variance partitioning of both models
VP = computeVariancePartitioning(m.1_evo, group = c(1,1,1,2,2,3,4,5,6), groupnames = c("C.Light","C.Nutrient","dLm","dRn","dSp","dWc"))
VP2 = computeVariancePartitioning(m.1_noevo, group = c(1,1,1,2,2), groupnames = c("C.Light","C.Nutrient"))

##### Figure 4a. Variance partitioning of evo model ####
par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 1.1, 11),oma=c(0,0,0,0))
plotVariancePartitioning(m.1_evo, VP = VP,cols=c("maroon4","gray47","darkorange","cyan3","firebrick1","mediumspringgreen","darkorchid4","midnightblue"),args.legend = list(x = "topright",inset=c(-0.5,0)),main=NULL)
p1 <- recordPlot()

##### Figure 4b. Variance partitioning of no_evo model ####
par(mar = c(5.1, 13, 1.1, 2.1),oma=c(0,0,0,0))
plotVariancePartitioning(m.1_noevo, VP = VP2,cols=c("maroon4","gray47","darkorchid4","midnightblue"),args.legend = list(x = "bottomleft",inset=c(-0.7,0)),main=NULL)
p2 <- recordPlot()

#################################
##### Figure 4c. Beta coefficients of the evo model with a 95% support level ####
par(mar = c(6, 4.1, 3, 2.1),oma=c(0,0,0,0))
#Have to rename to covaraites to have a more intuitive name
postBeta = getPostEstimate(m.1_evo, parName = "Beta")
plotBeta(m.1_evo, post = postBeta, param = "Support", supportLevel = 0.95)
p3 <- recordPlot()

###########################################################################################
##### Figure 4d. Predicted effect of evolution on species ####################################################################Influence of dLm on Sp

Gradient = constructGradient(m.1_evo,focalVariable = "dLm",
                             non.focalVariables = list("C.Light"=list(1),"Lightsq"=list(1),"C.Nutrients"=list(1),"Nutrientsq"=list(1),"dWc"=list(1),"dSp"=list(1),"dRn"=list(1)))

predY = predict(m.1_evo, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,ranLevels=Gradient$rLNew, expected=TRUE)

ma <- matrix(unlist(predY),ncol=1000)

da <- as.data.frame(t(exp(ma[41:60,])))


qa <- as.data.frame(apply(da,2,function(x) quantile(x,probs=c(0.025,0.975))))
me <- as.data.frame(apply(da,2,function(x) mean(x)))

Y_pred<-data.frame(q2.5=t(qa)[,1],q97.5=t(qa)[,2],mean=unname(me),dLm=Gradient$XDataNew$dLm)
Y_raw <- data.frame("Sp"=exp(m.1_evo$Y[,3]),"dLm"=m.1_evo$XData$dLm)

p4<-ggplot(Y_pred,aes(dLm))+
  geom_line(aes(y=mean))+
  geom_ribbon(aes(ymin=q2.5,ymax=q97.5),alpha=0.5,fill="lightblue")+
  geom_point(data=Y_raw,aes(x=dLm,y=Sp),color="grey")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none",strip.text = element_text(size=6)) + ylab("Count (Sp)")


#######################################
######Influence of dSp on Sp
Gradient = constructGradient(m.1_evo,focalVariable = "dSp",
                             non.focalVariables = list("C.Light"=list(1),"Lightsq"=list(1),"C.Nutrients"=list(1),"Nutrientsq"=list(1),"dWc"=list(1),"dLm"=list(1),"dRn"=list(1)))

predY = predict(m.1_evo, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,ranLevels=Gradient$rLNew, expected=TRUE)

ma <- matrix(unlist(predY),ncol=1000)

da <- as.data.frame(t(exp(ma[41:60,])))


qa <- as.data.frame(apply(da,2,function(x) quantile(x,probs=c(0.025,0.975))))
me <- as.data.frame(apply(da,2,function(x) mean(x)))

Y_pred<-data.frame(q2.5=t(qa)[,1],q97.5=t(qa)[,2],mean=unname(me),dSp=Gradient$XDataNew$dSp)
Y_raw <- data.frame("Sp"=exp(m.1_evo$Y[,3]),"dSp"=m.1_evo$XData$dSp)

p5<-ggplot(Y_pred,aes(dSp))+
  geom_line(aes(y=mean))+
  geom_ribbon(aes(ymin=q2.5,ymax=q97.5),alpha=0.5,fill="lightblue")+
  geom_point(data=Y_raw,aes(x=dSp,y=Sp),color="grey")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none",strip.text = element_text(size=6)) +ylab("Count (Sp)")

##################################
p1.g <- ggdraw(p1)
p2.g <- ggdraw(p2)
p3.g <- ggdraw(p3)
p4.g <- ggdraw(p4)
p5.g <- ggdraw(p5)

#Putting the plot together
d<-plot_grid(p1.g,p2.g,p3.g,plot_grid(p4.g,p5.g,nrow=2,align="v"),nrow=2,ncol=2,align="h",axis="b",labels = c('a','b','c','d'),label_size = 12,hjust=c(-3,-3,-10,0))

ggsave(file="./output/fig4.pdf",plot=d,width=297,height=210,units="mm")
