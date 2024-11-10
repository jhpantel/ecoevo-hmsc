## Code for Figure 4
## J.H. Pantel
## 09-11-2024

## libraries
library(ggplot2)
library(ggalluvial)
library(Hmsc)
library(vegan)
library(cowplot)

#### Figure 4a. Alpha diversity across d and h2 levels ####
## functions
CV <- function(x){
  a <- mean(x,na.rm=T)	#Get the average of all values
  b <- (x - a)^2					#Squared deviations from mean
  c <- sqrt((sum(b,na.rm=T))/(length(x[is.na(x) == F])-1)) #sum squared deviations, divide by sample N-1, take square root
  d <- c/a  #coefficient of variation
  return(CV=d)
}
standard_error <- function(x){
  st_err <- sd(x) / sqrt(length(x))	
  
  return(st_err)
}
## Read in and save relevant data
div_200 <- array(NA,dim=c(19,11,3,3),dimnames=list(c("d_zero","d_minus9","d_minus8","d_minus7","d_minus6","d_minus5","d_minus4","d_minus3","d_minus2","d_01","d_02","d_03","d_04","d_05","d_06","d_07","d_08","d_09","d_10"),c("h_0_","h_01_","h_02_","h_03_","h_04_","h_05_","h_06_","h_07_","h_08_","h_09_","h_10_"),c("alpha","gamma","beta"),c("mean","CV","SDE")))
d_lev <- c("d_zero","d_minus9","d_minus8","d_minus7","d_minus6","d_minus5","d_minus4","d_minus3","d_minus2","d_01","d_02","d_03","d_04","d_05","d_06","d_07","d_08","d_09","d_10")
h_lev <- c("h_0_","h_01_","h_02_","h_03_","h_04_","h_05_","h_06_","h_07_","h_08_","h_09_","h_10_")
## Inverse Simpson's diversity across additional heritability levels
for(i in 1:length(d_lev)){
  for(j in 1:length(h_lev)){
    ## Read in population size values
    result <- paste("./data/mc_v02/",h_lev[j],d_lev[i],"_res_v02.RData",sep="")
    tmp.env <- new.env()
    load(toString(result),envir=tmp.env)
    r <- get("N",pos=tmp.env) # alpha
    rm(tmp.env)
    r <- r[,,200]
    
    div_200[i,j,1,1] <- mean(diversity(r[rowSums(r) != 0,],index="invsimpson")) # mean
    div_200[i,j,1,2] <- CV(diversity(r[rowSums(r) != 0,],index="invsimpson")) # CV
    div_200[i,j,1,3] <- standard_error(diversity(r[rowSums(r) != 0,],index="invsimpson")) # SE
    
    # gamma
    div_200[i,j,2,1] <- if (is.null(dim(r[rowSums(r) != 0,])[1])) {NA} else if (dim(r[rowSums(r) != 0,])[1] < 2) {NA} else {diversity(colSums(r[rowSums(r) != 0,]),index="invsimpson")}
    
    # beta
    div_200[i,j,3,1] <- if (is.null(dim(r[rowSums(r) != 0,])[1])) {NA} else if (dim(r[rowSums(r) != 0,])[1] < 2) {NA} else {mean(diversity(colSums(r[rowSums(r) != 0,]),index="invsimpson") / diversity(r[rowSums(r) != 0,],index="invsimpson"))}
    
    div_200[i,j,3,2] <- if (is.null(dim(r[rowSums(r) != 0,])[1])) {NA} else if (dim(r[rowSums(r) != 0,])[1] < 2) {NA} else {CV(diversity(colSums(r[rowSums(r) != 0,]),index="invsimpson") / diversity(r[rowSums(r) != 0,],index="invsimpson"))}
    
    div_200[i,j,3,3] <- if (is.null(dim(r[rowSums(r) != 0,])[1])) {NA} else if (dim(r[rowSums(r) != 0,])[1] < 2) {NA} else {standard_error(diversity(colSums(r[rowSums(r) != 0,]),index="invsimpson") / diversity(r[rowSums(r) != 0,],index="invsimpson"))}
  }
}
## Plot
# Reshape for ggplot
mean.div_200 <- div_200[,,,1]
mean_div.gg <- reshape2::melt(mean.div_200)
se.div_200 <- div_200[,,,3]
se.div_200 <- reshape2::melt(se.div_200)
mean_div.gg$se <- se.div_200$value
colnames(mean_div.gg)[1:4] <- c("d_lev","h_lev","level","mean")
mean_div.gg$level <- as.factor(mean_div.gg$level)
h2_vals = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
mean_div.gg$h2 <- NA
for(i in 1:length(h2_vals)){
  mean_div.gg$h2[mean_div.gg$h_lev==unique(mean_div.gg$h_lev)[i]] <- h2_vals[i]
}

sub <- mean_div.gg[mean_div.gg$level=="gamma",]
sub2 <- sub[(sub$h2==0 | sub$h2==0.1 | sub$h2==0.3 | sub$h2==0.9),]
  
p0 <- ggplot(sub2,aes(x=d_lev,y=mean)) + geom_point(aes(colour=h2)) + scale_colour_gradient(low="lightgrey",high="forestgreen")+ geom_line(aes(group=interaction(level,h_lev)),col="black",linewidth=0.1) + coord_cartesian(ylim=c(0,15)) + theme_classic() + ylab(expression(gamma*" diversity")) + xlab("dispersal") + scale_x_discrete(labels=c(0,expression(10^-9),expression(10^-8),expression(10^-7),expression(10^-6),expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

#### Figure 4b. Stacked bar chart of mean VarPart results across d and h2 levels tested ####
## Enter VarPart data
data_bar <- data.frame(d=rep(NA,9),h2=rep(NA,9),Percentage=NA,level=rep(c("Env","deltaX","R(site)","R(time)"),times=9),col_col=rep(c("white","skyblue","darkgrey","orange"),times=9))

d_lev <- c("d_zero","d_minus9","d_minus8","d_minus7","d_minus6","d_minus5","d_minus4","d_minus3","d_minus2","d_01","d_02","d_03","d_04","d_05","d_06","d_07","d_08","d_09","d_10")
h_lev <- c("h_0_","h_01_","h_02_","h_03_","h_04_","h_05_","h_06_","h_07_","h_08_","h_09_","h_10_")
d_vals = c(0, 10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
h2_vals = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
## Get VarPart results across HMSC output
begin <- seq(1,36,4)
counter <- 1
for(v.d in c(1,4,8)){ # dlev
  for(v.h2 in c(2,4,10)){ # hlev
    ## Read in population size values
    result <- paste("./data/mc_v02/",h_lev[v.h2],d_lev[v.d],"_hmsc_v02.RData",sep="")
    load(result)
    s <- ncol(m.1$Y)
    VP <- computeVariancePartitioning(m.1, group = c(1,1,rep(2,s)), groupnames = c("Env","deltaX"))
    ng = dim(VP$vals)[1]
    leg = VP$groupnames
    for (r in 1:m.1$nr) {
      leg = c(leg, paste("Random: ", m.1$rLNames[r], sep = ""))
    }
    means = round(100 * rowMeans(VP$vals), 1)
    data_bar[begin[counter]:(begin[counter]+3),3] <- means
    data_bar[begin[counter]:(begin[counter]+3),1] <- d_vals[v.d]
    data_bar[begin[counter]:(begin[counter]+3),2] <- h2_vals[v.h2]
    counter <- counter + 1
  }
}

## h2=0.1
sub1 <- data_bar[data_bar$h2==0.1,]
sub1$d <- as.factor(sub1$d)
p1 <- ggplot(sub1,aes(y = Percentage, x = d, fill = as.character(level))) +
  geom_flow(aes(alluvium = level), alpha= .5, color = "white",
            curve_type = "xspline", 
            width = .5) +
  geom_col(width = .5,colour="black",linewidth=0.25) +
  scale_fill_manual(values = c("skyblue","white","darkgrey","orange")) +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),legend.position="none")
sub2 <- data_bar[data_bar$h2==0.3,]
sub2$d <- as.factor(sub2$d)
p2 <- ggplot(sub2,aes(y = Percentage, x = d, fill = as.character(level))) +
  geom_flow(aes(alluvium = level), alpha= .5, color = "white",
            curve_type = "xspline", 
            width = .5) +
  geom_col(width = .5,colour="black",linewidth=0.25) +
  scale_fill_manual(values = c("skyblue","white","darkgrey","orange")) +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),legend.position="none")
sub3 <- data_bar[data_bar$h2==0.9,]
sub3$d <- as.factor(sub3$d)
p3 <- ggplot(sub3,aes(y = Percentage, x = d, fill = as.character(level))) +
  geom_flow(aes(alluvium = level), alpha= .5, color = "white",
            curve_type = "xspline", 
            width = .5) +
  geom_col(width = .5,colour="black",linewidth=0.25) +
  scale_fill_manual(values = c("skyblue","white","darkgrey","orange")) +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),legend.position=c(0.8,0.5),legend.text = element_text(size=7),legend.title = element_blank())

#### Figure 4a. Conceptual MC figure ####
load(paste("./data/mc_v02/",h_lev[2],d_lev[8],"_res_v02.RData",sep=""))
mc <- data.frame(x=rep(NA,50*3),y=rep(NA,50*3),t=rep(c(0,10,200),each=50),div_a=rep(NA,50*3),xt=rep(NA,50*3))
# t0
mc[1:50,c(1,2)] <- xy
a  <- diversity(N[,,1],index="invsimpson")
a[is.finite(a)==FALSE] <- 0
mc[1:50,4] <- a # alpha diversity
b <- rowSums(xt[,,1]*N[,,1],na.rm=TRUE) / rowSums(N[,,1])
b[is.finite(b)==FALSE] <- 0
mc[1:50,5] <- b # community weighted mean trait value
# t10
mc[51:100,c(1,2)] <- xy
a  <- diversity(N[,,10],index="invsimpson")
a[is.finite(a)==FALSE] <- 0
mc[51:100,4] <- a # alpha diversity
b <- rowSums(xt[,,10]*N[,,10],na.rm=TRUE) / rowSums(N[,,10])
b[is.finite(b)==FALSE] <- 0
mc[51:100,5] <- b # community weighted mean trait value
# t200
mc[101:150,c(1,2)] <- xy
a  <- diversity(N[,,200],index="invsimpson")
a[is.finite(a)==FALSE] <- 0
mc[101:150,4] <- a # alpha diversity
b <- rowSums(xt[,,200]*N[,,200],na.rm=TRUE) / rowSums(N[,,200])
b[is.finite(b)==FALSE] <- 0
mc[101:150,5] <- b # community weighted mean trait value

c0 <- ggplot(mc,aes(x,y)) + geom_point(aes(color=xt,size=div_a)) + facet_wrap(~t) + scale_color_gradient(low = "lightgrey", high = "forestgreen") + theme(axis.text.x=element_blank(),axis.text.y=element_blank())

d <- plot_grid(plot_grid(c0,p0,ncol=1,rel_heights = c(1,1)),plot_grid(p1,p2,p3,nrow=3),ncol=2,rel_widths=c(2,1.2))

pdf("./output/fig4.pdf",paper="a4r")
d
dev.off()
