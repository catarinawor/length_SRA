#========================================
#Compare Om and SA when just run over one interation

#============================================

setwd("/Users/catarinawor/Documents/length_SRA/R")
source("read.admb.R")

require(reshape2)
require(tidyr)
require(ggplot2)



setwd("/Users/catarinawor/Documents/length_SRA/admb/")
est <- read.rep("SA/runone.rep")
om <- read.rep("OM/true_data_lsra.rep")


names(om)
names(est)


#plot parameter values
ompar<-c(om$Ro,om$reck)
omwt<-om$wt
omit<-om$it

estpar<-c(est$Ro,est$reck)
estwt<-est$wt
estit<-est$pit

estpar<-c(est$Ro,est$reck)

par(mfrow=c(2,2))
plot(ompar, pch=16, ylab="parameter values", xlab="Ro and kappa")
points(estpar,col="red")
plot(omwt, pch=16,  ylab="recruitment deviations",xlab="years")
lines(estwt, lwd=2, col="blue")
plot(omit, pch=16,  ylab="recruitment deviations",xlab="years")
lines(estit, lwd=2, col="blue")

#plot derived quantities
omdq<-c(om$reca,om$recb,om$phie)
estdq<-c(est$reca,est$recb,est$phie)

par(mfrow=c(1,3))
for(n in 1:3){
plot(omdq[n], pch=16)
points(estdq[n],col="red")}

par(mfrow=c(2,1))
plot(om$sbt, pch=16, main="Spawning biomass")
lines(est$sbt, lwd=2, col="blue")
plot((est$sbt-om$sbt)/om$sbt, pch=16, main="Spawning biomass error")






#plot annual U
par(mfrow=c(2,1))
plot(om$maxUy, pch=16)
lines(est$maxUy, lwd=2, col="blue")
plot((est$maxUy-om$maxUy)/om$maxUy, pch=16, main="umax deviations",ylim=c(-0.0020,0.))
abline(h=0.0,col="red",lwd=3)




#plot U at length


#plot Nat
dff<-as.data.frame(rbind(data.frame(source=rep("om",50),year=1:50,om$Nat),data.frame(source=rep("est",50),year=1:50,est$Nat)))

colnames(dff) <- c("source" ,"year",paste(1:ncol(om$Nat)))
summary(dff)
dim(dff)

?melt

mB  <- melt(dff,measure.vars=c(paste(1:ncol(om$Nat))))
head(mB)


dff<-data.frame(source=rep("bias",50),year=1:50,(est$Nat-om$Nat)/om$Nat)

colnames(dff) <- c("source" ,"year",paste(1:ncol(om$Nat)))
summary(dff)
dim(dff)


mB  <- melt(dff,measure.vars=c(paste(1:ncol(om$Nat))))
head(mB)
	# mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex","AgeErr"))
	# BroodYear <- mdf$Year-as.double(mdf$variable)
	# mdf <- cbind(mdf,BroodYear)
	# print(head(mdf,3))

	p <- ggplot(mB,aes((year),variable,size=abs(value)))
	p <- p + geom_point(aes(colour=as.factor(sign(value))),alpha=0.75)
	p <- p + ylab("deviations in N at age")
	p 
	#p <- p + geom_point(alpha=0.75,) 



#plot Ulength

dim(est$Ulength)
dfU<-as.data.frame(rbind(data.frame(source=rep("om",50),year=1:50,om$Ulength),data.frame(source=rep("est",50),year=1:50,est$Ulength)))

colnames(dfU) <- c("source" ,"year",paste(1:ncol(om$Ulength)))
summary(dfU)
dim(dfU)



mBU  <- melt(dfU,measure.vars=c(paste(1:ncol(om$Ulength))))
head(mBU)


dfU<-data.frame(source=rep("bias",50),year=1:50,(est$Ulength-om$Ulength)/om$Ulength)

colnames(dfU) <- c("source" ,"year",paste(1:ncol(om$Ulength)))
summary(dfU)
dim(dfU)


mBU  <- melt(dfU,measure.vars=c(paste(1:ncol(om$Ulength))))
head(mBU)
	# mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex","AgeErr"))
	# BroodYear <- mdf$Year-as.double(mdf$variable)
	# mdf <- cbind(mdf,BroodYear)
	# print(head(mdf,3))

	pU <- ggplot(mBU,aes((year),variable,size=abs(value)))
	pU <- pU + geom_point(aes(colour=as.factor(sign(value))),alpha=0.75)
	pU <- pU + ylab("deviations in U at length")
	pU 
	#p <- p + geom_point(alpha=0.75,) 



