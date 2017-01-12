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

estpar<-c(est$Ro,est$reck)
estwt<-est$wt

par(mfrow=c(1,2))
plot(ompar, pch=16, ylab="parameter values", xlab="Ro and kappa")
points(estpar,col="red")
plot(omwt, pch=16,  ylab="recruitment deviations",xlab="years")
lines(estwt, lwd=2, col="blue")

#plot derived quantities
omdq<-c(om$reca,om$recb,om$phie)
estdq<-c(est$reca,est$recb,est$phie)

par(mfrow=c(1,3))
for(n in 1:3){
plot(omdq[n], pch=16)
points(estdq[n],col="red")
	
}

par(mfrow=c(1,1))
plot(om$sbt, pch=16, title="Spawning biomass")
lines(est$sbt, lwd=2, col="blue")




#plot annual U
par(mfrow=c(1,1))
plot(om$maxUy, pch=16)
lines(est$maxUy, lwd=2, col="blue")


par(mfrow=c(1,1))
plot((est$maxUy-om$maxUy)/om$maxUy, pch=16, title="umax deviations")





#plot U at length


#plot Nat
dff<-as.data.frame(rbind(data.frame(source=rep("om",50),year=1:50,om$Nat),data.frame(source=rep("est",50),year=1:50,est$Nat)))

colnames(dff) <- c("source" ,"year",paste(1:ncol(om$Nat)))
summary(dff)
dim(dff)

?melt

mB  <- melt(dff,measure.vars=c(paste(1:ncol(om$Nat))))
head(mB)
	# mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex","AgeErr"))
	# BroodYear <- mdf$Year-as.double(mdf$variable)
	# mdf <- cbind(mdf,BroodYear)
	# print(head(mdf,3))

	p <- ggplot(mB,aes((year),variable,size=value))
	p <- p + geom_point(aes(colour=source, shape=source),alpha=0.75)
	p 
	#p <- p + geom_point(alpha=0.75,) 
	print(p)

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

	p <- ggplot(mB,aes((year),variable,size=value))
	p <- p + geom_point(aes(colour=as.factor(sign(value))),alpha=0.75)
	p <- p + ylab("deviations in U at age")
	p 
	#p <- p + geom_point(alpha=0.75,) 



#plot Nlt