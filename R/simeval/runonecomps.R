#========================================
#Compare Om and SA when just run over one interation

#============================================

setwd("/Users/catarinawor/Documents/length_SRA/R")
source("read.admb.R")

require(reshape2)
require(tidyr)
require(ggplot2)



setwd("/Users/catarinawor/Documents/length_SRA/admb/")
est <- read.rep("SA/length_sra.rep")
om <- read.rep("OM/true_data_lsra.rep")


names(om)
names(est)

plot(est$iyr, est$survB)
lines(est$iyr,est$q*est$predSurvB, lwd=2, col="blue")

length(om$wt)
length(est$wt)

est$Rinit
 
omit<-om$it

est$pit
est$P_al-om$P_al

#plot parameter values
ompar<-c(om$Ro,om$Rinit,om$reck)
omwt<-om$wt
om$rep_yr


estpar<-c(est$Ro,est$Rinit,est$reck)
estpar
estwt<-est$wt

length(estwt)
length(omwt)

estit<-est$q*est$predSurvB


exp(est$wt[1])

est$lvec
est$pvec


estpar-ompar
sum(estwt-omwt)

par(mfrow=c(1,1))
plot(ompar, pch=16, ylab="parameter values", xlab="Ro and kappa")
points(estpar,col="red")

#setwd("/Users/catarinawor/Documents/length_SRA/report")
#pdf("rec_it_devs.pdf")
#pdf("rec_it_devs_werr.pdf")
#pdf("rec_it_devs_SA.pdf")
par(mfrow=c(2,2))
plot(omwt, pch=16,  ylab="recruitment deviations",xlab="years")
lines(estwt, lwd=2, col="blue")
plot((omit), pch=16,  ylab="it deviations",xlab="years")
lines((estit), lwd=2, col="blue")
plot((estwt-omwt), pch=16,  ylab="recruitment bias",xlab="years")
abline(h=0,lw=2, col="red")
plot(((estit)-(omit)), pch=16,  ylab="it bias",xlab="years")
abline(h=0,lw=2, col="red")
#dev.off()


par(mfrow=c(2,2))
plot(omwt, pch=16,  ylab="recruitment deviations",xlab="years")
lines(estwt, lwd=2, col="blue")
plot((omit), pch=16,  ylab="it deviations",xlab="years")
lines((estit), lwd=2, col="blue")
plot((estwt-omwt)/omwt, pch=16,  ylab="recruitment deviations",xlab="years")
abline(h=0,lw=2, col="red")
plot(((estit)-(omit))/(omit), pch=16,  ylab="it deviations",xlab="years")
abline(h=0,lw=2, col="red")


sum(estit-omit)
	(estit-omit)/omit

sum(estit)

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
#pdf("umax_devs.pdf")
#pdf("umax_devs_werr.pdf")
#pdf("umax_devs_SA.pdf")
par(mfrow=c(2,1))
plot(om$maxUy, pch=16,ylim=range(c(om$maxUy,est$maxUy)),ylab="Umax")
lines(est$maxUy, lwd=2, col="blue")
plot((est$maxUy-om$maxUy)/om$maxUy, pch=16, ylab="umax deviations")
abline(h=0.0,col="red",lwd=3)
#dev.off()

par(mfrow=c(2,1))
plot(om$avgUy, pch=16,ylim=range(c(om$avgUy,est$avgUy)),ylab="Uavg")
lines(est$avgUy, lwd=2, col="blue")
plot((est$avgUy-om$avgUy)/om$avgUy, pch=16, ylab="uavg deviations")
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

	pU <- ggplot(mBU,aes(year,variable,size=abs(value)))
	pU <- pU + geom_point(aes(colour=as.factor(sign(value))),alpha=0.75)
	pU <- pU + ylab("deviations in U at length")
	pU 
#ggsave("Ulength_SA.pdf", plot = pU)
#ggsave("Ulength_werr.pdf", plot = pU)
#ggsave("Ulength.pdf", plot = pU)
	#p <- p + geom_point(alpha=0.75,) 


#plot NLt

dim(est$Nlt)
dim(om$Nlt)



dfnlt<-data.frame(source=rep("bias",50),year=1:50,(est$Nlt-om$Nlt)/om$Nlt)

colnames(dfnlt) <- c("source" ,"year",paste(1:ncol(om$Nlt)))


mBNlt  <- melt(dfnlt,measure.vars=c(paste(1:ncol(om$Nlt))))
head(mBNlt)
	# mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex","AgeErr"))
	# BroodYear <- mdf$Year-as.double(mdf$variable)
	# mdf <- cbind(mdf,BroodYear)
	# print(head(mdf,3))

	pNlt <- ggplot(mBNlt,aes((year),variable,size=abs(value)))
	pNlt <- pNlt + geom_point(aes(colour=as.factor(sign(value))),alpha=0.75)
	pNlt <- pNlt + title("deviations in U at length")
	pNlt 
