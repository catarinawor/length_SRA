## test GIT 
if (Sys.info()["nodename"] == "sager") setwd("~/Dropbox/MS_sra/perujmsra")
require(PBSmodelling)
source("read.admb.R")

obsdata = read.csv("peruobsdata.csv",sep=";",header=T)

minAgeVul = 2
maxAgeVul = 4

system("./perujmsra -nox -ind perujmsra.dat")
#system("./perujmsra -nox -ind perujmsra2.dat")

out = read.admb("perujmsra")
names(out)

# pdflabel = paste("parSet2",".pdf",sep="_")
 pdflabel = paste("parSet1",".pdf",sep="_")
 
pdf(file=pdflabel) 


par(mfcol=c(3,3),mar=c(4,4,1,1),oma=c(1,1,1,3), las=1, cex=0.7)

plot(out$yr,out$sbt/100000, xlab="Year", ylab="SSB /1e5", ylim=c(0,max(out$sbt/100000)),pch=19, type="b")
# plot(out$yr, true_ssb/true_ssb[1], xlab="Year", ylab="spawning depletion", ylim=c(0,1),pch=19)
# lines(out$yr,out$sbt/out$sbt[1], type="l", col="red", lwd=3)
#legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(19, -1), col=c(1,"red"), lwd=c(-1,3),bty="n")

plot(out$len,out$muUl ,pch=19, ylab="selectivity", xlim=c(range(out$len)), ylim=c(0,max(out$muUl)*1.1), xlab="length bin", type="b")

## Survey/Fishery time series
plot(out$iyr,out$yt/1e7, xlab="Year", ylab="Relative abundance / 1e7", ylim=c(0,max(out$yt/1e7)),pch=19)
idyt  <- which(out$yr %in% out$iyr)
lines(out$iyr, out$bt[idyt]*out$q/1e7, col="red", lwd=3)
legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(19, -1), col=c(1,"red"), lwd=c(-1,3),bty="n")

# plot(out$yr,obsdata$ct[11:44], xlab="Year", ylab="Catch / 1e5", ylim=c(0,max(obsdata$ct*1000)), pch=19, type="l")
# lines(out$yr, out$yit, col="red",lwd=3, type="l")
# legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(19, -1), col=c(1,"red"), lwd=c(-1,3),bty="n")

Uhat =  apply(out$Uage[,minAgeVul:maxAgeVul],1,mean)
plot(out$yr,Uhat, xlab="Year", ylab="Fishing Mortality", ylim=c(0,max(Uhat*1.5)), pch=19, type="l")

plot(out$yr, out$N[,1], xlab="Year", ylab="Recruits 1 year-old / 1e7",pch=19, type="b", ylim=c(0,max(out$N[,1]*1.1)))


## Stock recruitment
  div = 1e7
  st=seq(0, max(out$sbt), length=100)
  rt=out$kappa*out$Ro*st/(out$bo+(out$kappa-1)*st)
  plot(out$sbt/1e7,out$N[1:max(out$yr),1]/div,xlim=c(0,max(out$sbt)/div),ylim=c(0,max(out$N[1:max(out$yr),1]/div)), 
       xlab="Spawning biomass/1e7", ylab="Age-1 recruits / 1e7",pch=19)
  lines(st/div, rt/div, type="l",lwd=3,col="red")


## more test ##
Ufull = NULL
select_lt_bins = 10:40
for(y in 1:length(out$yr)) { 
  Ufull[y] = sum(out$Clt[y,select_lt_bins])/sum(out$Nlt[y,select_lt_bins])
}

plot(Ufull,type="b")

## vlt = (Clt/Nlt)/Ufull
vlt = matrix(NA,nrow=length(out$yr),ncol=length(out$len))
cn = out$Clt/out$Nlt
for(y in 1:length(out$yr)) { 
  vlt[y,] = cn[y,]/Ufull[y]
}

matplot(t(vlt/apply(vlt,1,max)),type="l",col="black",lwd=1)

## yield = (Nlt*vlt*wl)*Ufull
# awl = 2.38e-05 
# bwl = 2.7671 
# wl = awl*(out$len)^bwl
wlt = matrix(NA,nrow=length(out$yr),ncol=length(out$len))
for(y in 1:length(out$yr)) { 
  wlt[y,] = out$wl
}

ytl = (out$Nlt*vlt*out$wl)*Ufull
yieldt = apply(ytl,1,sum)
plot(out$yr,obsdata$ct[11:44], xlab="Year", ylab="Catch / 1e5", ylim=c(0,max(obsdata$ct)), pch=19)
lines(yieldt/1000,type="l",col="red",lwd=2)
legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(19, -1), col=c(1,"red"), lwd=c(-1,3),bty="n")


sum((out$Nlt[1,]*vlt[1,]*out$wl))*Ufull[1]

dev.off()


## plot length comps 
pdf(file="length_comps_PERU_1980_2013_parSet1.pdf") 
#pdf(file="length_comps_PERU_1980_2013_parSet2.pdf") 

ptrue_nat = prop.table(as.matrix(out$Clt), 1)

propltcomps <- ptrue_nat
yrs <- 1980:2013
nbins = dim(out$Nlt)[2]
yrslt <- dim(out$Nlt)[1]
## years with label
argyrs <- c(1992,2005,2013)
## label for bars
xlabel = 1:dim(out$Nlt)[2]
## set number of labels for each bar (1 means each bar has a label)
ibar = 2
labelname = "length bins (cm)"
phat <- prop.table(as.matrix(out$Nlt*out$Ulength), 1)

#x11()
par(mfcol=c(13,3),mar=c(0,2,1,0),oma=c(4,3,1,0), las=0.7)

for (j in 1:(yrslt)) {
  
  ## put xlabel on argyrs
  index <- which(yrs %in% argyrs)
  xarg <- rep(0,length(yrs))
  xarg[index] = yrs[index]
  ## index for xlabel in sequence
  isodd <- which(xlabel %% ibar != 0)
  xlabel[isodd] <- NA
  vsarg <- rep(0,length(yrs))
  
  if(xarg[j]!=0) arg = xlabel else arg = rep(NA,nbins-1)
  if(xarg[j]!=0) Xarg = labelname else Xarg = NA
  
  barplot(propltcomps[j,], main=yrs[1]-1+j , names.arg = F, cex.axis=0.8, ylim=c(0,max(propltcomps)))
  #          , legend.text=vsarg,args.legend = list(x = "topright", bty="n",bg = "n"))
   xx <- barplot(propltcomps[1,],plot=F)
  lines(xx, phat[j,], lwd=1.5, col="red")
  axis(1, xx, labels=arg ,col.axis=1,tick=F, las=2, line=-0.5, cex.axis=0.85)
  mtext(Xarg, 1, line=1.5, col="black", cex=0.75)
}

 dev.off()

