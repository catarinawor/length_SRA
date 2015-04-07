#Roberto smells funny
rm(list=ls()); gc(); #options("scipen"=100, "digits"=3)

seeds = seq(2+9, by=1, length=1000)
nsim = 1

# for(s in 1:nsim)
# {
#   if (Sys.info()["nodename"] =="sager") setwd("~/Dropbox/UBC_MSE")
#   
#   iseed = seeds[s]
#   tau = 0.1
#   use_hcr = "cr" 
#   use_ctl = "ctl_cr"
#   init_error = 0
#   rt_error = 1
#   obs_error = 1
#   Ftpattern= "updown"
# #   Ftpattern= "vpaF"
# #   Ftpattern= "F0"
# #   Ftpattern= "oneway_down"
# #   Ftpattern= "oneway"
# #   Ftpattern= "upconstant"
#   source("SEASM_16.R")
# }

### Commented out this section since we are not looking at VPA (at least for now)
##if (Sys.info()["nodename"] =="sager") setwd("~/Dropbox/MS_jurel/simVPA")
##
##datafilename = "Sim_SCAM_1ageold.csv"
##simJMdata =  read.csv(datafilename, header=T, sep=";", dec=".")
##minAgeVul = 6
##maxAgeVul = 12
##minage = 1; maxage = 12
##age = minage:maxage
##selName = "va"
##noTVsel = 1 # So, select va for yr 1
##newagecomp <- simJMdata[,paste0("a",age)]
##
##true_va = simJMdata[noTVsel, paste0(selName,age)]
##true_cal = simJMdata[,paste0("l",lbins)]
### true_cat = simJMdata[,paste0("truecat",paste0(".a",age))]
##true_ssb = simJMdata[,"true_ssb"]
##obs_yt   = simJMdata[,"true_allyt"]
##true_ct  = simJMdata[,"true_ct"]
##true_Ft2 = simJMdata[,"true_Ft"]
##true_rt = simJMdata[,"true_rt"]
##


#if (Sys.info()["nodename"] =="sager") setwd("~/Dropbox/MS_sra/simsra")
setwd("/Users/catarinawor/Documents/length_SRA/simeval/")


require(PBSmodelling)
source("read.admb.R")


Est_Tpl = "jmsra"
Sim_Tpl = "simsra"

system(paste('./',Sim_Tpl,' -ind ',Sim_Tpl,'.dat',sep=""), wait = TRUE)
input = read.rep("true_data_lsra.rep")
true_ct = input$true_ct
true_ut = input$true_ut
true_nat = input$true_Nat
true_cal = input$true_cal
true_rt = true_nat[,1]

system(paste('./',Est_Tpl,' -ind ',Est_Tpl,'.dat',sep=""), wait = TRUE)
out = read.admb(Est_Tpl)
names(out)
out$depletion

minAgeVul = 6
maxAgeVul = length(out$age)
minage = 1; maxage = 12
age = minage:maxage
nage = maxAgeVul

 #pdflabel = paste("simF",".pdf",sep="_")
 #pdf(file=pdflabel) 


 par(mfcol=c(3,3),mar=c(4,4,1,1),oma=c(1,1,1,3), las=1, cex=0.7)

#plot(out$yr,out$vt, xlab="Year", ylab="SSB /1e5", ylim=c(0,max(out$vbt)),pch=19)
#lines(out$yr,out$sbt, type="l", col="red", lwd=3)
#legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(19, -1), col=c(1,"red"), lwd=c(-1,3),bty="n")

plot(out$len,out$muUl ,pch=19, ylab="selectivity", xlim=c(range(out$len)), ylim=c(0,max(out$muUl)*1.1), xlab="length bin")


## Survey/Fishery time series
plot(out$yr,out$survB, xlab="Year", ylab="Relative abundance", ylim=c(0,max(out$survB)),pch=19)
idyt  <- which(out$yr %in% out$iyr)
lines(out$iyr, out$bt[idyt]*out$q, col="red", lwd=3)
legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(19, -1), col=c(1,"red"), lwd=c(-1,3),bty="n")

plot(out$yr,true_ct, xlab="Year", ylab="Catch", ylim=c(0,max(true_ct)), pch=19)
lines(out$yr, out$yield, col="red",lwd=3)
legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(19, -1), col=c(1,"red"), lwd=c(-1,3),bty="n")

Uhat =  apply(out$Uage[,minAgeVul:(nage-1)],1,mean)

plot(out$yr,true_ut, xlab="Year", ylab="Fishing Mortality", ylim=c(0,max(true_ut*1.5)), pch=19)
lines(out$yr, Uhat, lwd=3, col="red")
legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(19, -1), col=c(1,"red"), lwd=c(-1,3),bty="n")

plot(out$yr, true_rt, xlab="Year", ylab="Recruits 1 year-old / 1e7",pch=19, type="b", ylim=c(0,max(true_rt*1.1)))
lines(out$yr, out$N[,1], col="red", lwd=3)
legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(19, -1), col=c(1,"red"), lwd=c(-1,3),bty="n")

## Stock recruitment


  div = 1#1e7
  nyr = length(out$yr)
  st=seq(0, max(out$sbt), length=100)
  rt=(out$reca*st)/(1+out$recb*st)
  plot(out$sbt,out$N[1:nyr,1]/div,xlim=c(0,max(out$sbt)/div),ylim=c(0,max(out$N[1:nyr,1]/div)), 
       xlab="Spawning biomass/1e7", ylab="Age-1 recruits / 1e7",pch=19)
  lines(st/div, rt/div, type="l",lwd=3,col="red")



#pdf(file="length_comps_PERU_1980_2013.pdf") 

# ptrue_nat = prop.table(as.matrix(true_cal), 1)
# 
# propltcomps <- ptrue_nat
# yrs <- 1975:2012
# nbins = length(age)
# yrslt <- length(yrs)
# ## years with label
# argyrs <- c(1987,2000,2012)
# ## label for bars
# xlabel = 1:nage
# ## set number of labels for each bar (1 means each bar has a label)
# ibar = 1
# labelname = "age (years)"
# phat <- prop.table(as.matrix(out$Nlt), 1)
# 
# #x11()
# par(mfcol=c(13,3),mar=c(0,2,1,0),oma=c(4,3,1,0), las=0.7)
# 
# for (j in 1:(yrslt)) {
#   
#   ## put xlabel on argyrs
#   index <- which(yrs %in% argyrs)
#   xarg <- rep(0,length(yrs))
#   xarg[index] = yrs[index]
#   ## index for xlabel in sequence
#   isodd <- which(xlabel %% ibar != 0)
#   xlabel[isodd] <- NA
#   vsarg <- rep(0,length(yrs))
#   
#   if(xarg[j]!=0) arg = xlabel else arg = rep(NA,nbins-1)
#   if(xarg[j]!=0) Xarg = labelname else Xarg = NA
#   
#   barplot(propltcomps[j,], main=yrs[1]-1+j , names.arg = F, cex.axis=0.8, ylim=c(0,max(propltcomps)))
#   #          , legend.text=vsarg,args.legend = list(x = "topright", bty="n",bg = "n"))
#    xx <- barplot(propltcomps[1,],plot=F)
#   lines(xx, phat[j,], lwd=1.5, col="red")
#   axis(1, xx, labels=arg ,col.axis=1,tick=F, las=2, line=-0.5, cex.axis=0.85)
#   mtext(Xarg, 1, line=1.5, col="black", cex=0.75)
# }



# dev.off()


ptrue_nat = prop.table(as.matrix(true_cal), 1)

propltcomps <- ptrue_nat
yrs <- 1975:2013
nbins = dim(out$Nlt)[2]
yrslt <- dim(out$Nlt)[1]
## years with label
argyrs <- c(1987,2000,2012)
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





rm(list=ls()); 
if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
require(PBSmodelling)
source("read.admb.R")


Est_Tpl = "perujmsra"
Sim_Tpl = "simsra"

## Ft pattern 

Ft = c(0.009950166,0.043827127,0.076544905,0.108143165,0.138660214,0.168133048,0.196597397,0.224087769,0.250637491,0.276278750,
0.301042631,0.324959156,0.348057319,0.370365122,0.391909609,0.412716899,0.432812218,0.452219926,0.470963552,0.470963552,
0.452219926,0.432812218,0.412716899,0.391909609,0.370365122,0.348057319,0.324959156,0.301042631,0.276278750,0.250637491,
0.224087769,0.196597397,0.168133048,0.138660214,0.108143165,0.076544905,0.043827127,0.009950166)

## seed and store input and ouputs
seeds = seq(3, by=1, length=1000)
maxgrad_cr = NULL
hat_cr = NULL
ire_cr = NULL
re_cr = NULL
maxgrad_h = NULL
hat_h = NULL
ire_h = NULL
re_h = NULL
nsim = 50


for(s in 1:nsim) {

##  set seed and Ft  


saveSim = paste(Sim_Tpl,'.dat',sep="")
cat( file=saveSim, "## written:",date(),"\n", append=T )
cat( file=saveSim, "## seed\n", append=T )
cat( file=saveSim, seeds[s], "\n", append=T )
cat( file=saveSim, "## Ft\n", append=T )
cat( file=saveSim, Ft, "\n", append=T )

## run simulator
system(paste('./',Sim_Tpl,' -ind ',Sim_Tpl,'.dat',sep=""), wait = TRUE)
input = read.rep("true_data_lsra.rep")
names(input)
true_Ro = input$true_Ro
true_reck = input$true_reck
true_ct = input$true_ct
true_ut = input$true_ut

true_utend = input$true_ut[length(true_ut)]
true_nat = input$true_Nat
true_cal = input$true_cal
true_sbt = input$true_ssbt
true_depl = input$true_depl
true_rt = true_nat[,1]

## run estimator
system(paste('./',Est_Tpl,' -ind ',Est_Tpl,'.dat',sep=""), wait = TRUE)
out_cr = read.admb(Est_Tpl)
##names(out_cr)


## save sim-est outputs
maxgrad_cr <- rbind(maxgrad_cr, out_cr$fit$maxgrad)
##print(A$fit$std)
truePars = c(true_Ro,true_reck,true_depl,true_utend)
ihat_cr = c(out_cr$Ro,out_cr$kappa,out_cr$depletion,out_cr$maxUy[length(out_cr$maxUy)])
temp_ire_cr = (ihat_cr - truePars) / truePars
ire_cr <- rbind(ire_cr, temp_ire_cr)
hat_cr <- rbind(hat_cr, ihat_cr)

valid_maxgrad_cr = which(maxgrad_cr <= 0.0001)
valid_cr = which( hat_cr[,2] >= 2 & hat_cr[,2] <= true_reck*2)
valid_grad_cr = which( hat_cr[,2] >= 2 & hat_cr[,2] <= true_reck*2 & maxgrad_cr <= 0.0001)
if(s==nsim) { cat("# Valid Sim=", length(valid_grad_cr)) }

file.remove("simsra.dat")

}




plot_re = function(itheta,ivalid,h_cr,legend=T)  {
  
  colnames(itheta) = c("ro","kappa","Depletion","Uend")
  boxplot(itheta[ivalid,], ylim=c(-0.7,0.7), ylab="relative error", las=1)
  validSim = length(ivalid)
  abline(h=0)
  
  if(legend == T) {
    mtext(side=3,line=0,paste("steepness=",h_cr," Ft=","updowm"," simSigR=",0.6," simSS=",Nsample2," trueDepl=", round(trueDepl,2),"\n",
                              "simObsErr=",tau," Mt=",Mtpattern," Vul=",vapattern, " NSim=",validSim,"\n", "cvl=", l1cv," sigVul=", sigVul, 
                              " SimRhoR=",round(rhoR,2)," Fmsy=",round(Fmsy,3),"x",xF," true cr/h =",round(cr,3),"/",h ))
  }
  else
    mtext(side=3,line=0,paste("steepness=",h_cr," NSim=",validSim))
}



# pdflabel = paste("simF",Ftpattern,"iRE",".pdf",sep="_")

#pdf(file=pdflabel) 

par(mfcol=c(2,1),mar=c(4,1,1,1),oma=c(0,2,2.5,0), las=1)

plot_re(ire_cr,valid_cr,"CR",T)




