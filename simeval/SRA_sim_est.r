#===================================================================================


rm(list=ls()); 
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/length_SRA/simeval/")
require(PBSmodelling)
source("read.admb.R")


Est_Tpl = "jmsra"
Sim_Tpl = "simsra"



## seed and store input and ouputs
maxgrad_cr = NULL
hat_cr = NULL
ire_cr = NULL
re_cr = NULL
maxgrad_h = NULL
hat_h = NULL
ire_h = NULL
re_h = NULL
nsim = 1


for(s in 1:nsim) {

##  set seed and Ft  


#saveSim = paste(Sim_Tpl,'.dat',sep="")
#cat( file=saveSim, "## written:",date(),"\n", append=T )
#cat( file=saveSim, "## seed\n", append=T )
#cat( file=saveSim, seeds[s], "\n", append=T )
#cat( file=saveSim, "## Ft\n", append=T )
#cat( file=saveSim, Ft, "\n", append=T )

## run simulator
system(paste('./',Sim_Tpl,' -ind ',Sim_Tpl,'.dat',sep=""), wait = TRUE)
input = read.rep("true_data_lsra.rep")
#names(input)
true_Ro = input$true_Ro
true_reck = input$true_reck
true_ct = input$true_ct
true_ut = input$true_ut

true_utend = input$true_ut[length(true_ut)]
true_nat = input$true_Nat
true_cal = input$true_cal
true_sbt = input$true_sbt
true_depl = input$true_depl
true_rt = true_nat[,1]
true_q = input$true_q
## run estimator
system(paste('./',Est_Tpl,' -ind ',Est_Tpl,'.dat',sep=""), wait = TRUE)
out_cr = read.admb(Est_Tpl)
##names(out_cr)


## save sim-est outputs
maxgrad_cr <- rbind(maxgrad_cr, out_cr$fit$maxgrad)
##print(A$fit$std)
truePars = c(true_Ro,true_reck,true_depl[length(true_depl)],true_utend,true_q) #need to output q in SRasim
ihat_cr = c(out_cr$Ro,out_cr$kappa,out_cr$depletion,out_cr$maxUy[length(out_cr$maxUy)],out_cr$q)
temp_ire_cr = (ihat_cr - truePars) / truePars
ire_cr <- rbind(ire_cr, temp_ire_cr)

hat_cr <- rbind(hat_cr, ihat_cr)

valid_maxgrad_cr = which(maxgrad_cr <= 0.0001)
valid_cr = which( hat_cr[,2] >= 2 & hat_cr[,2] <= true_reck*2)
valid_grad_cr = which( hat_cr[,2] >= 2 & hat_cr[,2] <= true_reck*2 & maxgrad_cr <= 0.0001)
#if(s==nsim) { cat("# Valid Sim=", length(valid_grad_cr)) }

#file.remove("simsra.dat")

}



plot_re = function(itheta,ivalid,h_cr,legend=T)  {
  
  colnames(itheta) = c("ro","kappa","Depletion","Uend", "q")
  boxplot(itheta[ivalid,], ylim=c(-0.7,0.7), ylab="relative error", las=1)
  validSim = length(ivalid)
  abline(h=0)
  
  if(legend == T) {
    #mtext(side=3,line=0,paste("steepness=",h_cr," Ft=","updowm"," simSigR=",0.6," simSS=","na," trueDepl=", round(trueDepl,2),"\n",
    #                          "simObsErr=",tau," Mt=",Mtpattern," Vul=",vapattern, " NSim=",validSim,"\n", "cvl=", l1cv," sigVul=", sigVul, 
    #                          " SimRhoR=",round(rhoR,2)," Fmsy=",round(Fmsy,3),"x",xF," true cr/h =",round(cr,3),"/",h ))
  }
  else
    mtext(side=3,line=0,paste("steepness=",h_cr," NSim=",validSim))
}


# pdflabel = paste("simF",Ftpattern,"iRE",".pdf",sep="_")

#pdf(file=pdflabel) 

par(mfcol=c(1,1),mar=c(4,1,1,1),oma=c(0,2,2.5,0), las=1)
plot_re(ire_cr,valid_cr,"CR",T)


colnames(ire_cr) = c("ro","kappa","Depletion","Uend", "q")
barplot(ire_cr[valid_cr,], ylim=c(-0.7,0.7), ylab="relative error", las=1)
  


