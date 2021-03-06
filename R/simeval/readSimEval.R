
#==============================================
#Title: read_mse
#Author: Catarina Wor
#date: Oct 7th 2016
#Code to read OM and assessment results for  length_SRA
#
#==============================================
setwd("/Users/catarinawor/Documents/length_SRA/R")
source("read.admb.R")


readOutput <- function(dir)
{
	setwd(dir)
	assmt <- read.rep("SA/length_SRA.rep")
	par <- read.fit("SA/length_SRA")
	om <- read.rep("OM/true_data_lsra.rep")
	C <- list(SArep=assmt,SApar=par,OM=om)
	return( C );
}


	setwd("/Users/catarinawor/Documents/length_SRA/admb/")
	seed<- scan("seed.txt")	
	
	sims<-readOutput("/Users/catarinawor/Documents/length_SRA/admb/")
	
	#This will be the one with main paper results after the mistake has been corrected
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/base/")

	#SigR*2 after the mistake has been corrected 
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/base_sigR2/")
	
	#SigR*2 after the mistake has been corrected 
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/base_sigR2_q5/")
	
	#test with bias correction in q
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/test/")
	
	

	#This will be the one with biased Linf results after the mistake has been corrected
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/biasedlinf/")

	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/biasedlinfNew/")

	
	
	#no Rinit
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinit/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinitAEPtwo/")


	#no Rinit no k prior
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinit_nokprior/")

	#no Rinit no k prior no clt error and it err
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinit_nokprior_noCltErr/")
	
	#no Rinit no k and q prior, no clt error and no it er
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinit_noprior_noobserr/")
	
	#no Rinit no clt error 
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinit_noCltErr/")
	setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinit_noCltErrAEP/")
	
	#fast grow
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/fastgrowAEP/")
	

	#no Rinit, all error,  and no k and q prior done
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinit_noprior/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinit_nopriorAEP/")
		
	#no Rinit high clt error 
	# setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinithighcvIt/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinithighclterrAEP/")
	
	

	#no Rinit wrong k prior
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinitwrongkprior/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinitwrongkpriorAEP/")
	

	#no Rinit high survey cv
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinithighcvIt/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/noRinithighcvItAEP/")

	#no Rinit biased Linf
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/biasedlinf_noRinit/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/biasedlinf_noRinitAEP/")


	#FINAL
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/fspr/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/biasedlinf/")

		

	file.name <- paste("base",seed,".Rdata",sep="")
	
	save(sims,file=file.name)



##argsim = paste("./lagrangian_OM") 
##argest = paste("./lagrangian_est") 
###
##for(i in 1:5){
##	setwd("/Users/catarinawor/Documents/Lagrangian")
##	system(argsim)
##	system(argest)
##	seed<-scan("seed.txt")	
##	file.name <- paste("simest",seed,".Rdata",sep="")
##	sims<-readOutput("/Users/catarinawor/Documents/Lagrangian/")
##	setwd("/Users/catarinawor/Documents/Lagrangian/SimResult")
##	save(sims,file=file.name)
##}





#argsim = paste("./lagrangian_OM") 
#argest = paste("./lagrangian_est") 
#
#for(i in 1:2){
#
#	system(argsim)
# 	system(argest)
#	seed<-scan("seed.txt")	
#	file.name <- paste("simest",seed,".Rdata",sep="")
#	sims<-readOutput("/Users/catarinawor/Documents/Lagrangian/")
#	setwd("/Users/catarinawor/Documents/Lagrangian/SimResult")
#	save(sims,file=file.name)
#}
#
#	





#run.Simulation=function(N=1){
#	true_pars<-matrix(NA, nrow=N)
#	est_pars <-matrix(NA, nrow=N)
#	
#	argsim = paste("./lagrangian_OM") 
#	argest = paste("./lagrangian_est") 
#	
#	for(i in 1:N){
#		system(argsim)
#		system(argest)
#		sim = read.rep("lagrangian_OM.rep")
#		est = read.rep("lagrangian_est.rep")
#
#		true_pars[i,] <- c(sim$"mo", exp(sim$"log_tau_c"),sim$"maxPos50",sim$"maxPossd",sim$"cvPos")  
#		est_pars[i,] <- c(est$"mo",exp(est$"log_tau_c"),est$"maxPos50",est$"maxPossd",est$"cvPos")
#    
#    }
#}



#sim = read.rep("lagrangian_OM.rep")
#est = read.rep("lagrangian_est.rep")
#
#nomes <- names(sim)
