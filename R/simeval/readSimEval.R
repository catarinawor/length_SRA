
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
	
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/base/")
	#Rinint  is the one in the results
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/Rinit/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/qpr/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/new_bc/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/vulpen/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/vulpen10/")
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/vulpen2/")
	setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/vulpen6/")
	
	###setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/no_bc/")
	
	#setwd("/Users/catarinawor/Documents/length_SRA/R/simResult/cann/")


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
