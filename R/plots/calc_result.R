#==============================================
#Title: calc_PM
#Author: Catarina Wor
#date: Oct 21th 2016
#Read results from MSE runs and calculate performanace measures
#
#==============================================

#TODO
#
#low Clt error tau=0.05
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/lowCltErr"


# final formulation of estimation model
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/fspr"

# final formulation of estimation model - with Pla correction
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/base"

# final formulation of estimation model - with Pla correction  and SigR*2 qpriorsd .2
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/base_sigR2"

# final formulation of estimation model - with Pla correction, SigR*2 qpriorsd .5
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/base_sigR2_q5"

# correction on q formulation of estimation model - with Pla correction
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/test"



#appendix stuff Linf misspecification
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/biasedlinf"
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/biasedlinfNew"



#===================================================
#post revision directories

#base case
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/noRinit"
DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/noRinitAEPtwo"

#wrong prior in kappa - 10
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/noRinitwrongkprior"


#high clt err
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/noRinithighclterr"
DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/noRinithighclterrAEP"



#no Rinit no k prior,  no clt error and it err
#DIR<-"/Users/catarinawor/Documents/length_SRA/R/simResult/noRinit_noprior_noobserr/"
	

#no Rinit no k prior, 
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/noRinit_nokprior"


#no Rinit high survey cv
#DIR<-"/Users/catarinawor/Documents/length_SRA/R/simResult/noRinithighcvIt/"


#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/biasedlinf_noRinit"
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/biasedlinf_noRinitAEP"



#fast grow
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/fastgrow"

plotlib<-"/Users/catarinawor/Documents/Length_SRA/R/plots"
#.RFILES     <- list.files(.LIB,pattern="\\.[Rr]$")

SIMSdat<-list()

Rfiles <- list.files(DIR,pattern="\\.Rdata",full.name=TRUE)

plotfiles <- list.files(plotlib,pattern="plot",full.name=TRUE)
	
for(i in 1:length(Rfiles)){
	

	load(Rfiles[i])
	
	SIMSdat[[i]]<-sims


}



for(p in 1:length(plotfiles)){
	source(plotfiles[p])
}

length(SIMSdat)
M<-SIMSdat


plot_params_deriv( SIMSdat ,Rinit=FALSE, sv=TRUE, nome="Figure1")	
#plot_params( SIMSdat ,Rinit=TRUE, sv=F, nome="")
#plot_params( SIMSdat ,Rinit=F,sv=F, nome="")
plot_params_publ( SIMSdat ,Rinit=R)
plot_derivQuant( SIMSdat, sv=F, nome="fix" )
plot_derivQuant_publ( SIMSdat )
plot_SBiomass( SIMSdat , sv=F, nome="novp")
plot_udevs( SIMSdat )
plot_Recdevs( SIMSdat )
#plot_udevs_all( SIMSdat )
#plot_LLvals( SIMSdat )
plot_Sel( SIMSdat ,  sv=TRUE, nome="fspr" )
plot_Ult_Clt( SIMSdat , sv=T, nome="fspr")
#plot_expRecdevs( SIMSdat )
#plotLengthComps( SIMSdat )
plot_ugone(SIMSdat )	


names(SIMSdat[[i]])
names(SIMSdat[[i]]$SArep)
SIMSdat[[i]]$SArep$sigVul