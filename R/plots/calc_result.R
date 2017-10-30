#==============================================
#Title: calc_PM
#Author: Catarina Wor
#date: Oct 21th 2016
#Read results from MSE runs and calculate performanace measures
#
#==============================================

#TODO
#

#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/linf50priork"

#Vulpen=4
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/base2"

#Ro and Rinit estimate
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/Rinit"

#Ro and Rinit estimate eith narrower prior on kappa
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/Rinit_prior"

#same as above but with ssvul=1000
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/vulpen1000"

#no vulpen
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/novulpen"

#use prior in q and set sigR=2.0 in the last phase
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/wt2qprior"
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/smooth_sel"
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/smooth_sel_novpen"
DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/fspr"


#just Ro
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/linf50"

plotlib<-"/Users/catarinawor/Documents/Length_SRA/R/plots"
#.RFILES     <- list.files(.LIB,pattern="\\.[Rr]$")

SIMSdat<-list()

Rfiles <- list.files(DIR,pattern="\\.Rdata",full.name=TRUE)

plotfiles <- list.files(plotlib,pattern="plot",full.name=TRUE)
	
for(i in 1:length(Rfiles)){
	

	load(Rfiles[i])
	
	SIMSdat[[i]]<-sims


}


#load graphing routinesfor
for(p in 1:length(plotfiles)){
	source(plotfiles[p])
}

length(SIMSdat)
	
plot_params( SIMSdat ,Rinit=TRUE, sv=F, nome="")
#plot_params( SIMSdat ,Rinit=F,sv=F, nome="")
plot_params_publ( SIMSdat ,Rinit=TRUE)
plot_derivQuant( SIMSdat, sv=F, nome="novp" )
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