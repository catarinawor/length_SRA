#==============================================
#Title: calc_PM
#Author: Catarina Wor
#date: Oct 21th 2016
#Read results from MSE runs and calculate performanace measures
#
#==============================================

#TODO
#
 
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/base"


#this is the one it's currently "oficial"
DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/Rinit"
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/qpr"
##DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/no_bc"
#
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/new_bc"
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/vulpen"

#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/cann"


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
	
plot_params( SIMSdat ,Rinit=TRUE)
plot_params( SIMSdat ,Rinit=FALSE)
plot_params_publ( SIMSdat ,Rinit=TRUE)
plot_derivQuant( SIMSdat )
plot_derivQuant_publ( SIMSdat )
plot_SBiomass( SIMSdat )
plot_udevs( SIMSdat )
plot_Recdevs( SIMSdat )
#plot_udevs_all( SIMSdat )
#plot_LLvals( SIMSdat )
plot_Sel( SIMSdat )
#plot_expRecdevs( SIMSdat )
#plotLengthComps( SIMSdat )
plot_ugone(SIMSdat )	




