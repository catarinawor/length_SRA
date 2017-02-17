#==============================================
#Title: calc_PM
#Author: Catarina Wor
#date: Oct 21th 2016
#Read results form MSE runs and calculate performanace measures
#
#==============================================

#TODO
#

 
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/base"
DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/Rinit"


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
plot_derivQuant( SIMSdat )
plot_SBiomass( SIMSdat )
plot_udevs( SIMSdat )
plot_Recdevs( SIMSdat )
plot_LLvals( SIMSdat )
#plot_Sel( SIMSdat )
#plot_expRecdevs( SIMSdat )
#plotLengthComps( SIMSdat )	




