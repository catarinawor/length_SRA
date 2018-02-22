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


#appendix stuff Linf misspecification
#DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/biasedlinf"



#===================================================
#post revision directories
DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/noRinit"


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