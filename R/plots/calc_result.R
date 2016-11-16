#==============================================
#Title: calc_PM
#Author: Catarina Wor
#date: Oct 21th 2016
#Read results form MSE runs and calculate performanace measures
#
#==============================================

#TODO
#

 
DIR<-"/Users/catarinawor/Documents/Length_SRA/R/simResult/base"

plotlib<-"/Users/catarinawor/Documents/Length_SRA/R/plots"
#.RFILES     <- list.files(.LIB,pattern="\\.[Rr]$")

SIMSdat<-list()

Rfiles <- list.files(DIR,pattern="\\.Rdata",full.name=TRUE)

plotfiles <- list.files(plotlib,pattern="plot",full.name=TRUE)
	
for(i in 1:length(Rfiles)){
	

	load(Rfiles[i])
	
	SIMSdat[[i]]<-sims

}


#load graphing groutines
source(plotfiles)



length(SIMSdat)
	
plot_params( SIMSdat )
plot_derivQuant( SIMSdat )


	




