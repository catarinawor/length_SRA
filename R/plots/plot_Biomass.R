# ==============================================
#Title:plot_mainParams
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot true biomass fro lagrangian model.extract
# code adapted from iSCAM stuff (Martell et al)
# ==============================================
#for testing delete, when done



source("/Users/catarinawor/Documents/Length_SRA/R/plots/calc_quantile.R")
source("/Users/catarinawor/Documents/Length_SRA/R/plots/readscn.R")

#M<-SIMSdat

require(reshape2)
require(tidyr)
require(ggplot2)


plot_SBiomass <- function( M )
{
	cat("plot_SBiomass")


	n <- length( M )

	scn<-read_scnnames()

	cib<-NULL
	
	for(i in 1:n){

		bio <- data.frame(v=c(M[[i]]$OM$"sbt"[(M[[i]]$OM$rep_yr):(M[[i]]$OM$eyr)]))
		bio2<-data.frame(t(bio),scenario=scn[M[[i]]$OM$scnNumber])
		cib <- rbind(cib,bio2)
	}

	summary(cib)


		median<-NULL
		low<-NULL
		high<-NULL
		scenario<-NULL
		year<-NULL

		

	for(m in 1:length(scn)){
		qq<-apply(cib[cib$scenario==scn[m],-c(ncol(cib))],2,calc_quantile)
		median<-c(median,qq["50%",])
		low<-c(low,qq["2.5%",])
		high<-c(high,qq["97.5%",])
		scenario<-c(scenario,rep(scn[m],length(qq["50%",])))
		year<-c(year,(M[[1]]$OM$rep_yr):M[[i]]$OM$eyr)
	}

	

	fdf<-data.frame (Median=median,Low=low, High=high,Year=year, scenario=scenario)
	summary(fdf)
	
	p <- ggplot(fdf,aes(x=Year,y=Median)) + geom_line()
	p <- p + geom_ribbon(aes(ymax=High, ymin=Low),alpha=0.2)
	p <- p + labs(x="Year",y="Total Biomass")
	p <- p + ylim(min(fdf$Low),max(fdf$High))
	p <- p + theme_bw(11)
	p <- p + facet_wrap(~scenario)
	print(p)


	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("Sbiomass.pdf", plot=p)
	
}