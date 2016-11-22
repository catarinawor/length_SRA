# ==============================================
#Title:plot_mainParams
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot true biomass fro lagrangian model.extract
# code adapted from iSCAM stuff (Martell et al)
# ==============================================
#for testing delete, when done



source("/Users/catarinawor/Documents/Length_SRA/R/plots/calc_quantile.R")
#M<-SIMSdat

require(reshape2)
require(tidyr)
require(ggplot2)


plot_SBiomass <- function( M )
{
	cat("plot_SBiomass")


	n <- length( M )

	cib <- data.frame(M[[1]]$OM$"sbt"[(M[[1]]$OM$rep_yr):(M[[1]]$OM$eyr)])

	for(i in 2:n){

		bio <- data.frame(M[[i]]$OM$"sbt"[(M[[i]]$OM$rep_yr):(M[[i]]$OM$eyr)])
		cib <- cbind(cib,bio)
	}

	qq<-apply(cib,1,calc_quantile)

	median<-qq["50%",]
	low<-qq["2.5%",]
	high<-qq["97.5%",]

	fdf<-data.frame (Median=median,Low=low, High=high,Year=(M[[1]]$OM$rep_yr):M[[i]]$OM$eyr)

	
	p <- ggplot(fdf,aes(x=Year,y=Median)) + geom_line()
	p <- p + geom_ribbon(aes(ymax=High, ymin=Low),alpha=0.2)
	p <- p + labs(x="Year",y="Total Biomass")
	p <- p + ylim(min(fdf$Low),max(fdf$High))
	p <- p + theme_bw(11)
	print(p)


	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("Sbiomass.pdf", plot=p)
	
}