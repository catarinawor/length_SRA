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


plot_SBiomass <- function( M , sv=F, nome="")
{
	cat("plot_SBiomass")


	n <- length( M )

	scn<-read_scnnames()

	cib<-NULL

	cip<-NULL

	names(M[[i]]$SArep)
	
	for(i in 1:n){

		bio <- data.frame(v=c(M[[i]]$OM$"sbt"[(M[[i]]$OM$rep_yr):(M[[i]]$OM$eyr)]))
		bio2<-data.frame(t(bio),scenario=scn[M[[i]]$OM$scnNumber])
		cib <- rbind(cib,bio2)

		bip <- data.frame(v=c(M[[i]]$SArep$"sbt"))
		bip2<-data.frame(t(bip),scenario=scn[M[[i]]$OM$scnNumber])
		cip <- rbind(cip,bip2)
	}




		median<-NULL
		low<-NULL
		high<-NULL
		scenario<-NULL
		year<-NULL

		medianp<-NULL
		lowp<-NULL
		highp<-NULL
		scenariop<-NULL
		yearp<-NULL

		obs1<-NULL
		obs2<-NULL
		obs3<-NULL

		

	for(m in 1:length(scn)){
		qq<-apply(cib[cib$scenario==scn[m],-c(ncol(cib))],2,calc_quantile)
		median<-c(median,qq["50%",])
		low<-c(low,qq["2.5%",])
		high<-c(high,qq["97.5%",])
		scenario<-c(scenario,rep(scn[m],length(qq["50%",])))
		year<-c(year,(M[[i]]$OM$rep_yr):M[[i]]$OM$eyr)

		qqp<-apply(cip[cip$scenario==scn[m],-c(ncol(cip))],2,calc_quantile)
		medianp<-c(medianp,qqp["50%",])
		
		obs1<-c(obs1,unlist((cip[cip$scenario==scn[m],-c(ncol(cip))])[12,],use.names = FALSE))
		obs2<-c(obs2,unlist((cip[cip$scenario==scn[m],-c(ncol(cip))])[90,],use.names = FALSE))
		obs3<-c(obs3,unlist((cip[cip$scenario==scn[m],-c(ncol(cip))])[180,],use.names = FALSE))
		lowp<-c(lowp,qqp["2.5%",])
		highp<-c(highp,qqp["97.5%",])
		scenariop<-c(scenariop,rep(scn[m],length(qqp["50%",])))
		yearp<-c(yearp,M[[i]]$SArep$yr)
	}



	fdf<-data.frame(Median=c(median,rep(medianp,4)),Low=c(low,rep(lowp,4)), High=c(high,rep(highp,4)),
					Obs=c(rep(NA,length=length(median)),rep(NA,length=length(medianp)),as.vector(obs1),obs2,obs3),
					obsn=c(rep(NA,length=length(median)),rep(NA,length=length(medianp)),rep(c("1","2","3"),each=length(obs1))),	
					Year=c(year,rep(yearp,4)), scenario=c(scenario,rep(scenariop,4)), 
					type=rep(c("simulated",rep("estimated",4)),each=length(scenariop)))
	summary(fdf)
	
	p <- ggplot(fdf,aes(x=Year,y=Median,color=type,fill=type)) + geom_line(size=2)
	p <- p + geom_ribbon(aes(ymax=High, ymin=Low),alpha=0.2)
	#p <- p + geom_line(aes(x=Year, y=Obs,linetype=(obsn)),show_guide = FALSE)
	p <- p + labs(x="Year",y="Total Biomass")
	p <- p + theme_bw(16)
	p <- p + facet_wrap(~scenario,scale="free")
	p <- p + scale_colour_grey(start = 0., end = 0.7,labels = c("simulated", "estimated"))
	p <- p + scale_fill_grey(start = 0., end = 0.7,labels = c("simulated", "estimated"))
	p <- p + theme(axis.text = element_text(face="bold", size=12),
  			axis.title = element_text(face="bold", size=12),
  			strip.text = element_text(face="bold", size=15))
	print(p)





	if(sv==TRUE){
		setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
		ggsave(paste(nome,"SBiomass.pdf",sep=""), plot=p)
		
	}	
	
}


plot_SBiomass_pub <- function( M )
{
	cat("plot_SBiomass_pub")


	n <- length( M )

	scn<-read_scnnames()

	cib<-NULL

	cip<-NULL

	names(M[[i]]$SArep)
	
	for(i in 1:n){

		bio <- data.frame(v=c(M[[i]]$OM$"sbt"[(M[[i]]$OM$rep_yr):(M[[i]]$OM$eyr)]))
		bio2<-data.frame(t(bio),scenario=scn[M[[i]]$OM$scnNumber])
		cib <- rbind(cib,bio2)

		bip <- data.frame(v=c(M[[i]]$SArep$"sbt"))
		bip2<-data.frame(t(bip),scenario=scn[M[[i]]$OM$scnNumber])
		cip <- rbind(cip,bip2)
	}




		median<-NULL
		low<-NULL
		high<-NULL
		scenario<-NULL
		year<-NULL

		medianp<-NULL
		lowp<-NULL
		highp<-NULL
		scenariop<-NULL
		yearp<-NULL

		

	for(m in 1:length(scn)){
		qq<-apply(cib[cib$scenario==scn[m],-c(ncol(cib))],2,calc_quantile)
		median<-c(median,qq["50%",])
		low<-c(low,qq["2.5%",])
		high<-c(high,qq["97.5%",])
		scenario<-c(scenario,rep(scn[m],length(qq["50%",])))
		year<-c(year,(M[[i]]$OM$rep_yr):M[[i]]$OM$eyr)

		qqp<-apply(cip[cip$scenario==scn[m],-c(ncol(cip))],2,calc_quantile)
		medianp<-c(medianp,qqp["50%",])
		lowp<-c(lowp,qqp["2.5%",])
		highp<-c(highp,qqp["97.5%",])
		scenariop<-c(scenariop,rep(scn[m],length(qqp["50%",])))
		yearp<-c(yearp,M[[i]]$SArep$yr)
	}

	

	fdf<-data.frame (Median=c(median,medianp),Low=c(low,lowp), High=c(high,highp),Year=c(year,yearp), scenario=c(scenario,scenariop), type=rep(c("simulated","estimated"),each=length(scenariop)))
	summary(fdf)
	
	p <- ggplot(fdf,aes(x=Year,y=Median,color=type,fill=type)) + geom_line()
	p <- p + geom_ribbon(aes(ymax=High, ymin=Low),alpha=0.2)
	p <- p + labs(x="Year",y="Total Biomass")
	p <- p + theme_bw(11)
	p <- p + facet_wrap(~scenario,scale="free")
	print(p)


	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("SBiomass.pdf", plot=p)
	
}

