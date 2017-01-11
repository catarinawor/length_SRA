# ==============================================
#Title:plot_mainParams
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot deviation from true U
# code adapted from iSCAM stuff (Martell et al)
# ==============================================
#for testing delete, when done

source("/Users/catarinawor/Documents/Length_SRA/R/plots/calc_quantile.R")
source("/Users/catarinawor/Documents/Length_SRA/R/plots/readscn.R")

require(reshape2)
require(tidyr)
require(ggplot2)


plotLengthComps <- function( M )
{
	#n <- length(M)
	cat("plotLengthComps\n")
	

	#M<-SIMSdat

	st<- M[[1]]$OM$rep_yr-M[[1]]$OM$syr+1

	

	mdf <- NULL
	
	scn<-read_scnnames()


	n<-length(M)
	allB<-NULL


		

		
		
		#B   <- lapply(1:length(id),getDF)
		
	for( i in 1:n ){

		getDF <- function()
		{
			
			df <- data.frame((M[[i]]$SArep$Ulength-M[[i]]$OM$Ulength[st:nrow(M[[i]]$SArep$Ulength),])/M[[i]]$OM$Ulength[st:nrow(M[[i]]$SArep$Ulength),])
			df <- data.frame(scenario=scn[M[[i]]$OM$scnNumber],year=M[[i]]$SArep$yr,df)
			len <- M[[i]]$SArep$len
			colnames(df) <- c("Scenario","Year",paste(len))
			
			return(df)
		}

		B<-getDF()
		
		
	
		allB <- rbind(allB,B)
	}

	anos<-unique(allB$Year)

	head(allB)
	dim(allB)

	median<-NULL

	for(ss in 1:length(scn)){
		for(y in 1:length(anos)){

		tmpy<-allB[allB$Scenario==scn[ss]&allB$Year==anos[y],-c(1,2)]

		qq<-apply(tmpy,2,calc_quantile)

		median<-rbind(median,cbind(data.frame(scenario=scn[ss],year=anos[y],level="median"),t(qq["50%",])))
		#low<-rbind(low,cbind(data.frame(scenario=scn[ss],year=anos[y],level="low"),t(qq["2.5%",])))
		#high<-rbind(high,cbind(data.frame(scenario=scn[ss],year=anos[y],level="high"),t(qq["97.5%",])))
		
	}
}
	#totB<-rbind(median,low,high)
	
	mB  <- melt(median,id.vars=c("scenario","year","level"))
	# mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex","AgeErr"))
	# BroodYear <- mdf$Year-as.double(mdf$variable)
	# mdf <- cbind(mdf,BroodYear)
	# print(head(mdf,3))

	p <- ggplot(mB,aes((year),variable,size=value))
	p <- p + geom_point(aes(colour=as.factor(sign(value))),alpha=0.75) 
	#p <- p + geom_point(alpha=0.75,) 

	p <- p + scale_size_area(max_size=5)
	p <- p + labs(x="Year",y="length",size="bias")
	p <- p + facet_wrap(~scenario,scales="free")
	#p <- p + scale_colour_discrete(guide="none")
	p <- p + theme_bw(11)
	print(p)

	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("UlengthMedian.pdf", plot=p)
}