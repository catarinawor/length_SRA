# ==============================================
#Title:plot_Udevs
#Author: Catarina Wor
#date: November 2016
#Function to plot deviation from true U
# code adapted from iSCAM stuff (Martell et al)
# ==============================================
#for testing delete, when done



#source("/Users/catarinawor/Documents/Lagrangian/R/read_mse/Rplots/calc_quantile.R")
#M<-SIMSdat

require(reshape2)
require(tidyr)
require(ggplot2)

source("/Users/catarinawor/Documents/Length_SRA/R/plots/readscn.R")


plot_udevs <- function( M )
{
	cat("plot_Udevs")

	n <- length( M )
	mdf <- NULL
	
	scn<-read_scnnames()



	conv_n<-numeric(length=length(scn))
	
	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){
			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1

			
			est<-c(M[[i]]$SArep$avgUy)

			st<-M[[i]]$OM$rep_yr-M[[i]]$OM$syr+1

			true<-c(M[[i]]$OM$avgUy[st:length(M[[i]]$OM$avgUy)])
			


			bias<- (est- true) / true
			
			#df <- data.frame(Ro=bias[1], Rinit=bias[2], kappa=bias[3],Linf=bias[4],k=bias[5],to=bias[6],cvl=bias[7])

			df <- data.frame(bias=bias,  true =true,est=est,year=M[[i]]$SArep$yr,scenario=scn[M[[i]]$OM$scnNumber])
			mdf <- rbind(mdf,df)

			
		}
	}

	summary(mdf)


	
	#df2<-melt(mdf,variable.name = "parameter")


	
	p <- ggplot(mdf) 
	p <- p + geom_boxplot(aes(x=as.factor(year),y=bias))
	p <- p + geom_hline(yintercept=0, color="darkred", size=1.2, alpha=0.3)
	p <- p + labs(x="year",y=" U Bias")
	p <- p + theme_bw(11) 
	p <- p + coord_cartesian(ylim=c(-0.5, 0.5))
	p <- p + annotate("text" , x = 1, y = 0.4, label = paste("n = ",conv_n))
	p <- p + facet_wrap(~scenario)
	print(p)
##
##	##setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	##ggsave("U_bias.pdf", plot=p)
	
}

