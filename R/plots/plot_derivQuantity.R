# ==============================================
#Title:plot_mainParams
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot true biomass fro lagrangian model.extract
# code adapted from iSCAM stuff (Martell et al)
# ==============================================
#for testing delete, when done

source("/Users/catarinawor/Documents/Length_SRA/R/plots/readscn.R")


#source("/Users/catarinawor/Documents/Lagrangian/R/read_mse/Rplots/calc_quantile.R")
#M<-SIMSdat

require(reshape2)
require(tidyr)
require(ggplot2)


plot_derivQuant <- function( M )
{
	cat("plot_derivQuant")

	n <- length( M )
	scn<-read_scnnames()


	mdf <- NULL

	conv_n= 0

	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){

			conv_n <-  conv_n + 1


			est<-c(M[[i]]$SArep$depletion[length(M[[i]]$SArep$depletion)], 
				M[[i]]$SArep$q,
				M[[i]]$SArep$maxUy[length(M[[i]]$SArep$maxUy)])

			true<-c(M[[i]]$OM$depl[length(M[[i]]$OM$depl)], 
				M[[i]]$OM$q,
				M[[i]]$OM$ut[length(M[[i]]$OM$ut)])

			bias<- (est- true) / true

			df <- data.frame(Depletion=bias[1],q=bias[2],Uend=bias[3],scenario=scn[M[[i]]$OM$scnNumber])
			mdf <- rbind(mdf,df)
		}
	}

	
	df2<-melt(mdf,variable.name = "parameter",id="scenario")


	
	p <- ggplot(df2) 
	p <- p + geom_boxplot(aes(x=parameter,y=value, fill=parameter))
	p <- p + geom_hline(yintercept=0, color="darkred", size=1.2, alpha=0.3)
	p <- p + labs(x="Parameter",y="Bias")
	p <- p + ylim(-0.5, 0.5)
	p <- p + theme_bw(11)
	p <- p + facet_wrap(~scenario)
	print(p)

	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("derivQuant.pdf", plot=p)
	
}