# ==============================================
#Title:plot_mainParams
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot true biomass fro lagrangian model.extract
# code adapted from iSCAM stuff (Martell et al)
# ==============================================
#for testing delete, when done



#source("/Users/catarinawor/Documents/Lagrangian/R/read_mse/Rplots/calc_quantile.R")
#M<-SIMSdat

require(reshape2)
require(tidyr)
require(ggplot2)


plot_params <- function( M )
{
	cat("plot_params")

	n <- length( M )
	mdf <- NULL


	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-03){

			est<-c(M[[i]]$SArep$Ro ,
				M[[i]]$SArep$reck,
				M[[i]]$SArep$depletion[length(M[[i]]$SArep$depletion)], 
				M[[i]]$SArep$q,
				M[[i]]$SArep$maxUy[length(M[[i]]$SArep$maxUy)])

			true<-c(M[[i]]$OM$true_Ro ,
				M[[i]]$OM$true_reck,
				M[[i]]$OM$true_depl[length(M[[i]]$OM$true_depl)], 
				M[[i]]$OM$true_q,
				M[[i]]$OM$true_ut[length(M[[i]]$OM$true_ut)])

			bias<- (est- true) / true

			df <- data.frame(Ro=bias[1] ,kappa=bias[2],Depletion=bias[3],q=bias[4],Uend=bias[5])
			mdf <- rbind(mdf,df)
		}
	}

	
	df2<-melt(mdf,variable.name = "parameter")


	
	p <- ggplot(df2) 
	p <- p + geom_boxplot(aes(x=parameter,y=value, fill=parameter))
	p <- p + geom_hline(yintercept=0, color="darkred", size=1.2, alpha=0.3)
	p <- p + labs(x="Parameter",y="Bias")
	p <- p + theme_bw(11)
	print(p)

	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("main_params.pdf", plot=p)
	
}