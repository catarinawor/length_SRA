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


plot_LLvals <- function( M )
{
	cat("plot_LL")

	scn<-read_scnnames()


	n <- length( M )
	mdf <- NULL
	adf <- NULL
	


	conv_n<-numeric(length=length(scn))

	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){
			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1


			LLrec<-M[[i]]$SArep$pvec
			LLit<-M[[i]]$SArep$lvec
			
			#df <- data.frame(Ro=bias[1], Rinit=bias[2], kappa=bias[3],Linf=bias[4],k=bias[5],to=bias[6],cvl=bias[7])

			df <- data.frame(RecLL=LLrec,  ItLL=LLit, scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			#df <- data.frame(Ro=bias[1], Rinit=bias[2], kappa=bias[3],Linf=bias[4],k=bias[5],to=bias[6],cvl=bias[7], scenario=scn[M[[i]]$OM$scnNumber])

			mdf <- rbind(mdf,df)

			
		}
	}


	
	df2<-melt(mdf,variable.name = "LL",id=c("scenario","scnnumber"))

	df2$converge<-conv_n[df2$scnnumber]

	
	p <- ggplot(df2) 
	p <- p + geom_histogram(aes(x=value, fill=LL))
	p <- p + theme_bw(11) 
	p <- p + facet_wrap(scenario~LL,scales = "free",)
	p <- p + geom_text(data=df2, aes(x=1.0, y=0.44, label=converge), parse=TRUE)
	print(p)



	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("LL_vals.pdf", plot=p)
	
}