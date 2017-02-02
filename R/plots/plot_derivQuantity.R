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

	conv_n<-numeric(length=length(scn))

	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){

			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1


			est<-c(M[[i]]$SArep$depletion[length(M[[i]]$SArep$depletion)], 
				M[[i]]$SArep$q,
				M[[i]]$SArep$avgUy[length(M[[i]]$SArep$avgUy)])

			true<-c(M[[i]]$OM$depl[length(M[[i]]$OM$depl)], 
				M[[i]]$OM$q,
				M[[i]]$OM$avgUy[length(M[[i]]$OM$avgUy)])

			bias<- (est- true) / true

			df <- data.frame(Depletion=bias[1],q=bias[2],Uend=bias[3],scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			mdf <- rbind(mdf,df)
		}
	}

	
	df2<-melt(mdf,variable.name = "parameter",id=c("scenario","scnnumber"))

	df2$converge<-conv_n[df2$scnnumber]



	
	p <- ggplot(df2) 
	p <- p + geom_boxplot(aes(x=parameter,y=value, fill=parameter))
	p <- p + geom_hline(yintercept=0, color="darkred", size=1.2, alpha=0.3)
	p <- p + labs(x="Parameter",y="Bias")
	p <- p + ylim(-0.5, 0.5)
	p <- p + theme_bw(11)
	p <- p + facet_wrap(~scenario)
	p <- p + geom_text(data=df2, aes(x=1.0, y=0.44, label=converge), parse=TRUE)
	#p <- p + annotate("text" , x = 1.2, y = 0.4, label = paste("n = ",conv_n))
	print(p)

	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("derivQuant.pdf", plot=p)
	
}