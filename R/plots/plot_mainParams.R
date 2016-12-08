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


plot_params <- function( M )
{
	cat("plot_params")

	scn<-read_scnnames()


	n <- length( M )
	mdf <- NULL
	adf <- NULL
	


	conv_n<-numeric(length=length(scn))

	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){
			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1


			est<-c(M[[i]]$SArep$Ro,
				M[[i]]$SArep$Rinit,
				M[[i]]$SArep$reck,
				M[[i]]$SArep$Linf,
				M[[i]]$SArep$k,
				M[[i]]$SArep$to,
				M[[i]]$SArep$cvl)
				#M[[i]]$SArep$maxUy[length(M[[i]]$SArep$maxUy)])

			true<-c(M[[i]]$OM$Ro,
				M[[i]]$OM$Rinit,
				M[[i]]$OM$reck,
				M[[i]]$OM$Linf,
				M[[i]]$OM$k,
				M[[i]]$OM$to,
				M[[i]]$OM$cvl)
				#M[[i]]$OM$true_depl[length(M[[i]]$OM$true_depl)], 
				#M[[i]]$OM$true_q,
				#M[[i]]$OM$true_ut[length(M[[i]]$OM$true_ut)])

			
				bias<- (est- true) / true
			
			#df <- data.frame(Ro=bias[1], Rinit=bias[2], kappa=bias[3],Linf=bias[4],k=bias[5],to=bias[6],cvl=bias[7])

			#df <- data.frame(Ro=bias[1],  kappa=bias[2],Linf=bias[3],k=bias[4],to=bias[5],cvl=bias[6], scenario=scn[M[[i]]$OM$scnNumber])
			df <- data.frame(Ro=bias[1], Rinit=bias[2], kappa=bias[3],Linf=bias[4],k=bias[5],to=bias[6],cvl=bias[7], scenario=scn[M[[i]]$OM$scnNumber])

			mdf <- rbind(mdf,df)

			af <- data.frame(true = true, est = est, param=c("Ro", "Rinit","kappa","Linf","k","to","cvl"))
			#af <- data.frame(true = true, est = est, param=c("Ro", "kappa","Linf","k","to","cvl"), scenario=scn[M[[i]]$OM$scnNumber])
			
			adf <- rbind(adf,af)
		}
	}


	
	df2<-melt(mdf,variable.name = "parameter",id="scenario")


	
	p <- ggplot(df2) 
	p <- p + geom_boxplot(aes(x=parameter,y=value, fill=parameter))
	p <- p + geom_hline(yintercept=0, color="darkred", size=1.2, alpha=0.3)
	p <- p + labs(x="Parameter",y="Bias")
	p <- p + theme_bw(11) 
	p <- p + ylim(-0.5, 0.5)
	p <- p + facet_wrap(~scenario)
	p <- p + annotate("text" , x = 1.2, y = 0.4, label = paste("n = ",conv_n))
	print(p)

	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("main_params.pdf", plot=p)
	
}