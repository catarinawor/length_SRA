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


plot_params <- function( M , Rinit=T )
{
	cat("plot_params")

	scn<-read_scnnames()

	n <- length( M )
	mdf <- NULL
	adf <- NULL
	


	conv_n<-numeric(length=length(scn))

	

		for(i in 1:n){


		if(M[[i]]$SApar$maxgrad<1.0e-03){
			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1

			if(Rinit==TRUE){
				est<-c(M[[i]]$SArep$Ro,
					M[[i]]$SArep$Rinit,
					M[[i]]$SArep$reck)
				
				true<-c(M[[i]]$OM$Ro,
					M[[i]]$OM$Rinit,
					M[[i]]$OM$reck)
				
			
				bias<- (est- true) / true
			
				#tv <- data.frame(Ro=est[1], Rinit=est[2], kappa=est[3], scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)

				df <- data.frame(Ro=bias[1],  Rinit=bias[2], kappa=bias[3], scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
				#df <- data.frame(Ro=bias[1],  kappa=bias[2], scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)

				mdf <- rbind(mdf,df)

				af <- data.frame(true = true, est = est, param=c("Ro", "Rinit" , "kappa"), scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
				#af <- data.frame(true = true, est = est, param=c("Ro",  "kappa"), scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			

				adf <- rbind(adf,af)
			}else{
				est<-c(M[[i]]$SArep$Ro,
					#M[[i]]$SArep$Rinit,
					M[[i]]$SArep$reck)
				
				true<-c(M[[i]]$OM$Ro,
					#M[[i]]$OM$Rinit,
					M[[i]]$OM$reck)
				
			
					bias<- (est- true) / true
			
				#tv <- data.frame(Ro=est[1], Rinit=est[2], kappa=est[3], scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)

				#df <- data.frame(Ro=bias[1],  Rinit=bias[2], kappa=bias[3], scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
				df <- data.frame(Ro=bias[1],  kappa=bias[2], scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)

				mdf <- rbind(mdf,df)

				#af <- data.frame(true = true, est = est, param=c("Ro", "Rinit" , "kappa"), scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
				af <- data.frame(true = true, est = est, param=c("Ro",  "kappa"), scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			

				adf <- rbind(adf,af)
			}
		}
	}
	
	#adf[adf$param=="kappa",]
	#adf[adf$param=="Ro",]
	#adf[adf$param=="Rinit",]#

	#adf[adf$param=="kappa"&adf$scnnumber==1,]
	


	df2<-melt(mdf,variable.name = "parameter",id=c("scenario","scnnumber"))

	df2$converge<-conv_n[df2$scnnumber]

	
	p <- ggplot(df2) 
	p <- p + geom_boxplot(aes(x=scenario,y=value, fill=parameter))
	p <- p + geom_hline(yintercept=0, color="darkred", size=1.2, alpha=0.3)
	p <- p + labs(x="Parameter",y="Bias")
	p <- p + theme_bw(11) 
	p <- p + theme(axis.text = element_text(face="bold", size=12),
  axis.text.x= element_text(angle=45,hjust = 1),
  axis.title = element_text(face="bold", size=12))
	p <- p + coord_cartesian(ylim=c(-1, 1))
	p <- p + facet_wrap(~parameter)
	p <- p + geom_text(data=df2, aes(x=scnnumber, y=0.44, label=converge), parse=TRUE)
	print(p)



	#setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	#ggsave("main_params.pdf", plot=p)
	
}


plot_params_publ <- function( M , Rinit=T )
{
	cat("plot_params_publ")

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
				M[[i]]$SArep$reck)
				
			true<-c(M[[i]]$OM$Ro,
				M[[i]]$OM$Rinit,
				M[[i]]$OM$reck)
				
			
				bias<- (est- true) / true
			
			#df <- data.frame(Ro=bias[1], Rinit=bias[2], kappa=bias[3],Linf=bias[4],k=bias[5],to=bias[6],cvl=bias[7])

			df <- data.frame(Ro=bias[1],  Rinit=bias[2], kappa=bias[3], scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			#df <- data.frame(Ro=bias[1], Rinit=bias[2], kappa=bias[3],Linf=bias[4],k=bias[5],to=bias[6],cvl=bias[7], scenario=scn[M[[i]]$OM$scnNumber])

			mdf <- rbind(mdf,df)

			#af <- data.frame(true = true, est = est, param=c("Ro", "Rinit","kappa","Linf","k","to","cvl"))
			af <- data.frame(true = true, est = est, param=c("Ro", "Rinit" , "kappa"), scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			
			adf <- rbind(adf,af)
		}
	}

	
	df2<-melt(mdf,variable.name = "parameter",id=c("scenario","scnnumber"))

	#levels(df2$parameter)<-expression(c(R_[0],R_[init],kappa ))
	df2$converge<-conv_n[df2$scnnumber]
	df2$valuep<-df2$value*100
	summary(df2) 



	p <- ggplot(df2) 
	p <- p + geom_boxplot(aes(x=scenario,y=valuep))+coord_flip()
	p <- p + geom_hline(yintercept=0, color="black", size=1.2, alpha=0.3)
	p <- p + labs(x="Scenario",y="% Relative Error")
	p <- p + theme_bw(12) 
	p <- p + facet_wrap(~parameter,ncol = 1,labeller = label_parsed)
	p <- p + theme(axis.text = element_text(face="bold", size=12),
  axis.title = element_text(face="bold", size=12),
  strip.text = element_text(face="bold", size=16))
	print(p)




	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("main_params_publ.pdf", plot=p)

	#ggplot_build(p)$

	names(ggplot_build(p)$data[[1]])

	ggplot_build(p)$data[[1]]$middle[ggplot_build(p)$data[[1]]$PANEL==1]
	ggplot_build(p)$data[[1]]$middle[ggplot_build(p)$data[[1]]$PANEL==2]
	ggplot_build(p)$data[[1]]$middle[ggplot_build(p)$data[[1]]$PANEL==3]
}