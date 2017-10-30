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


plot_derivQuant <- function( M, sv=F, nome="" )
{
	cat("plot_derivQuant")

	n <- length( M )
	scn<-read_scnnames()


	mdf <- NULL
	cmdf <- NULL


	conv_n<-numeric(length=length(scn))

	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){

			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1

			names(M[[i]]$SArep)

			est<-c(M[[i]]$SArep$depletion[length(M[[i]]$SArep$depletion)],
				M[[i]]$SArep$ytarget[length(M[[i]]$SArep$ytarget)],
				M[[i]]$SArep$utarget[length(M[[i]]$SArep$utarget)],
				M[[i]]$SArep$q)
				#M[[i]]$SArep$avgUy[length(M[[i]]$SArep$avgUy)])

			true<-c(M[[i]]$OM$depl[length(M[[i]]$OM$depl)],
			 	M[[i]]$OM$ytarget[length(M[[i]]$OM$ytarget)],
				M[[i]]$OM$utarget[length(M[[i]]$OM$utarget)],
				M[[i]]$OM$q)
				#M[[i]]$OM$avgUy[length(M[[i]]$OM$avgUy)])

			bias<- (est- true) / true

			cdf <- data.frame(true_msy=true[2],est_msy=est[2],true_umsy=true[3],est_umsy=est[3],scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			cmdf <- rbind(cmdf,cdf)
			
			udf <- data.frame(true_umsy=true[2],est_msy=est[2],true_umsy=true[3],est_umsy=est[3],scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			umdf <- rbind(cmdf,cdf)
			

			df <- data.frame(Depletion=bias[1],Yield_target=bias[2],U_target=bias[3],q=bias[4],scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			mdf <- rbind(mdf,df)
		}
	}


	head(cmdf)

	
	df2<-melt(mdf,variable.name = "parameter",id=c("scenario","scnnumber"))

	df2$converge<-conv_n[df2$scnnumber]

	head(df2)


	nrow(df2[df2$scenario=="VC"&df2$parameter=="umsy"&df2$value>1,])
	
	p <- ggplot(df2) 
	p <- p + geom_boxplot(aes(x=scenario,y=value, fill=parameter))
	p <- p + geom_hline(yintercept=0, color="darkred", size=1.2, alpha=0.3)
	p <- p + labs(x="Parameter",y="Bias")
	p <- p + coord_cartesian(ylim=c(-1.0, 1.0))
	p <- p + theme_bw(11)
	p <- p +  theme(axis.text = element_text(face="bold", size=12),
  axis.text.x= element_text(angle=45,hjust = 1),
  axis.title = element_text(face="bold", size=12))
	p <- p + facet_wrap(~parameter)
	p <- p + geom_text(data=df2, aes(x=scnnumber, y=0.44, label=converge), parse=TRUE)
	#p <- p + annotate("text" , x = 1.2, y = 0.4, label = paste("n = ",conv_n))
	print(p)

	if(sv==TRUE){
		setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
		setwd("/Users/catarinawor/Documents/Length_SRA/report")
		ggsave(paste(nome,"derivQuant.pdf",sep=""), plot=p)
	
	}
	
}



plot_derivQuant_publ <- function( M )
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
				M[[i]]$SArep$ytarget[length(M[[i]]$SArep$ytarget)],
				M[[i]]$SArep$utarget[length(M[[i]]$SArep$utarget)],
				M[[i]]$SArep$q)
				

			true<-c(M[[i]]$OM$depl[length(M[[i]]$OM$depl)],
			 	M[[i]]$OM$ytarget[length(M[[i]]$OM$ytarget)],
				M[[i]]$OM$utarget[length(M[[i]]$OM$utarget)],
				M[[i]]$OM$q)
				#M[[i]]$OM$avgUy[length(M[[i]]$OM$avgUy)])

			bias<- (est- true) / true

			df <- data.frame(Depletion=bias[1],MSY=bias[2],UMSY=bias[3],q=bias[4],scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			mdf <- rbind(mdf,df)
		}
	}

	
	df2<-melt(mdf,variable.name = "parameter",id=c("scenario","scnnumber"))
	summary(df2)

	levels(df2$parameter)<-c("Depletion",expression("Yield"["target"]), expression("U"["target"]),"q")
	df2$valuep<-df2$value*100
	df2$converge<-conv_n[df2$scnnumber]

	df2$scenario<-factor(df2$scenario,levels = rev(levels(df2$scenario)),ordered = TRUE)



	
	p <- ggplot(df2) 
	p <- p + geom_boxplot(aes(x=scenario,y=value))+coord_flip(ylim=c(-1, 1))
	p <- p + geom_hline(yintercept=0, color="black", size=1.2, alpha=0.3)
	p <- p + labs(x="Scenario",y="Relative Proportional Error")
	p <- p + facet_wrap(~parameter,ncol = 1,labeller = label_parsed)
	p <- p + theme_bw(12)
	p <- p + theme(axis.text = element_text(face="bold", size=12),
  axis.text.x= element_text(angle=0,hjust = .5),
  axis.title = element_text(face="bold", size=12),
  strip.text = element_text(face="bold", size=16))
	print(p)

	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("derivQuant_publ.pdf", plot=p)

	ggplot_build(p)$data

	ggplot_build(p)$data[[1]]$middle[ggplot_build(p)$data[[1]]$PANEL==1]
	ggplot_build(p)$data[[1]]$middle[ggplot_build(p)$data[[1]]$PANEL==2]
	ggplot_build(p)$data[[1]]$middle[ggplot_build(p)$data[[1]]$PANEL==3]
	ggplot_build(p)$data[[1]]$middle[ggplot_build(p)$data[[1]]$PANEL==4]
	
}

