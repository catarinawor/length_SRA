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

			
			#est<-c(M[[i]]$SArep$avgUy)
			est<-apply(M[[i]]$SArep$Ulength[,1:(ncol(M[[i]]$SArep$Ulength)-5)],1, mean)
			#est<-c(M[[i]]$SArep$maxUy)


			st<-M[[i]]$OM$rep_yr-M[[i]]$OM$syr+1


			true<-apply(M[[i]]$OM$Ulength[st:length(M[[i]]$OM$avgUy),1:(ncol(M[[i]]$OM$Ulength)-5)],1, mean)
			#true<-c(M[[i]]$OM$avgUy[st:length(M[[i]]$OM$avgUy)])
			#true<-c(M[[i]]$OM$maxUy[st:length(M[[i]]$OM$maxUy)])
			
			names(M[[i]]$OM)

			bias<- (est- true) / true
			
			#df <- data.frame(Ro=bias[1], Rinit=bias[2], kappa=bias[3],Linf=bias[4],k=bias[5],to=bias[6],cvl=bias[7])

			df <- data.frame(bias=bias,  true =true,est=est,year=M[[i]]$SArep$yr,scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			mdf <- rbind(mdf,df)

			
		}
	}

	mdf$converge=conv_n[mdf$scnnumber]
	summary(mdf)


	
	#df2<-melt(mdf,variable.name = "parameter")


	
	p <- ggplot(mdf) 
	p <- p + geom_boxplot(aes(x=as.factor(year),y=bias))
	p <- p + geom_hline(yintercept=0, color="darkred", size=1.2, alpha=0.3)
	p <- p + labs(x="year",y=" avgU bias")
	p <- p + theme_bw(11) 
	p <- p + coord_cartesian(ylim=c(-0.5, 0.5))
	#p <- p + annotate("text" , x = 1, y = 0.4, label = paste("n = ",conv_n))
	p <- p + facet_wrap(~scenario)
	#p <- p + geom_text(data=mdf, aes(x=1.2, y=0.48, label=converge), parse=TRUE)
	print(p)

	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave("U_bias.pdf", plot=p)
	
}

calcugone <- function(U){

	#U<-M[[i]]$SArep$Ulength

	if(sum(U>1.0001)>0.0){
		return(1.)

	}else{
		return(0.)
	}

}


plot_ugone <- function( M )
{
	cat("plot_Udevs")

	n <- length( M )
	mdf <- NULL
	
	scn<-read_scnnames()



	conv_n<-numeric(length=length(scn))
	Ugo<-NULL
	


	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){
			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1

			tUg<-calcugone(M[[i]]$SArep$Ulength)

			df<-data.frame(Uerr=tUg,scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)

			Ugo<-rbind(Ugo,df)
			
		}
	}

	summary(Ugo)

	Uone<-aggregate(Ugo$Uerr,list(Ugo$scenario),function(x){sum(x)/length(x)})

	return(Uone)
	
}



plot_udevs_all <- function( M )
{
	cat("plot_Udevs_all")

	n <- length( M )
	mdf <- NULL
	
	scn<-read_scnnames()



	conv_n<-numeric(length=length(scn))
	
	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){
			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1

			
			st<-M[[i]]$OM$rep_yr-M[[i]]$OM$syr+1
			#est<-c(M[[i]]$SArep$avgUy)
			est<-M[[i]]$SArep$Ulength

			true<-M[[i]]$OM$Ulength[(M[[i]]$OM$rep_yr):M[[i]]$OM$eyr,]



			bias<- (est- true) / true

		
			#df <- data.frame(Ro=bias[1], Rinit=bias[2], kappa=bias[3],Linf=bias[4],k=bias[5],to=bias[6],cvl=bias[7])

			df <- data.frame(bias=melt(bias)$value,year=melt(bias)$Var1,len=melt(bias)$Var2,scenario=scn[M[[i]]$OM$scnNumber],scnnumber=M[[i]]$OM$scnNumber)
			mdf <- rbind(mdf,df)

			
		}
	}

	mdf$converge=conv_n[mdf$scnnumber]
	summary(mdf)


	
	#df2<-melt(mdf,variable.name = "parameter")

	for(an in 1:length(unique(mdf$len))){

		mmdf<-mdf[mdf$len==unique(mdf$len)[an],]
		summary(mmdf)
	
	p <- ggplot(mmdf) 
	p <- p + geom_boxplot(aes(x=scenario,y=(bias),fill=scenario))
	p <- p + geom_hline(yintercept=0, color="darkred", size=1.2, alpha=0.3)
	p <- p + labs(x="year",y=paste("len is ", unique(mmdf$len)), title="bias in Ulength")
	p <- p + theme_bw(11) 
	p <- p + ylim(-0.5, 0.5)
	#p <- p + annotate("text" , x = 1, y = 0.4, label = paste("n = ",conv_n))
	p <- p + facet_wrap(~year)
	#p <- p + geom_text(data=mdf, aes(x=1.2, y=0.48, label=converge), parse=TRUE)
	print(p)
	setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	ggsave(paste("U_bias_len",unique(mmdf$len),sep=""), plot=p)

   }

	
	
}