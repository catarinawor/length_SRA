
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


plot_Sel <- function( M )
{
	cat("plot_Sel")

	n <- length( M )
	
	scn<-read_scnnames()


	conv_n<-numeric(length=length(scn))

	cio<-NULL

	cip<-NULL
	
	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){
			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1

			estyrs<-M[[i]]$OM$rep_yr:M[[i]]$OM$eyr
			
			names(M[[i]]$SArep)

			nlen<-length(M[[i]]$SArep$len)

			umaxes_est<-apply(M[[i]]$SArep$Ulength[,(nlen-5):nlen],1,mean)
			sels_est<-(M[[i]]$SArep$Ulength)/umaxes_est

			umaxes_om<-apply(M[[i]]$OM$Ulength[estyrs,(nlen-5):nlen],1,mean)
			sels_OM<-(M[[i]]$OM$Ulength[estyrs,])/umaxes_om

			selom <- data.frame(sel=c(sels_OM),len=rep(1:ncol(sels_OM),each=length(estyrs)),yr=rep(estyrs,ncol(sels_OM)),type="OM", scenario=scn[M[[i]]$OM$scnNumber], scnNumber=scn[M[[i]]$OM$scnNumber])
			selest <- data.frame(sel=c(sels_est),len=rep(1:ncol(sels_est),each=length(estyrs)),yr=rep(estyrs,ncol(sels_est)),type="EST",scenario=scn[M[[i]]$OM$scnNumber],scnNumber=scn[M[[i]]$OM$scnNumber])

			
			cio <- rbind(cio,selom)
			cip <- rbind(cip,selest)

		}
	}

	##df<-rbind(cio,cip)
	##p <- ggplot(df,aes(x=as.factor(len),y=sel,color=type,fill=type)) 
	##p <- p + geom_boxplot()
	##p <- p + facet_wrap(yr~scenario,scale="free")
	##p <- p + theme_bw(11)
	##p 

		omd<-NULL
		esd<-NULL


	for(sc in 1:length(scn)){
		for(ll in 1:ncol(sels_OM)){
			for(y in 1:length(estyrs)){

			cioy<-c(calc_quantile(cio$sel[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]]))
			cipy<-c(calc_quantile(cip$sel[cip$yr==estyrs[y]&cip$len==ll&cio$scenario==scn[sc]]))

			co<-data.frame(median=cioy[3],low=cioy[1],high=cioy[5] , ll=ll,year=estyrs[y], type="om", scenario=cio$scenario[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]])
			ce<-data.frame(median=cipy[3],low=cipy[1],high=cipy[5] , ll=ll,year=estyrs[y], type="est", scenario=cip$scenario[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]])
			
			omd<-rbind(omd,co)
			esd<-rbind(esd,ce)
			}
		}

	}


		df2<-rbind(omd,esd)
		

		
		for(sc in 1:length(scn)){
			df3<-df2[df2$scenario==scn[sc],]
			summary(df3)

			p2 <- ggplot(df3,aes(x=(ll),y=median,color=type,fill=type)) 
			p2 <- p2 + geom_line()
			p2 <- p2 + geom_ribbon(aes(ymax=high, ymin=low),alpha=0.2)
			p2 <- p2 + theme_bw(11)
			p2 <- p2 + facet_wrap(~year,scale="free")
			p2 <- p2 + ggtitle(paste(scn[sc]))
			print(p2)

			setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
			ggsave(paste("sel_",scn[sc],".pdf",sep=""), plot=p2)
		}

		
		

		#p <- p + labs(x="Year",y="Total Biomass")
		#p <- p + ylim(min(fdf$Low),max(fdf$High))
	
	
}


plot_Sel_pub <- function( M )
{
	cat("plot_Sel")

	n <- length( M )
	
	scn<-read_scnnames()


	conv_n<-numeric(length=length(scn))

	cio<-NULL

	cip<-NULL
	
	for(i in 1:n){

		if(M[[i]]$SApar$maxgrad<1.0e-04){
			conv_n[M[[i]]$OM$scnNumber] <-  conv_n[M[[i]]$OM$scnNumber] + 1

			estyrs<-M[[i]]$OM$rep_yr:M[[i]]$OM$eyr
			
			names(M[[i]]$SArep)

			umaxes_est<-apply(M[[i]]$SArep$Ulength,1,mean)
			sels_est<-(M[[i]]$SArep$Ulength)/umaxes_est

			umaxes_om<-apply(M[[i]]$OM$Ulength[estyrs,],1,mean)
			sels_OM<-(M[[i]]$OM$Ulength[estyrs,])/umaxes_om

			selom <- data.frame(sel=c(sels_OM),len=rep(1:ncol(sels_OM),each=length(estyrs)),yr=rep(estyrs,ncol(sels_OM)),type="OM", scenario=scn[M[[i]]$OM$scnNumber], scnNumber=scn[M[[i]]$OM$scnNumber])
			selest <- data.frame(sel=c(sels_est),len=rep(1:ncol(sels_est),each=length(estyrs)),yr=rep(estyrs,ncol(sels_est)),type="EST",scenario=scn[M[[i]]$OM$scnNumber],scnNumber=scn[M[[i]]$OM$scnNumber])

			
			cio <- rbind(cio,selom)
			cip <- rbind(cip,selest)

		}
	}

	##df<-rbind(cio,cip)
	##p <- ggplot(df,aes(x=as.factor(len),y=sel,color=type,fill=type)) 
	##p <- p + geom_boxplot()
	##p <- p + facet_wrap(yr~scenario,scale="free")
	##p <- p + theme_bw(11)
	##p 

		omd<-NULL
		esd<-NULL


	for(sc in 1:length(scn)){
		for(ll in 1:ncol(sels_OM)){
			for(y in 1:length(estyrs)){

			cioy<-c(calc_quantile(cio$sel[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]]))
			cipy<-c(calc_quantile(cip$sel[cip$yr==estyrs[y]&cip$len==ll&cio$scenario==scn[sc]]))

			co<-data.frame(median=cioy[3],low=cioy[1],high=cioy[5] , ll=ll,year=estyrs[y], type="om", scenario=cio$scenario[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]])
			ce<-data.frame(median=cipy[3],low=cipy[1],high=cipy[5] , ll=ll,year=estyrs[y], type="est", scenario=cip$scenario[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]])
			
			omd<-rbind(omd,co)
			esd<-rbind(esd,ce)
	}}

		}


		df2<-rbind(omd,esd)

		summary(df2)
		repyr<-c(21,30,40,50)
		

		df22<-df2[df2$year==21|df2$year==30|df2$year==40|df2$year==50,]		
		df22$year<-as.factor(df22$year)
		summary(df22)

		p2 <- ggplot(df22,aes(x=ll,y=median,color=type,fill=type)) 
			p2 <- p2 + geom_line()
			p2 <- p2 + geom_ribbon(aes(ymax=high, ymin=low),alpha=0.2)
			p2 <- p2 + facet_grid(scenario~year,scale="free", labeller = label_both)
			p2 <- p2 + labs(x="Length",y="Selectivity")
			p2 <- p2 + theme_bw(12) 
			p2 <- p2 + scale_colour_grey(start = 0.2, end = 0.8,labels = c("simulated", "estimated"))
			p2 <- p2 + scale_fill_grey(start = 0.2, end = 0.8,labels = c("simulated", "estimated"))
			p2 <- p2 + theme(axis.text = element_text(face="bold", size=12),
  			axis.title = element_text(face="bold", size=12),
  			strip.text = element_text(face="bold", size=15))
  			p2 <- p2 + guides(fill = guide_legend(title = NULL),color=guide_legend(title = NULL))
			print(p2)
		

			setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
			ggsave("sel_publ.pdf", plot=p2)
		}

		
		

		#p <- p + labs(x="Year",y="Total Biomass")
		#p <- p + ylim(min(fdf$Low),max(fdf$High))
	
	


