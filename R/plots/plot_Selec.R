
# ==============================================
#Title:plot_Udevs
#Author: Catarina Wor
#date: November 2016
#Function to plot deviation from true U
# code adapted from iSCAM stuff (Martell et al)
# ==============================================
#for testing delete, when done



source("/Users/catarinawor/Documents/Lagrangian/R/read_mse/Rplots/calc_quantile.R")
#M<-SIMSdat

require(reshape2)
require(tidyr)
require(ggplot2)

source("/Users/catarinawor/Documents/Length_SRA/R/plots/readscn.R")


plot_Sel <- function( M, sv=FALSE, nome=""){
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

			umaxes_est<-apply(M[[i]]$SArep$Ulength,1,mean)
			sels_est<-((M[[i]]$SArep$Ulength)/umaxes_est)


			umaxes_om<-apply(M[[i]]$OM$Ulength[estyrs,],1,mean)
			sels_OM<-((M[[i]]$OM$Ulength[estyrs,])/umaxes_om)

			selom <- data.frame(sel=c(sels_OM),len=rep(1:ncol(sels_OM),each=length(estyrs)),yr=rep(estyrs,ncol(sels_OM)),type="OM", scenario=scn[M[[i]]$OM$scnNumber], scnNumber=scn[M[[i]]$OM$scnNumber])
			selest <- data.frame(sel=c(sels_est),len=rep(1:ncol(sels_est),each=length(estyrs)),yr=rep(estyrs,ncol(sels_est)),type="EST",scenario=scn[M[[i]]$OM$scnNumber],scnNumber=scn[M[[i]]$OM$scnNumber])

			
			cio <- rbind(cio,selom)
			cip <- rbind(cip,selest)

		}
	}

	#summary(cio)
	#df<-rbind(cio,cip)
	#for(sc in 1:length(scn)){
	#	df3<-df[df$scenario==scn[sc],]
	#	cio3<-cio[cio$scenario==scn[sc],]
	#	
	#	limo<-	cio3[1:(length(unique(cio3$yr))*length(unique(cio3$len))),]
	#
	#	p <- ggplot(df3,aes(x=as.factor(len),y=sel,color=type)) 
	#	p <- p + geom_boxplot(outlier.shape = NA)
	#	p <- p + facet_wrap(~yr,scale="free")
	#	p <- p + geom_line(data=limo,aes(y=sel,x=len), color="black")
	#	p <- p + ggtitle(paste(scn[sc]))
	#	p <- p + labs(x="Length",y="U")
	#	p <- p + theme_bw(12) 
	#	p <- p + scale_colour_grey(start = 0.1, end = 0.6,labels = c("simulated", "estimated"))
	#			
	#	print(p) 
	#
	#	if(sv==TRUE){
	#		setwd("/Users/catarinawor/Documents/Length_SRA/report")
	#		ggsave(paste(nome,"sel_",scn[sc],".pdf",sep=""), plot=p)
	#		
	#	}
	#
	#}
		
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
		
		summary(df2)
		
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
			ggsave(paste(nome,"_sel_",scn[sc],".pdf",sep=""), plot=p2)
		}

		
		

		#p <- p + labs(x="Year",y="Total Biomass")
		#p <- p + ylim(min(fdf$Low),max(fdf$High))
	
	
}


plot_Sel_pub <- function( M ){
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

			umaxes_est<-apply(M[[i]]$SArep$Ulength,1,mean)
			sels_est<-(M[[i]]$SArep$Ulength)/umaxes_est

			umaxes_om<-apply(M[[i]]$OM$Ulength[estyrs,],1,mean)
			sels_OM<-(M[[i]]$OM$Ulength[estyrs,])/umaxes_om

			selom <- data.frame(sel=c(sels_OM),len=rep(M[[i]][[1]]$len,each=length(estyrs)),yr=rep(estyrs,ncol(sels_OM)),type="OM", scenario=scn[M[[i]]$OM$scnNumber], scnNumber=scn[M[[i]]$OM$scnNumber], seed=i)
			selest <- data.frame(sel=c(sels_est),len=rep(M[[i]][[1]]$len,each=length(estyrs)),yr=rep(estyrs,ncol(sels_est)),type="EST",scenario=scn[M[[i]]$OM$scnNumber],scnNumber=scn[M[[i]]$OM$scnNumber],seed=i)

			
			cio <- rbind(cio,selom)
			cip <- rbind(cip,selest)

		}
	}

		omd<-NULL
		esd<-NULL

		mylen=unique(cio$len)

	for(sc in 1:length(scn)){
		for(ll in 1:ncol(sels_OM)){
			for(y in 1:length(estyrs)){

				summary(cio)

			cioy<-c(calc_quantile(cio$sel[cio$yr==estyrs[y]&cio$len==mylen[ll]&cio$scenario==scn[sc]]))
			cipy<-c(calc_quantile(cip$sel[cip$yr==estyrs[y]&cip$len==mylen[ll]&cio$scenario==scn[sc]]))

			obs1<-cip$sel[cio$yr==estyrs[y]&cio$len==mylen[ll]&cio$scenario==scn[sc]][12]
			obs2<-cip$sel[cio$yr==estyrs[y]&cio$len==mylen[ll]&cio$scenario==scn[sc]][90]
			obs3<-cip$sel[cio$yr==estyrs[y]&cio$len==mylen[ll]&cio$scenario==scn[sc]][180]

			co<-data.frame(median=rep(cioy[3],3),low=rep(cioy[1],3),high=rep(cioy[5],3),
						 obs=c(obs1,obs2,obs3), obsn=1:3, ll=rep(mylen[ll],3),year=rep(estyrs[y],3),
						  type="om", scenario=rep(unique(cio$scenario[cio$yr==estyrs[y]&cio$len==mylen[ll]&cio$scenario==scn[sc]]),3))
			

			ce<-data.frame(median=rep(cipy[3],3),low=rep(cipy[1],3),high=rep(cipy[5],3),
						 obs=c(obs1,obs2,obs3),obsn=1:3, ll=rep(mylen[ll],3),year=rep(estyrs[y],3),
						  type="est", scenario=rep(unique(cip$scenario[cio$yr==estyrs[y]&cio$len==mylen[ll]&cio$scenario==scn[sc]]),3))
			
			omd<-rbind(omd,co)
			esd<-rbind(esd,ce)
			}
		}
	}


		df2<-rbind(omd,esd)

		summary(df2)
		repyr<-c(21,30,40,50)
		#repyr<-c(25,35,45,49)
		

		df22<-df2[df2$year==21|df2$year==30|df2$year==40|df2$year==50,]		
		df22$year<-as.factor(df22$year)
		summary(df22)

		p2 <- ggplot(df22,aes(x=ll,y=median,color=type,fill=type)) 
			p2 <- p2 + geom_line()
			p2 <- p2 + geom_ribbon(aes(ymax=high, ymin=low),alpha=0.2)
			p2 <- p2 + geom_line(aes(x=ll,y=obs,linetype=as.factor(obsn)),show.legend = FALSE)
			p2 <- p2 + facet_grid(scenario~year, labeller = label_both)
			p2 <- p2 + labs(x="Length",y="Selectivity")
			p2 <- p2 + theme_bw(12) + scale_linetype_manual(values=c("dashed","twodash", "dotted"))
			p2 <- p2 + scale_colour_grey(start = 0.1, end = 0.6,labels = c("simulated", "estimated"))
			p2 <- p2 + scale_fill_grey(start = 0.1, end = 0.6,labels = c("simulated", "estimated"))
			p2 <- p2 + theme(axis.text = element_text(face="bold", size=12),
  			axis.title = element_text(face="bold", size=12),
  			strip.text = element_text(face="bold", size=15))
  			p2 <- p2 + guides(fill = guide_legend(title = NULL),color=guide_legend(title = NULL))
			print(p2)



		
		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
			ggsave("sel_publ.pdf", plot=p2, width = 14, height = 10)
		}
	

		
}		

		#p <- p + labs(x="Year",y="Total Biomass")
		#p <- p + ylim(min(fdf$Low),max(fdf$High))
	


plot_Ult_Clt <- function( M ,sv=F, nome="")
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

			umaxes_est<-apply(M[[i]]$SArep$Ulength,1,mean)
			sels_est<-(M[[i]]$SArep$Ulength)/umaxes_est

			cmaxes_est<-apply(M[[i]]$SArep$Clt,1,mean)
			cls_est<-(M[[i]]$SArep$Clt)/cmaxes_est

			#umaxes_om<-apply(M[[i]]$OM$Ulength[estyrs,],1,mean)
			#sels_OM<-(M[[i]]$OM$Ulength[estyrs,])/umaxes_om

			#selom <- data.frame(sel=c(sels_OM),len=rep(1:ncol(sels_OM),each=length(estyrs)),yr=rep(estyrs,ncol(sels_OM)),type="OM", scenario=scn[M[[i]]$OM$scnNumber], scnNumber=scn[M[[i]]$OM$scnNumber])
			clest <- data.frame(sel=c(cls_est),len=rep(1:ncol(cls_est),each=length(estyrs)),yr=rep(estyrs,ncol(cls_est)),type="Cls", scenario=scn[M[[i]]$OM$scnNumber], scnNumber=scn[M[[i]]$OM$scnNumber])
			selest <- data.frame(sel=c(sels_est),len=rep(1:ncol(sels_est),each=length(estyrs)),yr=rep(estyrs,ncol(sels_est)),type="Uls",scenario=scn[M[[i]]$OM$scnNumber],scnNumber=scn[M[[i]]$OM$scnNumber])

			
			cio <- rbind(cio,clest)
			cip <- rbind(cip,selest)

		}
	}

	#summary(cio)
	df<-rbind(cio,cip)
	for(sc in 1:length(scn)){
		df3<-df[df$scenario==scn[sc],]
		cio3<-cio[cio$scenario==scn[sc],]
	
		#limo<-	cio3[1:(length(unique(cio3$yr))*length(unique(cio3$len))),]
	
		p <- ggplot(df3,aes(x=as.factor(len),y=sel,color=type)) 
		p <- p + geom_boxplot(outlier.shape = NA)
		p <- p + facet_wrap(~yr,scale="free")
		p <- p + theme_bw(11)
		#p <- p + geom_line(data=limo,aes(y=sel,x=len), color="black")
		p <- p + ggtitle(paste(scn[sc]))
		p <- p + labs(x="Length",y=" ")
		p <- p + theme_bw(12) 
		p <- p + scale_colour_grey(start = 0.1, end = 0.6)
				
		print(p) 

		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Length_SRA/report")
			ggsave(paste(nome,"clt_ult_",scn[sc],".pdf",sep=""), plot=p)
		}
	
	}
		


	#for(sc in 1:length(scn)){
	#	for(ll in 1:ncol(sels_OM)){
	#		for(y in 1:length(estyrs)){

	#		cioy<-c(calc_quantile(cio$sel[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]]))
	#		cipy<-c(calc_quantile(cip$sel[cip$yr==estyrs[y]&cip$len==ll&cio$scenario==scn[sc]]))

	#		co<-data.frame(median=cioy[3],low=cioy[1],high=cioy[5] , ll=ll,year=estyrs[y], type="om", scenario=cio$scenario[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]])
	#		ce<-data.frame(median=cipy[3],low=cipy[1],high=cipy[5] , ll=ll,year=estyrs[y], type="est", scenario=cip$scenario[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]])
			
	#		omd<-rbind(omd,co)
	#		esd<-rbind(esd,ce)
	#		}
	#	}

	#}


	#	df2<-rbind(omd,esd)
		

		
	#	for(sc in 1:length(scn)){
	#		df3<-df2[df2$scenario==scn[sc],]
	
	#		summary(df3)

	#		p2 <- ggplot(df3,aes(x=(ll),y=median,color=type,fill=type)) 
	#		p2 <- p2 + geom_line()
	#		p2 <- p2 + geom_ribbon(aes(ymax=high, ymin=low),alpha=0.2)
	#		p2 <- p2 + theme_bw(11)
	#		p2 <- p2 + facet_wrap(~year,scale="free")
	#		p2 <- p2 + ggtitle(paste(scn[sc]))
	#		print(p2)

	#		setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	#		ggsave(paste("sel_",scn[sc],".pdf",sep=""), plot=p2)
	#	}

		
		

		#p <- p + labs(x="Year",y="Total Biomass")
		#p <- p + ylim(min(fdf$Low),max(fdf$High))
	
	
}






plot_USel <- function( M )
{
	cat("plot_USel")

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
			sels_est<-(M[[i]]$SArep$Ulength)#/umaxes_est

			umaxes_om<-apply(M[[i]]$OM$Ulength[estyrs,(nlen-5):nlen],1,mean)
			sels_OM<-(M[[i]]$OM$Ulength[estyrs,])#/umaxes_om

			selom <- data.frame(sel=c(sels_OM),len=rep(M[[i]][[1]]$len,each=length(estyrs)),yr=rep(estyrs,ncol(sels_OM)),type="OM", scenario=scn[M[[i]]$OM$scnNumber], scnNumber=scn[M[[i]]$OM$scnNumber])
			selest <- data.frame(sel=c(sels_est),len=rep(M[[i]][[1]]$len,each=length(estyrs)),yr=rep(estyrs,ncol(sels_est)),type="EST",scenario=scn[M[[i]]$OM$scnNumber],scnNumber=scn[M[[i]]$OM$scnNumber])

			
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

			#setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
			#ggsave(paste( scn[sc],"_sel",".pdf",sep=""), plot=p2)
		}

		
		

		#p <- p + labs(x="Year",y="Total Biomass")
		#p <- p + ylim(min(fdf$Low),max(fdf$High))
	
	
}





plot_Sel_biasLinf <- function( M, sv=FALSE, nome=""){
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

			umaxes_est<-apply(M[[i]]$SArep$Ulength,1,mean)
			sels_est<-(M[[i]]$SArep$Ulength)#/umaxes_est

			umaxes_om<-apply(M[[i]]$OM$Ulength[estyrs,],1,mean)
			sels_OM<-(M[[i]]$OM$Ulength[estyrs,])#/umaxes_om
			ncol(sels_est)

			selom <- data.frame(sel=c(sels_OM),len=rep(M[[i]]$SArep$len,each=length(estyrs)),yr=rep(estyrs,ncol(sels_OM)),type="OM", scenario=scn[M[[i]]$OM$scnNumber], scnNumber=scn[M[[i]]$OM$scnNumber])
			selest <- data.frame(sel=c(sels_est),len=rep(M[[i]]$SArep$len,each=length(estyrs)),yr=rep(estyrs,ncol(sels_est)),type="EST",scenario=scn[M[[i]]$OM$scnNumber],scnNumber=scn[M[[i]]$OM$scnNumber])

			
			cio <- rbind(cio,selom)
			cip <- rbind(cip,selest)

		}
	}

	#summary(cio_A)
	df<-rbind(cio,cip)

	head(df)

	df_A<-df[df$yr>46,]
	cio_A<-cio[cio$yr>46,]

	df_A$len<-as.factor(df_A$len)


	limo<-aggregate(cio_A$sel,by=list(cio_A$len,cio_A$scenario,cio_A$yr),median)
	limo$scenario<-	limo$Group.2
	limo$yr<-	limo$Group.3
	limo$len<-	as.numeric(as.factor(limo$Group.1))
	limo$sel<-	limo$x

	#limo$len[limo$scenario=="plus10"]<-	limo$len[limo$scenario=="plus10"]
	

	summary(df_A[df_A$scenario=="minus10",])

	p <- ggplot(df_A,aes(x=(len),y=sel,color=type)) 
		p <- p + geom_boxplot(outlier.shape = NA)
		p <- p + coord_cartesian(ylim=c(0,0.4))
		p <- p + facet_wrap(scenario~yr)
		p <- p + geom_line(data=limo,aes(y=sel,x=((len)), color="black"))
		p <- p + labs(x="Length",y=expression("U"["length"]))
		p <- p + theme_bw(16) 
		p <- p +scale_x_discrete(breaks=unique(df_A$len)[seq(1,33,2)],
			labels=unique(df_A$len)[seq(1,33,2)])
		p <- p + scale_colour_manual(values = c("black", "gray70","white"),labels = c("simulated","estimated",""))
				
		print(p) 



		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
			ggsave(paste(nome,"sel_ALinf.pdf",sep=""), plot=p)
			
		}
	
	}
		




	#for(sc in 1:length(scn)){
	#	for(ll in 1:ncol(sels_OM)){
	#		for(y in 1:length(estyrs)){

	#		cioy<-c(calc_quantile(cio$sel[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]]))
	#		cipy<-c(calc_quantile(cip$sel[cip$yr==estyrs[y]&cip$len==ll&cio$scenario==scn[sc]]))

	#		co<-data.frame(median=cioy[3],low=cioy[1],high=cioy[5] , ll=ll,year=estyrs[y], type="om", scenario=cio$scenario[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]])
	#		ce<-data.frame(median=cipy[3],low=cipy[1],high=cipy[5] , ll=ll,year=estyrs[y], type="est", scenario=cip$scenario[cio$yr==estyrs[y]&cio$len==ll&cio$scenario==scn[sc]])
			
	#		omd<-rbind(omd,co)
	#		esd<-rbind(esd,ce)
	#		}
	#	}

	#}


	#	df2<-rbind(omd,esd)
		

		
	#	for(sc in 1:length(scn)){
	#		df3<-df2[df2$scenario==scn[sc],]
	
	#		summary(df3)

	#		p2 <- ggplot(df3,aes(x=(ll),y=median,color=type,fill=type)) 
	#		p2 <- p2 + geom_line()
	#		p2 <- p2 + geom_ribbon(aes(ymax=high, ymin=low),alpha=0.2)
	#		p2 <- p2 + theme_bw(11)
	#		p2 <- p2 + facet_wrap(~year,scale="free")
	#		p2 <- p2 + ggtitle(paste(scn[sc]))
	#		print(p2)

	#		setwd("/Users/catarinawor/Documents/Length_SRA/R/plots/figs")
	#		ggsave(paste("sel_",scn[sc],".pdf",sep=""), plot=p2)
	#	}

		
		

		#p <- p + labs(x="Year",y="Total Biomass")
		#p <- p + ylim(min(fdf$Low),max(fdf$High))
	
	

	