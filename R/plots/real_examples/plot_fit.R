# ==============================================
#Title:plot_fit
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot fit to real data examples 
# ==============================================


require(reshape2)
require(tidyr)
require(ggplot2)
library(xtable)


#development stage
source("/Users/catarinawor/Documents/length_SRA/R/read.admb.R")

example_dir<-"/Users/catarinawor/Documents/length_SRA/examples/hake"

setwd(example_dir)
hake<-read.rep("length_SRA.rep")
 names( hake)

example_dir<-"/Users/catarinawor/Documents/length_SRA/examples/jack_mackerel"
setwd(example_dir)
jm<-read.rep("length_SRA.rep")
 names(jm)

setwd("/Users/catarinawor/Documents/length_SRA/report")


#plot of fit to time series


hk_it<- data.frame(year=hake$iyr,species=rep("hake",length(hake$iyr)),
		observed=hake$survB, predicted= hake$q*hake$predSurvB)
jm_it<- data.frame(year=jm$iyr,species=rep("jack mackerel",length(jm$iyr)),
		observed=jm$survB, predicted= jm$q*jm$predSurvB)

it_df<-rbind(hk_it,jm_it)

setwd("/Users/catarinawor/Documents/length_SRA/R/plots/figs")


pit<-ggplot(it_df)
pit<-pit+geom_line(aes(x=year,y=predicted),size = 1.2)
pit<-pit+geom_point(aes(x=year,y=observed),shape = 21, stroke = 1.2)
pit<-pit+facet_wrap(~species,scales="free")
pit<-pit+ ylab("Index of abundance")
pit<-pit+ theme_bw(16)
pit
ggsave("It_fit.pdf", plot=p)


#pdf("fit_survey.pdf")
par(mfrow=c(1,2))
plot(hake$iyr,hake$survB, pch=16)
lines(hake$iyr,hake$q*hake$predSurvB,lwd=2)
plot(jm$iyr,jm$survB,pch=16)
lines(jm$iyr,jm$q*jm$predSurvB,lwd=2)
#dev.off()
hake$yr

#===========================================================
##plot selectivities
dim(hake$Ulength)


apply(hake$Ulength,1,mean)

hakeselec<-(hake$Ulength)/hake$maxUy
hake$len

df<-melt(hakeselec,variable.name=c("year", "length"))
summary(df)


df$value[df$value==0]<-NA

df$len<-hake$len[as.numeric(df$Var2)]
df$year<-hake$yr[as.numeric(df$Var1)]

df1<-df[df$len<63,]

p2 <- ggplot(df,aes(x=len,y=value)) 
			p2 <- p2 + geom_line()
			p2 <- p2 + theme_bw(11)
			p2 <- p2 + facet_wrap(~year,scales="free") 
			print(p2)

			p2 <- p2 + coord_cartesian(ylim=c(0,7))
			print(p2)

#ggsave("selec_hake.pdf", plot=p2)
#jack jack_mackerel




jmselec<-(jm$Ulength)/jm$avgUy

df_jm<-melt(jmselec)
summary(df_jm)
df_jm$value[df_jm$value==0]<-NA

df_jm$len<-jm$len[as.numeric(df_jm$Var2)]
df_jm$year<-jm$yr[as.numeric(df_jm$Var1)]

#df1<-df[df$len<63,]

p_jm <- ggplot(df_jm,aes(x=len,y=value)) 
			p_jm <- p_jm + geom_line()
			p_jm <- p_jm + theme_bw(11)
			p_jm <- p_jm + facet_wrap(~year,scales="free")
			p_jm <- p_jm + coord_cartesian(ylim=c(0,7))
			
			print(p_jm)
#ggsave("selec_jm.pdf", plot=p_jm)

#===========================================================
##table of estimated parameters


params<-data.frame(species=c("hake","jack mackerel"),Ro=c(hake$Ro,jm$Ro),Rinit=c(hake$Rinit,jm$Rinit),
	reck=c(hake$reck,jm$reck),q=c(hake$q,jm$q))
print(xtable(params), type="latex", file="/Users/catarinawor/Documents/length_SRA/report/paramstab.tex")



#===========================================================
##plot selectivities

hk_msy<-data.frame(year=rep(hake$yr,2),value=c(hake$umsy,hake$msy),variable=rep(c("Umsy","msy"),each=length(hake$msy)),species=rep("hake",length(hake$msy)*2))
jm_msy<-data.frame(year=rep(jm$yr,2),value=c(jm$umsy,jm$msy),variable=rep(c("Umsy","msy"),each=length(jm$msy)), species=rep("jack mackerel",length(jm$msy)*2))

msy_df<-rbind(hk_msy,jm_msy)



pm<-ggplot(hk_msy)
pm<-pm+geom_line(aes(x=year,y=value))
pm<-pm+geom_point(aes(x=year,y=value))
pm<-pm+facet_wrap(~ variable, scales="free")
pm



pm<-ggplot(msy_df)
pm<-pm+geom_line(aes(x=year,y=value))
pm<-pm+geom_point(aes(x=year,y=value))
pm<-pm+facet_wrap(species~ variable, scales="free")
pm

pm<-pm+facet_wrap(~variable+species, scales="free")

?facet_wrap


	msy=c(hake$msy,jm$msy),umsy=c(hake$umsy,jm$umsy))
#pdf("msy_umsy_t.pdf")
par(mfrow=c(2,2))
plot(hake$yr,hake$umsy, type="b")
plot(hake$yr,hake$msy, type="b")
plot(jm$yr,jm$umsy, type="b")
plot(jm$yr,jm$msy, type="b")
#dev.off()

#============================
#random stuff
melt
length(hake$avgUy)

plot(example$sbt)

example$Rinit
example$Ro

ndf<-melt(example$Nat)


ndf$year<-example$yr[ndf$Var1]
ndf$age<-rep(example$age,each=length(example$yr))

summary(ndf)

p <- ggplot(ndf,aes((year),age,size=value))
	p <- p + geom_point(aes(colour=as.factor(sign(value))),alpha=0.75) 
	#p <- p + geom_point(alpha=0.75,) 

	p <- p + scale_size_area(max_size=5)
	p <- p + labs(x="Year",y="length",size="bias")
	p <- p + facet_wrap(~scenario,scales="free")
	#p <- p + scale_colour_discrete(guide="none")
	p <- p + theme_bw(11)
	print(p)


#====================================================================
#random stuff
names(hake)
hake$yield/hake$msy

apply(hake$Uage,1,max)/hake$umsy

hk_msy<-data.frame(year=rep(hake$yr,2),value=c(hake$umsy,hake$msy,hake$umsy/hake$umsy),variable=rep(c("Umsy","msy"),each=length(hake$msy)),species=rep("hake",length(hake$msy)*2))
jm_msy<-data.frame(year=rep(jm$yr,2),value=c(jm$umsy,jm$msy),variable=rep(c("Umsy","msy"),each=length(jm$msy)), species=rep("jack mackerel",length(jm$msy)*2))


apply(hake$Uage,1,max)





