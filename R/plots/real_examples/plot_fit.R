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
library(cowplot)


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


?merge

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
#ggsave("It_fit.pdf", plot=p)


#pdf("fit_survey.pdf")
par(mfrow=c(1,2))
plot(hake$iyr,hake$survB, pch=16)
lines(hake$iyr,hake$q*hake$predSurvB,lwd=2)
plot(jm$iyr,jm$survB,pch=16)
lines(jm$iyr,jm$q*jm$predSurvB,lwd=2)
#dev.off()
hake$yr


names(hake)

aly<-data.frame(year=hake$yr,predicted=hake$psurvB*hake$q)
jsu<-data.frame(year=hake$iyr,observed=hake$survB)

hk_all<-merge(aly,jsu, all.x=T)
hk_all<-hk_all[hk_all$year>1994,]



aly<-data.frame(year=hake$yr,predicted=hake$psurvB*hake$q)
jsu<-data.frame(year=hake$iyr,observed=hake$survB)

hk_all<-merge(aly,jsu, all.x=T)
hk_all<-hk_all[hk_all$year>1994,]
hk_all$species<-"Pacific hake"

alyj<-data.frame(year=jm$yr,predicted=jm$psurvB*jm$q)
jsuj<-data.frame(year=jm$iyr,observed=jm$survB)

jm_all<-merge(alyj,jsuj, all.x=T)
jm_all<-jm_all[jm_all$year>1985,]
jm_all$species<-"jack mackerel"

it_all<-rbind(hk_all,jm_all)

pit<-ggplot(it_all)
pit<-pit+geom_line(aes(x=year,y=predicted),size = 1.2)
pit<-pit+geom_point(aes(x=year,y=observed),shape = 21, stroke = 1.2)
pit<-pit+facet_wrap(~species,scales="free")
pit<-pit+ ylab("Index of abundance")
pit<-pit+ theme_bw(16)
pit

hake$predSurvB



#===========================================================

jm_bt<- data.frame(year=1:length(jm$bt),species=rep("jack mackerel",length(jm$bt)),
		biomass=jm$bt)

pbt<-ggplot(jm_bt)
pbt<-pbt+geom_line(aes(x=year,y=biomass),size = 1.2)
pbt
#===========================================================
##plot selectivities
names(hake)

dim(hake$Ulength)


apply(hake$Ulength,1,mean)

hakeselec<-(hake$Ulength)/hake$avgUy
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
			#p2 <- p2 + coord_cartesian(xlim=c(35,55))
			p2 <- p2 + facet_wrap(~year,scales="free",dir="v") 
			print(p2)

			#p2 <- p2 + coord_cartesian(ylim=c(0,7))
			#print(p2)

#ggsave("selec_hake.pdf", plot=p2)
#jack jack_mackerel




jmselec<-(jm$Ulength)/jm$maxUy

df_jm<-melt(jmselec)
summary(df_jm)
df_jm$value[df_jm$value==0]<-NA

df_jm$len<-jm$len[as.numeric(df_jm$Var2)]
df_jm$year<-jm$yr[as.numeric(df_jm$Var1)]

#df1<-df[df$len<63,]

p_jm <- ggplot(df_jm,aes(x=len,y=value)) 
			p_jm <- p_jm + geom_line()
			p_jm <- p_jm + theme_bw(11)
			p_jm <- p_jm + facet_wrap(~year,scales="free",dir="v")
			#p_jm <- p_jm + coord_cartesian(ylim=c(0,7))
			
			print(p_jm)
#ggsave("selec_jm.pdf", plot=p_jm)

#===========================================================
##table of estimated parameters


params<-data.frame(species=c("hake","jack mackerel"),Ro=c(hake$Ro,jm$Ro),Rinit=c(hake$Rinit,jm$Rinit),
	reck=c(hake$reck,jm$reck),q=c(hake$q,jm$q))
print(xtable(params), type="latex", file="/Users/catarinawor/Documents/length_SRA/report/paramstab.tex")



#===========================================================
##MSY



names(hake)
hake$yield/hake$msy  
hake$maxUy/hake$umsy

jm$yield/jm$msy

jm$maxUy/jm$umsy


hk_msy<-data.frame(year=rep(hake$yr,2),value=c(hake$umsy,hake$msy),variable=rep(c("Umsy","msy"),each=length(hake$msy)),species=rep("Pacific hake",length(hake$msy)*2))
hk_Umsy<-data.frame(year=rep(hake$yr,2),value=c(hake$umsy),variable=rep(c("Umsy","msy"),each=length(hake$msy)),species=rep("Pacific hake",length(hake$msy)))

jm_msy<-data.frame(year=rep(jm$yr,2),value=c(jm$umsy,jm$msy),variable=rep(c("Umsy","msy"),each=length(jm$msy)), species=rep("jack mackerel",length(jm$msy)*2))
jm_Umsy<-data.frame(year=rep(jm$yr,2),value=c(jm$umsy,jm$msy),variable=rep(c("Umsy","msy"),each=length(jm$msy)), species=rep("jack mackerel",length(jm$msy)*2))


msyall<-rbind(hk_msy,jm_msy)

setwd("/Users/catarinawor/Documents/length_SRA/R/plots/figs")


pit<-ggplot(it_all)
pit<-pit+geom_line(aes(x=year,y=predicted),size = 1.2)
pit<-pit+geom_point(aes(x=year,y=observed),shape = 21, stroke = 1.2)
pit<-pit+facet_wrap(~species,scales="free")
pit<-pit+ ylab("Index of abundance")
pit<-pit+ theme_bw(16) + xlab(" ")
pit <- pit +coord_cartesian(xlim = c(1975,2013))
pit



pit<-ggplot(it_all)
pit<-pit+geom_line(aes(x=year,y=predicted),size = 1.2)
pit<-pit+geom_point(aes(x=year,y=observed),shape = 21, stroke = 1.2)
pit<-pit+facet_wrap(~species,scales="free")
pit<-pit+ ylab("Index of abundance")
pit<-pit+ theme_bw(16) + xlab(" ")
pit <- pit +coord_cartesian(xlim = c(1975,2013))
pit





pm<-ggplot(msyall[msyall$variable=="Umsy",])
pm<-pm+geom_line(aes(x=year,y=value))
pm<-pm+geom_point(aes(x=year,y=value))
pm<-pm+  ylab("MSY")
pm<-pm+  xlab(" ") + theme_bw(16)
pm<-pm+facet_wrap(~ species, scales="free")
pm

umsyall<-msyall[msyall$variable=="Umsy",]

pu<-ggplot(umsyall)
pu<-pu+geom_line(aes(x=year,y=value))
pu<-pu+geom_point(aes(x=year,y=value))
pu<-pu+  ylab(expression(paste("U"["MSY"])))
pu<-pu+  xlab("Year ")+ theme_bw(16)
pu<-pu+facet_wrap(~ species, scales="free")
pu


plot_grid(pith,pitj,pm,pu, ncol=1)

theme_new <- theme_set(theme_bw(16) )
theme_new <- theme_update(plot.title = element_text(hjust = 0.5,face="bold"))


pith<-ggplot(hk_it)
pith<-pith+geom_line(aes(x=year,y=predicted),size = 1.2)
pith<-pith+geom_point(aes(x=year,y=observed),shape = 21, stroke = 1.2)
pith<-pith+ ylab("Index of abundance")
pith<-pith+ theme_new() + xlab(" ") 
pith<-pith+ ggtitle("Pacific hake")
pith<-pith +coord_cartesian(xlim = c(1975,2013))
pith

pitj<-ggplot(jm_it)
pitj<-pitj+geom_line(aes(x=year,y=predicted),size = 1.2)
pitj<-pitj+geom_point(aes(x=year,y=observed),shape = 21, stroke = 1.2)
pitj<-pitj+ ylab(" ") 
pitj<-pitj+ ggtitle("jack mackerel")
pitj<-pitj+ theme_new()  + xlab(" ")
pitj<-pitj +coord_cartesian(xlim = c(1980,2013))
pitj


msyhk<-hk_msy[hk_msy$variable=="msy",]
umsyhk<-hk_msy[hk_msy$variable=="Umsy",]

msyjm<-jm_msy[jm_msy$variable=="msy",]
umsyjm<-jm_msy[jm_msy$variable=="Umsy",]


pmh<-ggplot(msyhk)
pmh<-pmh+geom_line(aes(x=year,y=value))
pmh<-pmh+geom_point(aes(x=year,y=value))
pmh<-pmh+  ylab("MSY") 
pmh<-pmh+  xlab(" ") + theme_new
pmh


puh<-ggplot(umsyhk)
puh<-puh+geom_line(aes(x=year,y=value))
puh<-puh+geom_point(aes(x=year,y=value))
puh<-puh+  ylab(expression(paste("U"["MSY"])))
puh<-puh+  xlab("Year ")+ theme_new
puh


pmj<-ggplot(msyjm)
pmj<-pmj+geom_line(aes(x=year,y=value))
pmj<-pmj+geom_point(aes(x=year,y=value))
pmj<-pmj+  ylab(" ")
pmj<-pmj+  xlab(" ") + theme_new
pmj


puj<-ggplot(umsyjm)
puj<-puj+geom_line(aes(x=year,y=value))
puj<-puj+geom_point(aes(x=year,y=value))
puj<-puj+  ylab(" ")
puj<-puj+  xlab("Year ")+ theme_new
puj

setwd("/Users/catarinawor/Documents/length_SRA/R/plots/figs")

plot_grid(pith,pitj,pmh,pmj,puh,puj, ncol=2, align = 'h')
ggsave("real_examples.pdf")


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


apply(hake$Uage,1,max)



#==========================================
#selectivities at age

##plot selectivities
names(hake)

dim(hake$Uage)


apply(hake$Uage,1,mean)

hakeselec<-(hake$Uage)/apply(hake$Uage,1,max)

df<-melt(hakeselec,variable.name=c("year", "age"))
summary(df)



df$age<-hake$age[as.numeric(df$Var2)]
df$year<-hake$yr[as.numeric(df$Var1)]


p2 <- ggplot(df,aes(x=age,y=value)) 
			p2 <- p2 + geom_line()
			p2 <- p2 + theme_bw(11)
			#p2 <- p2 + coord_cartesian(xlim=c(35,55))
			p2 <- p2 + facet_wrap(~year,scales="free",dir="v") 
			print(p2)

			#p2 <- p2 + coord_cartesian(ylim=c(0,7))
			#print(p2)

#ggsave("selec_hake.pdf", plot=p2)
#jack jack_mackerel




jmselec<-(jm$Uage)/apply(jm$Uage,1,mean)

df_jm<-melt(jmselec)
summary(df_jm)

df_jm$age<-jm$age[as.numeric(df_jm$Var2)]
df_jm$year<-jm$yr[as.numeric(df_jm$Var1)]

#df1<-df[df$len<63,]

p_jm <- ggplot(df_jm,aes(x=age,y=value)) 
			p_jm <- p_jm + geom_line()
			p_jm <- p_jm + theme_bw(11)
			p_jm <- p_jm + facet_wrap(~year,scales="free")
			#p_jm <- p_jm + coord_cartesian(ylim=c(0,7))
			
			print(p_jm)
#ggsave("selec_jm.pdf", plot=p_jm)

