# ==============================================
#Title:plot_fit
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot fit to real data examples 
# ==============================================


require(reshape2)
require(tidyr)
require(ggplot2)


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


par(mfrow=c(1,2))
plot(hake$iyr,hake$survB, pch=16)
lines(hake$iyr,hake$q*hake$predSurvB,lwd=2)
plot(jm$iyr,jm$survB,pch=16)
lines(jm$iyr,jm$q*jm$predSurvB,lwd=2)




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



