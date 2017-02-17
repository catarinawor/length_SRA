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

example_dir<-"/Users/catarinawor/Documents/length_SRA/examples/hake"


source("/Users/catarinawor/Documents/length_SRA/R/read.admb.R")


setwd(example_dir)
example<-read.rep("length_SRA.rep")

names(example)

sv<-logical(length=length(example$yr))
o<-1
for(i in 1:length(example$yr)){
	if(example$yr[i]==example$iyr[o]){
		sv[i]<-TRUE
		o<-o+1
	}
}


plot(example$iyr,example$survB)
lines(example$iyr,example$psurvB[sv])




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



