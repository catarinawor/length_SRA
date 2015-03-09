#====================================================================================
#length based models - simulation of data
#June 2012

#====================================================================================

#run simulations!

setwd("/Users/catarinawor/Documents/length_SRA/simulation-evaluation")
#source("read.admb.R")

run.Simulation=function(N=1)
{
	for(i in 1:N){
	arg = paste("./sra_sim") 
	system(arg)
    #print(arg)
    }
}

run.Simulation()

SRA=read.table("SRA_bias.txt")

SRApar=read.table("SRA_param.txt")




lab=c("convergence","seed",rep("Recbias",38+12-1),"biasRo","biaskappa",rep("biasUlength_fin",77))
labpar=c("convergence","seed",rep("E_Rec",(38+12-1)),rep("T_Rec",(38+12-1)),"E_Ro","T_Ro","E_kappa",
	"T_kappa",rep("E_Ulength_fin",77),rep("T_Ulength_fin",77),rep("E_Umax",38),rep("T_Umax",38))


convSRA<-SRA[SRA[,1]<1e-1,]
convSRApar<-SRApar[SRApar[,1]<1e-1,]

dim(convSRA)
dim(convSRApar)



par(mfrow=c(2,2))
boxplot(convSRA[,lab=="Recbias"], main="RecDev")
abline(h=0,lwd=3,col="red")
boxplot(convSRA[,lab=="biasRo"], main="Ro")
abline(h=0,lwd=3,col="red")
boxplot(convSRA[,lab=="biaskappa"], main="kappa")
abline(h=0,lwd=3,col="red")
boxplot(convSRA[,lab=="biasUlength_fin"], main="Ul(eyr)")
abline(h=0,lwd=3,col="red")



par(mfrow=c(2,2))
boxplot(convSRA[,lab=="Recbias"]*100, main="RecDev",ylim=c(-100,1000))
abline(h=0,lwd=3,col="red")
boxplot(convSRA[,lab=="biasRo"]*100, main="Ro",ylim=c(-100,1200))
abline(h=0,lwd=3,col="red")
boxplot(convSRA[,lab=="biaskappa"]*100, main="kappa",ylim=c(-1000,40000000000000))
abline(h=0,lwd=3,col="red")
boxplot(convSRA[,lab=="biasUlength_fin"]*100, main="Ul(eyr)",ylim=c(-100,5e30))
abline(h=0,lwd=3,col="red")




convSRApar[labpar=="T_kappa"]
convSRApar[labpar=="E_kappa"]




convSRApar[labpar=="T_Ro"]
convSRApar[labpar=="E_Ro"]

#plot(1:50,as.vector(convSRApar[labpar=="T_Ulength_fin"]))

plot(1:34,as.vector(convSRApar[labpar=="T_Umax"][1,]),ylim=c(0,2))
lines(1:34,as.vector(convSRApar[labpar=="E_Umax"][1,]),lwd=2)


