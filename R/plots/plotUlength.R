# ==============================================
#Title:plot_mainParams
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot deviation from true U
# code adapted from iSCAM stuff (Martell et al)
# ==============================================
#for testing delete, when done

source("/Users/catarinawor/Documents/Length_SRA/R/plots/readscn.R")



plotAgeComps <- function( M )
{
	#n <- length(M)
	cat(".plotAgeComps\n")
	

	#M<-SIMSdat[[802]]


	names(M)
	names(M$SArep)
	names(M$OM)

	st<- M$OM$rep_yr-M$OM$syr+1

	dim(M$SArep$Clt)
	dim(M$OM$Clt[st:nrow(M$SArep$Clt),])

	

	mdf <- NULL
	
	scn<-read_scnnames()




	#for( i in 1:n ){


		

		getDF <- function()
		{
			#ix <- M$SArep$Clt
			
			df <- data.frame((M$SArep$Ulength-M$OM$Ulength[st:nrow(M$SArep$Ulength),])/M$OM$Ulength[st:nrow(M$SArep$Ulength),])
			df <- data.frame(scenario=scn[M$OM$scnNumber],year=M$SArep$yr,df)
			len <- M$SArep$len
			colnames(df) <- c("Scenario","Year",paste(len))
			
			return(df)
		}
		
		#B   <- lapply(1:length(id),getDF)
		B<-getDF()


		# A   <- data.frame(M[[i]]$d3_A)
		# # Ensure proportions are being plotted.
		# A[,-1:-6] <- A[,-1:-6]/rowSums(A[,-1:-6],na.rm=TRUE)
		# age <- seq(min(M[[i]]$n_A_sage),max(M[[i]]$n_A_nage))
		# # year gear area group sex
		# A   <- data.frame(Model=names(M)[i],A)
		# colnames(A) <- c("Model","Year","Gear","Area","Group","Sex","AgeErr",paste(age))
		# mdf <- rbind(mdf,A)

	#}
	mB  <- melt(B,id.vars=c("Scenario","Year"))
	

	# mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex","AgeErr"))
	# BroodYear <- mdf$Year-as.double(mdf$variable)
	# mdf <- cbind(mdf,BroodYear)
	# print(head(mdf,3))

	p <- ggplot(mB,aes((Year),variable,size=value))
	p <- p + geom_point(alpha=0.75,aes(colour=as.factor(sign(value)))) 
	p <- p + scale_size_area(max_size=5)
	p <- p + labs(x="Year",y="length",size="bias")
	p <- p + facet_wrap(~Scenario,scales="free")
	#p <- p + scale_colour_discrete(guide="none")
	p <- p + theme_bw(11)
	print(p)
}