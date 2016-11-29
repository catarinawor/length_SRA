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
	n <- length(M)
	cat(".plotAgeComps\n")
	
	mdf <- NULL
	
	scn<-read_scnnames()


	for( i in 1:n ){


		names(M[[i]]$SApar)

		getDF <- function(x)
		{
			ix <- id[x]
			
			df <- data.frame(M[[i]][ix])
			df <- data.frame(Model=names(M)[i],df)
			age <- seq(M[[i]]$n_A_sage[x],M[[i]]$n_A_nage[x])
			colnames(df) <- c("Model","Year","Gear","Area","Group","Sex","AgeErr",paste(age))
			
			return(df)
		}
		
		B   <- lapply(1:length(id),getDF)



		# A   <- data.frame(M[[i]]$d3_A)
		# # Ensure proportions are being plotted.
		# A[,-1:-6] <- A[,-1:-6]/rowSums(A[,-1:-6],na.rm=TRUE)
		# age <- seq(min(M[[i]]$n_A_sage),max(M[[i]]$n_A_nage))
		# # year gear area group sex
		# A   <- data.frame(Model=names(M)[i],A)
		# colnames(A) <- c("Model","Year","Gear","Area","Group","Sex","AgeErr",paste(age))
		# mdf <- rbind(mdf,A)

	}
	mB  <- melt(B,id.vars=c("Model","Year","Gear","Area","Group","Sex","AgeErr"))
	BroodYear <- mB$Year-as.double(mB$variable)
	mB  <- cbind(mB,BroodYear)

	# mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex","AgeErr"))
	# BroodYear <- mdf$Year-as.double(mdf$variable)
	# mdf <- cbind(mdf,BroodYear)
	# print(head(mdf,3))

	p <- ggplot(mB,aes((Year),variable,size=value))
	p <- p + geom_point(alpha=0.75,aes(colour=factor(BroodYear))) 
	p <- p + scale_size_area(max_size=5)
	p <- p + labs(x="Year",y="Age",size="Count")
	p <- p + facet_wrap(~L1+Model+Gear+AgeErr+Sex,scales="free")
	p <- p + scale_colour_discrete(guide="none")
	print(p + .THEME + theme(legend.position="top"))
}