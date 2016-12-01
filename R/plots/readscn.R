



read_scnnames <- function( )
{
	setwd("/Users/catarinawor/Documents/length_SRA/admb/OM")
	a<-scan("scenarios.txt", what="character", comment.char="#")
	return(a)
}