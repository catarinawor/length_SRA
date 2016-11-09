## Makefile for running simulation evaluation
## I don't think it works on pc! 
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis

OM=simsra
SIMDAT=simsra.dat
LAST := 10
NUMBERS := $(shell seq 1 ${LAST})


        

all: run 

run: $(OM) $(compi)
	-./$(OM) -ind $(SIMDAT) $(ARG)

compi: $(OM)
	admb $(OM) 
	


clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor 
