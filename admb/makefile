## Makefile for running simulation evaluation
## I don't think it works on pc! 
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis

DIREXEC= /Users/catarinawor/Documents/length_SRA/admb/SA
DIROM =  /Users/catarinawor/Documents/length_SRA/admb/OM
OM= simsra
OMDAT= simsra.dat
EXEC= SA/length_sra
DAT= SA/length_sra.dat
RDATA ='source(file.path("../","R/simeval/","readSimEval.R"))'


LAST := 100
NUMBERS := $(shell seq 1 ${LAST})


omi: 
	$(MAKE) --directory $(DIROM) 

sai:
	- $(MAKE) --directory $(DIREXEC) 
        
readRdat: 
	@echo $(RDATA) | R --vanilla --slave


all: omi sai readRdat

simeval: $(foreach var,$(NUMBERS), $(MAKE) all ;)


clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor 