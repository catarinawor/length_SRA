EXEC=length_sra
DAT=length_sra.dat
CTL=
ARG=
MCFLAG=-mcmc 10000 -mcsave 100 -nosdmcmc
NR=4
OM=OM/simsra
SIMDAT=OM/simsra.dat
LAST := 10
NUMBERS := $(shell seq 1 ${LAST})


        

all: $(EXEC)
	./$(EXEC) -ind $(DAT) $(ARG)

compi: $(EXEC)
	admb $(EXEC) 


run: ./$(EXEC) -ind $(DAT) $(ARG)

mcmc: $(EXEC) $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval


mceval: $(EXEC)
	cp $(CTL).psv $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval

simeval: $(foreach var,$(NUMBERS),./$(OM);\
	./$(EXEC) -ind $(DAT) $(ARG);)
	


clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor 
