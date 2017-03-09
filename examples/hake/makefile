
EXEC=length_SRA
DAT=length_SRA.dat
LOCDAT=hake.dat
LOCCTL=hake.ctl

DIST= /Users/catarinawor/Documents/length_SRA/admb/SA/length_SRA


all: $(EXEC) $(EXEC).rep

$(EXEC): $(DIST)
	cp $(DIST) $@

$(EXEC).rep: $(DIST) $(LOCDAT) $(LOCCTL)
	./$(EXEC) -ind $(DAT) 

clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor 
