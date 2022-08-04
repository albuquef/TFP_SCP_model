SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------
#CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio129/cplex/
#CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio129/concert/]

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio201/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio201/concert

CPCPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio201/cpoptimizer
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------
CCC = g++ 
# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
CCOPT = -m64 -g -O -fPIC -fexceptions -DNDEBUG -DIL_STD 
# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
CPLEXBINDIR = $(CPLEXDIR)/bin/$(BINDIST)

CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) 
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -ldl #-framework CoreFoundation -framework IOKit 


CPCPLEXBINDIR = $(CPCPLEXDIR)/bin/$(BINDIST)
CPCPLEXLIBDIR = $(CPCPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPCCLNDIRS  = -L$(CPCPLEXLIBDIR) -L$(CPCONCERTLIBDIR) 

# ------------------------------------------------------------
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR = $(CPLEXDIR)/include

CPCPLEXINCDIR = $(CPCPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
# ------------------------------------------------------------
SRCDIR = ./src
# ------------------------------------------------------------

all: make run 
run: run_mip 
run2: run_mip_decomp 

run_cycle: run_mip_cycle

# ------------------------------------------------------------
fegs: fegs.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(SRCDIR)/fegs $(SRCDIR)/fegs.o $(CCLNFLAGS)

fegs_decomp: fegs_decomp.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(SRCDIR)/fegs_decomp $(SRCDIR)/fegs_decomp.o $(CCLNFLAGS)

fegs_cycle: fegs_cycle.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(SRCDIR)/fegs_cycle $(SRCDIR)/fegs_cycle.o $(CCLNFLAGS)
# ------------------------------------------------------------
fegs.o: $(SRCDIR)/fegs.cpp
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/fegs.cpp -o $(SRCDIR)/fegs.o

fegs_decomp.o: $(SRCDIR)/fegs_decomp.cpp
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/fegs_decomp.cpp -o $(SRCDIR)/fegs_decomp.o

fegs_cycle.o: $(SRCDIR)/fegs_cycle.cpp
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/fegs_cycle.cpp -o $(SRCDIR)/fegs_cycle.o
# ------------------------------------------------------------
run_mip: fegs 
	$(SRCDIR)/fegs

run_mip_decomp: fegs_decomp 
	$(SRCDIR)/fegs_decomp

run_mip_cycle: fegs_cycle 
	$(SRCDIR)/fegs_cycle
# ------------------------------------------------------------
clean:
	rm -rf *.o $(SRCDIR)/*.o *~ $(SRCDIR)/fegs $(SRCDIR)/fegs_decomp $(SRCDIR)/fegs_cycle  #excluir o binario fegsCP
