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
#CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio129/concert/

# CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio201/cplex
# CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio201/concert
# CPCPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio201/cpoptimizer


CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio2211/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio2211/concert
CPCPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio2211/cpoptimizer

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
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -ldl  #-framework CoreFoundation -framework IOKit 


CPCPLEXBINDIR = $(CPCPLEXDIR)/bin/$(BINDIST)
CPCPLEXLIBDIR = $(CPCPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPCCLNDIRS  = -L$(CPCPLEXLIBDIR) -L$(CPCONCERTLIBDIR) 

# ------------------------------------------------------------
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR = $(CPLEXDIR)/include

CPCPLEXINCDIR = $(CPCPLEXDIR)/includerun_SPC

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
# ------------------------------------------------------------
SRCDIR = ./src
# ------------------------------------------------------------
DEFAULT = arg1
DEPS = teste.o main.o
OBJ = teste.o main.o

# all: run1 run2 run3
test: main
# runTest: main
# all: run1 run2 run3 
# run1: run_SPC 
# run2: run_SPC_simplified 
# run3: run_SIPC 


OBJECTS = main.o TFP_SCP.o graph.o matching.o TFP_SCP_simp.o

main: $(OBJECTS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) $(SRCDIR)/main.o $(SRCDIR)/TFP_SCP.o $(SRCDIR)/graph.o $(SRCDIR)/matching.o $(SRCDIR)/TFP_SCP_simp.o  $(CCLNFLAGS) -o $(SRCDIR)/main 

main.o: $(SRCDIR)/main.cpp
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/main.cpp -o $(SRCDIR)/main.o -w 

TFP_SCP.o: $(SRCDIR)/TFP_SCP.cpp $(SRCDIR)/TFP_SCP.h
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/TFP_SCP.cpp -o $(SRCDIR)/TFP_SCP.o -w 

graph.o: $(SRCDIR)/graph.cpp $(SRCDIR)/graph.h
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/graph.cpp -o $(SRCDIR)/graph.o -w 

matching.o: $(SRCDIR)/matching.cpp $(SRCDIR)/matching.h
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/matching.cpp -o $(SRCDIR)/matching.o -w 

TFP_SCP_simp.o: $(SRCDIR)/TFP_SCP_simp.cpp $(SRCDIR)/TFP_SCP_simp.h
	$(CCC) -c $(CCFLAGS) $(SRCDIR)/TFP_SCP_simp.cpp -o $(SRCDIR)/TFP_SCP_simp.o -w

# deafult values 
# PROB?=TFP_SCP  # ?= assignment operator. It only assigns the value if it is not defined at runtime.
PROB?=TFP_SCP_simp 
NUM_VERT?=14
SEC?=MTZ
METHOD?=MinMatching
GCLASS?= 1
CTYPE?= 1
INSTYPE?=random
GTYPE?=undirected 
VALID_INEQ?=null

test: main
	$(SRCDIR)/main -prob $(PROB) -n $(NUM_VERT) -gclass $(GCLASS) -ctype $(CTYPE) -itype $(INSTYPE) -gtype $(GTYPE) -sec $(SEC) -method $(METHOD) -validIneq $(VALID_INEQ)
# args respect the order



# ------------------------------------------------------------
# make run3 ARGS+="teste1" ARGS+="teste2"
# ------------------------------------------------------------

# ------------------------------------------------------------
# main: main.o
# 	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(SRCDIR)/main $(SRCDIR)/main.o $(CCLNFLAGS)

# TFP_SPC: TFP_SPC.o
# 	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(SRCDIR)/TFP_SPC $(SRCDIR)/TFP_SPC.o $(CCLNFLAGS)

# TFP_SPC_simplified: TFP_SPC_simplified.o
# 	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(SRCDIR)/TFP_SPC_simplified $(SRCDIR)/TFP_SPC_simplified.o $(CCLNFLAGS)

# TFP_SIPC: TFP_SIPC.o
# 	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(SRCDIR)/TFP_SIPC $(SRCDIR)/TFP_SIPC.o $(CCLNFLAGS) -w
# # ------------------------------------------------------------
# main.o: $(SRCDIR)/main.cpp
# 	$(CCC) -c $(CCFLAGS) $(SRCDIR)/main.cpp -o $(SRCDIR)/main.o -w 

# TFP_SPC.o: $(SRCDIR)/TFP_SPC.cpp
# 	$(CCC) -c $(CCFLAGS) $(SRCDIR)/TFP_SPC.cpp -o $(SRCDIR)/TFP_SPC.o -w

# TFP_SPC_simplified.o: $(SRCDIR)/TFP_SPC_simplified.cpp
# 	$(CCC) -c $(CCFLAGS) $(SRCDIR)/TFP_SPC_simplified.cpp -o $(SRCDIR)/TFP_SPC_simplified.o -w

# TFP_SIPC.o: $(SRCDIR)/TFP_SIPC.cpp
# 	$(CCC) -c $(CCFLAGS) $(SRCDIR)/TFP_SIPC.cpp -o $(SRCDIR)/TFP_SIPC.o -w

# -w remove all warnings 
# ------------------------------------------------------------
# test: main
# 	$(SRCDIR)/main $(ARGS)

# run_SPC: TFP_SPC 
# 	$(SRCDIR)/TFP_SPC $(ARGS)

# run_SPC_simplified: TFP_SPC_simplified 
# 	$(SRCDIR)/TFP_SPC_simplified $(ARGS)

# run_SIPC: TFP_SIPC 
# 	$(SRCDIR)/TFP_SIPC $(ARGS)

# runTest: main
# 	$(SRCDIR)/main $(METHOD) $(NUM_VERT) $(GRAPH_CLASS) $(CLASS_TYPE) $(INSTANCE_TYPE)

# ------------------------------------------------------------
# make run3 ARGS+="teste1" ARGS+="teste2"
# ------------------------------------------------------------
clean:
	rm -rf *.o $(SRCDIR)/*.o *~ $(SRCDIR)/output $(SRCDIR)/main $(SRCDIR)/TFP_SPC $(SRCDIR)/TFP_SPC_simplified $(SRCDIR)/TFP_SIPC  #delete binaries
 