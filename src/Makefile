SYSTEM     = x86-64_sles10_4.1
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /home/lapo/cplex/cplex125/cplex
CONCERTDIR    = /home/lapo/cplex/cplex125/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -pthread

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include


CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ 
	/bin/rm $(CX_EX) 

# ------------------------------------------------------------

# ------------------------------------------------------------
PRINCOBJ = main.o BranchAndCut.o ColumnGeneration.o Dinic.o Edge.o Instance.o Kruskals.o MinCut.o PathRelinking.o TreeGlover.o Utils.o VariableFixing.o
MAINHDRS = BranchAndCut.h ColumnGeneration.h Instance.h 
BCHDRS = BranchAndCut.h ColumnGeneration.h Utils.h
CGHDRS = ColumnGeneration.h Dinic.h Edge.h Instance.h Kruskals.h MinCut.h PathRelinking.h TypedefsAndDefines.h TreeGlover.h Utils.h VariableFixing.h
DINICHDRS = Dinic.h TypedefsAndDefines.h 
EDGEHDRS = Edge.h 
INSTHDRS = Instance.h 
KRUSKHDRS = Kruskals.h Edge.h Instance.h TypedefsAndDefines.h Utils.h 
MCHDRS = MinCut.h Dinic.h Instance.h TypedefsAndDefines.h Utils.h
PRHDRS = PathRelinking.h Edge.h Instance.h TypedefsAndDefines.h 
TGHDRS = TreeGlover.h Edge.h Instance.h TypedefsAndDefines.h
UTILSHDRS = Utils.h Edge.h Instance.h TypedefsAndDefines.h
VFHDRS = VariableFixing.h Edge.h Instance.h Kruskals.h TypedefsAndDefines.h Utils.h
# ------------------------------------------------------------
all: cgbc

cgbc: $(PRINCOBJ)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o cgbc $(PRINCOBJ) $(CCLNFLAGS)

main.o: main.cpp $(MAINHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c main.cpp $(CCLNFLAGS)	

BranchAndCut.o: BranchAndCut.cpp $(BCHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c BranchAndCut.cpp $(CCLNFLAGS)

ColumnGeneration.o: ColumnGeneration.cpp $(CGHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c ColumnGeneration.cpp $(CCLNFLAGS)

Dinic.o: Dinic.cpp $(DINICHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c Dinic.cpp $(CCLNFLAGS)

Edge.o: Edge.cpp $(EDGEHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c Edge.cpp $(CCLNFLAGS)

Instance.o: Instance.cpp $(INSTHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c Instance.cpp $(CCLNFLAGS)

Kruskals.o: Kruskals.cpp $(KRUSKHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c Kruskals.cpp $(CCLNFLAGS)

MinCut.o:MinCut.cpp $(MCHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c MinCut.cpp $(CCLNFLAGS)

PathRelinking.o: PathRelinking.cpp $(PRHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c PathRelinking.cpp $(CCLNFLAGS)

TreeGlover.o: TreeGlover.cpp $(TGHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c TreeGlover.cpp $(CCLNFLAGS)

Utils.o: Utils.cpp $(UTILSHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -c Utils.cpp $(CCLNFLAGS)

VariableFixing.o: VariableFixing.cpp $(VFHDRS)
	$(CCC) $(CCFLAGS) $(CCLNDIRS)  -c VariableFixing.cpp $(CCLNFLAGS)
