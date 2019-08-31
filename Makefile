CONCERTDIR = /opt/ibm/ILOG/CPLEX_Studio128/concert
CPLEXDIR = /opt/ibm/ILOG/CPLEX_Studio128/cplex
CPOPTIMIZERDIR = /opt/ibm/ILOG/CPLEX_Studio128/cpoptimizer

CFLAGSCPX = -DIL_STD -O -DNDEBUG -I$(CPOPTIMIZERDIR)/include -I$(CONCERTDIR)/include -I$(CPLEXDIR)/include  -lconcert -lilocplex -lcplex -fPIC -fstrict-aliasing -pedantic -Wall -fexceptions -frounding-math -Wno-long-long -m64 -DILOUSEMT -D_REENTRANT -DILM_REENTRANT

LDFLAGSCPX = -L$(CPOPTIMIZERDIR)/lib/x86-64_linux/static_pic -lcp -L$(CPLEXDIR)/lib/x86-64_linux/static_pic -L$(CONCERTDIR)/lib/x86-64_linux/static_pic -lconcert -lilocplex -lcplex -lpthread -lm -ldl

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = $(CFLAGSCPX)

# CHANGEME: Additional libraries
ADDLIBS = $(LDFLAGSCPX)

# C Compiler command
CC = gcc

# C Compiler options
CFLAGS = -O3 -pipe -DNDEBUG -Wimplicit -Wparentheses -Wsequence-point -Wreturn-type -Wcast-qual -Wall -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD

# additional C Compiler options for linking
CLINKFLAGS =  -Wl,--rpath -Wl,/opt/Ipopt-optimized-3.12.10/build/lib

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = -I`$(CYGPATH_W) /opt/Ipopt-optimized-3.12.10/build/include/coin`  $(ADDINCFLAGS)

# Linker flags
LIBS = -L/opt/Ipopt-optimized-3.12.10/build/lib -lipopt  -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl  -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/opt/intel/compilers_and_libraries_2018.3.222/linux/tbb/lib/intel64_lin/gcc4.7 -L/opt/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64_lin -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl  -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/opt/intel/compilers_and_libraries_2018.3.222/linux/tbb/lib/intel64_lin/gcc4.7 -L/opt/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64_lin -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64_lin -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath -lm  -ldl -lstdc++ -lm

2stpUpperSmoothing.so: 
	g++ $(CFLAGS) $(INCL) $(CLINKFLAGS) -shared -fPIC 2stpUpperSmoothing.cpp main.cpp -o 2stpUpperSmoothing.so $(LIBS) $(ADDLIBS) -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -lgfortran -fopenmp -lpthread -lmkl_mc -lmkl_def -lm -ldl 

2stpUpperSmoothing: 
	g++ $(CFLAGS) $(INCL) $(CLINKFLAGS) 2stpUpperSmoothing.cpp  main.cpp -o 2stpUpperSmoothing -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -lgfortran -fopenmp -lpthread -lm -lmkl_mc -lmkl_def -ldl $(LIBS) $(ADDLIBS)
