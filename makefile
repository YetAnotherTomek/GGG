#------------------------------#
# GGG MAKEFILE                 #
#                              #
# TR Sokolowski, Mar 2010      #
#------------------------------#

# Some variables used in this file, set to empty by default
# The values are set later for every compilation rule specifically
REMEMBER =
OPT64=

# Source and object files
TOOLS = ./Tools/
HF = $(TOOLS)mersenne_twister.hpp $(TOOLS)array2D.hpp $(TOOLS)tree_node.hpp $(TOOLS)binary_tree.hpp $(TOOLS)random_variable.hpp $(TOOLS)histo1D.hpp $(TOOLS)histo2D.hpp parameters.hpp reaction.hpp observable.hpp nucleus.hpp cortex.hpp simulation.hpp
OF = random_variable.o histo1D.o histo2D.o binary_tree.o parameters.o reaction.o observable.o nucleus.o cortex.o simulation.o

# Compiler and executable name
EXE = gillespie.X$(OPT64)
CFLAGS = $(EXTRAFLAGS)
EXTRAFLAGS = -O3 #-g3
LDFLAGS=
COMP = g++ $(CFLAGS) $(LDFLAGS)

# Archive file 'made' remembering the last compilation
MADE = ${strip $(shell cat made)}

ifeq ($(MADE),x64)
	OPT64 = 64
endif

normal : OPT64 =
normal : REMEMBER = normal
normal : $(EXE)

x64: CFLAGS = $(EXTRAFLAGS) -m64
x64: OPT64 = 64
x64: REMEMBER = x64
x64 : $(EXE)

debug: CFLAGS = $(EXTRAFLAGS) -g3
debug : OPT64 =
debug: REMEMBER = x64
debug: $(EXE)

debug64: CFLAGS = $(EXTRAFLAGS) -g3 -m64
debug64 : OPT64 = 64
debug64: REMEMBER = debug64
debug64: $(EXE)


$(EXE) : $(OF) gillespie.o
	$(COMP) -o $(EXE) gillespie.o $(OF) -lm
	@echo "  Made $(EXE)"
	$(shell echo $(REMEMBER) > made)
random_variable.o : $(TOOLS)random_variable.cpp $(TOOLS)random_variable.hpp makefile
	$(COMP) -c $(TOOLS)random_variable.cpp 
histo1D.o : $(TOOLS)histo1D.cpp $(TOOLS)histo1D.hpp makefile
	$(COMP) -c $(TOOLS)histo1D.cpp 
histo2D.o : $(TOOLS)histo2D.cpp $(TOOLS)histo2D.hpp $(TOOLS)array2D.hpp makefile
	$(COMP) -c $(TOOLS)histo2D.cpp 
binary_tree.o : $(TOOLS)binary_tree.cpp $(TOOLS)binary_tree.hpp $(TOOLS)tree_node.hpp makefile
	$(COMP) -c $(TOOLS)binary_tree.cpp 
parameters.o : parameters.cpp parameters.hpp observable.hpp observable.cpp makefile
	$(COMP) -c parameters.cpp 
reaction.o : reaction.hpp
	$(COMP) -c reaction.cpp 
observable.o : observable.hpp
	$(COMP) -c observable.cpp 
nucleus.o : nucleus.cpp nucleus.hpp reaction.hpp observable.hpp parameters.hpp $(TOOLS)array2D.hpp $(TOOLS)mersenne_twister.hpp makefile
	$(COMP) -c nucleus.cpp 
cortex.o : cortex.cpp cortex.hpp nucleus.hpp parameters.hpp $(TOOLS)binary_tree.hpp $(TOOLS)mersenne_twister.hpp makefile
	$(COMP) -c cortex.cpp 
simulation.o : simulation.cpp simulation.hpp cortex.hpp parameters.hpp $(TOOLS)random_variable.hpp $(TOOLS)histo1D.hpp $(TOOLS)mersenne_twister.hpp makefile
	$(COMP) -c simulation.cpp
gillespie.o : gillespie.cpp parameters.hpp simulation.hpp makefile
	$(COMP) -c gillespie.cpp

clean :
	rm $(OF)
	$(shell echo "clean" > made)
clear :
	rm $(OF)
	$(shell echo "clear" > made)
virgin :
	rm -v *.X *.X64 *.o
	$(shell echo "virgin" > made)
nothing: ;
