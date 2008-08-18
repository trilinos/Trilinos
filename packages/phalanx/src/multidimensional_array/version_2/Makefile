#####################################################################
## Example Makefile for a user application that does not use autoconf
## - Uses lapack concrete instations for group and vector 
## - Must use gnu-make (gmake) if the "shell" command is invoked
#####################################################################

CXX = g++
CXXFLAGS = -g -O2 -I.

#LIB_SOURCES := $(wildcard *.cpp) $(wildcard *.hpp)
LIB_SOURCES := $(wildcard *.cpp)
#LIB_SOURCES := main.cpp XPR_Variable.cpp XPR_Topology_Generic.cpp)
LIB_OBJECTS := $(addsuffix .o, $(basename $(LIB_SOURCES)))

main.exe: $(LIB_OBJECTS)
	$(CXX) $(CXXFLAGS) -o main.exe $(LIB_OBJECTS) $(LINK_FLAGS)

##main.o: main.cpp
##	$(TEUCHOS_CXX) $(CXXFLAGS) -c main.cpp

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o main.exe *~
