include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_ARCH)

MY_COMPILER_FLAGS = -DHAVE_CONFIG_H $(CXXFLAGS) -c -g\
	            -I$(TRILINOS_LIB)/include/
		    
MY_LINKER_FLAGS = $(LDFLAGS) $(TEST_C_OBJ) \
	-L$(TRILINOS_LIB)/lib/ \
	-laztecoo -lml -lepetra $(LIB_PATHS) $(ARCH_LIBS) -llapack -lblas

METIS_DIR = -L/home/msala/lib/Linux
METIS_INC = -I/home/msala/include
METIS_LIB = -lmetis-4.0

all: ex1

ex1: ex1.cpp
	$(CXX)     ex1.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex1.o   $(MY_LINKER_FLAGS)    -o ex1.exe
	/bin/rm -r ex1.o

ex2: ex2.cpp
	$(CXX)     ex2.cpp $(MY_COMPILER_FLAGS) $(METIS_INC)
	$(LINKER)  ex2.o   $(MY_LINKER_FLAGS)  \
	$(METIS_DIR)  $(METIS_LIB) -o ex2.exe
	/bin/rm -r ex2.o

clean:
	/bin/rm -f *.o ex?.exe
