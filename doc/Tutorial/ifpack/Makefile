include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_ARCH)

MY_COMPILER_FLAGS = -DHAVE_CONFIG_H $(CXXFLAGS) -c -g\
	            -I$(TRILINOS_LIB)/include/
		    
MY_LINKER_FLAGS = $(LDFLAGS) $(TEST_C_OBJ) \
	-L$(TRILINOS_LIB)/lib/ \
	-lifpack -laztecoo -lepetra $(LIB_PATHS) $(ARCH_LIBS) -llapack -lblas

all: ex1

ex1: ex1.cpp
	$(CXX)     ex1.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex1.o   $(MY_LINKER_FLAGS)    -o ex1.exe
	/bin/rm -r ex1.o

clean:
	/bin/rm -f *.o ex?.exe
