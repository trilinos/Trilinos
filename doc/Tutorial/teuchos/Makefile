include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_ARCH)

MY_COMPILER_FLAGS = -DHAVE_CONFIG_H $(CXXFLAGS) -c -g\
	            -I$(TRILINOS_LIB)/include/ \
		    -I/home/msala/Trilinos/packages/amesos/src/

MY_LINKER_FLAGS = $(LDFLAGS) $(TEST_C_OBJ) \
	-L$(TRILINOS_LIB)/lib/ \
	-lteuchos -llapack -lblas -lfrtbegin -lg2c

all: ex1 ex2 ex3 ex4

ex1: ex1.cpp
	$(CXX)     ex1.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex1.o   $(MY_LINKER_FLAGS)    -o ex1.exe
	/bin/rm -r ex1.o

ex2: ex2.cpp
	$(CXX)     ex2.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex2.o   $(MY_LINKER_FLAGS)    -o ex2.exe
	/bin/rm -r ex2.o

ex3: ex3.cpp
	$(CXX)     ex3.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex3.o   $(MY_LINKER_FLAGS)    -o ex3.exe
	/bin/rm -r ex3.o

ex4: ex4.cpp
	$(CXX)     ex4.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex4.o   $(MY_LINKER_FLAGS)    -o ex4.exe
	/bin/rm -r ex4.o

ex5: ex5.cpp
	$(CXX)     ex5.cpp $(MY_COMPILER_FLAGS) 
	$(LINKER)  ex5.o   $(MY_LINKER_FLAGS)    -o ex5.exe
	/bin/rm -r ex5.o

clean:
	/bin/rm -f *.o ex?.exe
