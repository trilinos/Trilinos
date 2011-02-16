#===============================================================================
#Compile with Trilinos
#===============================================================================

default: all

include ./Makefile.export.Trilinos_install

MY_LIB_DIR=/home/srajama/proj/Trilinos/Trilinos_HyperLU/tc_build/lib
RANLIB=/usr/bin/ranlib

OBJ = hyperlu_factor.o hyperlu_util.o Ifpack_HyperLU.o

all: library exec driver

library: $(OBJ) 
	$(Trilinos_AR) cr libhyperlu.a $(OBJ)
	$(RANLIB) libhyperlu.a

exec: $(OBJ) hyperlu.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -g $(OBJ) hyperlu.o $(MY_LIB_DIR)/libamesos.a $(MY_LIB_DIR)/libepetraext.a $(MY_LIB_DIR)/libepetra.a $(MY_LIB_DIR)/libteuchos.a  $(MY_LIB_DIR)/libisorropia.a  $(MY_LIB_DIR)/libzoltan.a $(MY_LIB_DIR)/libaztecoo.a $(Trilinos_TPL_LIBRARIES) -lgfortran -o hyperlu_factor.exe

driver: $(OBJ) hyperlu_driver.o
#$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -g $(OBJ) hyperlu_driver.o $(MY_LIB_DIR)/libml.a $(MY_LIB_DIR)/libifpack.a $(MY_LIB_DIR)/libamesos.a $(MY_LIB_DIR)/libepetraext.a $(MY_LIB_DIR)/libaztecoo.a $(MY_LIB_DIR)/libepetra.a $(MY_LIB_DIR)/libteuchos.a  $(MY_LIB_DIR)/libisorropia.a  $(MY_LIB_DIR)/libzoltan.a $(Trilinos_TPL_LIBRARIES) -lgfortran -o hyperlu_driver.exe
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -g $(OBJ) hyperlu_driver.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -lgfortran -o hyperlu_driver.exe

prober: hyperlu_prober.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -g hyperlu_prober.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -lgfortran -o hyperlu_prober.exe

clean:
	- rm a.out *.o *.a hyperlu_driver.exe hyperlu_factor.exe

hyperlu_factor.o: ./hyperlu_factor.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c -g $< -o $@

hyperlu_util.o: ./hyperlu_util.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c -g $< -o $@

Ifpack_HyperLU.o: ./Ifpack_HyperLU.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c -g $< -o $@

hyperlu.o: ./hyperlu.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c -g $< -o $@

hyperlu_driver.o: ./hyperlu_driver.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c -g $< -o $@

hyperlu_prober.o: ./hyperlu_prober.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c -g $< -o $@


#-------------------------------------------------------------------------------
