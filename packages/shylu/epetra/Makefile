#===============================================================================
#Compile with Trilinos
#===============================================================================

default: all

include ./Makefile.export.Trilinos_install

OBJ = hyperlu_factor.o hyperlu_util.o Ifpack_HyperLU.o hyperlu_probing_operator.o

all: library exec driver

library: $(OBJ) 
	$(Trilinos_AR) cr libhyperlu.a $(OBJ)

exec: $(OBJ) hyperlu.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) $(OBJ) hyperlu.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -lgfortran -o hyperlu_factor.exe

driver: $(OBJ) hyperlu_driver.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) $(OBJ) hyperlu_driver.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -lgfortran -o hyperlu_driver.exe

prober: hyperlu_prober.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) hyperlu_prober.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -lgfortran -o hyperlu_prober.exe

clean:
	- rm a.out *.o *.a hyperlu_driver.exe hyperlu_factor.exe

hyperlu_factor.o: ./hyperlu_factor.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

hyperlu_util.o: ./hyperlu_util.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

Ifpack_HyperLU.o: ./Ifpack_HyperLU.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

hyperlu.o: ./hyperlu.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

hyperlu_driver.o: ./hyperlu_driver.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

hyperlu_prober.o: ./hyperlu_prober.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

hyperlu_probing_operator.o: ./hyperlu_probing_operator.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

#-------------------------------------------------------------------------------
