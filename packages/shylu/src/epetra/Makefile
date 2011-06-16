#===============================================================================
#Compile with Trilinos
#===============================================================================

default: driver belos_driver

include ./Makefile.export.Trilinos_install

OBJ = shylu_factor.o shylu_util.o Ifpack_ShyLU.o shylu_probing_operator.o shylu_solve.o shylu_schur.o

all: driver belos_driver

library: $(OBJ) 
	$(Trilinos_AR) cr libshylu.a $(OBJ)

exec: $(OBJ) shylu.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) $(OBJ) shylu.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -o shylu_factor.exe

driver: $(OBJ) shylu_driver.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) $(OBJ) shylu_driver.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -lgfortran -o shylu_driver.exe

belos_driver: $(OBJ) shylu_belos_driver.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) $(OBJ) shylu_belos_driver.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -lgfortran -o shylu_belos_driver.exe

adriver: $(OBJ) amesos_driver.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) $(OBJ) amesos_driver.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -o amesos_driver.exe

prober: shylu_prober.o
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) shylu_prober.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -o shylu_prober.exe

#problem: AztecOO_problem.o
	#$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) AztecOO_problem.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -lgfortran -o AztecOO_problem.exe

#problem: Graph_problem.o
	#$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) Graph_problem.o $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) -lgfortran -o Graph_problem.exe

clean:
	- rm a.out *.o *.a shylu_driver.exe shylu_factor.exe shylu_belos_driver.exe

shylu_factor.o: ./shylu_factor.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

shylu_util.o: ./shylu_util.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

Ifpack_ShyLU.o: ./Ifpack_ShyLU.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

shylu.o: ./shylu.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

shylu_driver.o: ./shylu_driver.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

amesos_driver.o: ./amesos_driver.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

shylu_prober.o: ./shylu_prober.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

shylu_probing_operator.o: ./shylu_probing_operator.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

shylu_belos_driver.o: ./shylu_belos_driver.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

shylu_solve.o: ./shylu_solve.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

shylu_schur.o: ./shylu_schur.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

#Graph_problem.o: ./Graph_problem.cpp
	#$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

#AztecOO_problem.o: ./AztecOO_problem.cpp
	#$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_LIBRARY_DIRS) -c $< -o $@

#-------------------------------------------------------------------------------
