# -*- mode: makefile -*-
include $(TRILINOS_MAKEFILE)


data_redist_1: data_redist_1.cpp
initializing_tpetra_with_standalone_mpi: initializing_tpetra_with_standalone_mpi.cpp
initializing_tpetra_with_teuchos_mpi: initializing_tpetra_with_teuchos_mpi.cpp
initializing_tpetra_with_teuchos_serial: initializing_tpetra_with_teuchos_serial.cpp
map_contiguous: map_contiguous.cpp
map_contiguous_and_uniform: map_contiguous_and_uniform.cpp
map_contiguous_no_global_num: map_contiguous_no_global_num.cpp
map_cyclic: map_cyclic.cpp
matrix_construct_heat2d_1: matrix_construct_heat2d_1.cpp
matrix_construct_heat2d_2: matrix_construct_heat2d_2.cpp
matrix_fill_1: matrix_fill_1.cpp
power_method_1: power_method_1.cpp
vector: vector.cpp


# The compile line is enourmous.  Order matters
# <compiler> <flags> <includes> <source files> <output name> <library paths> <libraries> <TPL library paths> <TPL libraries>
% :: %.cpp
	$(Trilinos_CXX_COMPILER) $(PY_CFLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) $(Trilinos_CXX_COMPILER_FLAGS) $(GALERI_FLAG) $(PYTHON_FLAG) $< -o $@.exe $(Trilinos_LIBRARY_DIRS) $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARIES) $(PY_LDFLAGS)

a: a.cpp
	$(Trilinos_CXX_COMPILER) $(Trilinos_CXX_COMPILER_FLAGS) $(PY_CFLAGS) a.cpp -o a.x $(PY_LDFLAGS)

clean:
	rm -rf *.x *.o *.dSYM
