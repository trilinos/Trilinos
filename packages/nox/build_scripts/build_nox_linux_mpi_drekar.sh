#!/bin/bash
rm -rf CMakeCache.txt CMakeFiles
cmake \
-D CMAKE_INSTALL_PREFIX="/home/rppawlo/JUNK" \
-D Trilinos_ENABLE_DEBUG=OFF \
-D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
-D AztecOO_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
-D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES=OFF \
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON \
-D Trilinos_ENABLE_Fortran:BOOL=OFF \
-D Trilinos_ENABLE_TESTS:BOOL=OFF \
-D Trilinos_ENABLE_NOX:BOOL=ON \
-D NOX_ENABLE_TESTS:BOOL=ON \
-D NOX_ENABLE_EXAMPLES:BOOL=ON \
-D ML_ENABLE_PETSC:BOOL=OFF \
-D Trilinos_ENABLE_PyTrilinos:BOOL=OFF \
-D Trilinos_ENABLE_Sundance:BOOL=OFF \
-D TPL_ENABLE_MPI:BOOL=ON \
-D MPI_BASE_DIR:PATH="/home/rppawlo/install/gnu4.8.2/mpich" \
-D TPL_ENABLE_PETSC:BOOL=OFF \
-D PETSC_INCLUDE_DIRS:PATH="/home/rppawlo/install/gnu4.8.2/petsc-3.5/include" \
-D PETSC_LIBRARY_DIRS:PATH="/home/rppawlo/install/gnu4.8.2/petsc-3.5/lib" \
-D TPL_BLAS_LIBRARIES:PATH="/home/rppawlo/install/gnu4.8.2/blas/libblas.a" \
-D TPL_LAPACK_LIBRARIES:PATH="/home/rppawlo/install/gnu4.8.2/lapack/liblapack.a" \
-D CMAKE_CXX_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.8.2/mpich/bin/mpicxx" \
-D CMAKE_C_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.8.2/mpich/bin/mpicc" \
-D CMAKE_Fortran_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.8.2/mpich/bin/mpifort" \
-D CMAKE_CXX_FLAGS:STRING="-g -Werror -ansi -pedantic -ftrapv -Wall -Wno-long-long" \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="-L/usr/lib -lgfortran -ldl" \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
-D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
-D CMAKE_SKIP_RULE_DEPENDENCY=ON \
-D CMAKE_BUILD_TYPE:STRING=DEBUG \
-D Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=OFF \
-D Trilinos_ENABLE_EXPORT_MAKEFILES:BOOL=OFF \
-D Trilinos_DEPS_XML_OUTPUT_FILE:FILEPATH="" \
 ../Trilinos

# -G Ninja \
