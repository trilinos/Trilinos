#!/usr/bin/tcsh
rm -rf CMakeCache.txt
cmake \
-D CMAKE_INSTALL_PREFIX="/home/rppawlo/JUNK" \
-D Trilinos_ENABLE_DEBUG=OFF \
-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
-D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON \
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON \
-D Trilinos_ENABLE_NOX:BOOL=ON \
-D Trilinos_ENABLE_TESTS:BOOL=OFF \
-D NOX_ENABLE_TESTS:BOOL=ON \
-D Trilinos_ENABLE_PyTrilinos:BOOL=OFF \
-D Trilinos_ENABLE_Sundance:BOOL=OFF \
-D TPL_ENABLE_MPI:BOOL=ON \
-D MPI_BASE_DIR:PATH="/home/rppawlo/install/gnu4.7.2/openmpi" \
-D TPL_BLAS_LIBRARIES:PATH="/home/rppawlo/install/gnu4.7.2/blas/libblas.a" \
-D TPL_LAPACK_LIBRARIES:PATH="/home/rppawlo/install/gnu4.7.2/lapack/liblapack.a" \
-D CMAKE_CXX_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.7.2/openmpi/bin/mpiCC" \
-D CMAKE_C_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.7.2/openmpi/bin/mpicc" \
-D CMAKE_Fortran_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.7.2/openmpi/bin/mpif77" \
-D CMAKE_CXX_FLAGS:STRING="-g -O3 -ansi -pedantic -ftrapv -Wall -Wno-long-long" \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="-L/usr/lib -lgfortran" \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
-D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
-D CMAKE_SKIP_RULE_DEPENDENCY=ON \
-D CMAKE_BUILD_TYPE:STRING=DEBUG \
-D Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=OFF \
-D Trilinos_ENABLE_EXPORT_MAKEFILES:BOOL=OFF \
-D Trilinos_DEPS_XML_OUTPUT_FILE:FILEPATH="" \
 ../Trilinos

