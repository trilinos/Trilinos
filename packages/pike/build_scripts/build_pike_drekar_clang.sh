#!/bin/bash
rm -rf CMakeCache.txt CMakeFiles
cmake \
-D CMAKE_INSTALL_PREFIX="/home/rppawlo/install/clang/pike" \
-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
-D Trilinos_ENABLE_TESTS:BOOL=OFF \
-D Trilinos_ENABLE_PikeBlackBox:BOOL=ON \
-D Trilinos_ENABLE_PikeImplicit:BOOL=OFF \
-D Teuchos_ENABLE_LONG_LONG_INT:BOOL=ON \
-D Trilinos_ENABLE_DEBUG=OFF \
-D Pike_ENABLE_DEBUG:BOOL=ON \
-D Pike_ENABLE_TESTS:BOOL=ON \
-D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=OFF \
-D Pike_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
-D TPL_BLAS_LIBRARIES:PATH="/home/rppawlo/install/gnu4.8.1/blas/libblas.a" \
-D TPL_LAPACK_LIBRARIES:PATH="/home/rppawlo/install/gnu4.8.1/lapack/liblapack.a" \
-D TPL_ENABLE_MPI:BOOL=ON \
-D MPI_BASE_DIR:PATH="/home/rppawlo/install/gnu4.8.1/openmpi_clang" \
-D MPI_EXEC_MAX_NUMPROCS:STRING=6 \
-D COMPILER_VERSION="Clang-3.3" \
-D CMAKE_CXX_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.8.1/openmpi_clang/bin/mpiCC" \
-D CMAKE_C_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.8.1/openmpi_clang/bin/mpicc" \
-D CMAKE_Fortran_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.8.1/openmpi_clang/bin/mpif77" \
-D Trilinos_ENABLE_Fortran:BOOL=ON \
-D CMAKE_CXX_FLAGS:STRING="-g -Wall" \
-D CMAKE_C_FLAGS:STRING="" \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
-D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
-D CMAKE_SKIP_RULE_DEPENDENCY=ON \
-D CMAKE_BUILD_TYPE:STRING=DEBUG \
-D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF \
-D Trilinos_ENABLE_CHECKED_STL:BOOL=OFF \
-D Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=OFF \
../Trilinos 
