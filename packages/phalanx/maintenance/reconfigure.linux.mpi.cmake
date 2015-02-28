#!/bin/bash
rm -rf CMakeCache.txt CMakeFiles
cmake \
-D CMAKE_INSTALL_PREFIX="/home/rppawlo/JUNK" \
-D Trilinos_ENABLE_DEBUG=OFF \
-D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
-D Trilinos_ENABLE_KokkosClassic:BOOL=OFF \
-D Trilinos_ENABLE_KokkosCore:BOOL=ON \
-D Trilinos_ENABLE_KokkosAlgorithms:BOOL=ON \
-D Trilinos_ENABLE_Xpetra:BOOL=OFF \
-D Trilinos_ENABLE_MueLu:BOOL=OFF \
-D Trilinos_ENABLE_Phalanx:BOOL=ON \
-D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
-D Trilinos_ENABLE_TESTS:BOOL=OFF \
-D Phalanx_KOKKOS_DEVICE_TYPE:STRING="OPENMP" \
-D Phalanx_INDEX_SIZE_TYPE:STRING="INT" \
-D Phalanx_ENABLE_DEBUG:BOOL=OFF \
-D Phalanx_ENABLE_TESTS:BOOL=ON \
-D Phalanx_ENABLE_EXAMPLES:BOOL=ON \
-D Phalanx_EXPLICIT_TEMPLATE_INSTANTIATION=ON \
-D HAVE_INTREPID_KOKKOSCORE:BOOL=ON \
-D TPL_ENABLE_MPI:BOOL=ON \
-D MPI_BASE_DIR:PATH="/home/rppawlo/install/gnu4.8.2/mpich" \
-D TPL_ENABLE_HWLOC:BOOL=ON \
-D HWLOC_INCLUDE_DIRS:FILEPATH="/home/rppawlo/install/gnu4.8.2/hwloc/include" \
-D HWLOC_LIBRARY_DIRS:FILEPATH="/home/rppawlo/install/gnu4.8.2/hwloc/lib" \
-D TPL_ENABLE_Boost:BOOL=ON \
-D Boost_INCLUDE_DIRS:FILEPATH="/home/rppawlo/install/gnu4.8.2/boost/include" \
-D TPL_BLAS_LIBRARIES:PATH="/home/rppawlo/install/gnu4.8.2/blas/libblas.a" \
-D TPL_LAPACK_LIBRARIES:PATH="/home/rppawlo/install/gnu4.8.2/lapack/liblapack.a" \
-D CMAKE_CXX_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.8.2/mpich/bin/mpicxx" \
-D CMAKE_C_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.8.2/mpich/bin/mpicc" \
-D CMAKE_Fortran_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.8.2/mpich/bin/mpifort" \
-D CMAKE_CXX_FLAGS:STRING="-g -ansi -Wno-long-long -Wno-unused-local-typedefs -Wno-maybe-uninitialized" \
-D CMAKE_C_FLAGS:STRING="" \
-D CMAKE_Fortran_FLAGS:STRING="" \
-D Trilinos_ENABLE_CXX11:BOOL=OFF \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
-D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
-D CMAKE_SKIP_RULE_DEPENDENCY=ON \
-D CMAKE_BUILD_TYPE:STRING=RELEASE \
-D Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=OFF \
-D Trilinos_ENABLE_EXPORT_MAKEFILES:BOOL=OFF \
-D Trilinos_DEPS_XML_OUTPUT_FILE:FILEPATH="" \
-D Trilinos_ENABLE_OpenMP:BOOL=ON \
-D Kokkos_ENABLE_OpenMP:BOOL=ON \
../Trilinos

##-D Trilinos_EXTRA_LINK_FLAGS:STRING="-L/usr/lib -lgfortran" \

##-D CMAKE_CXX_FLAGS:STRING="-g -O3 -ansi -pedantic -ftrapv -Wall -Wno-long-long -Wno-unused-local-typedefs" \

##      -D TPL_ENABLE_CppUnit:BOOL=ON \
##      -D CppUnit_INCLUDE_DIRS:FILEPATH="/home/rppawlo/junk/include" \
##      -D CppUnit_LIBRARY_DIRS:FILEPATH="/home/rppawlo/junk/lib" \

##-D TPL_ENABLE_TVMET:BOOL=OFF \
##-D TVMET_INCLUDE_DIRS:FILEPATH="/home/rppawlo/local/include" \

