#!/usr/bin/tcsh
rm CMakeCache.txt
cmake \
      -D CMAKE_INSTALL_PREFIX="/home/rppawlo/trilinos_install" \
      -D Trilinos_ENABLE_DEBUG=ON \
      -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
      -D Trilinos_ENABLE_Phalanx:BOOL=ON \
      -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
      -D Trilinos_ENABLE_TESTS:BOOL=OFF \
      -D Phalanx_ENABLE_TESTS:BOOL=ON \
      -D Phalanx_ENABLE_EXAMPLES:BOOL=ON \
      -D Phalanx_EXPLICIT_TEMPLATE_INSTANTIATION=ON \
      -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_BASE_DIR:PATH="/home/rppawlo/local" \
      -D MPIEXEC_MAX_NUMPROCS:STRING=4 \
      -D TPL_ENABLE_Boost:BOOL=ON \
      -D Boost_INCLUDE_DIRS:FILEPATH="/home/rppawlo/Libs/Boost/boost_1_46_1" \
      -D TPL_ENABLE_TVMET:BOOL=OFF \
      -D TVMET_INCLUDE_DIRS:FILEPATH="/home/rppawlo/local/include" \
      -D CMAKE_CXX_COMPILER:FILEPATH="/home/rppawlo/local/bin/mpiCC" \
      -D CMAKE_C_COMPILER:FILEPATH="/home/rppawlo/local/bin/mpicc" \
      -D CMAKE_Fortran_COMPILER:FILEPATH="/home/rppawlo/local/bin/mpif77" \
      -D CMAKE_CXX_FLAGS:STRING="-g -O0 -ansi -pedantic -ftrapv -Wall -Wno-long-long -Wno-strict-aliasing -DBOOST_NO_HASH" \
      -D Trilinos_EXTRA_LINK_FLAGS:STRING="-L/usr/lib -lgfortran" \
      -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
      -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
      -D CMAKE_SKIP_RULE_DEPENDENCY=ON \
      -D Trilinos_ENABLE_STRONG_CXX_COMPILE_WARNINGS=OFF \
      -D Trilinos_ENABLE_STRONG_C_COMPILE_WARNINGS=OFF \
      -D Trilinos_ENABLE_SHADOW_WARNINGS=OFF \
      -D CMAKE_BUILD_TYPE:STRING=NONE \
       ../Trilinos

##      -D CMAKE_EXE_LINKER_FLAGS:STRING="-L/usr/lib -lgfortran" \

##      -D TPL_ENABLE_CppUnit:BOOL=ON \
##      -D CppUnit_INCLUDE_DIRS:FILEPATH="/home/rppawlo/junk/include" \
##      -D CppUnit_LIBRARY_DIRS:FILEPATH="/home/rppawlo/junk/lib" \
##      -D TPL_ENABLE_ADOLC:BOOL=ON \
##      -D ADOLC_INCLUDE_DIRS:FILEPATH="/home/rppawlo/junk/include" \
##      -D ADOLC_LIBRARY_DIRS:FILEPATH="/home/rppawlo/junk/lib" \

##      -D CMAKE_CXX_FLAGS:STRING="-DNDEBUG -O3 -ansi -pedantic -ftrapv -Wall -Wno-long-long -Wno-strict-aliasing -DBOOST_NO_HASH" \
