#!/usr/bin/tcsh
cmake \
      -D CMAKE_INSTALL_PREFIX="/home/rppawlo/JUNK15" \
      -D Trilinos_ENABLE_DEBUG=OFF \
      -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
      -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=ON \
      -D NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK=ON \
      -D Trilinos_ENABLE_NOX:BOOL=ON \
      -D Trilinos_ENABLE_TESTS:BOOL=OFF \
      -D NOX_ENABLE_TESTS:BOOL=ON \
      -D Trilinos_ENABLE_PyTrilinos:BOOL=OFF \
      -D Trilinos_ENABLE_Sundance:BOOL=OFF \
      -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_BASE_DIR:PATH="/home/rppawlo/local" \
      -D MPIEXEC_MAX_NUMPROCS:STRING=4 \
      -D CMAKE_CXX_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.6.3/openmpi/bin/mpiCC" \
      -D CMAKE_C_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.6.3/openmpi/bin/mpicc" \
      -D CMAKE_Fortran_COMPILER:FILEPATH="/home/rppawlo/install/gnu4.6.3/openmpi/bin/mpif77" \
      -D CMAKE_CXX_FLAGS:STRING="-g -O3 -ansi -pedantic -ftrapv -Wall -Wno-long-long" \
      -D Trilinos_EXTRA_LINK_FLAGS:STRING="-L/usr/lib -lgfortran" \
      -D HAVE_STRING_H:BOOL=ON \
      -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
      -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
      -D CMAKE_SKIP_RULE_DEPENDENCY=ON \
      -D CMAKE_BUILD_TYPE:STRING=DEBUG \
       ../Trilinos

##      -D CMAKE_EXE_LINKER_FLAGS:STRING="-L/usr/lib -lgfortran" \
