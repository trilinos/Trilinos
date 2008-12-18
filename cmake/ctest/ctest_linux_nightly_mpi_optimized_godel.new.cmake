
#
# Invocation specific options needed by standard support code
#

SET(CMAKE_BUILD_TYPE DEBUG)

SET(MPI_BASE_DIR "/usr/lib64/openmpi/1.2.5-gcc")

SET(CTEST_INITIAL_CACHE
"

CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}

Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF

CMAKE_CXX_COMPILER:FILEPATH=${MPI_BASE_DIR}/bin/mpiCC
CMAKE_C_COMPILER:FILEPATH=${MPI_BASE_DIR}/bin/mpicc
CMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/gfortran

CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}:STRING=-O3 -fprofile-arcs -ftest-coverage
CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}:STRING=-O3 -fprofile-arcs -ftest-coverage
CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}:STRING=-O5

CMAKE_EXE_LINKER_FLAGS:STRING=-fprofile-arcs -ftest-coverage

MEMORYCHECK_COMMAND:FILEPATH=/usr/local/bin/valgrind

DART_TESTING_TIMEOUT:STRING=600

CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE

Trilinos_ENABLE_Teuchos:BOOL=ON
Trilinos_ENABLE_TESTS:BOOL=ON
Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON

TPL_ENABLE_Boost:BOOL=ON

TPL_ENABLE_MPI:BOOL=ON
MPIEXEC_MAX_NUMPROCS:STRING=4
MPI_BASE_DIR:PATH=${MPI_BASE_DIR}

TPL_ENABLE_ParMETIS:BOOL=ON
ParMETIS_LIBRARY_DIRS:PATH=/home/kddevin/code/ParMETIS3_1

TPL_ENABLE_Scotch:BOOL=ON
Scotch_INCLUDE_DIRS:PATH=/home/kddevin/code/scotch_5.1/include
Scotch_LIBRARY_DIRS:PATH=/home/kddevin/code/scotch_5.1/lib

EpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON
EpetraExt_BUILD_BDF:BOOL=ON

"
)

# ToDo: Put this back!
#Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON

SET(BUILD_DIR_NAME "MPI_OPT")


#
# Include the standard support stuff
#

INCLUDE(${CTEST_SCRIPT_DIRECTORY}/ctest_base.cmake)


#
# Override default options and set other options
#

SET(CTEST_BUILD_NAME "Linux-gcc-mpi-opt-new")
SET(DEFAULT_TEST_TYPE Nightly)
SET(CTEST_CMAKE_COMMAND /usr/local/bin/cmake)
SET(CTEST_BUILD_COMMAND "make -j8 -i")


#
# Run the build
#

DO_TRILINOS_CTEST()
