
#
# Set the options specific to this build
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME MPI_OPT_SHARED)


SET(CTEST_INITIAL_CACHE "

BUILD_SHARED_LIBS:BOOL=ON

CMAKE_CXX_FLAGS_${BUILD_TYPE}:STRING=-O3 -fprofile-arcs -ftest-coverage
CMAKE_C_FLAGS_${BUILD_TYPE}:STRING=-O3 -fprofile-arcs -ftest-coverage
CMAKE_Fortran_FLAGS_${BUILD_TYPE}:STRING=-O5 -fprofile-arcs -ftest-coverage

DART_TESTING_TIMEOUT:STRING=600
CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE

Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON
Trilinos_ENABLE_AztecOO:BOOL=OFF

Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF

Trilinos_ENABLE_TESTS:BOOL=ON

Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON

TPL_ENABLE_Boost:BOOL=ON

EpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON
EpetraExt_BUILD_BDF:BOOL=ON

")

# 2008/12/19: rabartl: AztecOO is disabled above because does not
# currently work in shared library mode (see bug 4288).

# 2009/01/09: rabartl: FEI is disabled because AztecOO is disabled and
# FEI does not compile when AztecOO is disabled.  Brent said that he
# would look into this.

#
# Read in the platform-independent and platform-dependent options
#

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestSupport.godel.cmake")
