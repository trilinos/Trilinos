
#
# Set the options specific to this build
#

SET(BUILD_NAME SERIAL)
SET(BUILD_TYPE DEBUG)
SET(CTEST_DO_COVERAGE_TESTING TRUE)
SET(CTEST_DO_MEMORY_TESTING TRUE)


SET(CTEST_INITIAL_CACHE "

CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++
CMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc
CMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/gfortran

CMAKE_CXX_FLAGS_${BUILD_TYPE}:STRING=-g -O0 -fprofile-arcs -ftest-coverage
CMAKE_C_FLAGS)${BUILD_TYPE}:STRING=-g -O0 -fprofile-arcs -ftest-coverage
CMAKE_Fortran_FLAGS_${BUILD_TYPE}:STRING=-g -O0 -fprofile-arcs -ftest-coverage

CMAKE_EXE_LINKER_FLAGS:STRING=-fprofile-arcs -ftest-coverage -lm

DART_TESTING_TIMEOUT:STRING=600
CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE

Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF

Trilinos_ENABLE_Teuchos:BOOL=ON
Trilinos_ENABLE_TESTS:BOOL=ON
Trilinos_ENABLE_DEBUG:BOOL=ON_
Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON

TPL_ENABLE_Boost:BOOL=ON

EpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON
EpetraExt_BUILD_BDF:BOOL=ON

")


# 2008/12/18: rabartl: ToDo: Put this back
#Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON


#
# Read in the platform-independent and platform-dependent options
#

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestSupport.godel.cmake")
