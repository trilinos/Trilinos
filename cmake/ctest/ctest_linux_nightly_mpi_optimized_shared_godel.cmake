
SET(CTEST_SOURCE_NAME Trilinos)
SET(CTEST_BINARY_NAME BUILD)
SET(TEST_TYPE nightly)
SET(BUILD_TYPE release-shared)
SET(EXTRA_BUILD_TYPE mpi)
SET(HOSTTYPE Linux) # Have to set this manually on this machine for some reason?
SET(CTEST_DASHBOARD_ROOT /home/rabartl/PROJECTS/dashboards/Trilinos/MPI_OPT_SHARED)
SET(CTEST_CMAKE_COMMAND /usr/local/bin/cmake)

# Options for Nightly builds

SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)
#SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

SET(CTEST_CVS_CHECKOUT
  "cvs -q -d :ext:software.sandia.gov:/space/CVS co ${CTEST_SOURCE_NAME}"
  )
SET (CTEST_CVS_COMMAND
  "cvs -q"
  )

SET(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
SET(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

SET(TEST_TYPE $ENV{CTEST_TEST_TYPE})
IF (NOT TEST_TYPE)
  SET(TEST_TYPE Nightly)
ENDIF()

SET(CTEST_COMMAND 
  "\"${CTEST_EXECUTABLE_NAME}\" -D ${TEST_TYPE}Start"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ${TEST_TYPE}Update"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ${TEST_TYPE}Configure"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ${TEST_TYPE}Build"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ${TEST_TYPE}Submit"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ${TEST_TYPE}Test"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ${TEST_TYPE}Submit -A \"${CTEST_BINARY_DIRECTORY}/CMakeCache.txt\;${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}\""
)
# 2008/12/12: rabartl: Above, I am doing a second configure with raw
# CMake in order to get around a problem that I have been seeing where
# after the first pass through configuration the ParMETIS and Scotch
# paths are not added correctly.  I have tried to debug this but it is
# very painful to do so and I would rather just avoid this for now.

SET(CTEST_INITIAL_CACHE "

BUILD_SHARED_LIBS:BOOL=ON

Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS:BOOL=OFF

CMAKE_CXX_COMPILER:FILEPATH=/usr/lib64/openmpi/1.2.5-gcc/bin/mpiCC
CMAKE_C_COMPILER:FILEPATH=/usr/lib64/openmpi/1.2.5-gcc/bin/mpicc
CMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/gfortran

CMAKE_CXX_FLAGS:STRING=-O3 -ansi -Wall -Wshadow -Wunused-variable -Wunused-function -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage -fexceptions
CMAKE_C_FLAGS:STRING=-O3 -Wall -fprofile-arcs -ftest-coverage -fexceptions
CMAKE_EXE_LINKER_FLAGS:STRING=-fprofile-arcs -ftest-coverage

MAKECOMMAND:STRING=gmake -j8 -i

MEMORYCHECK_COMMAND:FILEPATH=/usr/local/bin/valgrind

DART_TESTING_TIMEOUT:STRING=600

CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE

TPL_ENABLE_Boost:BOOL=ON

TPL_ENABLE_MPI:BOOL=ON
MPIEXEC_MAX_NUMPROCS:STRING=4
MPI_EXTRA_LIBRARY:FILEPATH=""
MPI_INCLUDE_PATH:FILEPATH=/usr/lib64/openmpi/1.2.5-gcc/bin
MPI_COMPILER:FILEPATH=/usr/lib64/openmpi/1.2.5-gcc/bin/mpiCC
MPI_EXECUTABLE:FILEPATH=/usr/lib64/openmpi/1.2.5-gcc/bin/mpiexec

Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON
Trilinos_ENABLE_AztecOO:BOOL=OFF

Trilinos_ENABLE_TESTS:BOOL=ON

Teuchos_ENABLE_COMPLEX:BOOL=ON
Teuchos_ENABLE_EXTENDED:BOOL=ON
Teuchos_ENABLE_BOOST:BOOL=ON
Teuchos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON

EpetraExt_BUILD_GRAPH_REORDERINGS:BOOL=ON
EpetraExt_BUILD_BDF:BOOL=ON

BUILDNAME:STRING=${HOSTTYPE}-${TEST_TYPE}-${EXTRA_BUILD_TYPE}-${BUILD_TYPE}

CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}

")
