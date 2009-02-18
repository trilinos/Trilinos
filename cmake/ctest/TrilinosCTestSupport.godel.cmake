
#
# Set platform-specific options needed by the platform independent
# code
#

SET( CTEST_CVS_COMMAND "cvs -q -z3" )
SET( CTEST_CMAKE_COMMAND /home/rabartl/install/bin/cmake )

#
# Read in the platform-independent options
#

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestSupport.cmake")


#
# Update or add to platform-specific options
#


SET(MPI_BASE_DIR "/usr/lib64/openmpi/1.2.7-gcc")

# Need to set up the library path for executables created to run currectly
SET(ENV{LD_LIBRARY_PATH} "${MPI_BASE_DIR}/lib:$ENV{LD_LIBRARY_PATH}")


SET(CTEST_INITIAL_CACHE
"
${CTEST_INITIAL_CACHE}
MEMORYCHECK_COMMAND:FILEPATH=/usr/bin/valgrind
MAKECOMMAND:STRING=make -j8 -i
"
)

IF (COMM_TYPE STREQUAL MPI)

  SET(CTEST_INITIAL_CACHE
"
${CTEST_INITIAL_CACHE}
CMAKE_CXX_COMPILER:FILEPATH=${MPI_BASE_DIR}/bin/mpiCC
CMAKE_C_COMPILER:FILEPATH=${MPI_BASE_DIR}/bin/mpicc
CMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/gfortran
CMAKE_EXE_LINKER_FLAGS:STRING=-fprofile-arcs -ftest-coverage
TPL_ENABLE_MPI:BOOL=ON
MPIEXEC_MAX_NUMPROCS:STRING=4
MPI_BASE_DIR:PATH=${MPI_BASE_DIR}
"
  )

ELSE()

  SET(CTEST_INITIAL_CACHE
"
${CTEST_INITIAL_CACHE}
CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++
CMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc
CMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/f77
CMAKE_EXE_LINKER_FLAGS:STRING=-fprofile-arcs -ftest-coverage
"
  )

ENDIF()
