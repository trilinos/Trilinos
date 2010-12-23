
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.godel.gcc.cmake")

#
# Test only Zoltan with the C compiler with MPI.  Here the only reason
# that we are testing only Zoltan is that is the only package in
# Trilinos that can build with only C.
#   Use 64 bit global IDs in Zoltan.
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME ZOLTAN_C_64)

#SET(CTEST_DO_COVERAGE_TESTING TRUE)
#SET(CTEST_DO_MEMORY_TESTING TRUE)

SET(Trilinos_PACKAGES Zoltan)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DMPI_EXEC_MAX_NUMPROCS:STRING=11"
  "-DTrilinos_ENABLE_CXX:BOOL=OFF"
  "-DTrilinos_ENABLE_ForTrilinos:BOOL=OFF"
  "-DTrilinos_ENABLE_Fortran:BOOL=OFF"
  "-DZoltan_ENABLE_ULONG_IDS:BOOL=ON"
  "-DTPL_ENABLE_ParMETIS:BOOL=ON"
  "-DParMETIS_LIBRARY_DIRS:PATH=/home/kddevin/code/ParMETIS3_1"
  "-DTPL_ENABLE_Scotch:BOOL=ON"
  "-DScotch_INCLUDE_DIRS:PATH=/home/kddevin/code/scotch_5.1/include"
  "-DScotch_LIBRARY_DIRS:PATH=/home/kddevin/code/scotch_5.1/lib"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
