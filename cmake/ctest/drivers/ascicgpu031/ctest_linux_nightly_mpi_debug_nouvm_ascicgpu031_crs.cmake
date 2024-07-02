
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.ascicgpu031.gcc-cuda.crs.cmake")

#
# Set the options specific to this build case
#

# The variable BUILD_DIR_NAME is based COMM_TYPE, BUILD_TYPE, and BUILD_NAME_DETAILS.
# Tribits creates the variable listed under "Build Name" by prepending the OS type and compiler
# details to BUILD_DIR_NAME.
SET(COMM_TYPE MPI)
SET(BUILD_TYPE DEBUG)
SET(BUILD_NAME_DETAILS NOUVM_NODEPRECATED_CRS)

SET(CTEST_PARALLEL_LEVEL 8)
SET(CTEST_TEST_TYPE Nightly)
SET(Trilinos_TRACK Nightly)    # Set the CDash track
SET(CTEST_TEST_TIMEOUT 900)

SET(Trilinos_PACKAGES Amesos2 Belos Tpetra Ifpack2 MueLu Xpetra Zoltan2)

SET(EXTRA_CONFIGURE_OPTIONS
  ### TPLS ###
  "-DTPL_ENABLE_SuperLU:BOOL=OFF"
  "-DTPL_ENABLE_HWLOC:BOOL=OFF"
 
  ### PACKAGES CONFIGURATION ###
  "-DTpetra_INST_INT_INT:BOOL=OFF"
  "-DTpetra_INST_INT_LONG_LONG:BOOL=ON"
  "-DTpetra_INST_COMPLEX_FLOAT:BOOL=OFF"

  "-DKokkos_ENABLE_CUDA_UVM:BOOL=OFF"
  "-DTpetra_ENABLE_CUDA_UVM:BOOL=OFF"

  "-DTpetra_ENABLE_DEPRECATED_CODE:BOOL=OFF"

)

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
