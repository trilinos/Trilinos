# Include this file in the driver that calls TRIBITS_CTEST_DRIVER() to set the
# list of packages used by a mini build of SPARC.
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/SPARCMiniTrilinosPackagesEnables.cmake")
SET(Trilinos_PACKAGES ${SPARC_MiniTrilinos_Package_Enables})
# NOTE: Above, it is harmless to have the package enables set in
# SPARCMiniTrilinosPackagesEnables.cmake in the outer ctest -S driver script.
# We just include that file here to avoid duplicate code.
