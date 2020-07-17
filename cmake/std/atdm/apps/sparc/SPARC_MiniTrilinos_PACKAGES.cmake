# Include this file in the driver that calls TRIBITS_CTEST_DRIVER() to set the
# list of packages used by a mini build of SPARC.
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/SPARCMiniTrilinosPackagesList.cmake")
SET(Trilinos_PACKAGES ${SPARC_MiniTrilinos_Packages})
