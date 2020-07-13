# Include this file in the driver that calls TRIBITS_CTEST_DRIVER() to set the
# list of packages used by SPARC.
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/SPARCTrilinosPackagesList.cmake")
SET(Trilinos_PACKAGES ${SPARC_Trilinos_Packages})
# NOTE: Above, even though SPARC_Trilinos_Packages may list a bunch of
# packages that we don't want to test, the file ATDMDisables.cmake has
# <Package>_ENABLE_TESTS set to OFF for all fo the packages that we don't want
# to test.
