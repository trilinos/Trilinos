
TRIBITS_SUBPACKAGE(MpiComm)

#
# mfh 15 Jan 2014: This subpackage _requires_ MPI, so it doesn't make
# sense to have a test for MPI.
#
# TRIBITS_ADD_OPTION_AND_DEFINE(
#   ${PACKAGE_NAME}_ENABLE_MPI
#   KOKKOSMPICOMM_HAVE_MPI
#   "Enable MPI support in the Kokkos MpiComm subpackage."
#   "${TPL_ENABLE_MPI}"
#   )

ADD_SUBDIRECTORY(src)

#
# mfh 15 Jan 2014: Uncomment this when we add a unit_test directory.
#
#TRIBITS_ADD_TEST_DIRECTORIES(unit_test)

TRIBITS_SUBPACKAGE_POSTPROCESS()

