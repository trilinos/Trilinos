

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/)
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../base)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(SOURCES "")

FILE(GLOB SOURCES_ALL *.cpp)

# Temporarily disable tests. non-zero return codes on mpi processes
# other than rank 0 are causing trouble with test harness. See #11530
# One of the broken tests is fei_UBase_set_shared_ids.cpp.  There
# could be others as well.

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    fei_ubase
#    SOURCES ${SOURCES}
#    ARGS "--teuchos-suppress-startup-banner"
#    COMM serial mpi
#    NUM_MPI_PROCS 1
#    TESTONLYLIBS fei_test_utils
#   STANDARD_PASS_OUTPUT
#  )

#TRIBITS_ADD_TEST(
#    fei_ubase
#    NAME fei_ubase
#    ARGS "--teuchos-suppress-startup-banner"
#    COMM serial mpi
#    NUM_MPI_PROCS 2
#   STANDARD_PASS_OUTPUT
#  )

#TRIBITS_ADD_TEST(
#    fei_ubase
#    NAME fei_ubase
#    ARGS "--teuchos-suppress-startup-banner"
#    COMM serial mpi
#    NUM_MPI_PROCS 3
#   STANDARD_PASS_OUTPUT
#  )

#TRIBITS_ADD_TEST(
#    fei_ubase
#    NAME fei_ubase
#    ARGS "--teuchos-suppress-startup-banner"
#    COMM serial mpi
#    NUM_MPI_PROCS 4
#   STANDARD_PASS_OUTPUT
#  )

