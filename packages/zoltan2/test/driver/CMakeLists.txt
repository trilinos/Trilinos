#TRIBITS_ADD_EXECUTABLE_AND_TEST(driver
#  SOURCES ${HEADERS} driver.cpp
#  NUM_MPI_PROCS 4
#  COMM mpi serial
#  ARGS "basic.xml"
#  PASS_REGULAR_EXPRESSION "PASS"
#  FAIL_REGULAR_EXPRESSION "FAIL"
#)
#
#TRIBITS_COPY_FILES_TO_BINARY_DIR(Zoltan2_copy_test_options
#  SOURCE_FILES basic.xml
#  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
#  NOEXEPREFIX
#)

ADD_SUBDIRECTORY(driverinputs)

TRIBITS_ADD_EXECUTABLE(
   test_driver
   SOURCES ${HEADERS} test_driver.cpp
   COMM serial mpi
)

#TRIBITS_ADD_EXECUTABLE(
#   test_crs
#   SOURCES ${HEADERS} tpetra-crs-hello.cpp
#   COMM serial mpi
#)

IF (${PACKAGE_NAME}_ENABLE_Experimental)
# MultiJagged Tests
TRIBITS_ADD_TEST(
    test_driver
    NAME coordinateDriverTest
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/input_072315.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

TRIBITS_ADD_TEST(
    test_driver
    NAME rcbDriverTest
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
    "./driverinputs/rcbTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

ENDIF()

TRIBITS_ADD_TEST(
    test_driver
    NAME mjDriverTest
    NUM_MPI_PROCS 3
    COMM serial mpi
    ARGS
    "./driverinputs/multiJaggedTest.xml"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)

