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

#set(HEADERS
#    Zoltan2_MeshCoordinateTest.hpp
#    Zoltan2_TestInterface.hpp
#    Zoltan2_Tests.hpp
#)

ADD_SUBDIRECTORY(driverinputs)

#KDD TRIBITS_ADD_EXECUTABLE(
#KDD   test_driver
#KDD   SOURCES ${HEADERS} test_driver.cpp
#KDD   COMM serial mpi
#KDD )
