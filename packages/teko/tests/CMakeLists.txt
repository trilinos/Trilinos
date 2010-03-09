INCLUDE(PackageLibraryMacros)
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS "")
SET(SOURCES "")

APPEND_GLOB( HEADERS 
  ./src/Epetra/*.hpp
  ./src/*.hpp
)

APPEND_GLOB( SOURCES 
  ./src/Epetra/*.cpp
  ./src/*.cpp
)

PACKAGE_ADD_EXECUTABLE(
  testdriver
  SOURCES
      ${HEADERS}
      ${SOURCES}
      ./Test_Utils.cpp
      ./Test_Utils.hpp
      ./testdriver.cpp
  COMM mpi serial
  )

PACKAGE_ADD_TEST(
   testdriver
   ARGS
   NUM_MPI_PROCS 4
   FAIL_REGULAR_EXPRESSION "Teko tests failed"
   )

PACKAGE_ADD_TEST(
   testdriver
   ARGS
   NUM_MPI_PROCS 1
   FAIL_REGULAR_EXPRESSION "Teko tests failed"
   )

PACKAGE_COPY_FILES_TO_BINARY_DIR(testdriver_copyfiles
  SOURCE_FILES lsc_B_2.mm
     lsc_Bt_2.mm
     lsc_F_2.mm
     lsc_Qu_2.mm
     lsc_exact_2.mm
     lsc_rhs.mm
     nsjac.mm
     nsjac_1.mm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
  EXEDEPS testdriver
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  IterativePreconditionerFactory_test
  SOURCES
    unit_tests/tIterativePreconditionerFactory.cpp
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
