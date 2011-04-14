INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu

IF (${PACKAGE_NAME}_ENABLE_Epetra) 
#  IF (${PACKAGE_NAME}_ENABLE_Galeri) #TODO
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      AggregationExample
      SOURCES AggregationExample.cpp
      COMM serial mpi
      DEPLIBS muelu muelu-cthulhu
      )
#  ENDIF()
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ImportExportTest
  SOURCES ImportExportTest.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-cthulhu
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ImportExportTest2
  SOURCES ImportExportTest2.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-cthulhu
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ImportExportTest3
  SOURCES ImportExportTest3.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-cthulhu
  )
