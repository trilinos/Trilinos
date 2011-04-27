INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

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
