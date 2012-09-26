INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

IF (${PACKAGE_NAME}_ENABLE_EpetraExt) 
  #  IF (${PACKAGE_NAME}_ENABLE_Galeri) #TODO
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    AggregationExample
    SOURCES AggregationExample.cpp
    COMM serial mpi
    DEPLIBS muelu xpetra
    )
  #  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Epetra)
    TRIBITS_ADD_EXECUTABLE(
      VisualizeAggregates
      SOURCES visualizeAggregates.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )
    
    TRIBITS_COPY_FILES_TO_BINARY_DIR(VisualizeAggregates_example_copyfiles
      SOURCE_FILES Condif2Mat.mat Condif2Rhs.mat plotaggs.py convhull.py nodes0.txt
      )
  ENDIF()

ENDIF()
