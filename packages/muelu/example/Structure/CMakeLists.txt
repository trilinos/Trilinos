INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

IF (${PACKAGE_NAME}_ENABLE_EpetraExt) 

  PACKAGE_ADD_EXECUTABLE(
    Structure2D
    SOURCES Structure2D_epetra.cpp
    COMM serial mpi
    DEPLIBS muelu muelu-xpetra
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(Structure2d_example_copyfiles
    SOURCE_FILES stru2d.mat stru2d.rhs stru2d.ns
    )
  
ENDIF()

#TODO: run it as a test?