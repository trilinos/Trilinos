INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

IF (${PACKAGE_NAME}_ENABLE_EpetraExt) 

  TRIBITS_ADD_EXECUTABLE(
    Structure2D
    SOURCES Structure2D_epetra.cpp
    COMM serial mpi
    DEPLIBS muelu xpetra
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Structure2d_example_copyfiles
    SOURCE_FILES stru2d_A.txt stru2d_b.txt stru2d_ns.txt
    )
  
ENDIF()

#TODO: run it as a test?