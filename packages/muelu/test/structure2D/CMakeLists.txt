INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

IF (${PACKAGE_NAME}_ENABLE_EpetraExt) 

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Structure2D_Test
  SOURCES Structure2D_epetra.cpp
  COMM serial mpi
  DEPLIBS muelu muelu-xpetra
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Structure2d_test_copyfiles
SOURCE_FILES stru2d_A.txt stru2d_b.txt stru2d_ns.txt
)
  
ENDIF()
