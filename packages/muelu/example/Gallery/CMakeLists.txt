INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_Tpetra)  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    MatrixDriver
    SOURCES MueLu_MatrixDriver.cpp
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra)  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Epetra1DLaplace
    SOURCES Epetra1DLaplace.cpp
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra)  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Tpetra1DLaplace
    SOURCES Tpetra1DLaplace.cpp
    )
ENDIF()

# Note: this example only use muelu/gallery (do not requiere MueLu or Xpetra)
# muelu/gallery requiere Tpetra (TODO: direct dependency to build muelu/gallery)