INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Galeri)  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MatrixDriver
    SOURCES MueLu_MatrixDriver.cpp
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Galeri)  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Epetra1DLaplace
    SOURCES Epetra1DLaplace.cpp MueLu_MemoryProfiler.cpp
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Galeri)  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra1DLaplace
    SOURCES Tpetra1DLaplace.cpp MueLu_MemoryProfiler.cpp
    )
ENDIF()

# Note: this example only use muelu/gallery (do not requiere MueLu or Xpetra)
# muelu/gallery requiere Tpetra (TODO: direct dependency to build muelu/gallery)