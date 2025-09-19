TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(MatrixDriver SOURCES MueLu_MatrixDriver.cpp)

IF(${PACKAGE_NAME}_ENABLE_Epetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(Epetra1DLaplace SOURCES Epetra1DLaplace.cpp MueLu_MemoryProfiler.cpp)
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(Tpetra1DLaplace SOURCES Tpetra1DLaplace.cpp MueLu_MemoryProfiler.cpp)

# Note: this example only use muelu/gallery (do not requiere MueLu or Xpetra)
# muelu/gallery requiere Tpetra (TODO: direct dependency to build muelu/gallery)
