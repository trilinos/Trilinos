INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_Tpetra)  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    MatrixDriver
    SOURCES MueLu_MatrixDriver.cpp
    )
ENDIF()

# Note: this example only use muelu/gallery (do not requiere MueLu or Cthulhu)
# muelu/gallery requiere Tpetra (TODO: direct dependency to build muelu/gallery)