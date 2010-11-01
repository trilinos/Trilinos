# MatrixDriver

IF (Cthulhu_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(CMake/tpetra)
ENDIF()

IF (Cthulhu_ENABLE_Epetra)
  ADD_SUBDIRECTORY(CMake/epetra)
ENDIF()


# MatrixDriverTests
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)

SET(MatrixDriverTests_SOURCES
  MueLu_MatrixDriverTests.cpp
  )

IF (Cthulhu_ENABLE_Tpetra)
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  MatrixDriverTests
  SOURCES ${MatrixDriverTests_SOURCES}
  )
ENDIF()