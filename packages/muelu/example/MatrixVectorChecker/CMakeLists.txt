### MatrixVectorCheckerExample:

IF (Cthulhu_ENABLE_Tpetra)
  ADD_SUBDIRECTORY(CMake/tpetra)
ENDIF()

IF (Cthulhu_ENABLE_Epetra)
  ADD_SUBDIRECTORY(CMake/epetra)
ENDIF()

### read_MatrixMarket:

INCLUDE(PackageAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# note: this example don't use cthulhu. Only need Epetra + Galeri.
IF (${PACKAGE_NAME}_ENABLE_Epetra) 
  IF (${PACKAGE_NAME}_ENABLE_Galeri)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      read_MatrixMarket
      SOURCES MueLu_read_MatrixMarket.cpp
      COMM serial
      )
  ENDIF()
ENDIF()

#PACKAGE_COPY_FILES_TO_BINARY_DIR(read_MatrixMarket_example_copyfiles
#  SOURCE_FILES rgen.m
#)
