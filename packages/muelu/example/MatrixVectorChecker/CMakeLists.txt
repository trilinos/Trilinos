
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
SET(simple_example_SOURCES
  MueLu_read_MatrixMarket.cpp
  TumiStuff.h
  )

#PACKAGE_COPY_FILES_TO_BINARY_DIR(read_MatrixMarket_example_copyfiles
#  SOURCE_FILES rgen.m
#)
  
PACKAGE_ADD_EXECUTABLE_AND_TEST(
  read_MatrixMarket
  SOURCES ${simple_example_SOURCES}
  COMM serial
  )
