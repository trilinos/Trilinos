
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  TriKota_Simple
  SOURCES 
  Simple_Main.cpp
  Simple_ModelEval.cpp
  Simple_ModelEval.hpp
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS -v
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(TriKotaSimpleCopyDakotaIn
  DEST_FILES   dakota.in
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS  TriKota_Simple 
  )
