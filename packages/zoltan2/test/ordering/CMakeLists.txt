INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Ordering1
  SOURCES useMueLuGallery.hpp ordering1.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple.mtx"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  OrderingAMD
  SOURCES useMueLuGallery.hpp orderingAMD.cpp
  NUM_MPI_PROCS 4
  COMM serial
  ARGS
   "--inputFile=simple.mtx"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ordering_copy_files
  SOURCE_FILES
    simple.mtx
)
