INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Partitioning1
  SOURCES useMueLuGallery.hpp partitioning1.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS 
   "--inputFile=simple.mtx"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(partitioning_copy_files
  SOURCE_FILES
    simple.mtx
)
