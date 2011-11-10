INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Partitioning1
  SOURCES useMueLuGallery.hpp partitioning1.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS 
   "--inputFile=simple.mtx"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(Zoltan2_partitioning_copy_files
  SOURCE_FILES
    simple.mtx
)
