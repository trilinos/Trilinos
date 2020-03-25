
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Ordering1
  SOURCES ordering1.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_TEST(
  Ordering1
  NAME Ordering1_Ashear
  COMM serial mpi
  ARGS
   "--inputFile=Ashear"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PROJECT_NAME}_ENABLE_Galeri)

  TRIBITS_ADD_TEST(
    Ordering1
    NAME Ordering1_natural
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
     "--order_method=natural"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    Ordering1
    NAME Ordering1_random
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
     "--order_method=random"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    Ordering1
    NAME Ordering1_rcm
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
     "--order_method=rcm"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

  TRIBITS_ADD_TEST(
    Ordering1
    NAME Ordering1_sorteddegree
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
     "--order_method=sorted_degree"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

# Always run AMD test 
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    OrderingAMD
    SOURCES orderingAMD.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS
     "--inputFile=simple"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

IF (${PACKAGE_NAME}_ENABLE_Scotch)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    OrderingScotch
    SOURCES orderingScotch.cpp
    COMM serial mpi
    ARGS
     "--inputFile=simple_ordering"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

##########################################
IF (${PACKAGE_NAME}_ENABLE_Experimental)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ND
  SOURCES nd.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  ARGS
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )


  TRIBITS_ADD_TEST(
    ND
    NAME ND_4
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS
     "--inputFile=simple --numParts=4"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )

ENDIF()
##########################################


# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_order_tests
  SOURCE_FILES
    simple.mtx
    simple_coord.mtx
    simple_ordering.mtx
    Ashear.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

