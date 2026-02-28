TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER ${PANZER_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_Matrix_Free
  SOURCES test_fe_matrix_free.cpp
  NUM_MPI_PROCS 2
  COMM mpi
  PASS_REGULAR_EXPRESSION "TEST PASSED"
)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/stratimikos-chebyshev.xml"
               "${CMAKE_CURRENT_BINARY_DIR}/stratimikos-chebyshev.xml" COPYONLY)
               
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/stratimikos-jacobi.xml"
               "${CMAKE_CURRENT_BINARY_DIR}/stratimikos-jacobi.xml" COPYONLY)

# NVR: the below is not working for me; hence the above.
# TRIBITS_COPY_FILES_TO_BINARY_DIR(mf_cp
#   SOURCE_FILES
#     stratimikos-chebyshev.xml
#     stratimikos-jacobi.xml
#   SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
#   DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
# )