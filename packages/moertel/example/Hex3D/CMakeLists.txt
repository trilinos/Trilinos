INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)

IF(HAVE_MOERTEL_MPI)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Hex_3D
  SOURCES Hex_3D.cpp
  ARGS 
  COMM mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )

ELSE()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Hex_3D
  SOURCES Hex_3D.cpp
  ARGS 
  COMM serial
  PASS_REGULAR_EXPRESSION "Test passed!"
  XHOST trilinos-test.sandia.gov trilinos-test2.sandia.gov zan.sandia.gov
  )

ENDIF()


TRIBITS_COPY_FILES_TO_BINARY_DIR(
   Hex_3DCopyFiles
   DEST_FILES Hex_3D.grid
   SOURCE_DIR ${PACKAGE_SOURCE_DIR}/example/Hex3D
   EXEDEPS Hex_3D
   )
