INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

IF (${PACKAGE_NAME}_ENABLE_Experimental)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    rcbTest
    SOURCES rcbTest.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()
  
IF (${PACKAGE_NAME}_ENABLE_Experimental)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pqJaggedTest
    SOURCES PQJaggedTest.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
    ARGS 
    "I=1.03 P=2,2 O=0 C=4 F=simple"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

#IF (${PACKAGE_NAME}_ENABLE_Experimental)
#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    taskMapper
#    SOURCES taskMappingTest.cpp
#    NUM_MPI_PROCS 1
#    COMM serial mpi
#    ARGS 
#    "PART=part.txt PROC=proc.txt"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#  )
#ENDIF()

#IF (${PACKAGE_NAME}_ENABLE_Experimental)
#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    taskMapper3
#    SOURCES taskMappingTest3.cpp
#    NUM_MPI_PROCS 1
#    COMM serial mpi
#    ARGS 
#    "PROC=machine.txt NX=8 NY=8 NZ=8"
#    PASS_REGULAR_EXPRESSION "PASS"
#    FAIL_REGULAR_EXPRESSION "FAIL"
#  )
#ENDIF()


IF (${PACKAGE_NAME}_ENABLE_Experimental)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    taskMapper2
    SOURCES taskMappingTest2.cpp
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS 
    "PART=part.txt PROC=proc.txt"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Experimental)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    parallelTaskMapper
    SOURCES parallelTaskMappingTest.cpp
    NUM_MPI_PROCS 1
    COMM serial mpi
    ARGS 
    "PART=part.txt"
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  blockTest
  SOURCES blockTest.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Partitioning1
  SOURCES partitioning1.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS 
   "--inputFile=simple"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Partitioning1_OneProc
  SOURCES partitioning1.cpp
  NUM_MPI_PROCS 4
  COMM serial mpi
  ARGS 
   "--inputFile=simple --no-distribute"
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAIL"
  )

# Command that copies files to the executable directory.
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_partition_tests
  SOURCE_FILES
    simple.mtx
    simple_coord.mtx
    USAir97.mtx
    USAir97_coord.mtx
    GeomGenParam.txt
    proc.txt
    part.txt
    machine.txt
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)

