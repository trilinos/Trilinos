
SET(HEADERS "")
SET(SOURCES "")


SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

# The performance test currently depends on KokkosCore for its timers.
# Otherwise, it does not use anything else in KokkosCore directly.
#
IF (Tpetra_INST_DOUBLE)
    MESSAGE(STATUS "Tpetra: Performance-CGSolve test ENABLED")

    # This test makes sense only for Scalar == double
    TRIBITS_ADD_EXECUTABLE( Performance-CGSolve
                            SOURCES ${SOURCES}
                            COMM mpi
                            CATEGORIES BASIC PERFORMANCE
    )

    TRIBITS_ADD_TEST(
      Performance-CGSolve
      NAME Performance_StrongScaling_CGSolve
      ARGS "--size=200"
      COMM mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    TRIBITS_ADD_TEST(
      Performance-CGSolve
      NAME Performance_StrongScaling_CGSolve
      ARGS "--size=200"
      COMM mpi
      NUM_MPI_PROCS 4
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    TRIBITS_ADD_TEST(
      Performance-CGSolve
      NAME Performance_StrongScaling_CGSolve
      ARGS "--size=200"
      COMM mpi
      NUM_MPI_PROCS 9
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    TRIBITS_ADD_TEST(

      Performance-CGSolve
      NAME Performance_StrongScaling_CGSolve
      ARGS "--size=200"
      COMM mpi
      NUM_MPI_PROCS 16
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    TRIBITS_ADD_TEST(
      Performance-CGSolve
      NAME Performance_StrongScaling_CGSolve
      ARGS "--size=200"
      COMM mpi
      NUM_MPI_PROCS 25
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

ENDIF()

#TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyXmlForTpetraPerfCgSolve
#  SOURCE_DIR ${Tpetra_SOURCE_DIR}/test/
#  SOURCE_FILES Tpetra_PerformanceTests.xml
#  EXEDEPS Performance-CGSolve
#)


