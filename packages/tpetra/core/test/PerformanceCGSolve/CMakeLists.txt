TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# The performance test currently depends on KokkosCore for its timers.
# Otherwise, it does not use anything else in KokkosCore directly.
#
IF (Tpetra_INST_DOUBLE)
    MESSAGE(STATUS "Tpetra: Performance-CGSolve test ENABLED")

    # This test makes sense only for Scalar == double
    TRIBITS_ADD_EXECUTABLE( Performance-CGSolve
                            SOURCES cg_solve_file.cpp
                            COMM mpi
                            CATEGORIES BASIC PERFORMANCE
    )


    set(COUNTER 1)
    set(MAX_COUNT 11)

    WHILE(COUNTER LESS MAX_COUNT)
      math(EXPR NP "${COUNTER} * ${COUNTER}")

      TRIBITS_ADD_TEST(
        Performance-CGSolve
        NAME Performance_StrongScaling_CGSolve
        ARGS "--size=200"
        COMM mpi
        NUM_MPI_PROCS ${NP}
        STANDARD_PASS_OUTPUT
        RUN_SERIAL
        CATEGORIES PERFORMANCE
      )

      math(EXPR COUNTER "${COUNTER} + 1")
    ENDWHILE()

    IF (Tpetra_ENABLE_CUDA)
      MESSAGE(STATUS "Tpetra: Performance-CGSolve CUDA_LAUNCH_BLOCKING test ENABLED")
      set(COUNTER 1)
      set(MAX_COUNT 11)

      WHILE(COUNTER LESS MAX_COUNT)
        math(EXPR NP "${COUNTER} * ${COUNTER}")

        TRIBITS_ADD_TEST(
          Performance-CGSolve
          NAME Performance_StrongScaling_CGSolve_CUDA_LAUNCH_BLOCKING
          ARGS "--size=200"
          COMM mpi
          NUM_MPI_PROCS ${NP}
          ENVIRONMENT CUDA_LAUNCH_BLOCKING=1
          STANDARD_PASS_OUTPUT
          RUN_SERIAL
          CATEGORIES PERFORMANCE
        )

        math(EXPR COUNTER "${COUNTER} + 1")
      ENDWHILE()

    ENDIF()

ENDIF()
