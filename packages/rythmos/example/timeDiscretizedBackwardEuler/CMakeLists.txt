
INCLUDE(MultilineSet)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../../../epetraext/example/model_evaluator/DiagonalTransient)

FILE(GLOB SOURCES *.cpp)

TRIBITS_ADD_EXECUTABLE(
  timeDiscretizedBackwardEuler
  SOURCES ${SOURCES}
  TESTONLYLIBS rythmos_test_models
  COMM serial mpi  
  )

ASSERT_DEFINED(Stratimikos_ENABLE_Amesos)
IF (Stratimikos_ENABLE_Amesos)

  SET(STANDARD_TEST_CONFIG
    COMM serial mpi  
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  MULTILINE_SET(ARGS
    " --params-file=timeDiscretizedBackwardEuler1.xml"
    " --verb-level=VERB_EXTREME"
    " --num-time-steps=10"
    " --final-time=1e-3"
    )
  TRIBITS_ADD_TEST(
    timeDiscretizedBackwardEuler
    NAME timeDiscretizedBackwardEuler_amesos
    ARGS ${ARGS}
    XHOSTTYPE AIX
    ${STANDARD_TEST_CONFIG}
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Example_timeDiscretizedBackwardEuler_Amesos_CopyFiles
    DEST_FILES
      timeDiscretizedBackwardEuler1.xml
    SOURCE_PREFIX "_"
    EXEDEPS timeDiscretizedBackwardEuler
    )

ENDIF()


