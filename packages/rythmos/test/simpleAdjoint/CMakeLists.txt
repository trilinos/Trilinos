
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(MultilineSet)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../../../epetraext/example/model_evaluator/DiagonalTransient)

FILE(GLOB SOURCES *.cpp)

TRIBITS_ADD_EXECUTABLE(
  simpleAdjoint
  SOURCES ${SOURCES}
  DEPLIBS rythmos_test_models
  COMM serial mpi  
  )

ASSERT_DEFINED(Stratimikos_ENABLE_Amesos)
IF (Stratimikos_ENABLE_Amesos)

  TRIBITS_ADD_TEST(
    simpleAdjoint
    NAME simpleAdjoint_amesos
    ARGS
      "--num-time-steps=10 --soln-verb-level=VERB_EXTREME"
      "--params-file=simpleAdjoint_n4.xml --num-time-steps=10 --soln-verb-level=VERB_EXTREME"
    STANDARD_PASS_OUTPUT
    COMM serial mpi  
    NUM_MPI_PROCS 1
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_simpleAdjoint_CopyFiles
    DEST_FILES
      simpleAdjoint_n4.xml
    SOURCE_PREFIX "_"
    EXEDEPS simpleAdjoint
    )

ENDIF()
