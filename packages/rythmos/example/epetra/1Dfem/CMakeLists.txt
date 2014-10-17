
INCLUDE(MultilineSet)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

FILE(GLOB SOURCES *.cpp *.C)

TRIBITS_ADD_EXECUTABLE(
  1DfemTransient
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
    " --verbose "
    " --echo-command-line "
    )

  TRIBITS_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_amesos_BE
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --echo-command-line "
    )

  TRIBITS_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_amesos_IBDF
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --useNOX "
    " --echo-command-line "
    )

  TRIBITS_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_amesos_IBDF_NOX
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

ENDIF()


ASSERT_DEFINED(Stratimikos_ENABLE_AztecOO)
IF (Stratimikos_ENABLE_AztecOO)

  SET(STANDARD_TEST_CONFIG
    COMM serial mpi  
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --extra-linear-solver-params-file=Extra_AztecOO_Params.xml"
    " --linear-solver-params-used-file=lowsf.aztecoo.used.xml "
    " --echo-command-line "
    )

  TRIBITS_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_aztecoo_BE
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --extra-linear-solver-params-file=Extra_AztecOO_Params.xml"
    " --linear-solver-params-used-file=lowsf.aztecoo.used.xml "
    " --echo-command-line "
    )

  TRIBITS_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_aztecoo_IBDF
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

#  MULTILINE_SET(ARGS
#    " --verbose "
#    " --method=IRK "
#    " --numsteps=5 "
#    " --extra-linear-solver-params-file=Extra_AztecOO_Params.xml"
#    " --linear-solver-params-used-file=lowsf.aztecoo.used.xml "
#    " --echo-command-line "
#    )
#
#  TRIBITS_ADD_TEST(
#    1DfemTransient
#    NAME 1DfemTransient_aztecoo_IRK
#    ARGS ${ARGS}
#    ${STANDARD_TEST_CONFIG}
#    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --useNOX "
    " --extra-linear-solver-params-file=Extra_AztecOO_Params.xml"
    " --linear-solver-params-used-file=lowsf.aztecoo.used.xml "
    " --echo-command-line "
    )

  TRIBITS_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_aztecoo_IBDF_NOX
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Example_1DfemTransient_AztecOO_CopyFiles
    DEST_FILES
      Extra_AztecOO_Params.xml
      Extra_Belos_Params.xml
    SOURCE_PREFIX "_"
    EXEDEPS 1DfemTransient
    )
ENDIF()


ASSERT_DEFINED(Stratimikos_ENABLE_Belos)
IF (Stratimikos_ENABLE_Belos)

  SET(STANDARD_TEST_CONFIG
    COMM serial mpi  
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --extra-linear-solver-params-file=Extra_Belos_Params.xml"
    "  --linear-solver-params-used-file=lowsf.belos.used.xml "
    )
  TRIBITS_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_belos_BE
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --extra-linear-solver-params-file=Extra_Belos_Params.xml"
    "  --linear-solver-params-used-file=lowsf.belos.used.xml "
    )
  TRIBITS_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_belos_IBDF
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --verbose "
    " --method=BDF "
    " --useNOX "
    " --extra-linear-solver-params-file=Extra_Belos_Params.xml"
    "  --linear-solver-params-used-file=lowsf.belos.used.xml "
    )
  TRIBITS_ADD_TEST(
    1DfemTransient
    NAME 1DfemTransient_belos_IBDF_NOX
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

ENDIF()
