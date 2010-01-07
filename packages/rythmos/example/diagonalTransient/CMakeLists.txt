
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)
INCLUDE(MultilineSet)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

FILE(GLOB SOURCES *.cpp)

PACKAGE_ADD_EXECUTABLE(
  diagonalTransient
  SOURCES ${SOURCES}
  DEPLIBS rythmos_test_models
  COMM serial mpi  
  )

ASSERT_DEFINED(Stratimikos_ENABLE_Amesos)
IF (Stratimikos_ENABLE_Amesos)

  SET(STANDARD_TEST_CONFIG
    COMM serial mpi  
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_default
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --extra-params-file=ExtraParams2.xml"
    " --final-time=1e-3 --num-time-steps=10 --verb-level=VERB_MEDIUM "
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_amesos_base
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )
    
  MULTILINE_SET(ARGS
    " --extra-params-file=ExtraParams1.xml"
    "  --final-time=1e-3 "
    "  --num-time-steps=10 "
    "  --verb-level=VERB_MEDIUM "
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_amesos_coeff_s_pi
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --params-file=diagonalTransient.amesos.fwd-sens-solve.n16.p4.xml"
    " --final-time=1e-4 --num-time-steps=10 --verb-level=VERB_LOW"
    " --echo-command-line"
    " --fwd-sens-solve"
    " --max-state-error=1e-6"
    " --max-sens-error=1e-6"
    " --dump-final-solutions"
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_amesos_fwd-sens_n16_p4
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --params-file=diagonalTransient.amesos.fwd-sens-solve.n1.p1.xml "
    " --final-time=1e-4 "
    " --num-time-steps=10 "
    " --verb-level=VERB_LOW "
    " --echo-command-line "
    " --fwd-sens-solve "
    " --max-state-error=1e-6 "
    " --max-sens-error=1e-6 "
    " --dump-final-solutions "
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_amesos_fwd-sens_n1_p1
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --params-file=diagonalTransient.amesos.fwd-sens-solve.n5.p1.xml"
    " --final-time=1e-4 --num-time-steps=10 --verb-level=VERB_LOW"
    " --echo-command-line"
    " --fwd-sens-solve"
    " --max-state-error=1e-6"
    " --max-sens-error=1e-6"
    " --dump-final-solutions"
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_amesos_fwd-sens_n5_p1
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --params-file=diagonalTransient.amesos.fwd-BDF-fixed.n1.xml"
    " --final-time=1e-3 --num-time-steps=10 --verb-level=VERB_EXTREME"
    " --echo-command-line"
    " --max-state-error=1e-6"
    " --dump-final-solutions"
    " --use-BDF"
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_amesos_fwd-BDF-fixed
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --params-file=diagonalTransient.amesos.fwd-BDF-variable.n1.xml"
    " --final-time=1e-4 --num-time-steps=-1 --verb-level=VERB_EXTREME"
    " --echo-command-line"
    " --max-state-error=5e-6"
    " --dump-final-solutions"
    " --use-BDF"
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_amesos_fwd-BDF-variable
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  MULTILINE_SET(ARGS
    " --params-file=diagonalTransient.amesos.fwd-sens-solve-BDF.n1.p1.xml"
    " --final-time=1e-4 --num-time-steps=10 --verb-level=VERB_LOW"
    " --echo-command-line"
    " --fwd-sens-solve"
    " --max-state-error=1e-6"
    " --max-restate-error=1e-8"
    " --max-sens-error=1e-6"
    " --dump-final-solutions"
    " --use-BDF"
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_amesos_fwd-sens-BDF
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(Example_diagonalTransient_Amesos_CopyFiles
    DEST_FILES
      ExtraParams1.xml
      diagonalTransient.amesos.fwd-sens-solve.n16.p4.xml
      diagonalTransient.amesos.fwd-sens-solve.n1.p1.xml
      diagonalTransient.amesos.fwd-sens-solve.n5.p1.xml
      diagonalTransient.amesos.fwd-BDF-fixed.n1.xml
      diagonalTransient.amesos.fwd-BDF-variable.n1.xml
      diagonalTransient.amesos.fwd-sens-solve-BDF.n1.p1.xml
    SOURCE_PREFIX "_"
    EXEDEPS diagonalTransient
    )

ENDIF()


ASSERT_DEFINED(Stratimikos_ENABLE_AztecOO)
IF (Stratimikos_ENABLE_AztecOO)

  SET(STANDARD_TEST_CONFIG
    COMM serial mpi  
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
    )

  MULTILINE_SET(ARGS
    " --params-file=diagonalTransient_aztecoo.xml"
    " --final-time=1e-4 --num-time-steps=-1 --verb-level=VERB_EXTREME"
    " --echo-command-line"
    " --max-state-error=5e-6"
    " --dump-final-solutions"
    " --use-BDF"
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_aztecoo
    ARGS ${ARGS}
    ${STANDARD_TEST_CONFIG}
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(Example_diagonalTransient_AztecOO_CopyFiles
    DEST_FILES
      ExtraParams2.xml
      diagonalTransient_aztecoo.xml
    SOURCE_PREFIX "_"
    EXEDEPS diagonalTransient
    )

ENDIF()


ASSERT_DEFINED(Stratimikos_ENABLE_Belos)
IF (Stratimikos_ENABLE_Belos)

  MULTILINE_SET(ARGS
    " --params-file=diagonalTransient.belos.fwd-BE.n2.p2.xml"
    " --final-time=1e-3 --num-time-steps=10 --verb-level=VERB_EXTREME"
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_belos_BE
    ARGS ${ARGS}
    COMM serial mpi  
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
    )

  MULTILINE_SET(ARGS
    " --params-file=diagonalTransient.belos.fwd-IRK.n2.p2.xml"
    " --use-IRK"
    )
  PACKAGE_ADD_TEST(
    diagonalTransient
    NAME diagonalTransient_belos_IRK
    ARGS ${ARGS}
    COMM serial mpi  
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(Example_diagonalTransient_Belos_CopyFiles
    DEST_FILES
      diagonalTransient.belos.fwd-BE.n2.p2.xml
      diagonalTransient.belos.fwd-IRK.n2.p2.xml
    SOURCE_PREFIX "_"
    EXEDEPS diagonalTransient
    )

ENDIF()
