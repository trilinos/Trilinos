

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1Dfem
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemPrePostOperator.H
    1DfemPrePostOperator.C
    1Dfem.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  IF(NOX_ENABLE_ML)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      1DfemAndersonAcceleration
      SOURCES
      1DfemInterface.H
      1DfemInterface.C
      1DfemPrePostOperator.H
      1DfemPrePostOperator.C
      1DfemAndersonAcceleration.C
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "Test passed!"
      )
  ENDIF(NOX_ENABLE_ML)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemScaled
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemPrePostOperator.H
    1DfemPrePostOperator.C
    1DfemScaled.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemScaledInexactTrustRegion
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemPrePostOperator.H
    1DfemPrePostOperator.C
    1DfemScaledInexactTrustRegion.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemTransient
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemTransient.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  IF(NOX_ENABLE_ML)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      1DfemMLPreconditioner
      SOURCES
      1DfemInterface.H
      1DfemInterface.C
      1DfemMLPreconditioner.C
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "Test passed!"
      )
  ENDIF(NOX_ENABLE_ML)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecSolverPrecReuse
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecSolverPrecReuse.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecSolverPrecRebuild
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecSolverPrecRebuild.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecSolverPrecRecompute
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecSolverPrecRecompute.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecPrecRebuild
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecPrecRebuild.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecPrecReuse
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecPrecReuse.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemIfpackReuse
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemIfpackReuse.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemIfpackRecompute
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemIfpackRecompute.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemIfpackRebuild
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemIfpackRebuild.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemLeanMatrixFree
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemPrePostOperator.H
    1DfemPrePostOperator.C
    1DfemLeanMatrixFree.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemNonlinearCG
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemNonlinearCG.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    1DfemTensorBased
    SOURCES
    1DfemInterface.H
    1DfemInterface.C
    1DfemTensorBased.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  #  Need rest of tests...

  IF(NOX_ENABLE_EpetraExt)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      1DfemWithModelEvaluator
      SOURCES
      1DfemPrePostOperator.H
      1DfemPrePostOperator.C
      1DfemModelEvaluator.H
      1DfemModelEvaluator.C
      1DfemWithModelEvaluator.C
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "Test passed!"
      )
  ENDIF()

  IF(NOX_ENABLE_STRATIMIKOS_EPETRA_STACK)
    IF(NOX_ENABLE_Belos)
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        1DfemStratimikos
        SOURCES
        1DfemInterface.H
        1DfemInterface.C
        1DfemStratimikos.C
        COMM serial mpi
        PASS_REGULAR_EXPRESSION "Test passed!"
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        1DfemStratimikosReuse
        SOURCES
        1DfemInterface.H
        1DfemInterface.C
        1DfemStratimikosReuse.C
        COMM serial mpi
        PASS_REGULAR_EXPRESSION "Test passed!"
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        1DfemStratimikosInexact
        SOURCES
        1DfemInterface.H
        1DfemInterface.C
        1DfemStratimikosInexact.C
        COMM serial mpi
        PASS_REGULAR_EXPRESSION "Test passed!"
        )
    ENDIF(NOX_ENABLE_Belos)

    IF(NOX_ENABLE_ML)
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        1DfemStratimikosPrec
        SOURCES
        1DfemInterface.H
        1DfemInterface.C
        1DfemStratimikosPrec.C
        COMM serial mpi
        PASS_REGULAR_EXPRESSION "Test passed!"
        )
    ENDIF(NOX_ENABLE_ML)
  ENDIF(NOX_ENABLE_STRATIMIKOS_EPETRA_STACK)

ENDIF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)
