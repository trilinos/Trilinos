
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemTransient
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemTransient.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemMLPreconditioner
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemMLPreconditioner.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecSolverPrecReuse
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecSolverPrecReuse.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecSolverPrecRebuild
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecSolverPrecRebuild.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecSolverPrecRecompute
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecSolverPrecRecompute.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecPrecRebuild
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecPrecRebuild.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemAztecPrecReuse
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemAztecPrecReuse.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemIfpackReuse
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C 
    1DfemIfpackReuse.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemIfpackRecompute
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemIfpackRecompute.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemIfpackRebuild
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemIfpackRebuild.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    1DfemNonlinearCG
    SOURCES 
    1DfemInterface.H
    1DfemInterface.C
    1DfemNonlinearCG.C
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Test passed!"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
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
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
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

  IF(NOX_ENABLE_Stratimikos)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      1DfemStratimikos
      SOURCES 
      1DfemInterface.H
      1DfemInterface.C
      1DfemStratimikos.C
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "Test passed!"
      )

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      1DfemStratimikosPrec
      SOURCES 
      1DfemInterface.H
      1DfemInterface.C
      1DfemStratimikosPrec.C
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "Test passed!"
      )
  ENDIF()
  
  
  IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_BELOS)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      1DfemBelos
      SOURCES 
      1DfemInterface.H
      1DfemInterface.C
      1DfemBelos.C
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "Test passed!"
      )
  ENDIF()

ENDIF()
