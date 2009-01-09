
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
