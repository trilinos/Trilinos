
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA AND NOX_ENABLE_LOCA)
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Tcubed_Continuation_Example
    SOURCES 
    Basis.H
    Basis.C
    FiniteElementProblem.H
    FiniteElementProblem.C
    Problem_Interface.H
    Problem_Interface.C
    Example.C
    COMM serial mpi
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LOCA_Tcubed_TurningPoint_Example
    SOURCES 
    Basis.H
    Basis.C
    FiniteElementProblem.H
    FiniteElementProblem.C
    Problem_Interface.H
    Problem_Interface.C
    ExampleTPContinuation.C
    COMM serial mpi
    ARGS -v
    PASS_REGULAR_EXPRESSION "All tests passed"
    )

  IF(TPL_ENABLE_MF)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      LOCA_Tcubed_MultiContinuation_Example
      SOURCES 
      Basis.H
      Basis.C
      FiniteElementProblem.H
      FiniteElementProblem.C
      Problem_Interface.H
      Problem_Interface.C
      ExampleMultiContinuation.C       
      COMM serial mpi
      ARGS -v
      PASS_REGULAR_EXPRESSION "All tests passed"
      )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(
      LOCATcubedMultiContinuationCopyFiles
      SOURCE_FILES MFresultsDX.cfg MFresultsDX.net
      EXEDEPS LOCA_Tcubed_MultiContinuation_Example
      )

    ENDIF()

ENDIF()
