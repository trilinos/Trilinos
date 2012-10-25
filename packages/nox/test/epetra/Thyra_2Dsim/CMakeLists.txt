
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../utils)

ASSERT_DEFINED(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA
  NOX_ENABLE_Epetra NOX_ENABLE_EpetraExt
  NOX_ENABLE_Stratimikos Stratimikos_ENABLE_AztecOO)

SET(UNIT_TEST_DRIVER ${TEUCHOS_STD_UNIT_TEST_MAIN})

IF(NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA AND
    NOX_ENABLE_THYRA_EPETRA_ADAPTERS AND
    NOX_ENABLE_Epetra AND NOX_ENABLE_EpetraExt AND
    NOX_ENABLE_Stratimikos AND Stratimikos_ENABLE_AztecOO)
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_NonlinearSolver_NOX_2Dsim
    SOURCES EpetraModelEval2DSim.H EpetraModelEval2DSim.C 
            Thyra_NonlinearSolver_NOX_2Dsim.C
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test passed!"
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_2Dsim
    SOURCES EpetraModelEval2DSim.H EpetraModelEval2DSim.C Thyra_2Dsim.C
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test passed!"
    )
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_RowSumScaling
    SOURCES EpetraModelEval2DSim.H EpetraModelEval2DSim.C Thyra_RowSumScaling.C
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test passed!"
    )
   
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Thyra_PseudoTransient
    SOURCES  ${UNIT_TEST_DRIVER} Thyra_PseudoTransient.cpp NOX_RosenbrockModelEvaluator.hpp NOX_RosenbrockModelEvaluator.cpp
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test passed!"
    )
  
ENDIF()
