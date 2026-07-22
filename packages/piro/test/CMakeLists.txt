
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (Piro_ENABLE_NOX)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ThyraSolverTpetra
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
    SOURCES
    Main_ThyraSolver_Tpetra.cpp
    MockModelEval_A_Tpetra.cpp
    MockModelEval_A_Tpetra.hpp
    MatrixBased_LOWS.cpp
    MatrixBased_LOWS.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADDED_EXE_TARGET_NAME_OUT ThyraSolver_NAME
    ADDED_TESTS_NAMES_OUT ThyraSolverTpetra_NAME
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
      SecondOrderIntegratorThyra
      SOURCES
      Main_SecondOrderIntegrator_Thyra.cpp
      MockModelEval_C_Tpetra.cpp
      MockModelEval_C_Tpetra.hpp
      Piro_Test_MockObserver.hpp
      MatrixBased_LOWS.hpp
      MatrixBased_LOWS.cpp
      NUM_MPI_PROCS 1
      ARGS -v
      PASS_REGULAR_EXPRESSION "TEST PASSED"
    )

  IF (Piro_ENABLE_Tempus)
    #IF (TPL_ENABLE_MPI)
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        AnalysisDriverTranstientROL_MSD
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        SOURCES
        Main_AnalysisDriver_Tempus.cpp
        MockModelEval_A_Tpetra.cpp
        MockModelEval_A_Tpetra.hpp
        MockModelEval_B_Tpetra.cpp
        MockModelEval_B_Tpetra.hpp
        MockModelEval_B_Tpetra_2_parameters.cpp
        MockModelEval_B_Tpetra_2_parameters.hpp
        MatrixBased_LOWS.cpp
        MatrixBased_LOWS.hpp
        MassSpringDamperModel.cpp
        MassSpringDamperModel.hpp
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_ReducedSpace_Transient_MSD.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
        ADDED_TESTS_NAMES_OUT AnalysisDriverTranstientROL_MSD_NAME 
      )
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        AnalysisDriverTranstient_IR_ROL_MSD
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        SOURCES
        Main_AnalysisDriver_Tempus.cpp
        MockModelEval_A_Tpetra.cpp
        MockModelEval_A_Tpetra.hpp
        MockModelEval_B_Tpetra.cpp
        MockModelEval_B_Tpetra.hpp
        MockModelEval_B_Tpetra_2_parameters.cpp
        MockModelEval_B_Tpetra_2_parameters.hpp
        MatrixBased_LOWS.cpp
        MatrixBased_LOWS.hpp
        MassSpringDamperModel.cpp
        MassSpringDamperModel.hpp
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_ReducedSpace_Transient_Integrated_Response_MSD.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
        ADDED_TESTS_NAMES_OUT AnalysisDriverTranstient_IR_ROL_MSD_NAME 
      )
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        AnalysisDriverTranstientTempus_MSD
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        SOURCES
        Main_AnalysisDriver_Tempus.cpp
        MockModelEval_A_Tpetra.cpp
        MockModelEval_A_Tpetra.hpp
        MockModelEval_B_Tpetra.cpp
        MockModelEval_B_Tpetra.hpp
        MockModelEval_B_Tpetra_2_parameters.cpp
        MockModelEval_B_Tpetra_2_parameters.hpp
        MatrixBased_LOWS.cpp
        MatrixBased_LOWS.hpp
        MassSpringDamperModel.cpp
        MassSpringDamperModel.hpp
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_Tempus_Transient_MSD.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
        ADDED_TESTS_NAMES_OUT AnalysisDriverTranstientTempus_MSD_NAME 
      )
    #ENDIF()
    #TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #  AnalysisDriverTempus
    #  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
    #  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
    #  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
    #  EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
    #  SOURCES
    #  Main_AnalysisDriver_Tempus.cpp
    #  MockModelEval_A_Tpetra.cpp
    #  MockModelEval_A_Tpetra.hpp
    #  MockModelEval_B_Tpetra.cpp
    #  MockModelEval_B_Tpetra.hpp
    #  MockModelEval_B_Tpetra_2_parameters.cpp
    #  MockModelEval_B_Tpetra_2_parameters.hpp
    #  MatrixBased_LOWS.cpp
    #  MatrixBased_LOWS.hpp
    #  MassSpringDamperModel.cpp
    #  MassSpringDamperModel.hpp
    #  NUM_MPI_PROCS 1-4
    #  ARGS -v
    #  PASS_REGULAR_EXPRESSION "TEST PASSED"
    #  ADDED_TESTS_NAMES_OUT AnalysisDriverTempus_NAME 
    #)
  ENDIF()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    AnalysisDriverTpetra
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
    EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
    SOURCES
    Main_AnalysisDriver_Tpetra.cpp
    MockModelEval_A_Tpetra.cpp
    MockModelEval_A_Tpetra.hpp
    MockModelEval_B_Tpetra.cpp
    MockModelEval_B_Tpetra.hpp
    MockModelEval_B_Tpetra_2_parameters.cpp
    MockModelEval_B_Tpetra_2_parameters.hpp
    MatrixBased_LOWS.cpp
    MatrixBased_LOWS.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADDED_TESTS_NAMES_OUT AnalysisDriverTpetra_NAME 
  )
  
  IF (ThyraSolverTpetra_NAME)
    SET(ThyraSolverTpetra_EXENAME ThyraSolverTpetra)
  ELSE()
    SET(ThyraSolverTpetra_EXENAME)
  ENDIF()

  IF (AnalysisDriverTpetra_NAME)
    SET(AnalysisDriverTpetra_EXENAME AnalysisDriverTpetra)
  ELSE()
    SET(AnalysisDriverTpetra_EXENAME)
  ENDIF()

ENDIF(Piro_ENABLE_NOX)

IF (Piro_ENABLE_Tempus)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolver_UnitTests_Tpetra
    SOURCES
    Piro_TempusSolver_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_MockObserver.hpp
    Piro_Test_ThyraSupport.hpp
    Piro_Test_WeakenedModelEvaluator.hpp
    MockModelEval_A_Tpetra.hpp
    MockModelEval_A_Tpetra.cpp
    MatrixBased_LOWS.cpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolverForwardOnly_UnitTests_Tpetra
    SOURCES
    Piro_TempusSolverForwardOnly_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    MockModelEval_A_Tpetra.hpp
    MockModelEval_A_Tpetra.cpp
    MockModelEval_B_Tpetra.hpp
    MockModelEval_B_Tpetra.cpp
    MatrixBased_LOWS.cpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolver_SensitivitySinCos_Combined_FSA_UnitTests
    SOURCES
    Piro_TempusSolver_SensitivitySinCos_Combined_FSA_UnitTests.cpp
    Piro_TempusSolver_SensitivitySinCos_FSA.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    SinCosModel.cpp
    SinCosModel.hpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolver_SensitivitySinCos_Combined_FSA_Tangent_UnitTests
    SOURCES
    Piro_TempusSolver_SensitivitySinCos_Combined_FSA_Tangent_UnitTests.cpp
    Piro_TempusSolver_SensitivitySinCos_FSA.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    SinCosModel.cpp
    SinCosModel.hpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolver_SensitivitySinCos_Staggered_FSA_UnitTests
    SOURCES
    Piro_TempusSolver_SensitivitySinCos_Staggered_FSA_UnitTests.cpp
    Piro_TempusSolver_SensitivitySinCos_FSA.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    SinCosModel.cpp
    SinCosModel.hpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolver_SensitivitySinCos_Staggered_FSA_Tangent_UnitTests
    SOURCES
    Piro_TempusSolver_SensitivitySinCos_Staggered_FSA_Tangent_UnitTests.cpp
    Piro_TempusSolver_SensitivitySinCos_FSA.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    SinCosModel.cpp
    SinCosModel.hpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolver_AdjointSensitivitySinCos_Explicit_UnitTests
    SOURCES
    Piro_TempusSolver_AdjointSensitivitySinCos_Explicit_UnitTests.cpp
    Piro_TempusSolver_AdjointSensitivitySinCosUnitTests.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    SinCosModel.cpp
    SinCosModel.hpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolver_AdjointSensitivitySinCos_Implicit_UnitTests
    SOURCES
    Piro_TempusSolver_AdjointSensitivitySinCos_Implicit_UnitTests.cpp
    Piro_TempusSolver_AdjointSensitivitySinCosUnitTests.hpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    SinCosModel.cpp
    SinCosModel.hpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

ENDIF (Piro_ENABLE_Tempus)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copyTestInputFiles
  DEST_FILES   
        input_Solve_LOCA_Thyra.xml
        input_Solve_NOX_ForwardSensitivities_MatrixFree_Thyra.xml
        input_Solve_NOX_ForwardSensitivities_Thyra.xml
        input_Solve_NOX_AdjointSensitivities_ImplicitAdjointME_Thyra.xml
        input_Solve_NOX_AdjointSensitivities_ExplicitAdjointME_Thyra.xml
        input_Analysis_ROL_Tempus_Transient_MSD.xml
        input_Analysis_ROL_ReducedSpace_LineSearch.xml
        input_Analysis_ROL_ReducedSpace_LineSearch_AdjointSensitivities_CheckGradients.xml
        input_Analysis_ROL_ReducedSpace_LineSearch_HessianBasedDotProduct.xml
        input_Analysis_ROL_ReducedSpace_Transient.xml
        input_Analysis_ROL_ReducedSpace_Transient_2_parameters.xml
        input_Analysis_ROL_ReducedSpace_Transient_MSD.xml
        input_Analysis_ROL_ReducedSpace_Transient_Integrated_Response_MSD.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_HessianBasedDotProduct.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_HessianBasedDotProduct_2_parameters.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_BoundConstrained_NOXSolver.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_BoundConstrained_ExplicitAdjointME_NOXSolver.xml
        input_Analysis_ROL_FullSpace_AugmentedLagrangian_BoundConstrained.xml
        input_Tempus_BackwardEuler_SinCos.xml 
        input_tempus_be_nox_solver.yaml
        input_Solve_TR_Thyra.xml
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS ${ThyraSolverTpetra_EXENAME} ${AnalysisDriverTpetra_EXENAME}
)

