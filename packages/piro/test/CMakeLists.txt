
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

      # Problems whose optimum is determined by the state part of the response (small
      # parameter penalty, inactive bounds), so that they check the transient state
      # derivatives, the time discretization and the final time used by the optimizer.
      TRIBITS_ADD_TEST(
        AnalysisDriverTranstientROL_MSD
        NAME AnalysisDriverTranstientROL_MSD_StateDriven
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_ReducedSpace_Transient_MSD_StateDriven.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
      )
      TRIBITS_ADD_TEST(
        AnalysisDriverTranstient_IR_ROL_MSD
        NAME AnalysisDriverTranstient_IR_ROL_MSD_StateDriven
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_ReducedSpace_Transient_Integrated_Response_MSD_StateDriven.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
      )
      TRIBITS_ADD_TEST(
        AnalysisDriverTranstientTempus_MSD
        NAME AnalysisDriverTranstientTempus_MSD_StateDriven
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_Tempus_Transient_MSD_StateDriven.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
      )
      # Same problem as AnalysisDriverTranstientROL_MSD_StateDriven without adjoint model evaluator:
      # the adjoint Jacobians are the transposes of the Jacobians.
      TRIBITS_ADD_TEST(
        AnalysisDriverTranstientROL_MSD
        NAME AnalysisDriverTranstientROL_MSD_StateDriven_ImplicitAdjoint
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_ReducedSpace_Transient_MSD_StateDriven_ImplicitAdjoint.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
      )
      # Two parameter vectors, ordered differently than in the model (indices 1 and 0)
      TRIBITS_ADD_TEST(
        AnalysisDriverTranstientROL_MSD
        NAME AnalysisDriverTranstientROL_2_parameters
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_ReducedSpace_Transient_2_parameters.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
      )
      TRIBITS_ADD_TEST(
        AnalysisDriverTranstientTempus_MSD
        NAME AnalysisDriverTranstientTempus_2_parameters
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_Tempus_Transient_2_parameters.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
      )
      # Time-integrated and time-sampled (sum over response times) objectives with the Tempus driver,
      # compared with the minimizers of the continuous problems, with SDIRK 2 Stage 2nd order (and
      # variable time steps landing on the response times). With Backward Euler, the O(dt) error of the
      # Tempus adjoint gradient, dominated for this model by the last adjoint step (at the initial time),
      # is too large near the optimum for the optimization to converge. A smaller enough dt for the test to 
      # converge would lead to test timeout. Backward Euler is tested with the ROL dynamic constraint,
      # whose gradient is the exact gradient of the discrete objective.
      TRIBITS_ADD_TEST(
        AnalysisDriverTranstientTempus_MSD
        NAME AnalysisDriverTranstientTempus_MSD_TimeIntegral_SDIRK
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_Tempus_Transient_MSD_TimeIntegral_SDIRK.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
      )
      TRIBITS_ADD_TEST(
        AnalysisDriverTranstientTempus_MSD
        NAME AnalysisDriverTranstientTempus_MSD_SumOverTimes_SDIRK_Variable
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_Tempus_Transient_MSD_SumOverTimes_SDIRK_Variable.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
      )
      # Sum over response times with the ROL dynamic-constraint approach
      TRIBITS_ADD_TEST(
        AnalysisDriverTranstientROL_MSD
        NAME AnalysisDriverTranstientROL_MSD_SumOverTimes
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_ReducedSpace_Transient_MSD_SumOverTimes.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
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
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    AnalysisDriverHDSA
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_HdsaLib
    EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
    SOURCES
    Main_AnalysisDriver_HDSA.cpp
    MockModelEval_H_Tpetra.cpp
    MockModelEval_H_Tpetra.hpp
    MatrixBased_LOWS.cpp
    MatrixBased_LOWS.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADDED_TESTS_NAMES_OUT AnalysisDriverHDSA_NAME 
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

  IF (AnalysisDriverHDASA_NAME)
    SET(AnalysisDriverHDSA_EXENAME AnalysisDriverHDSA)
  ELSE()
    SET(AnalysisDriverHDSA_EXENAME)
  ENDIF()

ENDIF(Piro_ENABLE_NOX)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ThyraProductME_Constraint_SimOpt_UnitTests
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
  EXCLUDE_IF_NOT_TRUE Stratimikos_ENABLE_Belos
  SOURCES
  Piro_ThyraProductME_Constraint_SimOpt_UnitTests.cpp
  ${TEUCHOS_STD_UNIT_TEST_MAIN}
  MockModelEval_A_Tpetra.cpp
  MockModelEval_A_Tpetra.hpp
  MatrixBased_LOWS.cpp
  MatrixBased_LOWS.hpp
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

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

  # The mass-spring-damper model has locally replicated maps: one process
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ThyraProductME_Tempus_FinalObjective_UnitTests
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
    EXCLUDE_IF_NOT_TRUE Stratimikos_ENABLE_Belos
    SOURCES
    Piro_ThyraProductME_Tempus_FinalObjective_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    MassSpringDamperModel.cpp
    MassSpringDamperModel.hpp
    MatrixBased_LOWS.cpp
    MatrixBased_LOWS.hpp
    NUM_MPI_PROCS 1
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
        input_Analysis_HDSA.xml
        input_Analysis_HDSA_continuation.xml
        input_Analysis_HDSA_OED.xml
        input_Analysis_ROL_Tempus_Transient_MSD.xml
        input_Analysis_ROL_Tempus_Transient_MSD_StateDriven.xml
        input_Analysis_ROL_Tempus_Transient_2_parameters.xml
        input_Analysis_ROL_Tempus_Transient_MSD_TimeIntegral_SDIRK.xml
        input_Analysis_ROL_Tempus_Transient_MSD_SumOverTimes_SDIRK_Variable.xml
        input_Analysis_ROL_ReducedSpace_LineSearch.xml
        input_Analysis_ROL_ReducedSpace_LineSearch_AdjointSensitivities_CheckGradients.xml
        input_Analysis_ROL_ReducedSpace_LineSearch_HessianBasedDotProduct.xml
        input_Analysis_ROL_ReducedSpace_Transient.xml
        input_Analysis_ROL_ReducedSpace_Transient_2_parameters.xml
        input_Analysis_ROL_ReducedSpace_Transient_MSD.xml
        input_Analysis_ROL_ReducedSpace_Transient_Integrated_Response_MSD.xml
        input_Analysis_ROL_ReducedSpace_Transient_Integrated_Response_MSD_StateDriven.xml
        input_Analysis_ROL_ReducedSpace_Transient_MSD_StateDriven.xml
        input_Analysis_ROL_ReducedSpace_Transient_MSD_StateDriven_ImplicitAdjoint.xml
        input_Analysis_ROL_ReducedSpace_Transient_MSD_SumOverTimes.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_HessianBasedDotProduct.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_HessianBasedDotProduct_2_parameters.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_BoundConstrained_NOXSolver.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_BoundConstrained_ExplicitAdjointME_NOXSolver.xml
        input_Analysis_ROL_FullSpace_AugmentedLagrangian_BoundConstrained.xml
        input_Analysis_ROL_FullSpace_AugmentedLagrangian.xml
        input_Analysis_ROL_ReducedSpace_LineSearch_BoundConstrained_HessianBasedDotProduct_Invalid.xml
        input_Tempus_BackwardEuler_SinCos.xml 
        input_tempus_be_nox_solver.yaml
        input_Solve_TR_Thyra.xml
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS ${ThyraSolverTpetra_EXENAME} ${AnalysisDriverTpetra_EXENAME}
)

