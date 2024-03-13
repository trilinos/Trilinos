TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SolutionHistory
  SOURCES Tempus_UnitTest_SolutionHistory.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeStepControl
  SOURCES Tempus_UnitTest_TimeStepControl.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeStepControlStrategyConstant
  SOURCES Tempus_UnitTest_TimeStepControlStrategyConstant.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeStepControlStrategyBasicVS
  SOURCES Tempus_UnitTest_TimeStepControlStrategyBasicVS.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeStepControlStrategyIntegralController
  SOURCES Tempus_UnitTest_TimeStepControlStrategyIntegralController.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeStepControlStrategyComposite
  SOURCES Tempus_UnitTest_TimeStepControlStrategyComposite.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_BackwardEuler
  SOURCES Tempus_UnitTest_BackwardEuler.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_BDF2 
  SOURCES Tempus_UnitTest_BDF2.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_BackwardEuler
  SOURCES Tempus_UnitTest_DIRK_BackwardEuler.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_1Stage1stOrderRadauIA
  SOURCES Tempus_UnitTest_DIRK_1Stage1stOrderRadauIA.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_1StageTheta
  SOURCES Tempus_UnitTest_DIRK_1StageTheta.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_2Stage2ndOrderLobattoIIIB
  SOURCES Tempus_UnitTest_DIRK_2Stage2ndOrderLobattoIIIB.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_General
  SOURCES Tempus_UnitTest_DIRK_General.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_EDIRK_2Stage3rdOrder
  SOURCES Tempus_UnitTest_EDIRK_2Stage3rdOrder.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_EDIRK_2StageTheta
  SOURCES Tempus_UnitTest_EDIRK_2StageTheta.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_EDIRK_TrapezoidalRule
  SOURCES Tempus_UnitTest_EDIRK_TrapezoidalRule.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_21Pair
  SOURCES Tempus_UnitTest_SDIRK_21Pair.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_2Stage2ndOrder
  SOURCES Tempus_UnitTest_SDIRK_2Stage2ndOrder.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_2Stage3rdOrder
  SOURCES Tempus_UnitTest_SDIRK_2Stage3rdOrder.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_3Stage2ndOrder
  SOURCES Tempus_UnitTest_SDIRK_3Stage2ndOrder.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_3Stage4thOrder
  SOURCES Tempus_UnitTest_SDIRK_3Stage4thOrder.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_5Stage4thOrder
  SOURCES Tempus_UnitTest_SDIRK_5Stage4thOrder.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_5Stage5thOrder
  SOURCES Tempus_UnitTest_SDIRK_5Stage5thOrder.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_ImplicitMidpoint
  SOURCES Tempus_UnitTest_SDIRK_ImplicitMidpoint.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SSPDIRK22
  SOURCES Tempus_UnitTest_SSPDIRK22.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SSPDIRK23
  SOURCES Tempus_UnitTest_SSPDIRK23.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SSPDIRK32
  SOURCES Tempus_UnitTest_SSPDIRK32.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SSPDIRK33
  SOURCES Tempus_UnitTest_SSPDIRK33.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_3_8Rule
  SOURCES Tempus_UnitTest_ERK_3_8Rule.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_3Stage3rdOrder
  SOURCES Tempus_UnitTest_ERK_3Stage3rdOrder.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_3Stage3rdOrderHeun
  SOURCES Tempus_UnitTest_ERK_3Stage3rdOrderHeun.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_3Stage3rdOrderTVD
  SOURCES Tempus_UnitTest_ERK_3Stage3rdOrderTVD.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_4Stage3rdOrderRunge
  SOURCES Tempus_UnitTest_ERK_4Stage3rdOrderRunge.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_4Stage4thOrder
  SOURCES Tempus_UnitTest_ERK_4Stage4thOrder.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_5Stage3rdOrderKandG
  SOURCES Tempus_UnitTest_ERK_5Stage3rdOrderKandG.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_BogackiShampine32
  SOURCES Tempus_UnitTest_ERK_BogackiShampine32.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_ForwardEuler
  SOURCES Tempus_UnitTest_ERK_ForwardEuler.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_General
  SOURCES Tempus_UnitTest_ERK_General.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_Merson45
  SOURCES Tempus_UnitTest_ERK_Merson45.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_Midpoint
  SOURCES Tempus_UnitTest_ERK_Midpoint.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_Trapezoidal
  SOURCES Tempus_UnitTest_ERK_Trapezoidal.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ForwardEuler
  SOURCES Tempus_UnitTest_ForwardEuler.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_HHTAlpha
  SOURCES Tempus_UnitTest_HHTAlpha.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_IMEX_RK
  SOURCES Tempus_UnitTest_IMEX_RK.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_IMEX_RK_Partition
  SOURCES Tempus_UnitTest_IMEX_RK_Partition.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Leapfrog
  SOURCES Tempus_UnitTest_Leapfrog.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_NewmarkExplicitAForm
  SOURCES Tempus_UnitTest_NewmarkExplicitAForm.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_NewmarkImplicitAForm
  SOURCES Tempus_UnitTest_NewmarkImplicitAForm.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_NewmarkImplicitDForm
  SOURCES Tempus_UnitTest_NewmarkImplicitDForm.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_OperatorSplit
  SOURCES Tempus_UnitTest_OperatorSplit.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Trapezoidal
  SOURCES Tempus_UnitTest_Trapezoidal.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Subcycling
  SOURCES Tempus_UnitTest_Subcycling.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeEventBase
  SOURCES Tempus_UnitTest_TimeEventBase.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeEventRange
  SOURCES Tempus_UnitTest_TimeEventRange.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeEventList
  SOURCES Tempus_UnitTest_TimeEventList.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeEventRangeIndex
  SOURCES Tempus_UnitTest_TimeEventRangeIndex.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeEventListIndex
  SOURCES Tempus_UnitTest_TimeEventListIndex.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeEventComposite
  SOURCES Tempus_UnitTest_TimeEventComposite.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK
  SOURCES Tempus_UnitTest_ERK.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_NumericalUtils
  SOURCES Tempus_UnitTest_NumericalUtils.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Stepper_ErrorNorm
  SOURCES Tempus_UnitTest_Stepper_ErrorNorm.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_IntegratorBasic
  SOURCES Tempus_UnitTest_IntegratorBasic.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTest_BackwardEuler_CopyFiles
  DEST_FILES Tempus_BackwardEuler_SinCos.xml
  EXEDEPS UnitTest_BackwardEuler
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTest_BDF2_CopyFiles
  DEST_FILES Tempus_BDF2_SinCos.xml
  EXEDEPS UnitTest_BDF2
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTest_DIRK_BackwardEuler_CopyFiles
  DEST_FILES Tempus_DIRK_VanDerPol.xml
  EXEDEPS UnitTest_DIRK_BackwardEuler
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTest_NewmarkExplicitAForm_CopyFiles
  DEST_FILES Tempus_NewmarkExplicitAForm_HarmonicOscillator_Damped.xml
  EXEDEPS UnitTest_NewmarkExplicitAForm
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTest_NewmarkImplicitAForm_CopyFiles
  DEST_FILES Tempus_NewmarkImplicitAForm_HarmonicOscillator_Damped_SecondOrder.xml
  EXEDEPS UnitTest_NewmarkImplicitAForm
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTest_NewmarkImplicitDForm_CopyFiles
  DEST_FILES Tempus_NewmarkImplicitDForm_HarmonicOscillator_Damped_SecondOrder.xml
  EXEDEPS UnitTest_NewmarkImplicitDForm
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(UnitTest_IntegratorBasic_CopyFiles
  DEST_FILES Tempus_IntegratorBasic_default.xml Tempus_IntegratorBasic_ref.xml
  EXEDEPS UnitTest_IntegratorBasic
  )
