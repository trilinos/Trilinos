INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_TimeStepControl
  SOURCES Tempus_UnitTest_TimeStepControl.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_BDF2 
  SOURCES Tempus_UnitTest_BDF2.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_BackwardEuler
  SOURCES Tempus_UnitTest_BackwardEuler.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_BackwardEuler
  SOURCES Tempus_UnitTest_DIRK_BackwardEuler.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_1Stage1stOrderRadauIA
  SOURCES Tempus_UnitTest_DIRK_1Stage1stOrderRadauIA.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_1StageTheta
  SOURCES Tempus_UnitTest_DIRK_1StageTheta.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_2Stage2ndOrderLobattoIIIB
  SOURCES Tempus_UnitTest_DIRK_2Stage2ndOrderLobattoIIIB.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_DIRK_General
  SOURCES Tempus_UnitTest_DIRK_General.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_EDIRK_2Stage3rdOrder
  SOURCES Tempus_UnitTest_EDIRK_2Stage3rdOrder.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_EDIRK_2StageTheta
  SOURCES Tempus_UnitTest_EDIRK_2StageTheta.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_EDIRK_TrapezoidalRule
  SOURCES Tempus_UnitTest_EDIRK_TrapezoidalRule.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_21Pair
  SOURCES Tempus_UnitTest_SDIRK_21Pair.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_2Stage2ndOrder
  SOURCES Tempus_UnitTest_SDIRK_2Stage2ndOrder.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_2Stage3rdOrder
  SOURCES Tempus_UnitTest_SDIRK_2Stage3rdOrder.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_3Stage4thOrder
  SOURCES Tempus_UnitTest_SDIRK_3Stage4thOrder.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_5Stage4thOrder
  SOURCES Tempus_UnitTest_SDIRK_5Stage4thOrder.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_5Stage5thOrder
  SOURCES Tempus_UnitTest_SDIRK_5Stage5thOrder.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_SDIRK_ImplicitMidpoint
  SOURCES Tempus_UnitTest_SDIRK_ImplicitMidpoint.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_3_8Rule
  SOURCES Tempus_UnitTest_ERK_3_8Rule.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_3Stage3rdOrder
  SOURCES Tempus_UnitTest_ERK_3Stage3rdOrder.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_3Stage3rdOrderHeun
  SOURCES Tempus_UnitTest_ERK_3Stage3rdOrderHeun.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_3Stage3rdOrderTVD
  SOURCES Tempus_UnitTest_ERK_3Stage3rdOrderTVD.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_4Stage3rdOrderRunge
  SOURCES Tempus_UnitTest_ERK_4Stage3rdOrderRunge.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_4Stage4thOrder
  SOURCES Tempus_UnitTest_ERK_4Stage4thOrder.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_5Stage3rdOrderKandG
  SOURCES Tempus_UnitTest_ERK_5Stage3rdOrderKandG.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_BogackiShampine32
  SOURCES Tempus_UnitTest_ERK_BogackiShampine32.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_ForwardEuler
  SOURCES Tempus_UnitTest_ERK_ForwardEuler.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_General
  SOURCES Tempus_UnitTest_ERK_General.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_Merson45
  SOURCES Tempus_UnitTest_ERK_Merson45.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_Midpoint
  SOURCES Tempus_UnitTest_ERK_Midpoint.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ERK_Trapezoidal
  SOURCES Tempus_UnitTest_ERK_Trapezoidal.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_ForwardEuler
  SOURCES Tempus_UnitTest_ForwardEuler.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_HHTAlpha
  SOURCES Tempus_UnitTest_HHTAlpha.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_IMEX_RK
  SOURCES Tempus_UnitTest_IMEX_RK.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_IMEX_RK_Partition
  SOURCES Tempus_UnitTest_IMEX_RK_Partition.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Leapfrog
  SOURCES Tempus_UnitTest_Leapfrog.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_NewmarkExplicitAForm
  SOURCES Tempus_UnitTest_NewmarkExplicitAForm.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_NewmarkImplicitAForm
  SOURCES Tempus_UnitTest_NewmarkImplicitAForm.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_NewmarkImplicitDForm
  SOURCES Tempus_UnitTest_NewmarkImplicitDForm.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_OperatorSplit
  SOURCES Tempus_UnitTest_OperatorSplit.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Trapezoidal
  SOURCES Tempus_UnitTest_Trapezoidal.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UnitTest_Subcycling
  SOURCES Tempus_UnitTest_Subcycling.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

