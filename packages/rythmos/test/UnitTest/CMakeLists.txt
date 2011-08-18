
INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../../../epetraext/example/model_evaluator/DiagonalTransient)


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_NOX)
IF (${PACKAGE_NAME}_ENABLE_NOX)
  ADD_DEFINITIONS(-DHAVE_RYTHMOS_NOX)
ENDIF()


#FILE(GLOB SOURCES *.cpp)

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#    UnitTest
#    SOURCES ${SOURCES}
#    DEPLIBS rythmos_test_models
#    COMM serial
#    STANDARD_PASS_OUTPUT
#    )

IF(${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      ConvergenceTestHelpers_UnitTest
      SOURCES Rythmos_ConvergenceTestHelpers_UnitTest.cpp Rythmos_UnitTest.cpp 
      DEPLIBS rythmos_test_models
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    CubicSplineInterpolator_UnitTest
    SOURCES Rythmos_CubicSplineInterpolator_UnitTest.cpp Rythmos_UnitTest.cpp 
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    XHOSTTYPE Darwin
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    DataStore_UnitTest
    SOURCES Rythmos_DataStore_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    DefaultIntegrator_UnitTest
    SOURCES Rythmos_DefaultIntegrator_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ExplicitRK_UnitTest
    SOURCES Rythmos_ExplicitRK_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    HermiteInterpolator_UnitTest
    SOURCES Rythmos_HermiteInterpolator_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ImplicitRK_UnitTest
    SOURCES Rythmos_ImplicitRK_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    COMM serial # MPI should be okay (but see below)!
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )
# 2009/09/02: rabartl: The above test times out in an MPI build on even one
# process!  That is a very bad sign and means that might point to a very bad
# problem with this code!

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    IntegratorBuilder_UnitTest
    SOURCES Rythmos_IntegratorBuilder_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    COMM serial  # MPI should be okay (but see below)!
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )
# 2009/09/02: rabartl: The above test times out in an MPI build on even one
# process!  That is a very bad sign and means that might point to a very bad
# problem with this code!

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    InterpolationBuffer_UnitTest
    SOURCES Rythmos_InterpolationBuffer_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    LinearInterpolator_UnitTest
    SOURCES Rythmos_LinearInterpolator_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Quadrature_UnitTest
    SOURCES Rythmos_Quadrature_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
    RKButcherTableau_UnitTest
    SOURCES Rythmos_RKButcherTableau_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
    SinCosModel_UnitTest
    SOURCES Rythmos_SinCosModel_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )


IF (${PACKAGE_NAME}_ENABLE_Sacado)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      VanderPolModel_UnitTest
      SOURCES Rythmos_VanderPolModel_UnitTest.cpp Rythmos_UnitTest.cpp
      DEPLIBS rythmos_test_models
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      AdjointModelEvaluator_UnitTest
      SOURCES Rythmos_AdjointModelEvaluator_UnitTest.cpp Rythmos_UnitTest.cpp
      DEPLIBS rythmos_test_models
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )

  IF(${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
        ForwardSensitivity_UnitTest
        SOURCES Rythmos_ForwardSensitivity_UnitTest.cpp Rythmos_UnitTest.cpp
        DEPLIBS rythmos_test_models
        NUM_MPI_PROCS 1
        STANDARD_PASS_OUTPUT
        )

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
        ImplicitBDF_UnitTest
        SOURCES Rythmos_ImplicitBDF_UnitTest.cpp Rythmos_UnitTest.cpp
        DEPLIBS rythmos_test_models
        NUM_MPI_PROCS 1
        STANDARD_PASS_OUTPUT
        )
  ENDIF()

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      BasicDiscreteAdjointStepperTester_UnitTest
      SOURCES Rythmos_BasicDiscreteAdjointStepperTester_UnitTest.cpp Rythmos_UnitTest.cpp
      DEPLIBS rythmos_test_models
      HOST None
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    StepperBuilder_UnitTest
    SOURCES Rythmos_StepperBuilder_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    XHOSTTYPE_TEST AIX
    STANDARD_PASS_OUTPUT
    )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
    StepperHelpers_UnitTest
    SOURCES Rythmos_StepperHelpers_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
    TimeRange_UnitTest
    SOURCES Rythmos_TimeRange_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    PointwiseInterpolationBufferAppender_UnitTest
    SOURCES Rythmos_PointwiseInterpolationBufferAppender_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)
IF (${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_DEBUG)
  IF (${PACKAGE_NAME}_ENABLE_DEBUG)
    SET(STEPPER_VALITATOR_ARGS ARGS)
  ELSE()
    SET(STEPPER_VALITATOR_ARGS ARGS --not-unit-test=Rythmos_StepperValidator_ImplicitRK_UnitTest)
  ENDIF()

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      StepperValidator_UnitTest
      SOURCES Rythmos_StepperValidator_UnitTest.cpp Rythmos_UnitTest.cpp
      DEPLIBS rythmos_test_models
      NUM_MPI_PROCS 1
      ${STEPPER_VALITATOR_ARGS}
      STANDARD_PASS_OUTPUT
      )

ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    GlobalErrorEstimator_UnitTest
    SOURCES Rythmos_GlobalErrorEstimator_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    Thyra_UnitTest
    SOURCES Rythmos_Thyra_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

IF(${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      ForwardEulerStepper_UnitTest
      SOURCES Rythmos_ForwardEuler_UnitTest.cpp Rythmos_UnitTest.cpp
      DEPLIBS rythmos_test_models
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      BackwardEulerStepper_UnitTest
      SOURCES Rythmos_BackwardEuler_UnitTest.cpp Rythmos_UnitTest.cpp
      DEPLIBS rythmos_test_models
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      ForwardSensitivityExplicitModelEvaluator_UnitTest
      SOURCES Rythmos_ForwardSensitivityExplicitModelEvaluator_UnitTest.cpp Rythmos_UnitTest.cpp
      DEPLIBS rythmos_test_models
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )
  
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      ForwardSensitivityImplicitModelEvaluator_UnitTest
      SOURCES Rythmos_ForwardSensitivityImplicitModelEvaluator_UnitTest.cpp Rythmos_UnitTest.cpp
      DEPLIBS rythmos_test_models
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    StateSerializerStrategy_UnitTest
    SOURCES Rythmos_StateSerializerStrategy_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    StackedStepper_UnitTest
    SOURCES Rythmos_StackedStepper_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    COMM serial
    STANDARD_PASS_OUTPUT
    )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ForwardSensitivityStepperTester_UnitTest
    SOURCES Rythmos_ForwardSensitivityStepperTester_UnitTest.cpp Rythmos_UnitTest.cpp
    DEPLIBS rythmos_test_models
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

IF(${PACKAGE_NAME}_ENABLE_ThyraEpetraExtAdapters)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
      TestModelBuilder_UnitTest
      SOURCES Rythmos_TestModelBuilder_UnitTest.cpp Rythmos_UnitTest.cpp Rythmos_UnitTestModels.cpp
      DEPLIBS rythmos_test_models
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )
ENDIF()
