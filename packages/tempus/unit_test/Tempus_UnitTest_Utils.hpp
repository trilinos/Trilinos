// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_UnitTest_Utils_hpp
#define Tempus_UnitTest_Utils_hpp

#include "Tempus_StepperFactory.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperRKModifierBase.hpp"
#include "Tempus_StepperRKObserverBase.hpp"
#include "Tempus_StepperRKModifierXBase.hpp"
#include "Tempus_StepperRKModifierDefault.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;

using Tempus::StepperFactory;

/** \brief Unit test utility for Stepper construction through StepperFactory.
 */
void testFactoryConstruction(std::string stepperType,
  const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model)
{
  RCP<StepperFactory<double> > sf = Teuchos::rcp(new StepperFactory<double>());

  // Test using stepperType
  // Passing in model.
  auto stepper = sf->createStepper(stepperType, model);
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  // With setting model.
  stepper = sf->createStepper(stepperType);
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Test using ParameterList.
  // Passing in model.
  auto stepperPL = rcp_const_cast<ParameterList>(stepper->getValidParameters());
  stepper = sf->createStepper(stepperPL, model);
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  // With setting model.
  stepper = sf->createStepper(stepperPL);
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
}


/** \brief Unit test utility for ExplicitRK Stepper construction and accessors.
 */
void testExplicitRKAccessorsFullConstruction(std::string stepperType)
{
  // Default construction.
  RCP<Tempus::StepperExplicitRK<double> > stepper;
  if      (stepperType == "RK Explicit 3 Stage 3rd order")          stepper = rcp(new Tempus::StepperERK_3Stage3rdOrder<double>());
  else if (stepperType == "RK Explicit 3 Stage 3rd order by Heun")  stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderHeun<double>());
  else if (stepperType == "RK Explicit 3 Stage 3rd order TVD")      stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderTVD<double>());
  else if (stepperType == "RK Explicit 3/8 Rule")                   stepper = rcp(new Tempus::StepperERK_3_8Rule<double>());
  else if (stepperType == "RK Explicit 4 Stage 3rd order by Runge") stepper = rcp(new Tempus::StepperERK_4Stage3rdOrderRunge<double>());
  else if (stepperType == "RK Explicit 4 Stage")                    stepper = rcp(new Tempus::StepperERK_4Stage4thOrder<double>());
  else if (stepperType == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray") stepper = rcp(new Tempus::StepperERK_5Stage3rdOrderKandG<double>());
  else if (stepperType == "Bogacki-Shampine 3(2) Pair")             stepper = rcp(new Tempus::StepperERK_BogackiShampine32<double>());
  else if (stepperType == "RK Forward Euler")                       stepper = rcp(new Tempus::StepperERK_ForwardEuler<double>());
  else if (stepperType == "Merson 4(5) Pair")                       stepper = rcp(new Tempus::StepperERK_Merson45<double>());
  else if (stepperType == "RK Explicit Midpoint")                   stepper = rcp(new Tempus::StepperERK_Midpoint<double>());
  else if (stepperType == "RK Explicit Trapezoidal")                stepper = rcp(new Tempus::StepperERK_Trapezoidal<double>());
  else TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error - unknown stepperType = "+stepperType);


  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto modifier = rcp(new Tempus::StepperRKModifierDefault<double>());

  bool useFSAL              = stepper->getUseFSALDefault();
  std::string ICConsistency = stepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheckDefault();
  bool useEmbedded          = stepper->getUseEmbeddedDefault();

  // Test the set functions.
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  auto obs    = rcp(new Tempus::StepperRKObserverComposite<double>());
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
  stepper->setAppAction(modifier);                     stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseEmbedded(useEmbedded);                stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Full argument list construction.
  if        (stepperType == "RK Explicit 3 Stage 3rd order") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_3Stage3rdOrder<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_3Stage3rdOrder<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "RK Explicit 3 Stage 3rd order by Heun") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderHeun<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderHeun<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "RK Explicit 3 Stage 3rd order TVD") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderTVD<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderTVD<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "RK Explicit 3/8 Rule") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_3_8Rule<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_3_8Rule<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "RK Explicit 4 Stage 3rd order by Runge") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_4Stage3rdOrderRunge<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_4Stage3rdOrderRunge<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "RK Explicit 4 Stage") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_4Stage4thOrder<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_4Stage4thOrder<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_5Stage3rdOrderKandG<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_5Stage3rdOrderKandG<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "Bogacki-Shampine 3(2) Pair") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_BogackiShampine32<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_BogackiShampine32<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "RK Forward Euler") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_ForwardEuler<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_ForwardEuler<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "Merson 4(5) Pair") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_Merson45<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_Merson45<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "RK Explicit Midpoint") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_Midpoint<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_Midpoint<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  } else if (stepperType == "RK Explicit Trapezoidal") {
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    stepper = rcp(new Tempus::StepperERK_Trapezoidal<double>(
      model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
#endif
    stepper = rcp(new Tempus::StepperERK_Trapezoidal<double>(
      model, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded, modifier));
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  }
  else TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error - unknown stepperType = "+stepperType);
}


/** \brief Unit test class for RK Stepper Modifier AppAction.
 */
class StepperRKModifierTest
  : virtual public Tempus::StepperRKModifierBase<double>
{
public:

  /// Constructor
  StepperRKModifierTest()
    : testBEGIN_STEP(false),
      testBEGIN_STAGE(false),
      testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),
      testBEFORE_EXPLICIT_EVAL(false),
      testEND_STAGE(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testType("")
  {}

  /// Destructor
  virtual ~StepperRKModifierTest(){}

  /// Modify RK Stepper at action location.
  virtual void modify(
    Teuchos::RCP<Tempus::SolutionHistory<double> > sh,
    Teuchos::RCP<Tempus::StepperRKBase<double> > stepper,
    const typename Tempus::StepperRKAppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperRKAppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperRKAppAction<double>::BEGIN_STAGE:
      {
        testBEGIN_STAGE = true;
        break;
      }
      case StepperRKAppAction<double>::BEFORE_SOLVE:
      {
        testBEFORE_SOLVE = true;
        testDt = sh->getWorkingState()->getTimeStep()/10.0;
        sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      case StepperRKAppAction<double>::AFTER_SOLVE:
      {
        testAFTER_SOLVE = true;
        static bool done_once = false;
        if (!done_once) {
          testType = stepper->getStepperType() + " - Modifier";
          done_once = true;
        }
        stepper->setStepperType(testType);
        break;
      }
      case StepperRKAppAction<double>::BEFORE_EXPLICIT_EVAL:
      {
        testBEFORE_EXPLICIT_EVAL = true;
        break;
      }
      case StepperRKAppAction<double>::END_STAGE:
      {
        testEND_STAGE = true;
        break;
      }
      case StepperRKAppAction<double>::END_STEP:
      {
        testEND_STEP = true;
        auto x = sh->getWorkingState()->getX();
        testWorkingValue = get_ele(*(x), 0);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

  bool testBEGIN_STEP;
  bool testBEGIN_STAGE;
  bool testBEFORE_SOLVE;
  bool testAFTER_SOLVE;
  bool testBEFORE_EXPLICIT_EVAL;
  bool testEND_STAGE;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testType;
};


/** \brief Unit test class for RK Stepper Observer AppAction.
 */
class StepperRKObserverTest
  : virtual public Tempus::StepperRKObserverBase<double>
{
public:

  /// Constructor
  StepperRKObserverTest()
    : testBEGIN_STEP(false),
      testBEGIN_STAGE(false),
      testBEFORE_SOLVE(false),
      testAFTER_SOLVE(false),
      testBEFORE_EXPLICIT_EVAL(false),
      testEND_STAGE(false),
      testEND_STEP(false),
      testCurrentValue(-0.99),
      testWorkingValue(-0.99),
      testDt(-1.5),
      testType("")
  {}

  /// Destructor
  virtual ~StepperRKObserverTest(){}

  /// Observe RK Stepper at action location.
  virtual void observe(
    Teuchos::RCP<const Tempus::SolutionHistory<double> > sh,
    Teuchos::RCP<const Tempus::StepperRKBase<double> > stepper,
    const typename Tempus::StepperRKAppAction<double>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperRKAppAction<double>::BEGIN_STEP:
      {
        testBEGIN_STEP = true;
        auto x = sh->getCurrentState()->getX();
        testCurrentValue = get_ele(*(x), 0);
        break;
      }
      case StepperRKAppAction<double>::BEGIN_STAGE:
      {
        testBEGIN_STAGE = true;
        break;
      }
      case StepperRKAppAction<double>::BEFORE_SOLVE:
      {
        testBEFORE_SOLVE = true;
        testDt = sh->getWorkingState()->getTimeStep()/10.0;
        sh->getWorkingState()->setTimeStep(testDt);
        break;
      }
      case StepperRKAppAction<double>::AFTER_SOLVE:
      {
        testAFTER_SOLVE = true;
        testType = stepper->getStepperType();
        break;
      }
      case StepperRKAppAction<double>::BEFORE_EXPLICIT_EVAL:
      {
        testBEFORE_EXPLICIT_EVAL = true;
        break;
      }
      case StepperRKAppAction<double>::END_STAGE:
      {
        testEND_STAGE = true;
        break;
      }
      case StepperRKAppAction<double>::END_STEP:
      {
        testEND_STEP = true;
        auto x = sh->getWorkingState()->getX();
        testWorkingValue = get_ele(*(x), 0);
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

  bool testBEGIN_STEP;
  bool testBEGIN_STAGE;
  bool testBEFORE_SOLVE;
  bool testAFTER_SOLVE;
  bool testBEFORE_EXPLICIT_EVAL;
  bool testEND_STAGE;
  bool testEND_STEP;
  double testCurrentValue;
  double testWorkingValue;
  double testDt;
  std::string testType;
};


class StepperRKModifierXTest
  : virtual public Tempus::StepperRKModifierXBase<double>
{
public:

  /// Constructor
  StepperRKModifierXTest()
    : testX_BEGIN_STEP(false),
      testX_BEGIN_STAGE(false),
      testX_BEFORE_SOLVE(false),
      testX_AFTER_SOLVE(false),
      testX_BEFORE_EXPLICIT_EVAL(false),
      testXDOT_END_STAGE(false),
      testX_END_STEP(false),
      testX(-0.99),
      testXDot(-0.99),
      testDt(-1.5),
      testTime(-1.5)
  {}

  /// Destructor
  virtual ~StepperRKModifierXTest(){}

  /// Observe RK Stepper at end of takeStep.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<double> > x,
    const double time, const double dt,
    const typename Tempus::StepperRKModifierXBase<double>::MODIFIER_TYPE modType)
  {
    switch(modType) {
      case StepperRKModifierXBase<double>::X_BEGIN_STEP:
      {
        testX_BEGIN_STEP = true;
        testX = get_ele(*(x), 0);
        break;
      }
      case StepperRKModifierXBase<double>::X_BEGIN_STAGE:
      {
        testX_BEGIN_STAGE = true;
        break;
      }
      case StepperRKModifierXBase<double>::X_BEFORE_SOLVE:
      {
        testX_BEFORE_SOLVE = true;
        testDt = dt;
        break;
      }
      case StepperRKModifierXBase<double>::X_AFTER_SOLVE:
      {
        testX_AFTER_SOLVE = true;
        testTime = time;
        break;
      }
      case StepperRKModifierXBase<double>::X_BEFORE_EXPLICIT_EVAL:
      {
        testX_BEFORE_EXPLICIT_EVAL = true;
        break;
      }
      case StepperRKModifierXBase<double>::XDOT_END_STAGE:
      {
        testXDOT_END_STAGE = true;
        testXDot = get_ele(*(x), 0);
        break;
      }
      case StepperRKModifierXBase<double>::X_END_STEP:
      {
        testX_END_STEP = true;
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

  bool testX_BEGIN_STEP;
  bool testX_BEGIN_STAGE;
  bool testX_BEFORE_SOLVE;
  bool testX_AFTER_SOLVE;
  bool testX_BEFORE_EXPLICIT_EVAL;
  bool testXDOT_END_STAGE;
  bool testX_END_STEP;
  double testX;
  double testXDot;
  double testDt;
  double testTime;
};


/** \brief Unit test utility for ExplicitRK Stepper AppAction.
 */
void testExplicitRKAppAction(std::string stepperType,
                             Teuchos::FancyOStream &out, bool &success)
{
  // Default construction.
  RCP<Tempus::StepperExplicitRK<double> > stepper;
  if      (stepperType == "RK Explicit 3 Stage 3rd order")          stepper = rcp(new Tempus::StepperERK_3Stage3rdOrder<double>());
  else if (stepperType == "RK Explicit 3 Stage 3rd order by Heun")  stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderHeun<double>());
  else if (stepperType == "RK Explicit 3 Stage 3rd order TVD")      stepper = rcp(new Tempus::StepperERK_3Stage3rdOrderTVD<double>());
  else if (stepperType == "RK Explicit 3/8 Rule")                   stepper = rcp(new Tempus::StepperERK_3_8Rule<double>());
  else if (stepperType == "RK Explicit 4 Stage 3rd order by Runge") stepper = rcp(new Tempus::StepperERK_4Stage3rdOrderRunge<double>());
  else if (stepperType == "RK Explicit 4 Stage")                    stepper = rcp(new Tempus::StepperERK_4Stage4thOrder<double>());
  else if (stepperType == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray") stepper = rcp(new Tempus::StepperERK_5Stage3rdOrderKandG<double>());
  else if (stepperType == "Bogacki-Shampine 3(2) Pair")             stepper = rcp(new Tempus::StepperERK_BogackiShampine32<double>());
  else if (stepperType == "RK Forward Euler")                       stepper = rcp(new Tempus::StepperERK_ForwardEuler<double>());
  else if (stepperType == "General ERK")                            stepper = rcp(new Tempus::StepperERK_General<double>());
  else if (stepperType == "Merson 4(5) Pair")                       stepper = rcp(new Tempus::StepperERK_Merson45<double>());
  else if (stepperType == "RK Explicit Midpoint")                   stepper = rcp(new Tempus::StepperERK_Midpoint<double>());
  else if (stepperType == "RK Explicit Trapezoidal")                stepper = rcp(new Tempus::StepperERK_Trapezoidal<double>());
  else TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error - unknown stepperType = "+stepperType);

  auto testTypeOrig = stepper->getStepperType();

  // Test Modifier.
  {
    Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::SinCosModel<double>());
    stepper->setModel(model);
    auto modifier = rcp(new StepperRKModifierTest());
    stepper->setAppAction(modifier);
    stepper->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    auto testType = testTypeOrig + " - Modifier";

    // Create a SolutionHistory.
    auto solutionHistory = Tempus::createSolutionHistoryME(model);

    // Take one time step.
    stepper->setInitialConditions(solutionHistory);
    solutionHistory->initWorkingState();
    double dt = 0.1;
    solutionHistory->getWorkingState()->setTimeStep(dt);
    stepper->takeStep(solutionHistory);

    TEST_COMPARE(modifier->testBEGIN_STEP, ==, true);
    TEST_COMPARE(modifier->testBEGIN_STAGE, ==, true);
    TEST_COMPARE(modifier->testBEFORE_SOLVE, ==, true);
    TEST_COMPARE(modifier->testAFTER_SOLVE, ==, true);
    TEST_COMPARE(modifier->testBEFORE_EXPLICIT_EVAL, ==, true);
    TEST_COMPARE(modifier->testEND_STAGE, ==, true);
    TEST_COMPARE(modifier->testEND_STEP, ==, true);

    auto x = solutionHistory->getCurrentState()->getX();
    TEST_FLOATING_EQUALITY(modifier->testCurrentValue,get_ele(*(x), 0),1.0e-15);
    x = solutionHistory->getWorkingState()->getX();
    TEST_FLOATING_EQUALITY(modifier->testWorkingValue,get_ele(*(x), 0),1.0e-15);
    auto Dt = solutionHistory->getWorkingState()->getTimeStep();
    TEST_FLOATING_EQUALITY(modifier->testDt, Dt, 1.0e-15);

    TEST_COMPARE(modifier->testType, ==, testType);
  }

  // Test Observer.
  {
    Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::SinCosModel<double>());
    stepper->setModel(model);
    auto observer = rcp(new StepperRKObserverTest());
    stepper->setAppAction(observer);
    stepper->setStepperType(testTypeOrig);
    stepper->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
    auto testType = testTypeOrig;

    // Create a SolutionHistory.
    auto solutionHistory = Tempus::createSolutionHistoryME(model);

    // Take one time step.
    stepper->setInitialConditions(solutionHistory);
    solutionHistory->initWorkingState();
    double dt = 0.1;
    solutionHistory->getWorkingState()->setTimeStep(dt);
    stepper->takeStep(solutionHistory);

    TEST_COMPARE(observer->testBEGIN_STEP, ==, true);
    TEST_COMPARE(observer->testBEGIN_STAGE, ==, true);
    TEST_COMPARE(observer->testBEFORE_SOLVE, ==, true);
    TEST_COMPARE(observer->testAFTER_SOLVE, ==, true);
    TEST_COMPARE(observer->testBEFORE_EXPLICIT_EVAL, ==, true);
    TEST_COMPARE(observer->testEND_STAGE, ==, true);
    TEST_COMPARE(observer->testEND_STEP, ==, true);

    auto x = solutionHistory->getCurrentState()->getX();
    TEST_FLOATING_EQUALITY(observer->testCurrentValue,get_ele(*(x), 0),1.0e-15);
    x = solutionHistory->getWorkingState()->getX();
    TEST_FLOATING_EQUALITY(observer->testWorkingValue,get_ele(*(x), 0),1.0e-15);
    auto Dt = solutionHistory->getWorkingState()->getTimeStep();
    TEST_FLOATING_EQUALITY(observer->testDt, Dt, 1.0e-15);

    TEST_COMPARE(observer->testType, ==, testType);
  }

  // Test ModifierX.
  {
    Teuchos::RCP<const Thyra::ModelEvaluator<double> > model =
      rcp(new Tempus_Test::SinCosModel<double>());
    stepper->setModel(model);
    auto modifierX = rcp(new StepperRKModifierXTest());
    stepper->setAppAction(modifierX);
    stepper->setStepperType(testTypeOrig);
    stepper->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

    // Create a SolutionHistory.
    auto solutionHistory = Tempus::createSolutionHistoryME(model);

    // Take one time step.
    stepper->setInitialConditions(solutionHistory);
    solutionHistory->initWorkingState();
    double dt = 0.1;
    solutionHistory->getWorkingState()->setTimeStep(dt);
    stepper->takeStep(solutionHistory);

    TEST_COMPARE(modifierX->testX_BEGIN_STEP, ==, true);
    TEST_COMPARE(modifierX->testX_BEGIN_STAGE, ==, true);
    TEST_COMPARE(modifierX->testX_BEFORE_SOLVE, ==, true);
    TEST_COMPARE(modifierX->testX_AFTER_SOLVE, ==, true);
    TEST_COMPARE(modifierX->testX_BEFORE_EXPLICIT_EVAL, ==, true);
    TEST_COMPARE(modifierX->testXDOT_END_STAGE, ==, true);
    TEST_COMPARE(modifierX->testX_END_STEP, ==, true);

    auto x = solutionHistory->getCurrentState()->getX();
    TEST_FLOATING_EQUALITY(modifierX->testX, get_ele(*(x), 0), 1.0e-15);
    // Temporary memory for xDot is not guarranteed to exist outside the Stepper
    auto xDot = stepper->getStepperXDot(solutionHistory->getWorkingState());
    TEST_FLOATING_EQUALITY(modifierX->testXDot, get_ele(*(xDot), 0),1.0e-15);
    auto Dt = solutionHistory->getWorkingState()->getTimeStep();
    TEST_FLOATING_EQUALITY(modifierX->testDt, Dt, 1.0e-15);

    auto time = solutionHistory->getWorkingState()->getTime();
    TEST_FLOATING_EQUALITY(modifierX->testTime, time, 1.0e-15);
  }
}


} // namespace Tempus_Test
#endif // Tempus_UnitTest_Utils_hpp
