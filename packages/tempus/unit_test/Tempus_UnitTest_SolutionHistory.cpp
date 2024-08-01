//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_SolutionHistory.hpp"
#include "Tempus_InterpolatorFactory.hpp"

#include "../TestModels/DahlquistTestModel.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, Default_Construction)
{
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.

  auto model    = rcp(new Tempus_Test::SinCosModel<double>());
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX<double>(icSolution);
  sh->addState(icState);
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  // Test the get functions (i.e., defaults).
  TEST_COMPARE(sh->getNumStates(), ==, 1);
  TEST_COMPARE(sh->getInterpolator()->order(), ==, 0);
  TEST_COMPARE(sh->getStorageType(), ==, Tempus::STORAGE_TYPE_UNDO);
  TEST_COMPARE(sh->getStorageLimit(), ==, 2);

  // Test the set functions.
  sh->setName("Testing");
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());
  sh->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());
  sh->setStorageTypeString("Static");
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());
  sh->setStorageLimit(99);
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  TEST_COMPARE(sh->getName(), ==, "Testing");
  TEST_COMPARE(sh->getStorageType(), ==, Tempus::STORAGE_TYPE_STATIC);
  TEST_COMPARE(sh->getStorageTypeString(), ==, "Static");
  TEST_COMPARE(sh->getStorageLimit(), ==, 99);

  TEST_FLOATING_EQUALITY(sh->minTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 0.0, 1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, Full_Construction)
{
  std::string name = "Unit Test";

  auto history  = rcp(new std::vector<RCP<Tempus::SolutionState<double> > >);
  auto model    = rcp(new Tempus_Test::SinCosModel<double>());
  auto inArgsIC = model->getNominalValues();
  auto icSolution =
      rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  Teuchos::RCP<Tempus::SolutionState<double> > icState =
      Tempus::createSolutionStateX<double>(icSolution);
  history->push_back(icState);

  auto interpolator = Tempus::InterpolatorFactory<double>::createInterpolator();

  Tempus::StorageType storageType = Tempus::STORAGE_TYPE_STATIC;
  int storageLimit                = 99;

  auto sh = rcp(new Tempus::SolutionHistory<double>(name, history, interpolator,
                                                    storageType, storageLimit));
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  TEST_COMPARE(sh->getNumStates(), ==, 1);
  TEST_COMPARE(sh->getInterpolator()->order(), ==, 0);
  TEST_COMPARE(sh->getStorageType(), ==, Tempus::STORAGE_TYPE_STATIC);
  TEST_COMPARE(sh->getStorageLimit(), ==, 99);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, Create_Construction)
{
  auto sh_temp          = rcp(new Tempus::SolutionHistory<double>());
  RCP<ParameterList> pl = sh_temp->getNonconstParameterList();

  pl->setName("Unit Test");
  pl->set<std::string>("Storage Type", "Static");
  pl->set<int>("Storage Limit", 99);

  pl->sublist("Interpolator").set("Interpolator Type", "Lagrange");
  pl->sublist("Interpolator").set("Order", 1);

  auto sh = Tempus::createSolutionHistoryPL<double>(pl);
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.

  TEST_COMPARE(sh->getNumStates(), ==, 0);
  TEST_COMPARE(sh->getInterpolator()->order(), ==, 1);
  TEST_COMPARE(sh->getStorageType(), ==, Tempus::STORAGE_TYPE_STATIC);
  TEST_COMPARE(sh->getStorageLimit(), ==, 99);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, addState_With_Keep_Newest)
{
  // Setup SolutionHistory for testing.
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.

  auto model    = rcp(new Tempus_Test::DahlquistTestModel<double>(-1.0, false));
  auto inArgsIC = model->getNominalValues();
  auto icSoln   = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto state0   = Tempus::createSolutionStateX<double>(icSoln);
  state0->setTime(0.0);
  state0->setIndex(0);

  sh->setStorageTypeString("Keep Newest");
  sh->addState(state0);
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  TEST_COMPARE(sh->getNumStates(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 0.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->minTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 0.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(0, false) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 0.0, 1.0e-14);

  // ---------------------------------------------------------------------------

  // Second State -- should replace first state.
  auto state1 = Tempus::createSolutionStateX<double>(icSoln);
  state1->setTime(1.0);
  state1->setIndex(1);
  sh->addState(state1);

  TEST_COMPARE(sh->getNumStates(), ==, 1);  // Only 1 state!
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 1.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->minTime(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 1.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(0, false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(1, false) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 1.0, 1.0e-14);

  // ---------------------------------------------------------------------------

  // Third State -- old state should not be added.
  auto state2 = Tempus::createSolutionStateX<double>(icSoln);
  state2->setTime(-1.0);
  state2->setIndex(-1);
  sh->addState(state2);

  // Still second state.
  TEST_COMPARE(sh->getNumStates(), ==, 1);  // Only 1 state!
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 1.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->minTime(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 1.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(0, false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(1, false) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 1.0, 1.0e-14);

  // ---------------------------------------------------------------------------
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, addState_With_Undo)
{
  // Setup SolutionHistory for testing.
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.

  auto model    = rcp(new Tempus_Test::DahlquistTestModel<double>(-1.0, false));
  auto inArgsIC = model->getNominalValues();
  auto icSoln   = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto state0   = Tempus::createSolutionStateX<double>(icSoln);
  state0->setTime(0.0);
  state0->setIndex(0);
  sh->addState(state0);
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  TEST_COMPARE(sh->getNumStates(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 0.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->minTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 0.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(0, false) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 0.0, 1.0e-14);

  // ---------------------------------------------------------------------------

  // Second State
  auto state1 = Tempus::createSolutionStateX<double>(icSoln);
  state1->setTime(1.0);
  state1->setIndex(1);
  sh->addState(state1);

  TEST_COMPARE(sh->getNumStates(), ==, 2);
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 1.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->minTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 1.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(0, false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(1, false) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 1);
  TEST_COMPARE(sh->getStateTimeIndexNM1()->getIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexNM1()->getTime(), 0.0, 1.0e-14);

  // ---------------------------------------------------------------------------

  // Third State -- should be added and first state dropped.
  auto state2 = Tempus::createSolutionStateX<double>(icSoln);
  state2->setTime(2.0);
  state2->setIndex(2);
  sh->addState(state2);

  TEST_COMPARE(sh->getNumStates(), ==, 2);  // Only 2 states!
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 2.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 2);
  TEST_FLOATING_EQUALITY(sh->minTime(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 2.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(1, false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(2, false) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 2);
  TEST_COMPARE(sh->getStateTimeIndexNM1()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexNM1()->getTime(), 1.0, 1.0e-14);

  // ---------------------------------------------------------------------------

  // Fourth State -- old state should not be added.
  sh->addState(state0);

  // Still third and second states.
  TEST_COMPARE(sh->getNumStates(), ==, 2);  // Only 2 states!
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 2.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 2);
  TEST_FLOATING_EQUALITY(sh->minTime(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 2.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1(false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2(false) == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(1, false) != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(2, false) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 2);
  TEST_COMPARE(sh->getStateTimeIndexNM1()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexNM1()->getTime(), 1.0, 1.0e-14);

  // ---------------------------------------------------------------------------
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, addState_With_Static)
{
  // Setup SolutionHistory for testing.
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.
  sh->setStorageTypeString("Static");
  sh->setStorageLimit(7);

  auto model    = rcp(new Tempus_Test::DahlquistTestModel<double>(-1.0, false));
  auto inArgsIC = model->getNominalValues();

  // Sequential insertion.
  // ---------------------------------------------------------------------------
  for (size_t i = 0; i < 13; ++i) {
    auto icSoln = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
    auto stateI = Tempus::createSolutionStateX<double>(icSoln);
    stateI->setTime(i * 0.9);
    stateI->setTimeStep(0.9);
    stateI->setIndex(i);
    sh->addState(stateI);
  }
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  sh->describe(out, Teuchos::VERB_MEDIUM);

  TEST_COMPARE(sh->getNumStates(), ==, 7);
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 10.8, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 12);
  TEST_FLOATING_EQUALITY(sh->minTime(), 5.4, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 10.8, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(7) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 12);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 10.8, 1.0e-14);

  TEST_COMPARE(sh->getStateTimeIndexNM1()->getIndex(), ==, 11);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexNM1()->getTime(), 9.9, 1.0e-14);

  TEST_COMPARE(sh->getStateTimeIndexNM2()->getIndex(), ==, 10);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexNM2()->getTime(), 9.0, 1.0e-14);

  TEST_COMPARE(sh->getStateTimeIndex(7)->getIndex(), ==, 7);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndex(7)->getTime(), 6.3, 1.0e-14);

  // ---------------------------------------------------------------------------

  // "Random" insertion.
  // ---------------------------------------------------------------------------
  sh->clear();
  for (size_t i = 0; i < 3; ++i) {
    auto icSoln = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
    auto stateI = Tempus::createSolutionStateX<double>(icSoln);
    stateI->setTime(2 * i * 0.9 + 6.3);
    stateI->setTimeStep(0.9);
    stateI->setIndex(2 * i + 7);
    sh->addState(stateI);
  }
  for (size_t i = 0; i < 4; ++i) {
    auto icSoln = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
    auto stateI = Tempus::createSolutionStateX<double>(icSoln);
    stateI->setTime(2 * i * 0.9 + 5.4);
    stateI->setTimeStep(0.9);
    stateI->setIndex(2 * i + 6);
    sh->addState(stateI);
  }
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  sh->describe(out, Teuchos::VERB_MEDIUM);

  TEST_COMPARE(sh->getNumStates(), ==, 7);
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 10.8, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 12);
  TEST_FLOATING_EQUALITY(sh->minTime(), 5.4, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 10.8, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(7) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 12);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 10.8, 1.0e-14);

  TEST_COMPARE(sh->getStateTimeIndexNM1()->getIndex(), ==, 11);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexNM1()->getTime(), 9.9, 1.0e-14);

  TEST_COMPARE(sh->getStateTimeIndexNM2()->getIndex(), ==, 10);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexNM2()->getTime(), 9.0, 1.0e-14);

  TEST_COMPARE(sh->getStateTimeIndex(7)->getIndex(), ==, 7);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndex(7)->getTime(), 6.3, 1.0e-14);

  // ---------------------------------------------------------------------------
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, removeState)
{
  // Setup SolutionHistory for testing.
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.
  sh->setStorageTypeString("Static");
  sh->setStorageLimit(7);

  auto model    = rcp(new Tempus_Test::DahlquistTestModel<double>(-1.0, false));
  auto inArgsIC = model->getNominalValues();

  // ---------------------------------------------------------------------------
  for (size_t i = 0; i < 13; ++i) {
    auto icSoln = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
    auto stateI = Tempus::createSolutionStateX<double>(icSoln);
    stateI->setTime(i * 0.9);
    stateI->setTimeStep(0.9);
    stateI->setIndex(i);
    sh->addState(stateI);
  }
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  sh->describe(out, Teuchos::VERB_MEDIUM);

  sh->removeState(sh->getStateTimeIndex(6));
  sh->removeState(7.2);
  sh->removeState(sh->getStateTimeIndex(10));
  sh->removeState(10.8);

  sh->describe(out, Teuchos::VERB_MEDIUM);

  TEST_COMPARE(sh->getNumStates(), ==, 3);
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 9.9, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 11);
  TEST_FLOATING_EQUALITY(sh->minTime(), 6.3, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 9.9, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexN() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM1() == Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndexNM2() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getStateTimeIndex(7) != Teuchos::null);

  TEST_COMPARE(sh->getStateTimeIndexN()->getIndex(), ==, 11);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexN()->getTime(), 9.9, 1.0e-14);

  TEST_COMPARE(sh->getStateTimeIndexNM2()->getIndex(), ==, 9);
  TEST_FLOATING_EQUALITY(sh->getStateTimeIndexNM2()->getTime(), 8.1, 1.0e-14);

  // ---------------------------------------------------------------------------
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, initWorkingState_Passing)
{
  // Setup SolutionHistory for testing.
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.

  auto model    = rcp(new Tempus_Test::DahlquistTestModel<double>(-1.0, false));
  auto inArgsIC = model->getNominalValues();
  auto icSoln   = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto state0   = Tempus::createSolutionStateX<double>(icSoln);
  state0->setTime(0.0);
  state0->setTimeStep(1.0);
  state0->setIndex(0);
  sh->addState(state0);
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  // State before initializing workingState
  // with sh->getWorkingState() == Teuchos::null, i.e., workingState
  // was promoted from last time step or initial condition.
  TEST_COMPARE(sh->getNumStates(), ==, 1);
  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getCurrentState() == sh->getStateTimeIndexN());
  TEUCHOS_ASSERT(sh->getWorkingState() == Teuchos::null);

  sh->initWorkingState();

  // State after initializing workingState.
  TEST_COMPARE(sh->getNumStates(), ==, 2);
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 0.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->minTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 1.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getCurrentState() == sh->getStateTimeIndexNM1());
  TEUCHOS_ASSERT(sh->getWorkingState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getWorkingState() == sh->getStateTimeIndexN());

  TEST_COMPARE(sh->getCurrentState()->getIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->getCurrentState()->getTime(), 0.0, 1.0e-14);

  TEST_COMPARE(sh->getWorkingState()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getWorkingState()->getTime(), 1.0, 1.0e-14);

  // ---------------------------------------------------------------------------
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, initWorkingState_Failing)
{
  // Setup SolutionHistory for testing.
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.

  auto model    = rcp(new Tempus_Test::DahlquistTestModel<double>(-1.0, false));
  auto inArgsIC = model->getNominalValues();
  auto icSoln   = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto state0   = Tempus::createSolutionStateX<double>(icSoln);
  state0->setTime(0.0);
  state0->setTimeStep(1.0);
  state0->setIndex(0);
  sh->addState(state0);
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());

  auto currentSoln = icSoln->clone_v();
  Thyra::V_S(currentSoln.ptr(), std::numeric_limits<double>::quiet_NaN());
  auto state1 = Tempus::createSolutionStateX<double>(currentSoln);
  state1->setTime(1.0);
  state1->setTimeStep(1.0);
  state1->setIndex(1);
  sh->addWorkingState(state1);
  sh->getWorkingState()->setSolutionStatus(Tempus::Status::FAILED);

  // State before initializing workingState
  // with sh->getWorkingState() != Teuchos::null, i.e., workingState
  // was NOT promoted from last time step or initial condition.
  TEST_COMPARE(sh->getNumStates(), ==, 2);
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 0.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->minTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 1.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getCurrentState() == sh->getStateTimeIndexNM1());
  TEUCHOS_ASSERT(get_ele(*(sh->getCurrentState()->getX()), 0) == 1.0);

  TEUCHOS_ASSERT(sh->getWorkingState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getWorkingState() == sh->getStateTimeIndexN());
  TEUCHOS_ASSERT(std::isnan(get_ele(*(sh->getWorkingState()->getX()), 0)));  // !!!

  TEST_COMPARE(sh->getCurrentState()->getIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->getCurrentState()->getTime(), 0.0, 1.0e-14);

  TEST_COMPARE(sh->getWorkingState()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getWorkingState()->getTime(), 1.0, 1.0e-14);

  sh->initWorkingState();

  // State after initializing workingState.
  // Should be unchanged except the workingState->getX() is reset
  // to currentState->getX().
  TEST_COMPARE(sh->getNumStates(), ==, 2);
  TEST_FLOATING_EQUALITY(sh->getCurrentTime(), 0.0, 1.0e-14);
  TEST_COMPARE(sh->getCurrentIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->minTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(sh->maxTime(), 1.0, 1.0e-14);

  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getCurrentState() == sh->getStateTimeIndexNM1());
  TEUCHOS_ASSERT(get_ele(*(sh->getCurrentState()->getX()), 0) == 1.0);

  TEUCHOS_ASSERT(sh->getWorkingState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getWorkingState() == sh->getStateTimeIndexN());
  TEUCHOS_ASSERT(get_ele(*(sh->getWorkingState()->getX()), 0) == 1.0);  // !!!

  TEST_COMPARE(sh->getCurrentState()->getIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->getCurrentState()->getTime(), 0.0, 1.0e-14);

  TEST_COMPARE(sh->getWorkingState()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getWorkingState()->getTime(), 1.0, 1.0e-14);

  // ---------------------------------------------------------------------------
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, promoteWorkingState_Passing)
{
  // Setup SolutionHistory for testing.
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.

  auto model    = rcp(new Tempus_Test::DahlquistTestModel<double>(-1.0, false));
  auto inArgsIC = model->getNominalValues();
  auto icSoln   = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto state0   = Tempus::createSolutionStateX<double>(icSoln);
  state0->setTime(0.0);
  state0->setTimeStep(1.0);
  state0->setIndex(0);
  sh->addState(state0);
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());
  sh->initWorkingState();
  sh->getWorkingState()->setSolutionStatus(Tempus::Status::PASSED);

  // State before promoting workingState.
  // with workingState PASSing.
  TEST_COMPARE(sh->getNumStates(), ==, 2);
  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getCurrentState() == sh->getStateTimeIndexNM1());
  TEUCHOS_ASSERT(sh->getWorkingState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getWorkingState() == sh->getStateTimeIndexN());

  TEST_COMPARE(sh->getCurrentState()->getIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->getCurrentState()->getTime(), 0.0, 1.0e-14);

  TEST_COMPARE(sh->getWorkingState()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getWorkingState()->getTime(), 1.0, 1.0e-14);

  sh->promoteWorkingState();

  // State after promoting workingState.
  TEST_COMPARE(sh->getNumStates(), ==, 2);
  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getCurrentState() == sh->getStateTimeIndexN());
  TEUCHOS_ASSERT(sh->getWorkingState() == Teuchos::null);

  // Current state is now at t=1.0!
  TEST_COMPARE(sh->getCurrentState()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getCurrentState()->getTime(), 1.0, 1.0e-14);

  // ---------------------------------------------------------------------------
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(SolutionHistory, promoteWorkingState_Failing)
{
  // Setup SolutionHistory for testing.
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  TEUCHOS_TEST_FOR_EXCEPT(sh->isInitialized());  // Should NOT be initialized
                                                 // as no State is set.

  auto model    = rcp(new Tempus_Test::DahlquistTestModel<double>(-1.0, false));
  auto inArgsIC = model->getNominalValues();
  auto icSoln   = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto state0   = Tempus::createSolutionStateX<double>(icSoln);
  state0->setTime(0.0);
  state0->setTimeStep(1.0);
  state0->setIndex(0);
  sh->addState(state0);
  sh->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!sh->isInitialized());
  sh->initWorkingState();
  sh->getWorkingState()->setSolutionStatus(Tempus::Status::FAILED);

  // State before promoting workingState.
  // with workingState FAILing.
  TEST_COMPARE(sh->getNumStates(), ==, 2);
  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getCurrentState() == sh->getStateTimeIndexNM1());
  TEUCHOS_ASSERT(sh->getWorkingState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getWorkingState() == sh->getStateTimeIndexN());

  TEST_COMPARE(sh->getCurrentState()->getIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->getCurrentState()->getTime(), 0.0, 1.0e-14);

  TEST_COMPARE(sh->getWorkingState()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getWorkingState()->getTime(), 1.0, 1.0e-14);

  sh->promoteWorkingState();

  // State after promoting workingState.
  // Should be unchanged as we are trying to redo FAILing time step.
  TEST_COMPARE(sh->getNumStates(), ==, 2);
  TEUCHOS_ASSERT(sh->getCurrentState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getCurrentState() == sh->getStateTimeIndexNM1());
  TEUCHOS_ASSERT(sh->getWorkingState() != Teuchos::null);
  TEUCHOS_ASSERT(sh->getWorkingState() == sh->getStateTimeIndexN());

  // Current state is still at t=0.0!
  TEST_COMPARE(sh->getCurrentState()->getIndex(), ==, 0);
  TEST_FLOATING_EQUALITY(sh->getCurrentState()->getTime(), 0.0, 1.0e-14);

  // Working state is still at t=1.0!
  TEST_COMPARE(sh->getWorkingState()->getIndex(), ==, 1);
  TEST_FLOATING_EQUALITY(sh->getWorkingState()->getTime(), 1.0, 1.0e-14);

  // ---------------------------------------------------------------------------
}

}  // namespace Tempus_Unit_Test
