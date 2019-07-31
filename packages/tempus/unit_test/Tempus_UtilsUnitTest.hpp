// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_UtilsUnitTest_hpp
#define Tempus_UtilsUnitTest_hpp

#include "Teuchos_UnitTestHarness.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/HarmonicOscillatorModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ImplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEXPart_ImplicitModel.hpp"

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::ParameterList;


void StepperInitializeBasic(RCP<const Thyra::ModelEvaluator<double> > model,
                            RCP<Tempus::Stepper<double> > stepper,
                            RCP<Tempus::StepperObserver<double> > obs,
                            Teuchos::FancyOStream &out, bool &success)
{
  // --- Set methods from Stepper class. ---

  // Model
  stepper->setModel(model);          TEST_ASSERT(!(stepper->isInitialized()));
  stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

  // Solver
  auto solverPL = stepper->defaultSolverParameters();
  stepper->setSolver(solverPL);      TEST_ASSERT(!(stepper->isInitialized()));
  stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

  // Observer
  if (obs != Teuchos::null) {
    stepper->setObserver(obs);       TEST_ASSERT(!(stepper->isInitialized()));
    stepper->initialize();           TEST_ASSERT(stepper->isInitialized());
  }

  // Initial Conditions
  auto sh = rcp(new Tempus::SolutionHistory<double>());
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  auto xIC       = rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto xDotIC    = xIC->clone_v();
  auto xDotDotIC = xIC->clone_v();
  auto icState = rcp(new Tempus::SolutionState<double>(xIC,xDotIC,xDotDotIC));
  sh->addState(icState);
  sh->initialize();

  stepper->setInitialConditions(sh); TEST_ASSERT(!(stepper->isInitialized()));
  stepper->initialize();             TEST_ASSERT(stepper->isInitialized());

  // --- Set methods from ParameterListAcceptor. ---

  // Parameter
  auto pl = rcp_const_cast<ParameterList>(stepper->getValidParameters());
  stepper->setParameterList(pl);     TEST_ASSERT(!(stepper->isInitialized()));
  stepper->initialize();             TEST_ASSERT(stepper->isInitialized());
}

} // namespace Tempus_Test

#endif // Tempus_UtilsUnitTest_hpp
