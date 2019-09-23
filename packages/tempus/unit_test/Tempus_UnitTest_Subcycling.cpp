// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_StepperFactory.hpp"
//#include "Tempus_UnitTest_Utils.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::StepperFactory;
using Tempus::StepperExplicitRK;

// Comment out any of the following tests to exclude from build/run.
#define CONSTRUCTION
//#define STEPPERFACTORY_CONSTRUCTION


#ifdef CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, Default_Construction)
{
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperSubcycling<double>());
  auto sf = Teuchos::rcp(new Tempus::StepperFactory<double>());
  auto stepperBE = sf->createStepperBackwardEuler(model, Teuchos::null);
  stepper->setSubcyclingStepper(stepperBE);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized_);

  // Default values for construction.
  auto obs    = rcp(new Tempus::StepperSubcyclingObserver<double>());
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSALDefault();
  std::string ICConsistency = stepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheckDefault();

  // Test the set functions
  stepper->setSolver(solver);                          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized_);
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized_);
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized_);
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized_);
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized_);


  // Full argument list construction.
  auto scIntegrator = Teuchos::rcp(new Tempus::IntegratorBasic<double>());
  RCP<ParameterList> tempusPL = scIntegrator->getTempusParameterList();
  { // Set default subcycling Stepper to Forward Euler.
    tempusPL->sublist("Default Integrator")
                 .set("Stepper Name", "Default Subcycling Stepper");
    RCP<ParameterList> stepperPL = Teuchos::parameterList();
    stepperPL->set("Stepper Type", "Forward Euler");
    tempusPL->set("Default Subcycling Stepper", *stepperPL);
  }
  scIntegrator->setParameterList(Teuchos::null);

  stepper = rcp(new Tempus::StepperSubcycling<double>(
    model, obs, scIntegrator, useFSAL, ICConsistency, ICConsistencyCheck));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized_);

}
#endif // CONSTRUCTION


#ifdef STEPPERFACTORY_CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Subcycling, StepperFactory_Construction)
{
  // Read params from .xml file
  auto pList = getParametersFromXmlFile(
                 "../test/Subcycling/Tempus_Subcycling_VanDerPol.xml");
  auto tempusPL  = sublist(pList, "Tempus", true);
  auto stepperPL = sublist(tempusPL, "Demo Stepper", true);

  auto explicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  auto implicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());
  std::vector<RCP<const Thyra::ModelEvaluator<double> > > models;
  models.push_back(explicitModel);
  models.push_back(implicitModel);


  auto sf = Teuchos::rcp(new StepperFactory<double>());

  // Test using ParameterList.
  // Passing in model.
  auto stepper = sf->createMultiSteppers(stepperPL, models);
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

}
#endif // STEPPERFACTORY_CONSTRUCTION


} // namespace Tempus_Test
