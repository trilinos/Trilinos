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
#include "Tempus_UnitTest_Utils.hpp"
#include "Tempus_StepperRKButcherTableau.hpp"

#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ImplicitModel.hpp"
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
#define STEPPERFACTORY_CONSTRUCTION


#ifdef CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(OperatorSplit, Default_Construction)
{
  auto explicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ExplicitModel<double>());
  auto implicitModel = rcp(new Tempus_Test::VanDerPol_IMEX_ImplicitModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperOperatorSplit<double>());
  auto sf = Teuchos::rcp(new Tempus::StepperFactory<double>());
  auto subStepper1 =
    sf->createStepperForwardEuler(explicitModel, Teuchos::null);
  auto subStepper2 =
    sf->createStepperBackwardEuler(implicitModel, Teuchos::null);
  stepper->addStepper(subStepper1);
  stepper->addStepper(subStepper2);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto obs    = rcp(new Tempus::StepperOperatorSplitObserver<double>());

  bool useFSAL              = stepper->getUseFSALDefault();
  std::string ICConsistency = stepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheckDefault();
  int order = 2;

  // Test the set functions.
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setOrder(order);                            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setOrderMin(order);                         stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setOrderMax(order);                         stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Full argument list construction.
  std::vector<RCP<const Thyra::ModelEvaluator<double> > > models;
  models.push_back(explicitModel);
  models.push_back(implicitModel);

  std::vector<Teuchos::RCP<Tempus::Stepper<double> > > subStepperList;
  subStepperList.push_back(subStepper1);
  subStepperList.push_back(subStepper2);

  stepper = rcp(new Tempus::StepperOperatorSplit<double>(
    models, subStepperList, obs, useFSAL, ICConsistency, ICConsistencyCheck,
    order, order, order));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

}
#endif // CONSTRUCTION


#ifdef STEPPERFACTORY_CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(OperatorSplit, StepperFactory_Construction)
{
  // Read params from .xml file
  auto pList = getParametersFromXmlFile(
                 "../test/OperatorSplit/Tempus_OperatorSplit_VanDerPol.xml");
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
