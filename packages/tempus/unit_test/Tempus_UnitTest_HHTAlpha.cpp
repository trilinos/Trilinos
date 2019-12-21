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

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestModels/HarmonicOscillatorModel.hpp"
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

// Comment out any of the following tests to exclude from build/run.
#define CONSTRUCTION
#define STEPPERFACTORY_CONSTRUCTION


#ifdef CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(HHTAlpha, Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperHHTAlpha<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Default values for construction.
  auto solver = rcp(new Thyra::NOXNonlinearSolver());
  solver->setParameterList(Tempus::defaultSolverParameters());

  bool useFSAL              = stepper->getUseFSALDefault();
  std::string ICConsistency = stepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheckDefault();
  bool zeroInitialGuess     = stepper->getZeroInitialGuess();
  std::string schemeName    = "Newmark Beta User Defined";
  double beta               = 0.25;
  double gamma              = 0.5;
  double alpha_f            = 0.0;
  double alpha_m            = 0.0;


  // Test the set functions.
  stepper->setSolver(solver);                          stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setZeroInitialGuess(zeroInitialGuess);      stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  stepper->setSchemeName(schemeName);                  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setBeta(beta);                              stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setGamma(gamma);                            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAlphaF(alpha_f);                         stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setAlphaM(alpha_m);                         stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Full argument list construction.
  stepper = rcp(new Tempus::StepperHHTAlpha<double>(
    model, Teuchos::null, solver, useFSAL,
    ICConsistency, ICConsistencyCheck, zeroInitialGuess,
    schemeName, beta, gamma, alpha_f, alpha_m));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

}
#endif // CONSTRUCTION


#ifdef STEPPERFACTORY_CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(HHTAlpha, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::HarmonicOscillatorModel<double>());
  testFactoryConstruction("HHT-Alpha", model);
}
#endif // STEPPERFACTORY_CONSTRUCTION


} // namespace Tempus_Test
