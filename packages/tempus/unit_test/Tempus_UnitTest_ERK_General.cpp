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

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
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
TEUCHOS_UNIT_TEST(ERK_General, Default_Construction)
{
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());

  // Default construction.
  auto stepper = rcp(new Tempus::StepperERK_General<double>());
  stepper->setModel(model);
  stepper->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  // Default values for construction.
  auto obs    = rcp(new Tempus::StepperRKObserverComposite<double>());

  bool useFSAL              = stepper->getUseFSALDefault();
  std::string ICConsistency = stepper->getICConsistencyDefault();
  bool ICConsistencyCheck   = stepper->getICConsistencyCheckDefault();
  bool useEmbedded          = stepper->getUseEmbeddedDefault();

  int NumStages = 4;
  Teuchos::SerialDenseMatrix<int,double> A(NumStages,NumStages);
  Teuchos::SerialDenseVector<int,double> b(NumStages);
  Teuchos::SerialDenseVector<int,double> c(NumStages);
  Teuchos::SerialDenseVector<int,double> bstar(0);

  // Fill A:
  A(0,0) = 0.0;  A(0,1) = 0.0;  A(0,2) = 0.0;  A(0,3) = 0.0;
  A(1,0) = 0.5;  A(1,1) = 0.0;  A(1,2) = 0.0;  A(1,3) = 0.0;
  A(2,0) = 0.0;  A(2,1) = 0.5;  A(2,2) = 0.0;  A(2,3) = 0.0;
  A(3,0) = 0.0;  A(3,1) = 0.0;  A(3,2) = 1.0;  A(3,3) = 0.0;

  // Fill b:
  b(0) = 1.0/6.0;  b(1) = 1.0/3.0;  b(2) = 1.0/3.0;  b(3) = 1.0/6.0;

  // fill c:
  c(0) = 0.0;  c(1) = 0.5;  c(2) = 0.5;  c(3) = 1.0;

  int order = 4;


  // Test the set functions.
  stepper->setObserver(obs);                           stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseFSAL(useFSAL);                        stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistency(ICConsistency);            stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setICConsistencyCheck(ICConsistencyCheck);  stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());
  stepper->setUseEmbedded(useEmbedded);                stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

  stepper->setTableau(A, b, c, order, order, order);   stepper->initialize();  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());


  // Full argument list construction.
  stepper = rcp(new Tempus::StepperERK_General<double>(
    model, obs, useFSAL, ICConsistency, ICConsistencyCheck, useEmbedded,
    A, b, c, order, order, order, bstar));
  TEUCHOS_TEST_FOR_EXCEPT(!stepper->isInitialized());

}
#endif // CONSTRUCTION


#ifdef STEPPERFACTORY_CONSTRUCTION
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ERK_General, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("General ERK", model);
}
#endif // STEPPERFACTORY_CONSTRUCTION


} // namespace Tempus_Test
