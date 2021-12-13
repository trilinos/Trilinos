// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_UnitTest_Utils_hpp
#define Tempus_UnitTest_Utils_hpp


#include "Tempus_config.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "NOX_Thyra.H"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_IntegratorBasic.hpp"

#include "../TestModels/SinCosModel.hpp"


namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;

using Thyra::get_ele;

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


} // namespace Tempus_Test
#endif // Tempus_UnitTest_Utils_hpp
