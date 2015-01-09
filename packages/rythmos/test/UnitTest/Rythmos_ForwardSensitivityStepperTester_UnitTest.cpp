//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_ForwardSensitivityStepperTester.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_IntegratorBuilder.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "../SinCos/SinCosModel.hpp"


namespace Rythmos {


TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepperTester, create )
{
  const RCP<ForwardSensitivityStepperTester<double> >
    fwdSensStepperTester = forwardSensitivityStepperTester<double>();
  TEST_ASSERT(fwdSensStepperTester != Teuchos::null);
}


TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepperTester, linear )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  // Create the basic model
 
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  modelPL->set("Implicit model formulation", true);
  modelPL->set("Accept model parameters", true);
  RCP<SinCosModel> model = sinCosModel();
  model->setParameterList(modelPL);

  // Set up the IntegratorBuilder

  RCP<ParameterList> ibPL = Teuchos::getParametersFromXmlString(
    "<ParameterList>"
    "  <ParameterList name=\"Stepper Settings\">"
    "    <ParameterList name=\"Stepper Selection\">"
    "      <Parameter name=\"Stepper Type\" type=\"string\" value=\"Backward Euler\"/>"
    "    </ParameterList>"
    "  </ParameterList>"
    "  <ParameterList name=\"Integration Control Strategy Selection\">"
    "    <Parameter name=\"Integration Control Strategy Type\" type=\"string\""
    "      value=\"Simple Integration Control Strategy\"/>"
    "    <ParameterList name=\"Simple Integration Control Strategy\">"
    "      <Parameter name=\"Take Variable Steps\" type=\"bool\" value=\"false\"/>"
    "      <Parameter name=\"Fixed dt\" type=\"double\" value=\"0.5\"/>"
    "    </ParameterList>"
    "  </ParameterList>"
    "</ParameterList>"
    );
  
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>(ibPL);

  // Create the actual integrator ready to go

  RCP<Thyra::NonlinearSolverBase<double> >
    nlSolver = timeStepNonlinearSolver<double>();

  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();

  RCP<IntegratorBase<double> >
    integrator = ib->create(model, ic, nlSolver);

  // Create the fwd sens stepper and an integrator

  RCP<ForwardSensitivityStepper<double> > stateAndSensStepper =
    forwardSensitivityStepper<double>();
  stateAndSensStepper->initializeSyncedSteppers(
    model, 0, ic, integrator->getNonconstStepper(), nlSolver
    );

  MEB::InArgs<double> state_and_sens_ic =
    createStateAndSensInitialCondition(*stateAndSensStepper, ic);
  stateAndSensStepper->setInitialCondition(state_and_sens_ic);
  
  RCP<IntegratorBase<double> >
    sensIntegrator = integrator->cloneIntegrator();

  sensIntegrator->setStepper(stateAndSensStepper, integrator->getFwdTimeRange().upper());
  
  // Test the stepper

  RCP<ParameterList> fsstPL = Teuchos::getParametersFromXmlString(
    "<ParameterList>"
    "  <ParameterList name=\"FD Calc\">"
    "    <Parameter name=\"FD Method\" type=\"string\" value=\"order-four-central\"/>"
    "    <Parameter name=\"FD Step Length\" type=\"double\" value=\"1e-3\"/>"
    "    <Parameter name=\"FD Step Select Type\" type=\"string\" value=\"Relative\"/>"
    "  </ParameterList>"
    "  <Parameter name=\"Error Tol\" type=\"double\" value=\"1e-10\"/>"
    "</ParameterList>"
    );
  const RCP<ForwardSensitivityStepperTester<double> >
    fwdSensStepperTester = forwardSensitivityStepperTester<double>(fsstPL);

  fwdSensStepperTester->setVerbLevel(Teuchos::VERB_EXTREME);
  fwdSensStepperTester->setOStream(Teuchos::rcpFromRef(out));

  TEST_ASSERT(fwdSensStepperTester->testForwardSens(sensIntegrator.ptr()));

}


} // namespace Rythmos
