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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_BasicDiscreteAdjointStepperTester.hpp"
#include "Rythmos_IntegratorBuilder.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_AdjointModelEvaluator.hpp"
#include "Thyra_LinearNonlinearSolver.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"

#include "../SinCos/SinCosModel.hpp"
#include "../VanderPol/VanderPolModel.hpp"


namespace Rythmos {


Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_MEDIUM;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::setVerbosityLevelOption(
    "verb-level", &verbLevel,
    "Verbosity level for objects.",
    &Teuchos::UnitTestRepository::getCLP() );
}


TEUCHOS_UNIT_TEST( BasicDiscreteAdjointStepperTester, rawNonlinearAdjoint )
{

  using Teuchos::outArg;
  using Teuchos::describe;
  using Teuchos::getParametersFromXmlString;
  typedef Thyra::ModelEvaluatorBase MEB;

  //
  out << "\nA) Create the nonlinear ME ...\n";
  //

  RCP<VanderPolModel> stateModel = vanderPolModel(
    getParametersFromXmlString(
      "<ParameterList>"
      "  <Parameter name=\"Implicit model formulation\" type=\"bool\" value=\"1\"/>"
      "</ParameterList>"
      )
    );

  //
  out << "\nB) Create the nonlinear solver ...\n";
  //

  RCP<TimeStepNonlinearSolver<double> > nlSolver = timeStepNonlinearSolver<double>(
    getParametersFromXmlString(
      "<ParameterList>"
      "  <Parameter name=\"Default Tol\" type=\"double\" value=\"1.0e-10\"/>"
      "  <Parameter name=\"Default Max Iters\" type=\"int\" value=\"20\"/>"
      "</ParameterList>"
      )
    );

  //
  out << "\nC) Create the integrator for the forward state problem ...\n";
  //

  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>(
    Teuchos::getParametersFromXmlString(
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
      "      <Parameter name=\"Fixed dt\" type=\"double\" value=\"0.5\"/>" // Gives 2 time steps!
      "    </ParameterList>"
      "  </ParameterList>"
      "  <ParameterList name=\"Interpolation Buffer Settings\">"
      "    <ParameterList name=\"Trailing Interpolation Buffer Selection\">"
      "      <Parameter name=\"Interpolation Buffer Type\" type=\"string\" value=\"Interpolation Buffer\"/>"
      "    </ParameterList>"
      "  </ParameterList>"
      "</ParameterList>"
      )
    );
  
  MEB::InArgs<double> ic = stateModel->getNominalValues();
  RCP<IntegratorBase<double> > integrator = ib->create(stateModel, ic, nlSolver);
  //integrator->setVerbLevel(Teuchos::VERB_EXTREME);

  // ToDo: Set the trailing IB to pick up the entire state solution!

  // 
  out << "\nD) Solve the basic forward problem ...\n";
  //

  const TimeRange<double> fwdTimeRange = integrator->getFwdTimeRange();
  const double t_final = fwdTimeRange.upper();
  RCP<const Thyra::VectorBase<double> > x_final, x_dot_final;
  get_fwd_x_and_x_dot( *integrator, t_final, outArg(x_final), outArg(x_dot_final) );

  out << "\nt_final = " << t_final << "\n";
  out << "\nx_final: " << *x_final;
  out << "\nx_dot_final: " << *x_dot_final;

  //
  out << "\nE) Create the basic adjoint model (no distributed response) ...\n";
  //

  RCP<AdjointModelEvaluator<double> > adjModel =
    adjointModelEvaluator<double>(
      stateModel, fwdTimeRange
      );
  adjModel->setFwdStateSolutionBuffer(integrator);

  //
  out << "\nF) Create a stepper and integrator for the adjoint ...\n";
  //
  
  RCP<Thyra::LinearNonlinearSolver<double> > adjTimeStepSolver =
    Thyra::linearNonlinearSolver<double>();
  RCP<Rythmos::StepperBase<double> > adjStepper =
    integrator->getStepper()->cloneStepperAlgorithm();

  //
  out << "\nG) Set up the initial condition for the adjoint at the final time ...\n";
  //
  
  const RCP<const Thyra::VectorSpaceBase<double> >
    f_space = stateModel->get_f_space();
  
  // lambda(t_final) = x_final
  const RCP<Thyra::VectorBase<double> > lambda_ic = createMember(f_space);
  V_V( lambda_ic.ptr(), *x_final );
  
  // lambda_dot(t_final,i) = 0.0
  const RCP<Thyra::VectorBase<double> > lambda_dot_ic = createMember(f_space);
  Thyra::V_S( lambda_dot_ic.ptr(), 0.0 );
  
  MEB::InArgs<double> adj_ic = adjModel->getNominalValues();
  adj_ic.set_x(lambda_ic);
  adj_ic.set_x_dot(lambda_dot_ic);
  out << "\nadj_ic: " << describe(adj_ic, Teuchos::VERB_EXTREME);

  RCP<Rythmos::IntegratorBase<double> > adjIntegrator =
    ib->create(adjModel, adj_ic, adjTimeStepSolver);
  
  //
  out << "\nH) Integrate the adjoint backwards in time (using backward time) ...\n";
  //
  
  adjStepper->setInitialCondition(adj_ic);
  adjIntegrator->setStepper(adjStepper, fwdTimeRange.length());
  
  const double adj_t_final = fwdTimeRange.length();
  RCP<const Thyra::VectorBase<double> > lambda_final, lambda_dot_final;
  get_fwd_x_and_x_dot( *adjIntegrator, adj_t_final,
    outArg(lambda_final), outArg(lambda_dot_final) );

  out << "\nadj_t_final = " << adj_t_final << "\n";
  out << "\nlambda_final: " << *lambda_final;
  out << "\nlambda_dot_final: " << *lambda_dot_final;

}


TEUCHOS_UNIT_TEST( BasicDiscreteAdjointStepperTester, create )
{
  const RCP<BasicDiscreteAdjointStepperTester<double> >
    adjStepperTester = basicDiscreteAdjointStepperTester<double>();
  TEST_ASSERT(adjStepperTester != Teuchos::null);
}


TEUCHOS_UNIT_TEST( BasicDiscreteAdjointStepperTester, linear )
{

  using Teuchos::getParametersFromXmlString;
  typedef Thyra::ModelEvaluatorBase MEB;

  //
  out << "\nA) Create the nonlinear ME ...\n";
  //

  RCP<VanderPolModel> stateModel = vanderPolModel(
    getParametersFromXmlString(
      "<ParameterList>"
      "  <Parameter name=\"Implicit model formulation\" type=\"bool\" value=\"1\"/>"
      "</ParameterList>"
      )
    );

  //
  out << "\nB) Create the nonlinear solver ...\n";
  //

  RCP<TimeStepNonlinearSolver<double> > nlSolver = timeStepNonlinearSolver<double>(
    getParametersFromXmlString(
      "<ParameterList>"
      "  <Parameter name=\"Default Tol\" type=\"double\" value=\"1.0e-10\"/>"
      "  <Parameter name=\"Default Max Iters\" type=\"int\" value=\"20\"/>"
      "</ParameterList>"
      )
    );

  //
  out << "\nC) Create the integrator for the forward state problem ...\n";
  //

  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>(
    getParametersFromXmlString(
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
      "      <Parameter name=\"Fixed dt\" type=\"double\" value=\"0.5\"/>" // Gives 2 time steps!
      "    </ParameterList>"
      "  </ParameterList>"
      "  <ParameterList name=\"Interpolation Buffer Settings\">"
      "    <ParameterList name=\"Trailing Interpolation Buffer Selection\">"
      "      <Parameter name=\"Interpolation Buffer Type\" type=\"string\" value=\"Interpolation Buffer\"/>"
      "    </ParameterList>"
      "  </ParameterList>"
      "</ParameterList>"
      )
    );
  
  MEB::InArgs<double> ic = stateModel->getNominalValues();
  RCP<IntegratorBase<double> > integrator = ib->create(stateModel, ic, nlSolver);
  //integrator->setVerbLevel(Teuchos::VERB_EXTREME);

  // ToDo: Set the trailing IB to pick up the entire state solution!

  const TimeRange<double> fwdTimeRange = integrator->getFwdTimeRange();

  //
  out << "\nD) Create the basic adjoint model (no distributed response) ...\n";
  //

  RCP<AdjointModelEvaluator<double> > adjModel =
    adjointModelEvaluator<double>(
      stateModel, fwdTimeRange
      );
  adjModel->setFwdStateSolutionBuffer(integrator);

  //
  out << "\nE) Test the adjoint stepper against forward sensitivities ...\n";
  //

  const RCP<BasicDiscreteAdjointStepperTester<double> > adjStepperTester =
    basicDiscreteAdjointStepperTester<double>();

  adjStepperTester->setOStream(Teuchos::fancyOStream(Teuchos::rcpFromRef(out)));
  adjStepperTester->setVerbLevel(verbLevel);

  const bool result = adjStepperTester->testAdjointStepper(*adjModel, integrator.ptr());
  if (!result)
    success = false;

}


} // namespace Rythmos
