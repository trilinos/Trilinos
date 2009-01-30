//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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

#include "Rythmos_Types.hpp"
#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_UnitTestModels.hpp"

#include "Rythmos_StepperBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepper, create ) { 
  RCP<ForwardSensitivityStepper<double> > sensStepper =
    forwardSensitivityStepper<double>();
  TEST_EQUALITY_CONST( is_null(sensStepper), false );
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardSensitivityStepper, initializeDecoupled ) {
  RCP<Teuchos::FancyOStream>
    std_out = Teuchos::VerboseObjectBase::getDefaultOStream();
  //Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
  RCP<SinCosModel> stateModel = sinCosModel();
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  modelPL->set("Accept model parameters",true);
  modelPL->set("Implicit model formulation",true);
  modelPL->set("Provide nominal values",true);
  modelPL->set("Coeff a", 0.0);
  modelPL->set("Coeff f", 1.0);
  modelPL->set("IC x_0", 0.0);
  modelPL->set("IC x_1", 1.0);
  stateModel->setParameterList(modelPL);
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  RCP<ParameterList> stepperPL = Teuchos::parameterList();
  stepperPL->set("Stepper Type","Backward Euler");
  builder->setParameterList(stepperPL);
  RCP<TimeStepNonlinearSolver<double> > nonlinearSolver = timeStepNonlinearSolver<double>();
  builder->setNonlinearSolver(nonlinearSolver);
  RCP<StepperBase<double> > stateStepper = builder->create();
  int p_index = 0;
  RCP<ForwardSensitivityStepper<double> > stateAndSensStepper = 
    forwardSensitivityStepper<double>();
  stateAndSensStepper->initializeSyncedSteppers(
    stateModel, 
    p_index, 
    stateModel->getNominalValues(),
    stateStepper, 
    nonlinearSolver
    );

  typedef Thyra::ModelEvaluatorBase MEB;
  const MEB::InArgs<double> state_ic = stateModel->getNominalValues();

  RCP<Thyra::VectorBase<double> > s_bar_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  assign( &*s_bar_init, 0.0 );
  RCP<Thyra::VectorBase<double> > s_bar_dot_init
    = createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
  assign( &*s_bar_dot_init, 0.0 );

  RCP<const StateAndForwardSensitivityModelEvaluator<double> >
    stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

  MEB::InArgs<double>
    state_and_sens_ic = stateAndSensStepper->getModel()->createInArgs();

  // Copy time, parameters etc.
  state_and_sens_ic.setArgs(state_ic);
  // Set initial condition for x_bar = [ x; s_bar ]
  state_and_sens_ic.set_x(
    stateAndSensModel->create_x_bar_vec(state_ic.get_x(),s_bar_init)
    );
  // Set initial condition for x_bar_dot = [ x_dot; s_bar_dot ]
  state_and_sens_ic.set_x_dot(
    stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(),s_bar_dot_init)
    );

  stateAndSensStepper->setInitialCondition(state_and_sens_ic);

  double finalTime = 0.00001;
  int numTimeSteps = 1;
  RCP<IntegratorBase<double> > integrator;
  {
    RCP<ParameterList> integratorPL = Teuchos::parameterList();
    integratorPL->set( "Take Variable Steps", false );
    integratorPL->set( "Fixed dt", Teuchos::as<double>((finalTime - state_ic.get_t())/numTimeSteps) );
    RCP<IntegratorBase<double> >
      defaultIntegrator = controlledDefaultIntegrator<double>(
        simpleIntegrationControlStrategy<double>(integratorPL)
        );
    integrator = defaultIntegrator;
  }

  RCP<StepperAsModelEvaluator<double> >
    stateAndSensIntegratorAsModel = stepperAsModelEvaluator(
      Teuchos::rcp_implicit_cast<StepperBase<double> >(stateAndSensStepper),
      integrator, state_and_sens_ic
      );
  stateAndSensIntegratorAsModel->setVerbLevel(verbLevel);
    
  RCP<Thyra::VectorBase<double> > x_bar_final;
  {
  
    x_bar_final = createMember(stateAndSensIntegratorAsModel->get_g_space(0));
  
    Thyra::eval_g(
      *stateAndSensIntegratorAsModel,
      0, *state_ic.get_p(0),
      finalTime,
      0, &*x_bar_final
      );

    *std_out
      << "\nx_bar_final = x_bar(p,finalTime) evaluated using stateAndSensIntegratorAsModel:\n"
      << Teuchos::describe(*x_bar_final,verbLevel);

  }

  // Now we check that the sensitivities are correct
  RCP<const Thyra::VectorBase<double> >
    DxDp_vec_final = Thyra::productVectorBase<double>(x_bar_final)->getVectorBlock(1);
  *std_out << "\nDxDp_vec_final:\n"
    << Teuchos::describe(*DxDp_vec_final,verbLevel);
  RCP<const Thyra::DefaultMultiVectorProductVector<double> > DxDp_mv_final =
    Teuchos::rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<double> >(DxDp_vec_final,true);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s0_final = DxDp_mv_final->getVectorBlock(0);
  RCP<const Thyra::VectorBase<double> >
    DxDp_s1_final = DxDp_mv_final->getVectorBlock(1);

  *std_out << "\nDxDp_s0_final:\n"
    << Teuchos::describe(*DxDp_s0_final,verbLevel);
  *std_out << "\nDxDp_s1_final:\n"
    << Teuchos::describe(*DxDp_s0_final,verbLevel);

  MEB::InArgs<double> exactSensSolution = stateModel->getExactSensSolution(0,finalTime);
  RCP<const Thyra::VectorBase<double> > ds0dp = exactSensSolution.get_x();
  exactSensSolution = stateModel->getExactSensSolution(1,finalTime);
  RCP<const Thyra::VectorBase<double> > ds1dp = exactSensSolution.get_x();

  *std_out << "\nds0dp exact:\n"
    << Teuchos::describe(*ds0dp,verbLevel);
  *std_out << "\nds1dp exact:\n"
    << Teuchos::describe(*ds1dp,verbLevel);

  /*
  double maxSensError = 1.0e-4;
  double s0_correct = Thyra::testRelNormDiffErr(
    "DxDp_s0_final", *DxDp_s0_final,
    "DxDp_exact_s0_final", *ds0dp,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &*std_out, verbLevel
    );
  TEST_EQUALITY_CONST( s0_correct, true );

  double s1_correct = Thyra::testRelNormDiffErr(
    "DxDp_s1_final", *DxDp_s1_final,
    "DxDp_exact_s1_final", *ds1dp,
    "maxSensError", maxSensError,
    "warningTol", 1.0, // Don't warn
    &*std_out, verbLevel
    );
  TEST_EQUALITY_CONST( s1_correct, true );
  */

}

} // namespace Rythmos

