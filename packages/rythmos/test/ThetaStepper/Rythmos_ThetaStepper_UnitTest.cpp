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
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Rythmos_Types.hpp"
#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_ThetaStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "../SinCos/SinCosModel.hpp"

#include "Rythmos_StepperBuilder.hpp"

namespace {
  const std::string StepperType_name = "Stepper Type";
}

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_ThetaStepper, createImplicitEuler) {
  // Model
  RCP<SinCosModel> model = sinCosModel(true);

  RCP<ParameterList> modelPL = 
    Teuchos::getParametersFromXmlFile("modelParams.xml");
  modelPL->validateParametersAndSetDefaults(*model->getValidParameters());

  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();

  // Solver
  RCP<TimeStepNonlinearSolver<double> > solver = 
    timeStepNonlinearSolver<double>();

  // test reading params from .xml file
  RCP<ParameterList> stepperParamList = 
    Teuchos::getParametersFromXmlFile("implicitEulerParams.xml");

  // Stepper
  RCP<ThetaStepper<double> > stepper = 
    thetaStepper<double>(model, solver, stepperParamList);
  TEST_ASSERT( !is_null(stepper) );

  stepper->setInitialCondition(ic);

  // for testing purposes only
  stepper->setVerbLevel(Teuchos::VERB_EXTREME);

  double dt = 1.0;
  double dt_taken = 0.0;
  TEST_NOTHROW( dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED) );
  TEST_EQUALITY_CONST( dt_taken, dt );
}

TEUCHOS_UNIT_TEST( Rythmos_ThetaStepper, createTrapezoid) {
  // Model
  RCP<SinCosModel> model = sinCosModel(true);

  RCP<ParameterList> modelPL = 
    Teuchos::getParametersFromXmlFile("modelParams.xml");
  modelPL->validateParametersAndSetDefaults(*model->getValidParameters());

  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();

  // Solver
  RCP<TimeStepNonlinearSolver<double> > solver = 
    timeStepNonlinearSolver<double>();

  // test generating internal param list
  RCP<ParameterList> stepperParamList = Teuchos::parameterList();
  ParameterList& pl = stepperParamList->sublist("Step Control Settings");
  pl.set("Theta Stepper Type", "Trapezoid");

  //RCP<ParameterList> stepperParamList = 
  //  Teuchos::getParametersFromXmlFile("trapezoidParams.xml");

  // Stepper
  RCP<ThetaStepper<double> > stepper = 
    thetaStepper<double>(model, solver, stepperParamList);
  TEST_ASSERT( !is_null(stepper) );

  stepper->setInitialCondition(ic);

  // for testing purposes only
  stepper->setVerbLevel(Teuchos::VERB_EXTREME);

  double dt = 1.0;
  double dt_taken = 0.0;
  TEST_NOTHROW( dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED) );
  TEST_EQUALITY_CONST( dt_taken, dt );

  // take extra steps to verify 2nd order method is working
  TEST_NOTHROW( dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED) );
  TEST_EQUALITY_CONST( dt_taken, dt );

  TEST_NOTHROW( dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED) );
  TEST_EQUALITY_CONST( dt_taken, dt );
}

TEUCHOS_UNIT_TEST( Rythmos_ThetaStepper, createThetaStepper ) {
  // Verify the builder operates correctly for ThetaStepper
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  {
    // Specify which stepper we want
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(StepperType_name, "Theta");
    // Specify a BackwardEuler setting
    RCP<ParameterList> tsSettings = Teuchos::sublist(pl,"Theta");
    RCP<ParameterList> vopl = Teuchos::sublist(tsSettings,"VerboseObject");
    vopl->set("Verbosity Level","none");
    builder->setParameterList(pl);
  }

  // Create the stepper
  RCP<StepperBase<double> > stepper = builder->create();
  TEST_EQUALITY( is_null(stepper), false );
  // Verify we got the correct stepper
  RCP<ThetaStepper<double> > thetaStepper = 
    Teuchos::rcp_dynamic_cast<ThetaStepper<double> >(stepper,false);
  TEST_EQUALITY( is_null(thetaStepper), false );
  // Verify appropriate settings have propagated into the stepper correctly
  Teuchos::EVerbosityLevel verbLevel = thetaStepper->getVerbLevel();
  TEST_EQUALITY( verbLevel, Teuchos::VERB_NONE );
}

} // namespace Rythmos

