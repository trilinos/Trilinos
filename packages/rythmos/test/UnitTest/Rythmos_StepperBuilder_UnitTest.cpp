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
#include "Teuchos_ParameterList.hpp"

#include "Rythmos_StepperBuilder.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, setParameterList ) {
  RCP<ParameterList> pl = Teuchos::parameterList();
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  TEST_NOTHROW(builder->setParameterList(pl));
  // Test that StepperBuilder validates its parameter list
  pl->set("Hello","World");
  TEST_THROW(builder->setParameterList(pl), std::logic_error);
}

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, getNonconstParameterList ) {
  StepperBuilder<double> builder;
  RCP<ParameterList> pl = builder.getNonconstParameterList();
  TEST_EQUALITY( is_null(pl), true );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, getValidParameters ) {
  StepperBuilder<double> builder;
  RCP<const ParameterList> pl = builder.getValidParameters();
  TEST_EQUALITY( is_null(pl), false );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, unsetParameterList ) {
  StepperBuilder<double> builder;
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set(StepperType_name, ForwardEuler_name);
  builder.setParameterList(pl);
  RCP<ParameterList> oldPL = builder.unsetParameterList();
  // Did I get my parameter list back?
  TEST_EQUALITY( oldPL.get(), pl.get() );
  // Does the builder now have a null parameter list?
  TEST_EQUALITY( is_null(builder.getParameterList()), true );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, createBEStepper ) {
  // Verify the builder operates correctly for Backward Euler
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  {
    // Specify which stepper we want
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(StepperType_name, BackwardEuler_name);
    // Specify a BackwardEuler setting
    RCP<ParameterList> beSettings = Teuchos::sublist(pl,BackwardEulerSettings_name);
    RCP<ParameterList> vopl = Teuchos::sublist(beSettings,"VerboseObject");
    vopl->set("Verbosity Level","none");
    builder->setParameterList(pl);
  }
  // Verify we get an exception when we try to create the BE stepper without a nonlinear solver
  TEST_THROW( builder->create(), std::logic_error );
  {
    // Specify a nonlinear solver
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
    builder->setNonlinearSolver(nlSolver);
  }
  // Create the stepper
  RCP<StepperBase<double> > stepper = builder->create();
  TEST_EQUALITY( is_null(stepper), false );
  // Verify we got the correct stepper
  RCP<BackwardEulerStepper<double> > beStepper = Teuchos::rcp_dynamic_cast<BackwardEulerStepper<double> >(stepper,false);
  TEST_EQUALITY( is_null(beStepper), false );
  // Verify it has the correct nonlinear solver
  RCP<const Thyra::NonlinearSolverBase<double> > nlSolver = beStepper->getSolver();
  RCP<const TimeStepNonlinearSolver<double> > tsnlSolver = Teuchos::rcp_dynamic_cast<const TimeStepNonlinearSolver<double> >(nlSolver,false);
  TEST_EQUALITY( is_null(tsnlSolver), false );
  // Verify appropriate settings have propagated into the stepper correctly
  Teuchos::EVerbosityLevel verbLevel = beStepper->getVerbLevel();
  TEST_EQUALITY( verbLevel, Teuchos::VERB_NONE );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, createIBDFStepper ) {
  // Verify the builder operates correctly for IBDF Stepper
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  {
    // Specify which stepper we want
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(StepperType_name, ImplicitBDF_name);
    // Specify a IBDF setting
    RCP<ParameterList> ibdfSettings = Teuchos::sublist(pl,ImplicitBDFSettings_name);
    RCP<ParameterList> vopl = Teuchos::sublist(ibdfSettings,"VerboseObject");
    vopl->set("Verbosity Level","none");
    builder->setParameterList(pl);
  }
  // Verify we get an exception when we try to create the IBDF stepper without a nonlinear solver
  TEST_THROW( builder->create(), std::logic_error );
  {
    // Specify a nonlinear solver
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
    builder->setNonlinearSolver(nlSolver);
  }
  // Create the stepper
  RCP<StepperBase<double> > stepper = builder->create();
  TEST_EQUALITY( is_null(stepper), false );
  // Verify we got the correct stepper
  RCP<ImplicitBDFStepper<double> > ibdfStepper = Teuchos::rcp_dynamic_cast<ImplicitBDFStepper<double> >(stepper,false);
  TEST_EQUALITY( is_null(ibdfStepper), false );
  // Verify it has the correct nonlinear solver
  RCP<const Thyra::NonlinearSolverBase<double> > nlSolver = ibdfStepper->getSolver();
  RCP<const TimeStepNonlinearSolver<double> > tsnlSolver = Teuchos::rcp_dynamic_cast<const TimeStepNonlinearSolver<double> >(nlSolver,false);
  TEST_EQUALITY( is_null(tsnlSolver), false );
  // Verify appropriate settings have propagated into the stepper correctly
  Teuchos::EVerbosityLevel verbLevel = ibdfStepper->getVerbLevel();
  TEST_EQUALITY( verbLevel, Teuchos::VERB_NONE );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, createFEStepper ) {
  // Verify the builder operates correctly for FE Stepper
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  {
    // Specify which stepper we want
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(StepperType_name, ForwardEuler_name);
    // Specify a FE setting
    RCP<ParameterList> feSettings = Teuchos::sublist(pl,ForwardEulerSettings_name);
    RCP<ParameterList> vopl = Teuchos::sublist(feSettings,"VerboseObject");
    vopl->set("Verbosity Level","none");
    builder->setParameterList(pl);
  }
  // Create the stepper
  RCP<StepperBase<double> > stepper = builder->create();
  TEST_EQUALITY( is_null(stepper), false );
  // Verify we got the correct stepper
  RCP<ForwardEulerStepper<double> > feStepper = Teuchos::rcp_dynamic_cast<ForwardEulerStepper<double> >(stepper,false);
  TEST_EQUALITY( is_null(feStepper), false );
  // Verify appropriate settings have propagated into the stepper correctly
  Teuchos::EVerbosityLevel verbLevel = feStepper->getVerbLevel();
  TEST_EQUALITY( verbLevel, Teuchos::VERB_NONE );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, createERKStepper ) {
  // Verify the builder operates correctly for ERK Stepper
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  {
    // Specify which stepper we want
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(StepperType_name, ExplicitRK_name);
    // Specify a ERK setting
    RCP<ParameterList> erkSettings = Teuchos::sublist(pl,ExplicitRKSettings_name);
    RCP<ParameterList> vopl = Teuchos::sublist(erkSettings,"VerboseObject");
    vopl->set("Verbosity Level","none");
    builder->setParameterList(pl);
  }
  // Create the stepper
  RCP<StepperBase<double> > stepper = builder->create();
  TEST_EQUALITY( is_null(stepper), false );
  // Verify we got the correct stepper
  RCP<ExplicitRKStepper<double> > erkStepper = Teuchos::rcp_dynamic_cast<ExplicitRKStepper<double> >(stepper,false);
  TEST_EQUALITY( is_null(erkStepper), false );
  // Verify appropriate settings have propagated into the stepper correctly
  Teuchos::EVerbosityLevel verbLevel = erkStepper->getVerbLevel();
  TEST_EQUALITY( verbLevel, Teuchos::VERB_NONE );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, createIRKStepper ) {
  // Verify the builder operates correctly for IRK Stepper
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  {
    // Specify which stepper we want
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(StepperType_name, ImplicitRK_name);
    // Specify a IRK setting
    RCP<ParameterList> irkSettings = Teuchos::sublist(pl,ImplicitRKSettings_name);
    RCP<ParameterList> vopl = Teuchos::sublist(irkSettings,"VerboseObject");
    vopl->set("Verbosity Level","none");
    builder->setParameterList(pl);
  }
  // Verify we get an exception when we try to create the IRK stepper without a nonlinear solver
  TEST_THROW( builder->create(), std::logic_error );
  {
    // Specify a nonlinear solver
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
    builder->setNonlinearSolver(nlSolver);
  }
  // Create the stepper
  RCP<StepperBase<double> > stepper = builder->create();
  TEST_EQUALITY( is_null(stepper), false );
  // Verify we got the correct stepper
  RCP<ImplicitRKStepper<double> > irkStepper = Teuchos::rcp_dynamic_cast<ImplicitRKStepper<double> >(stepper,false);
  TEST_EQUALITY( is_null(irkStepper), false );
  // Verify it has the correct nonlinear solver
  RCP<const Thyra::NonlinearSolverBase<double> > nlSolver = irkStepper->getSolver();
  RCP<const TimeStepNonlinearSolver<double> > tsnlSolver = Teuchos::rcp_dynamic_cast<const TimeStepNonlinearSolver<double> >(nlSolver,false);
  TEST_EQUALITY( is_null(tsnlSolver), false );
  // Verify appropriate settings have propagated into the stepper correctly
  Teuchos::EVerbosityLevel verbLevel = irkStepper->getVerbLevel();
  TEST_EQUALITY( verbLevel, Teuchos::VERB_NONE );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperBuilder, createETPStepper ) {
  // Verify the builder operates correctly for ETP Stepper
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  {
    // Specify which stepper we want
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(StepperType_name, ExplicitTP_name);
    // Specify a ETP setting
    RCP<ParameterList> etpSettings = Teuchos::sublist(pl,ExplicitTPSettings_name);
    RCP<ParameterList> vopl = Teuchos::sublist(etpSettings,"VerboseObject");
    vopl->set("Verbosity Level","none");
    builder->setParameterList(pl);
  }
  // Create the stepper
  RCP<StepperBase<double> > stepper = builder->create();
  TEST_EQUALITY( is_null(stepper), false );
  // Verify we got the correct stepper
  RCP<ExplicitTaylorPolynomialStepper<double> > etpStepper = Teuchos::rcp_dynamic_cast<ExplicitTaylorPolynomialStepper<double> >(stepper,false);
  TEST_EQUALITY( is_null(etpStepper), false );
  // Verify appropriate settings have propagated into the stepper correctly
  Teuchos::EVerbosityLevel verbLevel = etpStepper->getVerbLevel();
  TEST_EQUALITY( verbLevel, Teuchos::VERB_NONE );
}

} // namespace Rythmos 


