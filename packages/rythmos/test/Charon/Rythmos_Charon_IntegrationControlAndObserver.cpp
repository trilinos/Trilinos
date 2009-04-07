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

#include "Rythmos_Charon_IntegrationControlAndObserver.hpp"
#include "Rythmos_Charon_ImplicitBDFStepperErrWtVecCalc.hpp"
#include "Rythmos_Charon_ImplicitBDFStepperStepControl.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ImplicitBDFStepperErrWtVecCalc.hpp"
#include "Rythmos_ImplicitBDFStepperStepControl.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"

#include "Thyra_ModelEvaluator.hpp"

namespace RythmosCharon {

using Teuchos::is_null;

TEUCHOS_UNIT_TEST( Rythmos_Charon, create ) {
  RCP<CharonIntegrationControlAndObserver> icao = Teuchos::rcp(new CharonIntegrationControlAndObserver());
  TEST_ASSERT( !is_null(icao) );

  RCP<CharonImplicitBDFStepperErrWtVecCalc> ibsewvc = Teuchos::rcp(new CharonImplicitBDFStepperErrWtVecCalc());
  TEST_ASSERT( !is_null(ibsewvc) );

  RCP<CharonImplicitBDFStepperStepControl<double> > ibssc = Teuchos::rcp(new CharonImplicitBDFStepperStepControl<double>());
  TEST_ASSERT( !is_null(ibssc) );
}


TEUCHOS_UNIT_TEST( Rythmos_Charon, ImplicitBDF_RythmosControl ) {
  RCP<Rythmos::ImplicitBDFStepper<double> > stepper;
  stepper = Rythmos::implicitBDFStepper<double>();

  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

  RCP<Rythmos::ImplicitBDFStepperErrWtVecCalc<double> > errWtVecCalc = 
    Teuchos::rcp(new Rythmos::ImplicitBDFStepperErrWtVecCalc<double> ());
  errWtVecCalc->setParameterList(pl);
  RCP<Rythmos::ImplicitBDFStepperStepControl<double> > stepControl = 
    Teuchos::rcp(new Rythmos::ImplicitBDFStepperStepControl<double> ());
  stepControl->setErrWtVecCalc(errWtVecCalc);
  stepControl->setParameterList(pl);
  stepper->setStepControlStrategy(stepControl);
  TEST_ASSERT( !is_null(stepper) );
}

TEUCHOS_UNIT_TEST( Rythmos_Charon, ImplicitBDF_CharonControl ) {
  RCP<Rythmos::ImplicitBDFStepper<double> > stepper;
  stepper = Rythmos::implicitBDFStepper<double>();

  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

  RCP<RythmosCharon::CharonImplicitBDFStepperErrWtVecCalc> errWtVecCalc = 
    Teuchos::rcp(new RythmosCharon::CharonImplicitBDFStepperErrWtVecCalc());
  errWtVecCalc->setParameterList(pl);
  RCP<RythmosCharon::CharonImplicitBDFStepperStepControl<double> > stepControl = 
    Teuchos::rcp(new RythmosCharon::CharonImplicitBDFStepperStepControl<double>());
  stepControl->setErrWtVecCalc(errWtVecCalc);
  stepControl->setParameterList(pl);
  stepper->setStepControlStrategy(stepControl);
  TEST_ASSERT( !is_null(stepper) );
}

TEUCHOS_UNIT_TEST( Rythmos_Charon, BackwardEuler ) {
  RCP<Rythmos::BackwardEulerStepper<double> > stepper;
  stepper = Rythmos::backwardEulerStepper<double>();
  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  stepper->setParameterList(pl);
  TEST_ASSERT( !is_null(stepper) );
}

TEUCHOS_UNIT_TEST( Rythmos_Charon, integrator_BE ) {
  RCP<RythmosCharon::CharonIntegrationControlAndObserver> controlAndObserver = 
    Teuchos::rcp(new RythmosCharon::CharonIntegrationControlAndObserver());
  RCP<Rythmos::IntegratorBase<double> > integrator;
  {
    RCP<Rythmos::DefaultIntegrator<double> > dInt = 
      Rythmos::defaultIntegrator<double>(
          controlAndObserver,
          controlAndObserver
          );
    RCP<Rythmos::SimpleIntegrationControlStrategy<double> > intCont = 
      Rythmos::simpleIntegrationControlStrategy<double>();
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Take Variable Steps",false);
    pl->set("Fixed dt",0.1);
    intCont->setParameterList(pl);
    dInt->setIntegrationControlStrategy(intCont);
    integrator = dInt;
  }
  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->sublist("VerboseObject").set("Verbosity Level","none");
  integrator->setParameterList(pl);
  //integrator->setOStream(out);
  //integrator->setVerbLevel(Teuchos::VERB_MEDIUM);
  // Model
  RCP<Thyra::ModelEvaluator<double> > model;
  model = Rythmos::sinCosModel(true);
  // Solver
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
  nlSolver = Rythmos::timeStepNonlinearSolver<double>();
  // IC
  Thyra::ModelEvaluatorBase::InArgs<double> ic;
  ic = model->getNominalValues();
  // Stepper
  RCP<Rythmos::StepperBase<double> > stepper;
  {
    RCP<Rythmos::BackwardEulerStepper<double> > beStepper = Rythmos::backwardEulerStepper<double>();
    beStepper->setModel(model);
    beStepper->setSolver(nlSolver);
    beStepper->setInitialCondition(ic);
    RCP<Teuchos::ParameterList> bePL = Teuchos::parameterList();
    beStepper->setParameterList(bePL);
    stepper = beStepper;
  }

  double finalTime = 1.0;
  integrator->setStepper(stepper,finalTime);
  integrator->getFwdPoints(
      Teuchos::tuple<double>(finalTime),
      0, // x_vec
      0, // xdot_vec
      0  // accuracy_vec
      );
  TEST_ASSERT( !is_null(integrator) );
}

TEUCHOS_UNIT_TEST( Rythmos_Charon, integrator_IBDF ) {
  RCP<RythmosCharon::CharonIntegrationControlAndObserver> controlAndObserver = 
    Teuchos::rcp(new RythmosCharon::CharonIntegrationControlAndObserver());;
  RCP<Rythmos::IntegratorBase<double> > integrator; 
  {
    RCP<Rythmos::DefaultIntegrator<double> > dInt = 
      Rythmos::defaultIntegrator<double>(
          controlAndObserver,
          controlAndObserver
          );
    RCP<Rythmos::SimpleIntegrationControlStrategy<double> > intCont = 
      Rythmos::simpleIntegrationControlStrategy<double>();
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Take Variable Steps",false);
    pl->set("Fixed dt",0.1);
    intCont->setParameterList(pl);
    dInt->setIntegrationControlStrategy(intCont);
    integrator = dInt;
  }
  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->sublist("VerboseObject").set("Verbosity Level","none");
  integrator->setParameterList(pl);
  //integrator->setOStream(out);
  //integrator->setVerbLevel(Teuchos::VERB_MEDIUM);
  // Model
  RCP<Thyra::ModelEvaluator<double> > model;
  model = Rythmos::sinCosModel(true);
  // Solver
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
  nlSolver = Rythmos::timeStepNonlinearSolver<double>();
  // IC
  Thyra::ModelEvaluatorBase::InArgs<double> ic;
  ic = model->getNominalValues();
  // Stepper
  RCP<Rythmos::StepperBase<double> > stepper;
  {
    RCP<Rythmos::ImplicitBDFStepper<double> > ibdfStepper = Rythmos::implicitBDFStepper<double>();
    RCP<Teuchos::ParameterList> errWtPL = Teuchos::parameterList();
    RCP<Rythmos::ImplicitBDFStepperErrWtVecCalc<double> > errWtVecCalc = 
      Teuchos::rcp(new Rythmos::ImplicitBDFStepperErrWtVecCalc<double> ());
    errWtVecCalc->setParameterList(errWtPL);
    RCP<Teuchos::ParameterList> stepContPL = Teuchos::parameterList();
    RCP<Rythmos::ImplicitBDFStepperStepControl<double> > stepControl = 
      Teuchos::rcp(new Rythmos::ImplicitBDFStepperStepControl<double>());
    stepControl->setErrWtVecCalc(errWtVecCalc);
    stepControl->setParameterList(stepContPL);
    ibdfStepper->setStepControlStrategy(stepControl);
    ibdfStepper->setModel(model);
    ibdfStepper->setSolver(nlSolver);
    ibdfStepper->setInitialCondition(ic);
    RCP<Teuchos::ParameterList> ibdfPL = Teuchos::parameterList();
    ibdfStepper->setParameterList(ibdfPL);
    stepper = ibdfStepper;
  }

  double finalTime = 1.0;
  integrator->setStepper(stepper,finalTime);
  integrator->getFwdPoints(
      Teuchos::tuple<double>(finalTime),
      0, // x_vec
      0, // xdot_vec
      0  // accuracy_vec
      );
  TEST_ASSERT( !is_null(integrator) );
}

} // namespace RythmosCharon



