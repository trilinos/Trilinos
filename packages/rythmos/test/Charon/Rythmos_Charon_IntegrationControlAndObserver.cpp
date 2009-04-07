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

#ifdef HAVE_RYTHMOS_EXPERIMENTAL

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_Types.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

#include "Rythmos_Charon_IntegrationControlAndObserver.hpp"

namespace RythmosCharon {

TEUCHOS_UNIT_TEST( Rythmos_Charon, create ) {
  Teuchos::RCP<CharonIntegrationControlAndObserver> icao = Teuchos::rcp(new CharonIntegrationControlAndObserver());
  TEST_ASSERT( !is_null(icao) );

  Teuchos::RCP<CharonImplicitBDFStepperErrWtVecCalc> ibsewvc = Teuchos::rcp(new CharonImplicitBDFStepperErrWtVecCalc());
  TEST_ASSERT( !is_null(ibsewvc) );

  Teuchos::RCP<CharonImplicitBDFStepperStepControl> ibssc = Teuchos::rcp(new CharonImplicitBDFStepperStepControl());
  TEST_ASSERT( !is_null(ibssc) );
}


TEUCHOS_UNIT_TEST( Rythmos_Charon, ImplicitBDF_RythmosControl ) {
  Teuchos::RCP<Rythmos::SolverAcceptingStepperBase<double> > stepper;
  stepper = Rythmos::implicitBDFStepper<double>();

  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::ParameterList();

  Teuchos::RCP<Rythmos::ImplicitBDFStepperErrWtVecCalc> errWtVecCalc = 
    Teuchos::rcp(new Rythmos::ImplicitBDFStepperErrWtVecCalc());
  errWtVecCalc->setParameterList(pl);
  Teuchos::RCP<Rythmos::ImplicitBDFStepperStepControl> stepControl = 
    Teuchos::rcp(new Rythmos::ImplicitBDFStepperStepControl());
  stepControl->setErrWtVecCalc(errWtVecCalc);
  stepControl->setParameterList(pl);
  stepper->setStepControlStrategy(stepControl);
  TEST_ASSERT( !is_null(stepper) );
}
TEUCHOS_UNIT_TEST( Rythmos_Charon, ImplicitBDF_CharonControl ) {
  Teuchos::RCP<Rythmos::SolverAcceptingStepperBase<double> > stepper;
  stepper = Rythmos::implicitBDFStepper<double>();

  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::ParameterList();

  Teuchos::RCP<RythmosCharon::CharonImplicitBDFStepperErrWtVecCalc> errWtVecCalc = 
    Teuchos::rcp(new RythmosCharon::CharonImplicitBDFStepperErrWtVecCalc());
  errWtVecCalc->setParameterList(pl);
  Teuchos::RCP<RythmosCharon::CharonImplicitBDFStepperStepControl> stepControl = 
    Teuchos::rcp(new RythmosCharon::CharonImplicitBDFStepperStepControl());
  stepControl->setErrWtVecCalc(errWtVecCalc);
  stepControl->setParameterList(pl);
  stepper->setStepControlStrategy(stepControl);
  TEST_ASSERT( !is_null(stepper) );
}

TEUCHOS_UNIT_TEST( Rythmos_Charon, BackwardEuler ) {
  Teuchos::RCP<Rythmos::BackwardEulerStepper<double> > stepper;
  stepper = Rythmos::backwardEulerStepper<double>();
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::ParameterList();
  stepper->setParameterList(pl);
  TEST_ASSERT( !is_null(stepper) );
}

TEUCHOS_UNIT_TEST( Rythmos_Charon, integrator_BE ) {
  Teuchos::RCP<RythmosCharon::CharonIntegrationControlAndObserver> controlAndObserver = 
    Teuchos::rcp(new RythmosCharon::CharonIntegrationControlAndObserver());;
  Teuchos::RCP<Rythmos::IntegratorBase<double> > integrator = 
    Rythmos::defaultIntegrator<double>(
        controlAndObserver,
        controlAndObserver
        );
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  integrator->setParameterList(pl);
  integrator->setOStream(out);
  integrator->setVerbLevel(Teuchos::VERB_MEDIUM);
  // Model
  Teuchos::RCP<Thyra::ModelEvaluator> model;
  model = Rythmos::sinCosModel(true);
  // Solver
  Teuchos::RCP<Thyra::NonlinearSolverBase<double> nlSolver;
  nlSolver = Rythmos::timeStepNonlinearSolver<double>();
  // IC
  Thyra::ModelEvaluatorBase::InArgs<double> ic;
  ic = model->getNominalValues();
  // Stepper
  Teuchos::RCP<Rythmos::StepperBase<double> > stepper;
  stepper = Rythmos::backwardEulerStepper<double>();
  stepper->setModel(model);
  stepper->setSolver(nlSolver);
  stepper->setInitialCondition(ic);

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
  Teuchos::RCP<RythmosCharon::CharonIntegrationControlAndObserver> controlAndObserver = 
    Teuchos::rcp(new RythmosCharon::CharonIntegrationControlAndObserver());;
  Teuchos::RCP<Rythmos::IntegratorBase<double> > integrator = 
    Rythmos::defaultIntegrator<double>(
        controlAndObserver,
        controlAndObserver
        );
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  integrator->setParameterList(pl);
  integrator->setOStream(out);
  integrator->setVerbLevel(Teuchos::VERB_MEDIUM);
  // Model
  Teuchos::RCP<Thyra::ModelEvaluator> model;
  model = Rythmos::sinCosModel(true);
  // Solver
  Teuchos::RCP<Thyra::NonlinearSolverBase<double> nlSolver;
  nlSolver = Rythmos::timeStepNonlinearSolver<double>();
  // IC
  Thyra::ModelEvaluatorBase::InArgs<double> ic;
  ic = model->getNominalValues();
  // Stepper
  Teuchos::RCP<Rythmos::StepperBase<double> > stepper;
  stepper = Rythmos::implicitBDFStepper<double>();
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::ParameterList();
    Teuchos::RCP<Rythmos::ImplicitBDFStepperErrWtVecCalc> errWtVecCalc = 
      Teuchos::rcp(new Rythmos::ImplicitBDFStepperErrWtVecCalc());
    errWtVecCalc->setParameterList(pl);
    Teuchos::RCP<Rythmos::ImplicitBDFStepperStepControl> stepControl = 
      Teuchos::rcp(new Rythmos::ImplicitBDFStepperStepControl());
    stepControl->setErrWtVecCalc(errWtVecCalc);
    stepControl->setParameterList(pl);
    stepper->setStepControlStrategy(stepControl);
  }
  stepper->setModel(model);
  stepper->setSolver(nlSolver);
  stepper->setInitialCondition(ic);

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

#endif // HAVE_RYTHMOS_EXPERIMENTAL

