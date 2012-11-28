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

#include "Rythmos_IntegratorBuilder.hpp"
#include "Rythmos_IntegratorBuilder_Helpers.hpp"
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_IntegrationControlStrategyAcceptingIntegratorBase.hpp"
#include "Rythmos_ImplicitBDFStepperStepControl.hpp"
#include "Rythmos_ImplicitBDFStepperErrWtVecCalc.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_PointwiseInterpolationBufferAppender.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_HermiteInterpolator.hpp"
#include "Rythmos_CubicSplineInterpolator.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_UnitTestModels.hpp"
#include "Rythmos_RKButcherTableau.hpp"

#ifdef Rythmos_ENABLE_NOX
#  include "Thyra_NonlinearSolver_NOX.hpp"
#endif

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, construct ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  TEST_ASSERT( !is_null(ib) );
  TEST_NOTHROW( ib = Teuchos::null );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, constructFoolish ) {
  RCP<FoolishIntegrator> fi;
  TEST_NOTHROW( fi = rcp(new FoolishIntegrator()) );

  RCP<FoolishIntegrationControlStrategy> fics;
  TEST_NOTHROW( fics = rcp(new FoolishIntegrationControlStrategy()) );

  RCP<FoolishStepper> fs;
  TEST_NOTHROW( fs = rcp(new FoolishStepper()) );

  RCP<FoolishStepControlStrategy> fscs;
  TEST_NOTHROW( fscs = rcp(new FoolishStepControlStrategy()) );

  RCP<FoolishInterpolationBuffer> fib;
  TEST_NOTHROW( fib = rcp(new FoolishInterpolationBuffer()) );

  RCP<FoolishInterpolationBufferAppender> fiba;
  TEST_NOTHROW( fiba = rcp(new FoolishInterpolationBufferAppender()) );

  RCP<FoolishErrWtVecCalc> fewvc;
  TEST_NOTHROW( fewvc = rcp(new FoolishErrWtVecCalc()) );

  RCP<FoolishInterpolator> fInterp;
  TEST_NOTHROW( fInterp = rcp(new FoolishInterpolator()) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setIntegratorFactory ) {

  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setIntegratorFactory(
    Teuchos::abstractFactoryStd< IntegratorBase<double>, FoolishIntegrator >(),
    "Foolish Integrator");
  ib->setIntegratorFactory(
    Teuchos::abstractFactoryStd< IntegratorBase<double>,
                                 DefaultIntegrator<double> >(),
    "Other Default Integrator");

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Integrator Settings")
     .sublist("Integrator Selection")
     .set("Integrator Type","Foolish Integrator");
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Explicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type","Explicit 4 Stage");
  ib->setParameterList(pl);
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  TEST_NOTHROW( integrator = ib->create(model,ic,nlSolver) );
  RCP<FoolishIntegrator> fInt =
    Teuchos::rcp_dynamic_cast<FoolishIntegrator>(integrator, false);
  TEST_ASSERT( !is_null(fInt) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setIntegratorFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW( 
    ib->setIntegratorFactory(
      Teuchos::abstractFactoryStd< IntegratorBase<double>,
                                   FoolishIntegrator >(),
      "Default Integrator"),
    std::logic_error);
#else // TEUCHOS_DEBUG
  TEST_NOTHROW( 
    ib->setIntegratorFactory(
      Teuchos::abstractFactoryStd< IntegratorBase<double>,
                                   FoolishIntegrator >(),
      "Default Integrator") );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setIntegrationControlFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setIntegrationControlFactory(
    Teuchos::abstractFactoryStd< IntegrationControlStrategyBase<double>,
                                 FoolishIntegrationControlStrategy >(),
    "Foolish Integration Control");
  ib->setIntegrationControlFactory(
    Teuchos::abstractFactoryStd< IntegrationControlStrategyBase<double>,
                                 SimpleIntegrationControlStrategy<double> >(),
    "Other Simple Integration Control");
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Integration Control Strategy Selection")
     .set("Integration Control Strategy Type","Foolish Integration Control");
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Explicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type","Explicit 4 Stage");
  ib->setParameterList(pl);
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  TEST_NOTHROW( integrator = ib->create(model,ic,nlSolver) );
  RCP<IntegrationControlStrategyAcceptingIntegratorBase<double> > specialInt = 
    Teuchos::rcp_dynamic_cast<
      IntegrationControlStrategyAcceptingIntegratorBase<double> >(integrator,
                                                                  false);
  TEST_ASSERT( !is_null(specialInt) );
  RCP<const IntegrationControlStrategyBase<double> > ibControl = 
    specialInt->getIntegrationControlStrategy();
  TEST_ASSERT( !is_null(ibControl) );
  RCP<const FoolishIntegrationControlStrategy> specialIBControl = 
    Teuchos::rcp_dynamic_cast<
      const FoolishIntegrationControlStrategy>(ibControl,false);
  TEST_ASSERT( !is_null(specialIBControl) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setIntegrationControlFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW( 
    ib->setIntegrationControlFactory(
      Teuchos::abstractFactoryStd< IntegrationControlStrategyBase<double>,
                                   FoolishIntegrationControlStrategy >(),
      "Simple Integration Control Strategy"),
    std::logic_error);
#else // TEUCHOS_DEBUG
  TEST_NOTHROW( 
    ib->setIntegrationControlFactory(
      Teuchos::abstractFactoryStd< IntegrationControlStrategyBase<double>,
                                   FoolishIntegrationControlStrategy >(),
      "Simple Integration Control Strategy") );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setStepperBuilder ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<StepperBuilder<double> > sb = stepperBuilder<double>();
  sb->setStepperFactory(
    Teuchos::abstractFactoryStd< StepperBase<double>, FoolishStepper >(),
    "Foolish Stepper");
  sb->setStepperFactory(
    Teuchos::abstractFactoryStd< StepperBase<double>,
                                 BackwardEulerStepper<double> >(),
    "Other Backward Euler Stepper");
  ib->setStepperBuilder(sb);

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Foolish Stepper");
  ib->setParameterList(pl);
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  integrator = ib->create(model,ic,nlSolver);
  RCP<const StepperBase<double> > stepper = integrator->getStepper();
  TEST_ASSERT( !is_null(stepper) );
  RCP<const FoolishStepper> fStepper =
    Teuchos::rcp_dynamic_cast<const FoolishStepper>(stepper,false);
  TEST_ASSERT( !is_null(fStepper) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setStepControlFactory ) {

  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setStepControlFactory(
      Teuchos::abstractFactoryStd< StepControlStrategyBase<double>,
                                   FoolishStepControlStrategy >(),
      "Foolish Step Control");
  ib->setStepControlFactory(
      Teuchos::abstractFactoryStd< StepControlStrategyBase<double>,
                                   ImplicitBDFStepperStepControl<double> >(),
      "Other Implicit BDF Step Control");

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Implicit BDF");
  pl->sublist("Stepper Settings")
     .sublist("Step Control Settings")
     .sublist("Step Control Strategy Selection")
     .set("Step Control Strategy Type","Foolish Step Control");
  ib->setParameterList(pl);

  // Model:
  RCP<SinCosModel> model = sinCosModel(true);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();

  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  integrator = ib->create(model,ic,nlSolver);
  RCP<const StepperBase<double> > stepper = integrator->getStepper();
  TEST_ASSERT( !is_null(stepper) );
  RCP<const ImplicitBDFStepper<double> > ibdfStepper = 
    Teuchos::rcp_dynamic_cast<const ImplicitBDFStepper<double> >(stepper,false);
  TEST_ASSERT( !is_null(ibdfStepper) );
  RCP<const StepControlStrategyBase<double> > stepControl =
    ibdfStepper->getStepControlStrategy();
  TEST_ASSERT( !is_null(stepControl) );
  RCP<const FoolishStepControlStrategy> fStepControl =
   Teuchos::rcp_dynamic_cast<
      const FoolishStepControlStrategy>(stepControl,false);
  TEST_ASSERT( !is_null(fStepControl) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setStepControlFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW( 
      ib->setStepControlFactory(
        Teuchos::abstractFactoryStd< StepControlStrategyBase<double>,
                                     FoolishStepControlStrategy >(),
        "Implicit BDF Stepper Step Control Strategy"),
      std::logic_error);
#else // TEUCHOS_DEBUG
  TEST_NOTHROW( 
      ib->setStepControlFactory(
        Teuchos::abstractFactoryStd< StepControlStrategyBase<double>,
                                     FoolishStepControlStrategy >(),
        "Implicit BDF Stepper Step Control Strategy") );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setInterpolationBufferFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setInterpolationBufferFactory(
      Teuchos::abstractFactoryStd< InterpolationBufferBase<double>,
                                   FoolishInterpolationBuffer >(),
      "Foolish InterpolationBuffer");
  ib->setInterpolationBufferFactory(
      Teuchos::abstractFactoryStd< InterpolationBufferBase<double>,
                                   InterpolationBuffer<double> >(),
      "Other InterpolationBuffer");

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Explicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type","Explicit 4 Stage");
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .set("Interpolation Buffer Type","Foolish InterpolationBuffer");
  ib->setParameterList(pl);
  
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null 

  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  integrator = ib->create(model,ic,nlSolver);
  RCP<TrailingInterpolationBufferAcceptingIntegratorBase<double> >
    tibaIntegrator = Teuchos::rcp_dynamic_cast<
      TrailingInterpolationBufferAcceptingIntegratorBase<double> > (integrator,
                                                                    false);
  TEST_ASSERT( !is_null(tibaIntegrator) );
  RCP<const InterpolationBufferBase<double> > interpBuffer =
    tibaIntegrator->getTrailingInterpolationBuffer();
  TEST_ASSERT( !is_null(interpBuffer) );
  RCP<const FoolishInterpolationBuffer> fInterpBuffer = 
    Teuchos::rcp_dynamic_cast<
      const FoolishInterpolationBuffer>(interpBuffer,false);
  TEST_ASSERT( !is_null(fInterpBuffer) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder,
                   setInterpolationBufferFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW( 
      ib->setInterpolationBufferFactory(
        Teuchos::abstractFactoryStd< InterpolationBufferBase<double>,
                                     FoolishInterpolationBuffer >(),
        "Interpolation Buffer"),
      std::logic_error);
#else // TEUCHOS_DEBUG
  TEST_NOTHROW(
      ib->setInterpolationBufferFactory(
        Teuchos::abstractFactoryStd< InterpolationBufferBase<double>,
                                     FoolishInterpolationBuffer >(),
        "Interpolation Buffer") );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder,
                   setInterpolationBufferAppenderFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setInterpolationBufferAppenderFactory(
    Teuchos::abstractFactoryStd< InterpolationBufferAppenderBase<double>,
                                 FoolishInterpolationBufferAppender >(),
    "Foolish InterpolationBufferAppender");
  ib->setInterpolationBufferAppenderFactory(
    Teuchos::abstractFactoryStd<InterpolationBufferAppenderBase<double>,
                               PointwiseInterpolationBufferAppender<double> >(),
    "Other InterpolationBufferAppender");

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Explicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type","Explicit 4 Stage");
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Interpolation Buffer Appender Selection")
     .set("Interpolation Buffer Appender Type",
          "Foolish InterpolationBufferAppender");
  ib->setParameterList(pl);
  
  // Model:
  RCP<SinCosModel> model = sinCosModel(false);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null 

  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  integrator = ib->create(model,ic,nlSolver);
  RCP<InterpolationBufferAppenderAcceptingIntegratorBase<double> >
    ibaaIntegrator = Teuchos::rcp_dynamic_cast<
      InterpolationBufferAppenderAcceptingIntegratorBase<double> >(integrator,
                                                                   false);
  TEST_ASSERT( !is_null(ibaaIntegrator) );
  RCP<const InterpolationBufferAppenderBase<double> > interpBufferAppender =
    ibaaIntegrator->getInterpolationBufferAppender();
  TEST_ASSERT( !is_null(interpBufferAppender) );
  RCP<const FoolishInterpolationBufferAppender> fInterpBufferAppender = 
    Teuchos::rcp_dynamic_cast<
      const FoolishInterpolationBufferAppender>(interpBufferAppender,false);
  TEST_ASSERT( !is_null(fInterpBufferAppender) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder,
                   setInterpolationBufferAppenderFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW(
    ib->setInterpolationBufferAppenderFactory(
      Teuchos::abstractFactoryStd< InterpolationBufferAppenderBase<double>,
                                   FoolishInterpolationBufferAppender >(),
      "Pointwise Interpolation Buffer Appender"),
    std::logic_error);
#else // TEUCHOS_DEBUG
  TEST_NOTHROW(
    ib->setInterpolationBufferAppenderFactory(
      Teuchos::abstractFactoryStd< InterpolationBufferAppenderBase<double>,
                                   FoolishInterpolationBufferAppender >(),
      "Pointwise Interpolation Buffer Appender") );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setErrWtVecCalcFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setErrWtVecCalcFactory(
      Teuchos::abstractFactoryStd< ErrWtVecCalcBase<double>,
                                   FoolishErrWtVecCalc >(),
      "Foolish ErrWtVecCalc");
  ib->setErrWtVecCalcFactory(
      Teuchos::abstractFactoryStd< ErrWtVecCalcBase<double>,
                                   ImplicitBDFStepperErrWtVecCalc<double> >(),
      "Other ErrWtVecCalc");

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Implicit BDF");
  pl->sublist("Stepper Settings")
     .sublist("Step Control Settings")
     .sublist("Step Control Strategy Selection")
     .set("Step Control Strategy Type",
          "Implicit BDF Stepper Step Control Strategy");
  pl->sublist("Stepper Settings")
     .sublist("Step Control Settings")
     .sublist("Error Weight Vector Calculator Selection")
     .set("Error Weight Vector Calculator Type","Foolish ErrWtVecCalc");
  ib->setParameterList(pl);
  
  // Model:
  RCP<SinCosModel> model = sinCosModel(true);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();

  // Create Integrator
  RCP<IntegratorBase<double> > integrator;
  integrator = ib->create(model,ic,nlSolver);
  RCP<const StepperBase<double> > stepper = integrator->getStepper();
  RCP<const StepControlStrategyAcceptingStepperBase<double> > scsaStepper = 
    Teuchos::rcp_dynamic_cast<
      const StepControlStrategyAcceptingStepperBase<double> >(stepper,false);
  TEST_ASSERT( !is_null(scsaStepper) );
  RCP<const StepControlStrategyBase<double> > stepControl =
    scsaStepper->getStepControlStrategy();
  RCP<const ImplicitBDFStepperStepControl<double> > ibdfStepControl = 
    Teuchos::rcp_dynamic_cast<
      const ImplicitBDFStepperStepControl<double> >(stepControl,false);
  TEST_ASSERT( !is_null(ibdfStepControl) );
  RCP<const ErrWtVecCalcBase<double> > myErrWtVecCalc =
    ibdfStepControl->getErrWtVecCalc();
  RCP<const FoolishErrWtVecCalc> foolishErrWtVecCalc = 
    Teuchos::rcp_dynamic_cast<const FoolishErrWtVecCalc>(myErrWtVecCalc,false);
  TEST_ASSERT( !is_null(foolishErrWtVecCalc) );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setErrWtVecCalcFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW(
      ib->setErrWtVecCalcFactory(
        Teuchos::abstractFactoryStd< ErrWtVecCalcBase<double>,
                                     FoolishErrWtVecCalc >(),
        "Implicit BDF Stepper Error Weight Vector Calculator"),
      std::logic_error);
#else // TEUCHOS_DEBUG
  TEST_NOTHROW(
      ib->setErrWtVecCalcFactory(
        Teuchos::abstractFactoryStd< ErrWtVecCalcBase<double>,
                                     FoolishErrWtVecCalc >(),
        "Implicit BDF Stepper Error Weight Vector Calculator") );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setInterpolatorFactory ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  ib->setInterpolatorFactory(
      Teuchos::abstractFactoryStd< InterpolatorBase<double>,
                                   FoolishInterpolator >(),
      "Foolish Interpolator");
  ib->setInterpolatorFactory(
      Teuchos::abstractFactoryStd< InterpolatorBase<double>,
                                   LinearInterpolator<double> >(),
      "Other Interpolator");

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Integrator Settings")
     .sublist("Integrator Selection")
     .set("Integrator Type","Default Integrator");
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Backward Euler");
  pl->sublist("Stepper Settings")
     .sublist("Interpolator Selection")
     .set("Interpolator Type","Foolish Interpolator");
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .set("Interpolation Buffer Type","Interpolation Buffer");
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Interpolator Selection")
     .set("Interpolator Type","Cubic Spline Interpolator");
  ib->setParameterList(pl);
  
  // Model:
  RCP<SinCosModel> model = sinCosModel(true);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();

  // Create Integrator
  RCP<IntegratorBase<double> > integrator;

  // Create Stepper
  integrator = ib->create(model,ic,nlSolver);
  RCP<const StepperBase<double> > stepper = integrator->getStepper();
  RCP<const BackwardEulerStepper<double> > beStepper = 
    Teuchos::rcp_dynamic_cast<const BackwardEulerStepper<double> >(stepper,
                                                                   true);
  // First test that BackwardEuler Stepper got the Foolish Interpolator
  RCP<const InterpolatorBase<double> > interp = beStepper->getInterpolator();
  RCP<const FoolishInterpolator> fInterp =
    Teuchos::rcp_dynamic_cast<const FoolishInterpolator>(interp,false);
  TEST_ASSERT( !is_null(fInterp) );

  // Second test that InterpolationBuffer got the Cubic Spline Interpolator
  RCP<DefaultIntegrator<double> > dIntegrator = 
    Teuchos::rcp_dynamic_cast<DefaultIntegrator<double> >(integrator,true);
  RCP<const InterpolationBufferBase<double> > tInterp =
    dIntegrator->getTrailingInterpolationBuffer();
  RCP<const InterpolationBuffer<double> > interpBuffer = 
    Teuchos::rcp_dynamic_cast<const InterpolationBuffer<double> >(tInterp,true);
  RCP<const InterpolatorBase<double> > interpolator =
    interpBuffer->getInterpolator();
  RCP<const CubicSplineInterpolator<double> > cInterp = 
    Teuchos::rcp_dynamic_cast<
      const CubicSplineInterpolator<double> >(interpolator,false);
  TEST_ASSERT( !is_null(cInterp) );

}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setInterpolatorFactory_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
#ifdef TEUCHOS_DEBUG
  TEST_THROW(
      ib->setInterpolatorFactory(
        Teuchos::abstractFactoryStd< InterpolatorBase<double>,
                                     FoolishInterpolator >(),
        "Linear Interpolator"),
      std::logic_error);
#else // TEUCHOS_DEBUG
  TEST_NOTHROW(
      ib->setInterpolatorFactory(
        Teuchos::abstractFactoryStd< InterpolatorBase<double>,
                                     FoolishInterpolator >(),
        "Linear Interpolator") );
  TEST_THROW( ib->getValidParameters(), std::logic_error );
#endif // TEUCHOS_DEBUG
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setWFactoryObject ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  // First, we verify that if we don't set one, we don't get one.  :-)
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Integrator Settings")
     .sublist("Integrator Selection")
     .set("Integrator Type","Default Integrator");
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Implicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type",
          "Implicit 2 Stage 4th order Gauss");
  ib->setParameterList(pl);
  
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();
  {
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    RCP<const StepperBase<double> > stepper = integrator->getStepper();
    RCP<const ImplicitRKStepper<double> > irkStepper = 
      Teuchos::rcp_dynamic_cast<const ImplicitRKStepper<double> >(stepper,true);
    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > wFactory =
      irkStepper->get_W_factory();
    TEST_ASSERT( is_null(wFactory) );
  }

  // Second, we verify that if we set one, we get one.
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > wF =
    getWFactory<double>(Teuchos::parameterList());
  ib->setWFactoryObject(wF);
  {
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    RCP<const StepperBase<double> > stepper = integrator->getStepper();
    RCP<const ImplicitRKStepper<double> > irkStepper = 
      Teuchos::rcp_dynamic_cast<const ImplicitRKStepper<double> >(stepper,true);
    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > wFactory =
      irkStepper->get_W_factory();
    TEST_ASSERT( !is_null(wFactory) );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setWFactoryObject_bad ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > wFactory; // null
  TEST_THROW( ib->setWFactoryObject(wFactory), std::logic_error );
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setRKButcherTableauBuilder ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<RKButcherTableauBuilder<double> > rkbtB =
    rKButcherTableauBuilder<double>();

  rkbtB->setRKButcherTableauFactory(
      Teuchos::abstractFactoryStd< RKButcherTableauBase<double>,
                                   Implicit1Stage2ndOrderGauss_RKBT<double> >(),
      "Another One Stage RKBT");
  rkbtB->setRKButcherTableauFactory(
      Teuchos::abstractFactoryStd< RKButcherTableauBase<double>,
                                   Implicit1Stage2ndOrderGauss_RKBT<double> >(),
      "Yet Another One Stage RKBT");
  ib->setRKButcherTableauBuilder(rkbtB);

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*ib->getValidParameters());
  pl->sublist("Integrator Settings")
     .sublist("Integrator Selection")
     .set("Integrator Type","Default Integrator");
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Implicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type","Another One Stage RKBT");
  ib->setParameterList(pl);
  
  // Model:
  RCP<SinCosModel> model = sinCosModel(true);
  // IC:
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  // Nonlinear Solver:
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();

  // Create Integrator
  RCP<IntegratorBase<double> > integrator;

  // Create Stepper
  integrator = ib->create(model,ic,nlSolver);
  RCP<const StepperBase<double> > stepper = integrator->getStepper();
  RCP<const ImplicitRKStepper<double> > irkStepper = 
    Teuchos::rcp_dynamic_cast<const ImplicitRKStepper<double> >(stepper,true);
  RCP<const RKButcherTableauBase<double> > rkbt =
    irkStepper->getRKButcherTableau();
  RCP<const Implicit1Stage2ndOrderGauss_RKBT<double> > BE_rkbt = 
    Teuchos::rcp_dynamic_cast<
      const Implicit1Stage2ndOrderGauss_RKBT<double> >(rkbt,false);
  TEST_ASSERT( !is_null(BE_rkbt) );
}

/*
TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, explicitRKwithImplicitRKBT ) {
  // Should throw
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, implicitRKwithExplicitRKBT ) {
  // What happens in this case?
}
*/


TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_ExplicitRK ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(false);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Explicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type","Explicit 4 Stage");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  TEST_ASSERT( !is_null(integrator) );
  {
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_ASSERT( !is_null(stepper_out) );
    RCP<const ExplicitRKStepper<double> > erkStepper =
      Teuchos::rcp_dynamic_cast<const ExplicitRKStepper<double> >(stepper_out,
                                                                  false);
    TEST_ASSERT( !is_null(erkStepper) );
    RCP<const Thyra::ModelEvaluator<double> > model_out =
      stepper_out->getModel();
    TEST_ASSERT( !is_null(model_out) );
    RCP<const SinCosModel> sinCosModel =
      Teuchos::rcp_dynamic_cast<const SinCosModel>(model_out,false);
    TEST_ASSERT( !is_null(sinCosModel) );
    RCP<const RKButcherTableauBase<double> > rkbt =
      erkStepper->getRKButcherTableau();
    TEST_ASSERT( !is_null(rkbt) );
    RCP<const Explicit4Stage4thOrder_RKBT<double> > rkbt4 =
      Teuchos::rcp_dynamic_cast<
        const Explicit4Stage4thOrder_RKBT<double> >(rkbt,false);
    TEST_ASSERT( !is_null(rkbt4) );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_ForwardEuler ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(false);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Forward Euler");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  TEST_ASSERT( !is_null(integrator) );
  {
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_ASSERT( !is_null(stepper_out) );
    RCP<const ForwardEulerStepper<double> > feStepper =
      Teuchos::rcp_dynamic_cast<const ForwardEulerStepper<double> >(stepper_out,
                                                                    false);
    TEST_ASSERT( !is_null(feStepper) );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_ImplicitRK ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Implicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type","Backward Euler");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  TEST_ASSERT( !is_null(integrator) );
  {
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_ASSERT( !is_null(stepper_out) );
    RCP<const ImplicitRKStepper<double> > irkStepper =
      Teuchos::rcp_dynamic_cast<const ImplicitRKStepper<double> >(stepper_out,
                                                                  false);
    TEST_ASSERT( !is_null(irkStepper) );
    RCP<const Thyra::ModelEvaluator<double> > model_out =
      stepper_out->getModel();
    TEST_ASSERT( !is_null(model_out) );
    RCP<const SinCosModel> sinCosModel =
      Teuchos::rcp_dynamic_cast<const SinCosModel>(model_out,false);
    TEST_ASSERT( !is_null(sinCosModel) );
    RCP<const Thyra::NonlinearSolverBase<double> > nlSolver_out =
      irkStepper->getSolver();
    TEST_ASSERT( !is_null(nlSolver_out) );
    RCP<const TimeStepNonlinearSolver<double> > tsnlSolver =
      Teuchos::rcp_dynamic_cast<
        const TimeStepNonlinearSolver<double> >(nlSolver_out, false);
    TEST_ASSERT( !is_null(tsnlSolver) );
    RCP<const RKButcherTableauBase<double> > rkbt =
      irkStepper->getRKButcherTableau();
    TEST_ASSERT( !is_null(rkbt) );
    RCP<const BackwardEuler_RKBT<double> > rkbtBE =
      Teuchos::rcp_dynamic_cast<const BackwardEuler_RKBT<double> >(rkbt,false);
    TEST_ASSERT( !is_null(rkbtBE) );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_ImplicitBDF ) {
  {
    // Base test that all the pointers are of the right type
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    pl->sublist("Integration Control Strategy Selection")
       .set("Integration Control Strategy Type",
            "Simple Integration Control Strategy");
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
       .set("Interpolation Buffer Type","Interpolation Buffer");
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Interpolator Selection")
       .set("Interpolator Type","Cubic Spline Interpolator");
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","Implicit BDF");
    pl->sublist("Stepper Settings")
       .sublist("Step Control Settings")
       .sublist("Step Control Strategy Selection")
       .set("Step Control Strategy Type",
            "Implicit BDF Stepper Step Control Strategy");
    pl->sublist("Stepper Settings")
       .sublist("Step Control Settings")
       .sublist("Error Weight Vector Calculator Selection")
       .set("Error Weight Vector Calculator Type",
            "Implicit BDF Stepper Error Weight Vector Calculator");
    ib->setParameterList(pl);
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    TEST_ASSERT( !is_null(integrator) );
    
    // Stepper = Implicit BDF
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_ASSERT( !is_null(stepper_out) );
    RCP<const ImplicitBDFStepper<double> > ibdfStepper =
      Teuchos::rcp_dynamic_cast<const ImplicitBDFStepper<double> >(stepper_out,
                                                                   false);
    TEST_ASSERT( !is_null(ibdfStepper) );
    // Model = SinCosModel
    RCP<const Thyra::ModelEvaluator<double> > model_out = 
      stepper_out->getModel();
    TEST_ASSERT( !is_null(model_out) );
    RCP<const SinCosModel> sinCosModel =
      Teuchos::rcp_dynamic_cast<const SinCosModel>(model_out,false);
    TEST_ASSERT( !is_null(sinCosModel) );
    // Solver = TimeStepNonlinearSolver
    RCP<const Thyra::NonlinearSolverBase<double> > nlSolver_out =
      ibdfStepper->getSolver();
    TEST_ASSERT( !is_null(nlSolver_out) );
    RCP<const TimeStepNonlinearSolver<double> > tsnlSolver =
      Teuchos::rcp_dynamic_cast<
        const TimeStepNonlinearSolver<double> >(nlSolver_out, false);
    TEST_ASSERT( !is_null(tsnlSolver) );
    // Step Control = ImplicitBDFStepperStepControl
    RCP<const StepControlStrategyBase<double> > stepControl_out =
      ibdfStepper->getStepControlStrategy();
    TEST_ASSERT( !is_null(stepControl_out) );
    RCP<const ImplicitBDFStepperStepControl<double> > ibdfStepControl =
      Teuchos::rcp_dynamic_cast<
        const ImplicitBDFStepperStepControl<double> >(stepControl_out, false);
    TEST_ASSERT( !is_null(ibdfStepControl) );
    // ErrWtVecCalc = ImplicitBDFStepperErrWtVecCalc
    RCP<const ErrWtVecCalcBase<double> > myErrWtVecCalc =
      ibdfStepControl->getErrWtVecCalc();
    RCP<const ImplicitBDFStepperErrWtVecCalc<double> > myIBDFErrWtVecCalc = 
      Teuchos::rcp_dynamic_cast<
        const ImplicitBDFStepperErrWtVecCalc<double> >(myErrWtVecCalc);
    TEST_ASSERT( !is_null(myIBDFErrWtVecCalc) );
    // Integrator = DefaultIntegrator
    RCP<DefaultIntegrator<double> > defInt =
      Teuchos::rcp_dynamic_cast<DefaultIntegrator<double> >(integrator,false);
    TEST_ASSERT( !is_null(defInt) );
    // Integration Control = SimpleIntegrationControlStrategy
    RCP<const IntegrationControlStrategyBase<double> > integrationControl_out =
      defInt->getIntegrationControlStrategy();
    TEST_ASSERT( !is_null(integrationControl_out) );
    RCP<const SimpleIntegrationControlStrategy<double> > simpleIControl =
      Teuchos::rcp_dynamic_cast<
        const SimpleIntegrationControlStrategy<double> >(integrationControl_out,
                                                         false);
    TEST_ASSERT( !is_null(simpleIControl) );
    // TrailingInterpolationBuffer = InterpolationBuffer
    RCP<const InterpolationBufferBase<double> > interpolationBuffer =
      defInt->getTrailingInterpolationBuffer();
    TEST_ASSERT( !is_null(interpolationBuffer) );
    RCP<const InterpolationBuffer<double> > defaultInterpolationBuffer =
      Teuchos::rcp_dynamic_cast<
        const InterpolationBuffer<double> >(interpolationBuffer, false);
    TEST_ASSERT( !is_null(defaultInterpolationBuffer) );
    // Interpolator = CubicSplineInterpolator
    RCP<const InterpolatorBase<double> > interp =
      defaultInterpolationBuffer->getInterpolator();
    TEST_ASSERT( !is_null(interp) );
    RCP<const CubicSplineInterpolator<double> > cInterp = 
      Teuchos::rcp_dynamic_cast<const CubicSplineInterpolator<double> >(interp,
                                                                        false);
    TEST_ASSERT( !is_null(cInterp) );
  }
  {
    // Test that we can set Integration Control Strategy to None.
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","Implicit BDF");
    pl->sublist("Integration Control Strategy Selection")
       .set("Integration Control Strategy Type","None");
    ib->setParameterList(pl);
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    TEST_ASSERT( !is_null(integrator) );
    // Integrator = DefaultIntegrator
    RCP<DefaultIntegrator<double> > defInt =
      Teuchos::rcp_dynamic_cast<DefaultIntegrator<double> >(integrator,false);
    TEST_ASSERT( !is_null(defInt) );
    // Integration Control = None
    RCP<const IntegrationControlStrategyBase<double> > integrationControl_out =
      defInt->getIntegrationControlStrategy();
    TEST_ASSERT( is_null(integrationControl_out) );
  }
  {
    // Test that we can set Step Control Strategy to None
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","Implicit BDF");
    pl->sublist("Stepper Settings")
       .sublist("Step Control Settings")
       .sublist("Step Control Strategy Selection")
       .set("Step Control Strategy Type","None");
    ib->setParameterList(pl);
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    TEST_ASSERT( !is_null(integrator) );
    // Stepper = Implicit BDF
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_ASSERT( !is_null(stepper_out) );
    RCP<const ImplicitBDFStepper<double> > ibdfStepper =
      Teuchos::rcp_dynamic_cast<const ImplicitBDFStepper<double> >(stepper_out,
                                                                   false);
    TEST_ASSERT( !is_null(ibdfStepper) );
    // Step Control = None
    RCP<const StepControlStrategyBase<double> > stepControl_out =
      ibdfStepper->getStepControlStrategy();
    TEST_ASSERT( is_null(stepControl_out) );
  }
  {
    // Test that we can set TrailingInterpolationBuffer to None
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","Implicit BDF");
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
       .set("Interpolation Buffer Type","None");
    ib->setParameterList(pl);
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    TEST_ASSERT( !is_null(integrator) );
    // Trailing Interpolation Buffer = null
    RCP<TrailingInterpolationBufferAcceptingIntegratorBase<double> >
      tibaIntegrator = Teuchos::rcp_dynamic_cast<
        TrailingInterpolationBufferAcceptingIntegratorBase<double> >(integrator,
                                                                     false);
    TEST_ASSERT( !is_null(tibaIntegrator) );
    RCP<const InterpolationBufferBase<double> > trailingIB =
      tibaIntegrator->getTrailingInterpolationBuffer();
    TEST_ASSERT( is_null(trailingIB) );
  }
  {
    // Test that we can set the Interpolator to None
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->setParameters(*(ib->getValidParameters()));
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","Implicit BDF");
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Trailing Interpolation Buffer Selection")
       .set("Interpolation Buffer Type","Interpolation Buffer");
    pl->sublist("Interpolation Buffer Settings")
       .sublist("Interpolator Selection")
       .set("Interpolator Type","None");
    ib->setParameterList(pl);
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();
    RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
    TEST_ASSERT( !is_null(integrator) );
    // Trailing Interpolation Buffer = null
    RCP<TrailingInterpolationBufferAcceptingIntegratorBase<double> >
      tibaIntegrator = Teuchos::rcp_dynamic_cast<
        TrailingInterpolationBufferAcceptingIntegratorBase<double> >(integrator,
                                                                     false);
    TEST_ASSERT( !is_null(tibaIntegrator) );
    RCP<const InterpolationBufferBase<double> > trailingIB =
      tibaIntegrator->getTrailingInterpolationBuffer();
    TEST_ASSERT( !is_null(trailingIB) );
    RCP<const InterpolationBuffer<double> > defaultIB = 
      Teuchos::rcp_dynamic_cast<const InterpolationBuffer<double> >(trailingIB,
                                                                    false);
    TEST_ASSERT( !is_null(defaultIB) );
    RCP<const InterpolatorBase<double> > interp = defaultIB->getInterpolator();
    TEST_ASSERT( !is_null(interp) ); // Default interpolator is created
                                     // in this case
    RCP<const LinearInterpolator<double> > lInterp = 
     Teuchos::rcp_dynamic_cast<const LinearInterpolator<double> >(interp,false);
    TEST_ASSERT( !is_null(lInterp) );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_BackwardEuler ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Backward Euler");
  pl->sublist("Stepper Settings")
     .sublist("Interpolator Selection")
     .set("Interpolator Type","Cubic Spline Interpolator");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  TEST_ASSERT( !is_null(integrator) );
  {
    RCP<const StepperBase<double> > stepper_out = integrator->getStepper();
    TEST_ASSERT( !is_null(stepper_out) );
    RCP<const BackwardEulerStepper<double> > beStepper =
      Teuchos::rcp_dynamic_cast<
        const BackwardEulerStepper<double> >(stepper_out,false);
    TEST_ASSERT( !is_null(beStepper) );
    RCP<const Thyra::ModelEvaluator<double> > model_out =
      stepper_out->getModel();
    TEST_ASSERT( !is_null(model_out) );
    RCP<const SinCosModel> sinCosModel =
      Teuchos::rcp_dynamic_cast<const SinCosModel>(model_out,false);
    TEST_ASSERT( !is_null(sinCosModel) );
    RCP<const Thyra::NonlinearSolverBase<double> > nlSolver_out =
      beStepper->getSolver();
    TEST_ASSERT( !is_null(nlSolver_out) );
    RCP<const TimeStepNonlinearSolver<double> > tsnlSolver =
      Teuchos::rcp_dynamic_cast<
        const TimeStepNonlinearSolver<double> >(nlSolver_out, false);
    TEST_ASSERT( !is_null(tsnlSolver) );
    RCP<const InterpolatorBase<double> > interp = beStepper->getInterpolator();
    TEST_ASSERT( !is_null(interp) );
    RCP<const CubicSplineInterpolator<double> > cInterp =
      Teuchos::rcp_dynamic_cast<const CubicSplineInterpolator<double> >(interp,
                                                                        false);
    TEST_ASSERT( !is_null(cInterp) );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, fullyInitialized_ERK ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(false);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Explicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type","Explicit 4 Stage");
  pl->sublist("Integration Control Strategy Selection")
     .set("Integration Control Strategy Type",
          "Simple Integration Control Strategy");
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Take Variable Steps",false);
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Fixed dt",0.1);
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .set("Interpolation Buffer Type","Interpolation Buffer");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  Teuchos::Array<double> time_vec;
  time_vec.push_back(pl->sublist("Integrator Settings")
                        .get<double>("Final Time"));
  integrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  TEST_ASSERT( true ); 
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, fullyInitialized_FE ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(false);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Forward Euler");
  pl->sublist("Integration Control Strategy Selection")
     .set("Integration Control Strategy Type",
          "Simple Integration Control Strategy");
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Take Variable Steps",false);
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Fixed dt",0.1);
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .set("Interpolation Buffer Type","Interpolation Buffer");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver; // null
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  Teuchos::Array<double> time_vec;
  time_vec.push_back(pl->sublist("Integrator Settings")
                        .get<double>("Final Time"));
  integrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  TEST_ASSERT( true ); 
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, fullyInitialized_BE ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Backward Euler");
  pl->sublist("Integration Control Strategy Selection")
     .set("Integration Control Strategy Type",
          "Simple Integration Control Strategy");
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Take Variable Steps",false);
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Fixed dt",0.1);
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .set("Interpolation Buffer Type","Interpolation Buffer");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  Teuchos::Array<double> time_vec;
  time_vec.push_back(pl->sublist("Integrator Settings")
                        .get<double>("Final Time"));
  integrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  TEST_ASSERT( true ); 
}

#ifdef Rythmos_ENABLE_NOX

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, fullyInitialized_BE_NOX ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Backward Euler");
  pl->sublist("Integration Control Strategy Selection")
     .set("Integration Control Strategy Type",
          "Simple Integration Control Strategy");
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Take Variable Steps",false);
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Fixed dt",0.1);
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .set("Interpolation Buffer Type","Interpolation Buffer");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    rcp(new Thyra::NOXNonlinearSolver);
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  Teuchos::Array<double> time_vec;
  time_vec.push_back(pl->sublist("Integrator Settings")
                        .get<double>("Final Time"));
  integrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  TEST_ASSERT( true ); 
}

#endif // Rythmos_ENABLE_NOX

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, fullyInitialized_DIRK ) {
  // DIRK/SDIRK w/o WFactory
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  modelPL->sublist("Stratimikos").set("Linear Solver Type","AztecOO");
  modelPL->sublist("Stratimikos").set("Preconditioner Type","None");
  modelPL->sublist("DiagonalTransientModel").set("NumElements",2);
  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(modelPL);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Implicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .sublist("Implicit RK")
     .sublist("VerboseObject")
     .set("Verbosity Level","none");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type",
          "Singly Diagonal IRK 2 Stage 3rd order");
  pl->sublist("Integration Control Strategy Selection")
     .set("Integration Control Strategy Type",
          "Simple Integration Control Strategy");
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Take Variable Steps",false);
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Fixed dt",0.1);
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .set("Interpolation Buffer Type","Interpolation Buffer");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  Teuchos::Array<double> time_vec;
  time_vec.push_back(pl->sublist("Integrator Settings")
                        .get<double>("Final Time"));
  integrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  TEST_ASSERT( true ); 
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, fullyInitialized_IRK ) {
  // Dense RKBT w/ WFactory
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> modelPL = Teuchos::parameterList();
  modelPL->sublist("Stratimikos").set("Linear Solver Type","AztecOO");
  modelPL->sublist("Stratimikos").set("Preconditioner Type","None");
  modelPL->sublist("DiagonalTransientModel").set("NumElements",2);
  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(modelPL);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Implicit RK");
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .sublist("Implicit RK")
     .sublist("VerboseObject")
     .set("Verbosity Level","none");
  pl->sublist("Stepper Settings")
     .sublist("Runge Kutta Butcher Tableau Selection")
     .set("Runge Kutta Butcher Tableau Type",
          "Implicit 3 Stage 6th order Gauss");
  pl->sublist("Integration Control Strategy Selection")
     .set("Integration Control Strategy Type",
          "Simple Integration Control Strategy");
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Take Variable Steps",false);
  pl->sublist("Integration Control Strategy Selection")
     .sublist("Simple Integration Control Strategy")
     .set("Fixed dt",0.1);
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .set("Interpolation Buffer Type","Interpolation Buffer");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();
  ib->setWFactoryObject(getWFactory<double>(modelPL));
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  Teuchos::Array<double> time_vec;
  time_vec.push_back(pl->sublist("Integrator Settings")
                        .get<double>("Final Time"));
  integrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  TEST_ASSERT( true ); 
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, fullyInitialized_IBDF_minimal ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Implicit BDF");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  Teuchos::Array<double> time_vec;
  time_vec.push_back(pl->sublist("Integrator Settings")
                        .get<double>("Final Time"));
  integrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  TEST_ASSERT( true ); 
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, fullyInitialized_IBDF_all ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<SinCosModel> model = sinCosModel(true);
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setParameters(*(ib->getValidParameters()));
  pl->sublist("Integrator Settings")
     .sublist("Integrator Selection")
     .set("Integrator Type","Default Integrator");
  pl->sublist("Stepper Settings")
     .sublist("Stepper Selection")
     .set("Stepper Type","Implicit BDF");
  pl->sublist("Stepper Settings")
     .sublist("Step Control Settings")
     .sublist("Error Weight Vector Calculator Selection")
     .set("Error Weight Vector Calculator Type",
          "Implicit BDF Stepper Error Weight Vector Calculator");
  pl->sublist("Stepper Settings")
     .sublist("Step Control Settings")
     .sublist("Step Control Strategy Selection")
     .set("Step Control Strategy Type",
          "Implicit BDF Stepper Step Control Strategy");
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Trailing Interpolation Buffer Selection")
     .set("Interpolation Buffer Type","Interpolation Buffer");
  pl->sublist("Interpolation Buffer Settings")
     .sublist("Interpolator Selection")
     .set("Interpolator Type","Hermite Interpolator");
  //pl->sublist("Interpolation Buffer Settings")
  //   .sublist("Interpolator Selection")
  //   .set("Interpolator Type","Linear Interpolator");
  //pl->sublist("Interpolation Buffer Settings")
  //   .sublist("Interpolator Selection")
  //   .set("Interpolator Type","Cubic Spline Interpolator");
  pl->sublist("Integration Control Strategy Selection")
     .set("Integration Control Strategy Type",
          "Simple Integration Control Strategy");
  ib->setParameterList(pl);
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
    timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > integrator = ib->create(model,ic,nlSolver);
  Teuchos::Array<double> time_vec;
  time_vec.push_back(pl->sublist("Integrator Settings")
                        .get<double>("Final Time"));
  integrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  TEST_ASSERT( true ); 
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, create_invalid ) {
  {
    // Integrator == null
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Integrator Settings")
       .sublist("Integrator Selection")
       .set("Integrator Type","None");
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","Explicit RK");
    pl->sublist("Stepper Settings")
       .sublist("Runge Kutta Butcher Tableau Selection")
       .set("Runge Kutta Butcher Tableau Type","Explicit 4 Stage");
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model = sinCosModel(false);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
  {
    // Stepper == null
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","None");
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model = sinCosModel(false);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
  {
    // Model == null
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    RCP<ParameterList> pl = Teuchos::parameterList();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model; 
    Thyra::ModelEvaluatorBase::InArgs<double> ic;
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
  {
    // Implicit Stepper and Solver == null
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","Implicit BDF");
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
  {
    // RK Butcher Tableau Accepting Stepper and RKBT == null  (ERK)
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","Explicit RK");
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model = sinCosModel(false);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver;
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
  {
    // RK Butcher Tableau Accepting Stepper and RKBT == null (IRK)
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Stepper Settings")
       .sublist("Stepper Selection")
       .set("Stepper Type","Implicit RK");
    RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
    ib->setParameterList(pl);
    RCP<IntegratorBase<double> > integrator;
    RCP<SinCosModel> model = sinCosModel(true);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
    RCP<Thyra::NonlinearSolverBase<double> > nlSolver =
      timeStepNonlinearSolver<double>();
    TEST_THROW( integrator = ib->create(model,ic,nlSolver), std::logic_error );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, setParameterList ) {
  // does it validate the list?
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Hello","World");
  TEST_THROW(ib->setParameterList(pl), std::logic_error);
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, getValidParameters ) {
  // Create a hand-version of the Valid Parameter List 
  // "Foo Settings" => sublist of various options that IntegratorBuilder owns
  // "Foo Selection" => sublist of various options that ObjectBuilder owns
  RCP<ParameterList> validPL = Teuchos::parameterList();
  {
    ParameterList& integratorSettingsPL =
      validPL->sublist("Integrator Settings");
    {
      integratorSettingsPL.set("Final Time",1.0);
      integratorSettingsPL.set("Land On Final Time", true);
      integratorSettingsPL.sublist("Integrator Selection")
                          .disableRecursiveValidation();
    }
    validPL->sublist("Integration Control Strategy Selection")
                     .disableRecursiveValidation();
    ParameterList& stepperSettingsPL = validPL->sublist("Stepper Settings");
    {
      stepperSettingsPL.sublist("Stepper Selection")
                       .disableRecursiveValidation();
      ParameterList& stepControlSettingsPL =
        stepperSettingsPL.sublist("Step Control Settings");
      {
        stepControlSettingsPL.sublist("Step Control Strategy Selection")
                             .disableRecursiveValidation();
        stepControlSettingsPL
                           .sublist("Error Weight Vector Calculator Selection")
                           .disableRecursiveValidation();
      }
      stepperSettingsPL.sublist("Interpolator Selection")
                       .disableRecursiveValidation();
      //stepperSettingsPL.sublist("Nonlinear Solver Selection")
      //                 .disableRecursiveValidation();
      //stepperSettingsPL
      //              .sublist("Runge-Kutta Stepper Butcher Tableau Selection")
      //              .disableRecursiveValidation();
    }
    ParameterList& ibSettingsPL =
      validPL->sublist("Interpolation Buffer Settings");
    {
      ibSettingsPL.sublist("Trailing Interpolation Buffer Selection")
                  .disableRecursiveValidation();
      ibSettingsPL.sublist("Interpolation Buffer Appender Selection")
                  .disableRecursiveValidation();
      ibSettingsPL.sublist("Interpolator Selection")
                  .disableRecursiveValidation();
    }
    //ParameterList& observerSettingsPL =
    //  validPL->sublist("Integration Observer Settings");
    //{
    //  observerSettingsPL.sublist("Integration Observer Enabled")
    //                    .disableRecursiveValidation();
    //  observerSettingsPL.sublist("Integration Observer List")
    //                    .disableRecursiveValidation();
    //}
  }

  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<const ParameterList> pl = ib->getValidParameters();
  TEST_ASSERT( !is_null(pl) );

  //TEST_NOTHROW( validPL->validateParameters(*pl) );  
  validPL->validateParameters(*pl);  
}

TEUCHOS_UNIT_TEST( Rythmos_IntegratorBuilder, printParams ) {
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<const ParameterList> pl = ib->getValidParameters();
  std::cout << "Valid Parameter List for IntegratorBase:" << std::endl;
  std::ofstream fout("Rythmos_ParameterList.txt");
  pl->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
                                                       .indent(4));
  fout.close();
  TEST_ASSERT( true );
}

} // namespace Rythmos 



