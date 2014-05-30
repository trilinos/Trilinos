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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef Rythmos_INTEGRATOR_BUILDER_DEF_H
#define Rythmos_INTEGRATOR_BUILDER_DEF_H

// Rythmos classes:
#include "Rythmos_IntegratorBuilder_decl.hpp"
#include "Rythmos_IntegrationControlStrategyAcceptingIntegratorBase.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"

// Teuchos:
#include "Teuchos_as.hpp"

// Specific objects to seed the builder:
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_RampingIntegrationControlStrategy.hpp"
#include "Rythmos_FixedStepControlStrategy.hpp"
#include "Rythmos_SimpleStepControlStrategy.hpp"
#include "Rythmos_FirstOrderErrorStepControlStrategy.hpp"
#include "Rythmos_ImplicitBDFStepperStepControl.hpp"
#include "Rythmos_ImplicitBDFStepperRampingStepControl.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_PointwiseInterpolationBufferAppender.hpp"
//#include "Rythmos_SmartInterpolationBufferAppender.hpp"
#include "Rythmos_ImplicitBDFStepperErrWtVecCalc.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_HermiteInterpolator.hpp"
#include "Rythmos_CubicSplineInterpolator.hpp"

// Includes for the Forward Sensitivity Integrator Builder:
#include "Rythmos_ForwardSensitivityStepper.hpp"


namespace {

  // Valid ParameterList names:
  static std::string integratorSettings_name = "Integrator Settings";
  static std::string integratorSettings_docs =
    "These parameters are used directly in setting up the Integrator.";
  static std::string integratorSelection_name = "Integrator Selection";
  static std::string integratorSelection_docs =
    "Select the Integrator to be used.";
  static std::string integrationControlSelection_name =
    "Integration Control Strategy Selection";
  static std::string integrationControlSelection_docs =
    "Note that some settings conflict between step control and integration "
    "control.  In general, the integration control decides which steps will "
    "be fixed or variable, not the stepper.  When the integration control "
    "decides to take variable steps, the step control is then responsible "
    "for choosing appropriate step-sizes.";
  static std::string stepperSettings_name = "Stepper Settings";
  static std::string stepperSettings_docs =
    "This parameter list sets various parameters for the Stepper.";
  static std::string stepperSelection_name = "Stepper Selection";
  static std::string stepperSelection_docs =
    "Selects the Stepper for the time integration.  It should be that "
    "some time integrators can be accessed through different Steppers, "
    "e.g., Backward Euler can be obtained through the `Backward Euler', "
    "a first-order `Implicit BDF', or a one-stage `Implicit RK' Stepper."
    "Special note for `Implicit RK' Stepper:  If a fully implicit RK Butcher "
    "tableau is chosen, then the stepper will not be fully initialized "
    "unless a W factory object is set on the IntegratorBuilder through "
    "setWFactoryObject.";
  static std::string ForwardEulerStepper_name = "Forward Euler";
  static std::string ForwardEulerStepper_docs =
    "This is the basic Forward Euler method: x_n = x_{n-1} + dt*x_dot_{n-1}";
  static std::string BackwardEulerStepper_name = "Backward Euler";
  static std::string BackwardEulerStepper_docs =
    "This is the basic Backward Euler method: x_n = x_{n-1} + dt*x_dot_n";
  static std::string ImplicitBDFStepper_name = "Implicit BDF";
  static std::string ImplicitBDFStepper_docs =
    "This Stepper provides a re-implementation of the algorithms in "
    "the LLNL Sundials code IDA. This is an implicit BDF integrator "
    "for DAEs which uses variable step-sizes and variable-orders "
    "first through fourth.";
  static std::string rkButcherTableauSelection_name =
    "Runge Kutta Butcher Tableau Selection";
  static std::string rkButcherTableauSelection_docs =
    "Only the Explicit RK Stepper and the Implicit RK Stepper accept an "
    "RK Butcher Tableau.";
  static std::string ExplicitRKStepper_name = "Explicit RK";
  static std::string ExplicitRKStepper_docs =
    "This Stepper has many explicit time-integrators using Runge-Kutta "
    "formulation and the Butcher Tableau specification.  See `"
    +rkButcherTableauSelection_name+"' ParameterList for available options.";
  static std::string ImplicitRKStepper_name = "Implicit RK";
  static std::string ImplicitRKStepper_docs =
    "This Stepper has many implicit time-integrators using Runge-Kutta "
    "formulation and the Butcher Tableau specification.  See `"
    +rkButcherTableauSelection_name+"' ParameterList for available options.";
  static std::string stepControlSettings_name = "Step Control Settings";
  static std::string stepControlSettings_docs =
    "Not all step control strategies are compatible with each stepper.  "
    "If the strategy has the name of a stepper in its name, then it only "
    "works with that stepper.";
  static std::string stepControlSelection_name =
    "Step Control Strategy Selection";
  static std::string stepControlSelection_docs =
    "Used to select the Control Strategy for the stepper.";
  static std::string errWtVecSelection_name =
    "Error Weight Vector Calculator Selection";
  static std::string errWtVecSelection_docs =
    "Not all ErrWtVec calculators are compatible with each step control "
    "strategy.  If the calculator has the name of a stepper or another "
    "step control strategy in its name, then it only works with that step "
    "control strategy.";
  static std::string interpolationBufferSettings_name =
    "Interpolation Buffer Settings";
  static std::string interpolationBufferSettings_docs =
    "This parameter list sets various parameters for the InterpolationBuffer.";
  static std::string interpolationBufferSelection_name =
    "Trailing Interpolation Buffer Selection";
  static std::string interpolationBufferSelection_docs =
    "Used to select the Interpolation Buffer.";
  static std::string interpolationBufferAppenderSelection_name =
   "Interpolation Buffer Appender Selection";
  static std::string interpolationBufferAppenderSelection_docs =
   "Used to select the Interpolation Buffer Appender.";
  static std::string initialTime_name = "Initial Time";
  static int initialTime_default = 0; // Should be Scalar(0.0)
  static std::string initialTime_docs =
    "The initial time to start integration.";
  static std::string finalTime_name = "Final Time";
  static int finalTime_default = 1; // Should be Scalar(1.0)
  static std::string finalTime_docs = "The final time to end integration.";
  static std::string landOnFinalTime_name = "Land On Final Time";
  static bool landOnFinalTime_default = true;
  static std::string landOnFinalTime_docs =
    "Exactly land on the final time; do not step past final time and "
    "interpolate.";
  static std::string interpolatorSelection_name = "Interpolator Selection";
  static std::string interpolatorSelection_docs =
    "Choose the interpolator to use.";
  static std::string stepperInterpolatorSelection_docs =
    "Note all Steppers accept an interpolator.  Currently, only the "
    "BackwardEuler stepper does.";

  // Builder names:
  static std::string integratorBuilder_name = "Rythmos::Integrator";
  static std::string integratorBuilderType_name = "Integrator Type";
  static std::string integrationControlBuilder_name =
    "Rythmos::IntegrationControlStrategy";
  static std::string integrationControlBuilderType_name =
    "Integration Control Strategy Type";
  static std::string stepControlBuilder_name =
    "Rythmos::StepControlStrategy";
  static std::string stepControlBuilderType_name =
    "Step Control Strategy Type";
  static std::string interpolationBufferBuilder_name =
    "Rythmos::InterpolationBuffer";
  static std::string interpolationBufferBuilderType_name =
    "Interpolation Buffer Type";
  static std::string interpolationBufferAppenderBuilder_name =
    "Rythmos::InterpolationBufferAppender";
  static std::string interpolationBufferAppenderBuilderType_name =
    "Interpolation Buffer Appender Type";
  static std::string errWtVecCalcBuilder_name = "Rythmos::ErrWtVecCalc";
  static std::string errWtVecCalcBuilderType_name =
    "Error Weight Vector Calculator Type";
  static std::string interpolatorBuilder_name = "Rythmos::Interpolator";
  static std::string interpolatorBuilderType_name = "Interpolator Type";

  // Specific object names:
  static std::string defaultIntegrator_name = "Default Integrator";
  static std::string defaultIntegrator_docs =
    "This Integrator will accept an IntergationControlStrategy, and "
    "can have an IntegrationObserver.  The client can specify the "
    "maximum number of time steps allowed.  The Integrator will loop "
    "over the Stepper until it reaches the requested time. For each "
    "step, the step size will be determined through a couple "
    "mechanisms/filters.  If an Integration Control Strategy has "
    "been specified, the step size and the step type (fixed or "
    "variable) will be determined by it.  Otherwise the step size "
    "will be set to the maximum real value and the step type will "
    "be variable.  Next if the step size is beyond the final time "
    "and the `"+landOnFinalTime_name+"' is specified, the step size is "
    "adjusted to advance the state to the final time.  The Stepper "
    "is passed this step size and type to advance the state.  The "
    "DefaultIntegrator determines the step size and type taken by "
    "the Stepper, and if the step has failed.  If the "
    "IntegrationControlStrategy handles failures, it can suggest "
    "another step size and retry with the Stepper.  Otherwise, the "
    "Integrator will fall through with a failure.  With a successful "
    "step of the Stepper, the Integrator repeats the above until it "
    "reaches the requested time.  Multiple requested times can be "
    "passed to the Integrator.";
  static std::string simpleIntegrationControl_name =
    "Simple Integration Control Strategy";
  static std::string simpleIntegrationControl_docs =
    "This Integration Control Strategy is meant to be simple with "
    "very few parameters controlling it.  Basically the client can "
    "select fixed step type (the Stepper can only take the requested "
    "step size) or variable step type (the Stepper can adjust the step "
    "size to meet accuracy, order, or other criteria).  For fixed step "
    "type, the client can specify the step size and number of steps. "
    "For variable step type, the client can set the maximum step size "
    "allowable.";
  static std::string rampingIntegrationControl_name =
    "Ramping Integration Control Strategy";
  static std::string rampingIntegrationControl_docs =
    "This Integration Control Strategy is very similar to `"
    +simpleIntegrationControl_name+"' except for handling an initial "
    "constant-sized steps followed by a ramped-fixed-sized steps, "
    "and finally variable- or fixed-sized steps.  The client needs to "
    "additionally set the initial step size and the maximum number of "
    "step failures allowed.";
  static std::string rampErrIntegrationControl_name =
    "Ramp and Error Integration Control Strategy";
  static std::string fixedStepControl_name = "Fixed Step Control Strategy";
  static std::string fixedStepControl_docs =
    "This Step Control Strategy can be used for Steppers setup for "
    "variable step type (a stepper that can adjust its step size based "
    "on accuracy, order or other criteria), but would like to make fixed "
    "step sizes or used fixed step size as its default.\n";
  static std::string simpleStepControl_name = "Simple Step Control Strategy";
  static std::string simpleStepControl_docs =
    "This Step Control Strategy starts with the initial step size, "
    "and simply increases or decreases the step size by the "
    "appropriate factor which is based on the change in the "
    "solution relative to the specified relative and absolute "
    "tolerances (|dx| < r*|x| + a) and if solution status from the "
    "solver passes.  Additionally the step size is bounded by the "
    "miminum and maximum step size, and the stepper will fail if "
    "the step size fails more than the specified value.";
  static std::string implicitBDFStepControl_name =
    "Implicit BDF Stepper Step Control Strategy";
  static std::string implicitBDFStepControl_docs =
    "This Step Control Strategy is specifically for use with the `"
    +ImplicitBDFStepper_name+"' Stepper.  The parameters in this list "
    "and sublist are directly related to those available in SUNDAILS/IDA.  "
    "See Hindmarsh, `The PVODE and IDA Algorithms', 2000 for more details. ";
  static std::string implicitBDFStepperErrWtVecCalc_name =
    "Implicit BDF Stepper Error Weight Vector Calculator";
  static std::string implicitBDFStepperErrWtVecCalc_docs =
    "This Error Weight Vector Calculator is specifically for use with the `"
    +ImplicitBDFStepper_name+"' Stepper.";
  static std::string firstOrderErrorStepControl_name =
    "First Order Error Step Control Strategy";
  static std::string firstOrderErrorStepControl_docs =
    "This Step Control Strategy produces a step size based on a first-order "
    "predictor (Forward Euler) and a first-order solution (Backward Euler) by "
    "by using a weight norm of the difference between the predicted and "
    "solution.  See Gresho and Sani, `Incompressible Flow and the Finite "
    "Element Method', Vol. 1, 1998, p. 268.";
  static std::string implicitBDFRampingStepControl_name =
    "Implicit BDF Stepper Ramping Step Control Strategy";
  static std::string implicitBDFRampingStepControl_docs =
    "This Step Control Strategy is specifically for use with the `"
    +ImplicitBDFStepper_name+"' Stepper, and has a two-phase approach: "
    "constant step sizes and followed by variable step sizes.  The step "
    "size is adjusted based on the WRMS, see "
    +implicitBDFStepperErrWtVecCalc_name;
  static std::string defaultInterpolationBuffer_name = "Interpolation Buffer";
  static std::string defaultInterpolationBuffer_docs =
    "Sets parameters for the Interpolation Buffer.";
  static std::string pointwiseInterpolationBufferAppender_name =
    "Pointwise Interpolation Buffer Appender";
  static std::string pointwiseInterpolationBufferAppender_docs =
    "Appender that just transfers nodes without any regard for accuracy or "
    "order.";
//  static std::string smartInterpolationBufferAppender_name =
//    "Smart Interpolation Buffer Appender";
  static std::string linearInterpolator_name = "Linear Interpolator";
  static std::string linearInterpolator_docs =
    "This provides a simple linear interpolation between time nodes.";
  static std::string hermiteInterpolator_name = "Hermite Interpolator";
  static std::string hermiteInterpolator_docs =
    "This provides a piecewise cubic Hermite interpolation on each interval "
    "where the data is the solution and its time derivatives at the end "
    "points of the interval.  It will match 3rd degree polynomials exactly "
    "with both function values and derivatives.";
  static std::string cubicSplineInterpolator_name = "Cubic Spline Interpolator";
  static std::string cubicSplineInterpolator_docs =
    "This provides a cubic spline interpolation between time nodes.";

} // namespace


namespace Rythmos {


template<class Scalar>
IntegratorBuilder<Scalar>::IntegratorBuilder()
{
  this->initializeDefaults_();
}


template<class Scalar>
IntegratorBuilder<Scalar>::~IntegratorBuilder()
{
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setIntegratorFactory(
  const RCP<const Teuchos::AbstractFactory<IntegratorBase<Scalar> > > &integratorFactory,
  const std::string &integratorName
  )
{
  integratorBuilder_->setObjectFactory(integratorFactory, integratorName);
  validPL_ = Teuchos::null;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setIntegrationControlFactory(
  const RCP<const Teuchos::AbstractFactory<IntegrationControlStrategyBase<Scalar> > > &integrationControlFactory,
  const std::string &integrationControlName
  )
{
  integrationControlBuilder_->setObjectFactory(integrationControlFactory,
                                               integrationControlName);
  validPL_ = Teuchos::null;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setStepperBuilder(
    const RCP<StepperBuilder<Scalar> > &stepperBuilder
    )
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(stepperBuilder));
  stepperBuilder_ = stepperBuilder;
  validPL_ = Teuchos::null;
}


template<class Scalar>
RCP<StepperBuilder<Scalar> > IntegratorBuilder<Scalar>::getStepperBuilder()
{
  return stepperBuilder_;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setRKButcherTableauBuilder(
    const RCP<RKButcherTableauBuilder<Scalar> > & rkbtBuilder
    )
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(rkbtBuilder));
  rkbtBuilder_ = rkbtBuilder;
  validPL_ = Teuchos::null;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setStepControlFactory(
  const RCP<const Teuchos::AbstractFactory<StepControlStrategyBase<Scalar> > > &
    stepControlStrategyFactory,
  const std::string &stepControlName
  )
{
  stepControlBuilder_->setObjectFactory(stepControlStrategyFactory,
                                        stepControlName);
  validPL_ = Teuchos::null;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setInterpolationBufferFactory(
  const RCP<const Teuchos::AbstractFactory<InterpolationBufferBase<Scalar> > > &interpolationBufferFactory,
  const std::string &interpolationBufferName
  )
{
  interpolationBufferBuilder_->setObjectFactory(interpolationBufferFactory, interpolationBufferName);
  validPL_ = Teuchos::null;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setInterpolationBufferAppenderFactory(
  const RCP<const Teuchos::AbstractFactory<InterpolationBufferAppenderBase<Scalar> > > &interpolationBufferAppenderFactory,
  const std::string &interpolationBufferAppenderName
  )
{
  interpolationBufferAppenderBuilder_->setObjectFactory(
    interpolationBufferAppenderFactory, interpolationBufferAppenderName);
  validPL_ = Teuchos::null;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setErrWtVecCalcFactory(
    const RCP<const Teuchos::AbstractFactory<ErrWtVecCalcBase<Scalar> > > &
      errWtVecCalcFactory,
    const std::string &errWtVecCalcFactoryName
    )
{
  errWtVecCalcBuilder_->setObjectFactory(errWtVecCalcFactory,
                                         errWtVecCalcFactoryName);
  validPL_ = Teuchos::null;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setInterpolatorFactory(
    const RCP<const Teuchos::AbstractFactory<InterpolatorBase<Scalar> > > &
      interpolatorFactory,
    const std::string &interpolatorFactoryName
    )
{
  interpolatorBuilder_->setObjectFactory(interpolatorFactory,
                                         interpolatorFactoryName);
  validPL_ = Teuchos::null;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setWFactoryObject(
    const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &wFactoryObject
    )
{
  TEUCHOS_ASSERT( !is_null(wFactoryObject) );
  wFactoryObject_ = wFactoryObject;
}


template<class Scalar>
void IntegratorBuilder<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*this->getValidParameters());
  paramList_ = paramList;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
IntegratorBuilder<Scalar>::getValidParameters() const
{
  if (is_null(validPL_)) {
    RCP<ParameterList> pl = Teuchos::parameterList();

    // Integrator Settings
    ParameterList& integratorSettingsPL =
      pl->sublist(integratorSettings_name,false,integratorSettings_docs);
    {
      // Initial Time
      integratorSettingsPL.set(initialTime_name,
                               Teuchos::as<Scalar>(initialTime_default),
                               initialTime_docs);
      // Final Time
      integratorSettingsPL.set(finalTime_name,
                               Teuchos::as<Scalar>(finalTime_default),
                               finalTime_docs);
      // Land On Final Time
      integratorSettingsPL.set(landOnFinalTime_name,landOnFinalTime_default,
                               landOnFinalTime_docs);
      // Integrator Selection
      ParameterList& integratorSelectionPL =
        integratorSettingsPL.sublist(integratorSelection_name,false,
                                     integratorSelection_docs)
                            .disableRecursiveValidation();
        // Default Integrator
        integratorSelectionPL.sublist(defaultIntegrator_name,false,
                                      defaultIntegrator_docs)
                              .disableRecursiveValidation();
      integratorSelectionPL
        .setParameters(*(integratorBuilder_->getValidParameters()));
    }

    // Integration Control Selection
    ParameterList& integrationControlSelectionPL =
      pl->sublist(integrationControlSelection_name,false,
                  integrationControlSelection_docs)
         .disableRecursiveValidation();
      // Simple Integration Control Strategy
      integrationControlSelectionPL.sublist(simpleIntegrationControl_name,false,
                                            simpleIntegrationControl_docs)
                                   .disableRecursiveValidation();
      // Ramping Integration Control Strategy
      integrationControlSelectionPL.sublist(rampingIntegrationControl_name,
                                            false,
                                            rampingIntegrationControl_docs)
                                   .disableRecursiveValidation();
    integrationControlSelectionPL
      .setParameters(*(integrationControlBuilder_->getValidParameters()));

    // Stepper Settings
    ParameterList& stepperSettingsPL = pl->sublist(stepperSettings_name, false,
                                                   stepperSettings_docs);
    {
      // Stepper Selection
      ParameterList& stepperSelectionPL =
        stepperSettingsPL.sublist(stepperSelection_name, false,
                                  stepperSelection_docs)
                         .disableRecursiveValidation();
        // Forward Euler
        stepperSelectionPL.sublist(ForwardEulerStepper_name, false,
                                   ForwardEulerStepper_docs)
                          .disableRecursiveValidation();
        // Backward Euler
        stepperSelectionPL.sublist(BackwardEulerStepper_name, false,
                                   BackwardEulerStepper_docs)
                          .disableRecursiveValidation();
        // Implicit BDF
        stepperSelectionPL.sublist(ImplicitBDFStepper_name, false,
                                   ImplicitBDFStepper_docs)
                          .disableRecursiveValidation();
        // Explicit RK
        stepperSelectionPL.sublist(ExplicitRKStepper_name, false,
                                   ExplicitRKStepper_docs)
                          .disableRecursiveValidation();
        // Implicit RK
        stepperSelectionPL.sublist(ImplicitRKStepper_name, false,
                                   ImplicitRKStepper_docs)
                          .disableRecursiveValidation();
      stepperSelectionPL
        .setParameters(*(stepperBuilder_->getValidParameters()));
      // Step Control Settings
      ParameterList& stepControlSettingsPL =
        stepperSettingsPL.sublist(stepControlSettings_name, false,
                                  stepControlSettings_docs);
      {
        // Step Control Selection
        ParameterList& stepControlSelectionPL =
          stepControlSettingsPL.sublist(stepControlSelection_name,false,
                                        stepControlSelection_docs)
                               .disableRecursiveValidation();
          // Fixed Step Control Strategy
          stepControlSelectionPL.sublist(fixedStepControl_name,false,
                                         fixedStepControl_docs)
                                .disableRecursiveValidation();
          // Simple Step Control Strategy
          stepControlSelectionPL.sublist(simpleStepControl_name,false,
                                         simpleStepControl_docs)
                                .disableRecursiveValidation();
          // First Order Error Step Control Strategy
          stepControlSelectionPL.sublist(firstOrderErrorStepControl_name,false,
                                         firstOrderErrorStepControl_docs)
                                .disableRecursiveValidation();
          // Implicit BDF Stepper Step Control Strategy
          stepControlSelectionPL.sublist(implicitBDFStepControl_name,false,
                                         implicitBDFStepControl_docs)
                                .disableRecursiveValidation();
          // Implicit BDF Stepper Ramping Step Control Strategy
          stepControlSelectionPL.sublist(implicitBDFRampingStepControl_name,
                                         false,
                                         implicitBDFRampingStepControl_docs)
                                .disableRecursiveValidation();
        stepControlSelectionPL
          .setParameters(*(stepControlBuilder_->getValidParameters()));

        // ErrWtVec Selection
        ParameterList& errWtVecSelectionPL =
          stepControlSettingsPL.sublist(errWtVecSelection_name,false,
                                        errWtVecSelection_docs)
                               .disableRecursiveValidation();
          // Implicit BDF Stepper Error Weight Vector Calculator
          errWtVecSelectionPL.sublist(implicitBDFStepperErrWtVecCalc_name,
                                      false,
                                      implicitBDFStepperErrWtVecCalc_docs)
                             .disableRecursiveValidation();
        errWtVecSelectionPL
          .setParameters(*(errWtVecCalcBuilder_->getValidParameters()));
      }
      // Interpolator Selection
      ParameterList& interpolatorSelectionPL =
        stepperSettingsPL.sublist(interpolatorSelection_name,false,
                                  stepperInterpolatorSelection_docs)
                         .disableRecursiveValidation();
        // Linear Interpolator
        interpolatorSelectionPL.sublist(linearInterpolator_name, false,
                                        linearInterpolator_docs)
                               .disableRecursiveValidation();
        // Hermite Interpolator
        interpolatorSelectionPL.sublist(hermiteInterpolator_name, false,
                                        hermiteInterpolator_docs)
                               .disableRecursiveValidation();
        // Cubic Spline Interpolator
        interpolatorSelectionPL.sublist(cubicSplineInterpolator_name, false,
                                        cubicSplineInterpolator_docs)
                               .disableRecursiveValidation();
      interpolatorSelectionPL
        .setParameters(*(interpolatorBuilder_->getValidParameters()));

      // RKBT Selection
      ParameterList& rkbtSelectionPL =
        stepperSettingsPL.sublist(rkButcherTableauSelection_name,false,
                                  rkButcherTableauSelection_docs)
                         .disableRecursiveValidation();
      rkbtSelectionPL.setParameters(*(rkbtBuilder_->getValidParameters()));
      // Nonlinear Solver Selection (TODO)
    }

    // Interpolation Buffer Settings
    ParameterList& interpolationBufferSettingsPL =
      pl->sublist(interpolationBufferSettings_name,false,
                  interpolationBufferSettings_docs);
    {
      // Interpolation Buffer Selection
      ParameterList& interpolationBufferSelectionPL =
        interpolationBufferSettingsPL.sublist(interpolationBufferSelection_name,
                                              false,
                                              interpolationBufferSelection_docs)
                                     .disableRecursiveValidation();
        // Interpolation Buffer
        interpolationBufferSelectionPL
          .sublist(defaultInterpolationBuffer_name, false,
                   defaultInterpolationBuffer_docs)
          .disableRecursiveValidation();
      interpolationBufferSelectionPL
        .setParameters(*(interpolationBufferBuilder_->getValidParameters()));
      // Interpolation Buffer Appender Selection
      ParameterList& interpolationBufferAppenderSelectionPL =
        interpolationBufferSettingsPL
          .sublist(interpolationBufferAppenderSelection_name, false,
                   interpolationBufferAppenderSelection_docs)
          .disableRecursiveValidation();
        // Pointwise Interpolation Buffer Appender
        interpolationBufferAppenderSelectionPL
          .sublist(pointwiseInterpolationBufferAppender_name,false,
                   pointwiseInterpolationBufferAppender_docs)
          .disableRecursiveValidation();
      interpolationBufferAppenderSelectionPL
        .setParameters(*(interpolationBufferAppenderBuilder_->getValidParameters()));
      // Interpolator Selection
      ParameterList& interpolatorSelectionPL =
        interpolationBufferSettingsPL.sublist(interpolatorSelection_name,false,
                                              interpolatorSelection_docs)
                                     .disableRecursiveValidation();
        // Linear Interpolator
        interpolatorSelectionPL.sublist(linearInterpolator_name, false,
                                        linearInterpolator_docs)
                               .disableRecursiveValidation();
        // Hermite Interpolator
        interpolatorSelectionPL.sublist(hermiteInterpolator_name, false,
                                        hermiteInterpolator_docs)
                               .disableRecursiveValidation();
        // Cubic Spline Interpolator
        interpolatorSelectionPL.sublist(cubicSplineInterpolator_name, false,
                                        cubicSplineInterpolator_docs)
                               .disableRecursiveValidation();
      interpolatorSelectionPL
        .setParameters(*(interpolatorBuilder_->getValidParameters()));
    }

    // Integration Observer Settings

    validPL_ = pl;
  }
  return validPL_;
}


template<class Scalar>
RCP<ParameterList> IntegratorBuilder<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<ParameterList> IntegratorBuilder<Scalar>::unsetParameterList()
{
  RCP<ParameterList> pl = paramList_;
  paramList_ = Teuchos::null;
  return pl;
}


template<class Scalar>
RCP<const ParameterList> IntegratorBuilder<Scalar>::getParameterList() const
{
  return paramList_;
}


// Where should we throw exceptions?
// 1.  If the integrator comes back null (done)
// 2.  If the stepper comes back null (done)
// 3.  If model is null (done)
// 4.  If the stepper is implicit and nlSolver is null (done)
// 5.  If the stepper accepts an RKBT but "None" is selected (done)
//
// a.  Its okay if the integration control comes back null, the
//     IntegrationControlStrategyAcceptingIntegratorBase will deal with it
// b.  Its okay if the step control comes back null, the
//     StepControlStrategyAcceptingStepperBase will deal with it
template<class Scalar>
RCP<IntegratorBase<Scalar> >
IntegratorBuilder<Scalar>::create(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& initialCondition,
    const RCP<Thyra::NonlinearSolverBase<Scalar> >& nlSolver
    ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( is_null(model), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  The model passed in is null!"
      );
  TEUCHOS_TEST_FOR_EXCEPTION( is_null(paramList_), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  Please set a parameter "
      "list on this class before calling create."
      );
  RCP<ParameterList> integratorSettingsPL = sublist(paramList_,
                                                    integratorSettings_name);

  // Create the integrator first
  RCP<ParameterList> integratorSelectionPL = sublist(integratorSettingsPL,
                                                     integratorSelection_name);
  integratorBuilder_->setParameterList(integratorSelectionPL);
  RCP<IntegratorBase<Scalar> > integrator = integratorBuilder_->create();
  TEUCHOS_TEST_FOR_EXCEPTION( is_null(integrator), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  The integrator came back "
      "null from the ObjectBuilder!"
      );

  // Check for IntegrationControlStrategy and set it on the integrator
  RCP<IntegrationControlStrategyAcceptingIntegratorBase<Scalar> >
    icsaIntegrator =
      Teuchos::rcp_dynamic_cast<IntegrationControlStrategyAcceptingIntegratorBase<Scalar> >(integrator,false);
  if (!is_null(icsaIntegrator)) {
    RCP<ParameterList> integrationControlSelectionPL =
      sublist(paramList_,integrationControlSelection_name);
    integrationControlBuilder_->setParameterList(integrationControlSelectionPL);
    RCP<IntegrationControlStrategyBase<Scalar> > integrationControl =
      integrationControlBuilder_->create();
    if (!is_null(integrationControl)) {
      icsaIntegrator->setIntegrationControlStrategy(integrationControl);
    }
  }
  RCP<ParameterList> interpolationBufferSettingsPL =
    sublist(paramList_,interpolationBufferSettings_name);

  // Check for a trailing interpolation buffer and set it on the integrator
  RCP<TrailingInterpolationBufferAcceptingIntegratorBase<Scalar> >
    tibaIntegrator =
      Teuchos::rcp_dynamic_cast<TrailingInterpolationBufferAcceptingIntegratorBase<Scalar> >(integrator,false);
  if (!is_null(tibaIntegrator)) {
    RCP<ParameterList> interpolationBufferSelectionPL =
      sublist(interpolationBufferSettingsPL,interpolationBufferSelection_name);
    interpolationBufferBuilder_->setParameterList(interpolationBufferSelectionPL);
    RCP<InterpolationBufferBase<Scalar> > ib =
      interpolationBufferBuilder_->create();
    if (!is_null(ib)) {
      // Check for an interpolator
      RCP<InterpolatorAcceptingObjectBase<Scalar> > iaobIB =
        Teuchos::rcp_dynamic_cast<InterpolatorAcceptingObjectBase<Scalar> >(ib,false);
      if (!is_null(iaobIB)) {
        RCP<ParameterList> interpolatorSelectionPL =
          sublist(interpolationBufferSettingsPL,interpolatorSelection_name);
        interpolatorBuilder_->setParameterList(interpolatorSelectionPL);
        RCP<InterpolatorBase<Scalar> > interpolator =
          interpolatorBuilder_->create();
        if (!is_null(interpolator)) {
          iaobIB->setInterpolator(interpolator);
        }
      }
      tibaIntegrator->setTrailingInterpolationBuffer(ib);
    }
  }

  // Check for an InterpolationBufferAppender and set it on the integrator
  RCP<InterpolationBufferAppenderAcceptingIntegratorBase<Scalar> > ibaaIntegrator =
    Teuchos::rcp_dynamic_cast<InterpolationBufferAppenderAcceptingIntegratorBase<Scalar> >(integrator,false);
  if (!is_null(ibaaIntegrator)) {
    RCP<ParameterList> interpolationBufferAppenderSelectionPL =
      sublist(interpolationBufferSettingsPL,
              interpolationBufferAppenderSelection_name);
    interpolationBufferAppenderBuilder_->setParameterList(interpolationBufferAppenderSelectionPL);
    RCP<InterpolationBufferAppenderBase<Scalar> > interpolationBufferAppender =
       interpolationBufferAppenderBuilder_->create();
    if (!is_null(interpolationBufferAppender)) {
      ibaaIntegrator->setInterpolationBufferAppender(interpolationBufferAppender);
    }
  }
  RCP<ParameterList> stepperSettingsPL =
    sublist(paramList_,stepperSettings_name);

  // Create the Stepper
  RCP<ParameterList> stepperSelectionPL = sublist(stepperSettingsPL,
                                                  stepperSelection_name);
  stepperBuilder_->setParameterList(stepperSelectionPL);
  RCP<StepperBase<Scalar> > stepper = stepperBuilder_->create();
  TEUCHOS_TEST_FOR_EXCEPTION( is_null(stepper), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  The stepper came back "
      "null from the StepperBuilder!");

  // Create the Step Control
  RCP<ParameterList> stepControlSettingsPL =
    sublist(stepperSettingsPL,stepControlSettings_name);
  RCP<StepControlStrategyAcceptingStepperBase<Scalar> > scsaStepper =
    Teuchos::rcp_dynamic_cast<StepControlStrategyAcceptingStepperBase<Scalar> >(stepper,false);
  if (!is_null(scsaStepper)) {
    RCP<ParameterList> stepControlSelectionPL =
      sublist(stepControlSettingsPL,stepControlSelection_name);
    stepControlBuilder_->setParameterList(stepControlSelectionPL);
    RCP<StepControlStrategyBase<Scalar> > stepControl =
      stepControlBuilder_->create();
    if (!is_null(stepControl)) {
      // Create the ErrWtVecCalc
      RCP<ErrWtVecCalcAcceptingStepControlStrategyBase<Scalar> >
        ewvcaStepControl =
          Teuchos::rcp_dynamic_cast<ErrWtVecCalcAcceptingStepControlStrategyBase<Scalar> >(stepControl,false);
      if (!is_null(ewvcaStepControl)) {
        RCP<ParameterList> errWtVecSelectionPL =
          sublist(stepControlSettingsPL,errWtVecSelection_name);
        errWtVecCalcBuilder_->setParameterList(errWtVecSelectionPL);
        RCP<ErrWtVecCalcBase<Scalar> > errWtVecCalc =
          errWtVecCalcBuilder_->create();
        if (!is_null(errWtVecCalc)) {
          ewvcaStepControl->setErrWtVecCalc(errWtVecCalc);
        }
      }
      scsaStepper->setStepControlStrategy(stepControl);
    }
  }

  // Check for an Interpolator
  RCP<InterpolatorAcceptingObjectBase<Scalar> > iaobStepper =
    Teuchos::rcp_dynamic_cast<InterpolatorAcceptingObjectBase<Scalar> >(stepper,
                                                                        false);
  if (!is_null(iaobStepper)) {
    RCP<ParameterList> interpolatorSelectionPL =
      sublist(stepperSettingsPL,interpolatorSelection_name);
    interpolatorBuilder_->setParameterList(interpolatorSelectionPL);
    RCP<InterpolatorBase<Scalar> > interpolator =
      interpolatorBuilder_->create();
    if (!is_null(interpolator)) {
      iaobStepper->setInterpolator(interpolator);
    }
  }

  // Check for an RKBT Selection
  RCP<RKButcherTableauAcceptingStepperBase<Scalar> > rkbtaStepper =
    Teuchos::rcp_dynamic_cast<RKButcherTableauAcceptingStepperBase<Scalar> >(stepper,false);
  if (!is_null(rkbtaStepper)) {
    RCP<ParameterList> rkButcherTableauSelectionPL =
      sublist(stepperSettingsPL,rkButcherTableauSelection_name);
    rkbtBuilder_->setParameterList(rkButcherTableauSelectionPL);
    RCP<RKButcherTableauBase<Scalar> > rkbt = rkbtBuilder_->create();
    TEUCHOS_TEST_FOR_EXCEPTION( is_null(rkbt), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  The Stepper accepts a "
      "RK Butcher Tableau, but none were specified!"
      );
    rkbtaStepper->setRKButcherTableau(rkbt);
  }

  // Check for a W Factory
  RCP<ImplicitRKStepper<Scalar> > irkStepper =
    Teuchos::rcp_dynamic_cast<ImplicitRKStepper<Scalar> >(stepper,false);
  if (!is_null(irkStepper)) {
    if (!is_null(wFactoryObject_)) {
      irkStepper->set_W_factory(wFactoryObject_);
    }
  }

  // Check for Nonlinear Solver Selection (TODO)
  // Set model on stepper
  stepper->setModel(model);
  // Set initial condition on stepper
  stepper->setInitialCondition(initialCondition);
  // Set nonlinear solver on stepper
  RCP<SolverAcceptingStepperBase<Scalar> > saStepper =
    Teuchos::rcp_dynamic_cast<SolverAcceptingStepperBase<Scalar> >(stepper,
                                                                   false);
  if(!is_null(saStepper)) {
    TEUCHOS_TEST_FOR_EXCEPTION( is_null(nlSolver), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  The nonlinear solver passed "
      "in is null and the stepper is implicit!"
      );
    saStepper->setSolver(nlSolver);
  }
  Scalar finalTime = integratorSettingsPL->get<Scalar>(
    finalTime_name, Teuchos::as<Scalar>(finalTime_default));
  bool landOnFinalTime = integratorSettingsPL->get<bool>(
    landOnFinalTime_name, landOnFinalTime_default);
  integrator->setStepper(stepper,finalTime,landOnFinalTime);
  return integrator;
}

template<class Scalar>
void IntegratorBuilder<Scalar>::initializeDefaults_()
{

  using Teuchos::abstractFactoryStd;

  // Integrator
  integratorBuilder_ = Teuchos::objectBuilder<IntegratorBase<Scalar> >();
  integratorBuilder_->setObjectName(integratorBuilder_name);
  integratorBuilder_->setObjectTypeName(integratorBuilderType_name);
  integratorBuilder_->setObjectFactory(
      abstractFactoryStd< IntegratorBase<Scalar>, DefaultIntegrator<Scalar> >(),
      defaultIntegrator_name);

  // Integration Control Strategy
  integrationControlBuilder_ =
    Teuchos::objectBuilder<IntegrationControlStrategyBase<Scalar> >();
  integrationControlBuilder_->setObjectName(integrationControlBuilder_name);
  integrationControlBuilder_->setObjectTypeName(integrationControlBuilderType_name);
  integrationControlBuilder_->setObjectFactory(
      abstractFactoryStd< IntegrationControlStrategyBase<Scalar>,
                          SimpleIntegrationControlStrategy<Scalar> >(),
      simpleIntegrationControl_name);
  integrationControlBuilder_->setObjectFactory(
      abstractFactoryStd< IntegrationControlStrategyBase<Scalar>,
                          RampingIntegrationControlStrategy<Scalar> >(),
      rampingIntegrationControl_name);
  integrationControlBuilder_->setDefaultObject("None");

  // Stepper Builder
  stepperBuilder_ = stepperBuilder<Scalar>();

  // RKBT Builder
  rkbtBuilder_ = rKButcherTableauBuilder<Scalar>();

  // Step Control Strategy
  stepControlBuilder_ =
    Teuchos::objectBuilder<StepControlStrategyBase<Scalar> >();
  stepControlBuilder_->setObjectName(stepControlBuilder_name);
  stepControlBuilder_->setObjectTypeName(stepControlBuilderType_name);
  stepControlBuilder_->setObjectFactory(
      abstractFactoryStd< StepControlStrategyBase<Scalar>,
                          FixedStepControlStrategy<Scalar> >(),
      fixedStepControl_name);
  stepControlBuilder_->setObjectFactory(
      abstractFactoryStd< StepControlStrategyBase<Scalar>,
                          SimpleStepControlStrategy<Scalar> >(),
      simpleStepControl_name);
  stepControlBuilder_->setObjectFactory(
      abstractFactoryStd< StepControlStrategyBase<Scalar>,
                          FirstOrderErrorStepControlStrategy<Scalar> >(),
      firstOrderErrorStepControl_name);
  stepControlBuilder_->setObjectFactory(
      abstractFactoryStd< StepControlStrategyBase<Scalar>,
                          ImplicitBDFStepperStepControl<Scalar> >(),
      implicitBDFStepControl_name);
  stepControlBuilder_->setObjectFactory(
      abstractFactoryStd< StepControlStrategyBase<Scalar>,
                          ImplicitBDFStepperRampingStepControl<Scalar> >(),
      implicitBDFRampingStepControl_name);
  stepControlBuilder_->setDefaultObject("None");

  // Trailing Interpolation Buffer
  interpolationBufferBuilder_ =
    Teuchos::objectBuilder<InterpolationBufferBase<Scalar> >();
  interpolationBufferBuilder_->setObjectName(interpolationBufferBuilder_name);
  interpolationBufferBuilder_->setObjectTypeName(
    interpolationBufferBuilderType_name);
  interpolationBufferBuilder_->setObjectFactory(
      abstractFactoryStd< InterpolationBufferBase<Scalar>,
                          InterpolationBuffer<Scalar> >(),
      defaultInterpolationBuffer_name);
  interpolationBufferBuilder_->setDefaultObject("None");

  // Interpolation Buffer Appender
  interpolationBufferAppenderBuilder_ =
    Teuchos::objectBuilder<InterpolationBufferAppenderBase<Scalar> >();
  interpolationBufferAppenderBuilder_->setObjectName(
    interpolationBufferAppenderBuilder_name);
  interpolationBufferAppenderBuilder_->setObjectTypeName(
    interpolationBufferAppenderBuilderType_name);
//  interpolationBufferAppenderBuilder_->setObjectFactory(
//      abstractFactoryStd< InterpolationBufferAppenderBase<Scalar>,
//                          SmartInterpolationBufferAppender<Scalar> >(),
//      smartInterpolationBufferAppender_name);
  interpolationBufferAppenderBuilder_->setObjectFactory(
      abstractFactoryStd< InterpolationBufferAppenderBase<Scalar>,
                          PointwiseInterpolationBufferAppender<Scalar> >(),
      pointwiseInterpolationBufferAppender_name
      );
  interpolationBufferAppenderBuilder_->setDefaultObject("None");

  // ErrWtVecCalc
  errWtVecCalcBuilder_ = Teuchos::objectBuilder<ErrWtVecCalcBase<Scalar> >();
  errWtVecCalcBuilder_->setObjectName(errWtVecCalcBuilder_name);
  errWtVecCalcBuilder_->setObjectTypeName(errWtVecCalcBuilderType_name);
  errWtVecCalcBuilder_->setObjectFactory(
      abstractFactoryStd< ErrWtVecCalcBase<Scalar>,
                          ImplicitBDFStepperErrWtVecCalc<Scalar> >(),
      implicitBDFStepperErrWtVecCalc_name);
  errWtVecCalcBuilder_->setDefaultObject("None");

  // Interpolator
  interpolatorBuilder_ = Teuchos::objectBuilder<InterpolatorBase<Scalar> >();
  interpolatorBuilder_->setObjectName(interpolatorBuilder_name);
  interpolatorBuilder_->setObjectTypeName(interpolatorBuilderType_name);
  interpolatorBuilder_->setObjectFactory(
      abstractFactoryStd< InterpolatorBase<Scalar>,
                          LinearInterpolator<Scalar> >(),
      linearInterpolator_name);
  interpolatorBuilder_->setObjectFactory(
      abstractFactoryStd< InterpolatorBase<Scalar>,
                          HermiteInterpolator<Scalar> >(),
      hermiteInterpolator_name);
  interpolatorBuilder_->setObjectFactory(
      abstractFactoryStd< InterpolatorBase<Scalar>,
                          CubicSplineInterpolator<Scalar> >(),
      cubicSplineInterpolator_name);
  interpolatorBuilder_->setDefaultObject("None");

}


} // namespace Rythmos


template<class Scalar>
Teuchos::RCP<Rythmos::IntegratorBuilder<Scalar> >
Rythmos::integratorBuilder()
{
  return rcp(new IntegratorBuilder<Scalar>);
}


template<class Scalar>
Teuchos::RCP<Rythmos::IntegratorBuilder<Scalar> >
Rythmos::integratorBuilder(const RCP<ParameterList> &paramList)
{
  const RCP<IntegratorBuilder<Scalar> > ib = integratorBuilder<Scalar>();
  ib->setParameterList(paramList);
  return ib;
}

template<class Scalar>
Teuchos::RCP<Rythmos::IntegratorBase<Scalar> > Rythmos::createForwardSensitivityIntegrator(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    const int& p_index,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& model_ic,
    const RCP<Thyra::NonlinearSolverBase<Scalar> >& nlSolver,
    const RCP<Teuchos::ParameterList>& integratorBuilderPL
    )
{
  RCP<IntegratorBuilder<Scalar> > ib = integratorBuilder<Scalar>(integratorBuilderPL);
  RCP<IntegratorBase<Scalar> > sensIntegrator = ib->create(model,model_ic,nlSolver);
  RCP<ForwardSensitivityStepper<Scalar> > stateAndSensStepper =
    forwardSensitivityStepper<Scalar>();
  stateAndSensStepper->initializeSyncedSteppers(
    model, p_index, model_ic, sensIntegrator->getNonconstStepper(), nlSolver
    );
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> state_and_sens_ic =
    createStateAndSensInitialCondition(*stateAndSensStepper, model_ic);
  stateAndSensStepper->setInitialCondition(state_and_sens_ic);
  sensIntegrator->setStepper(stateAndSensStepper, sensIntegrator->getFwdTimeRange().upper());
  return sensIntegrator;
}



//
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_INTEGRATOR_BUILDER_INSTANT(SCALAR) \
  \
  template class IntegratorBuilder< SCALAR >; \
  \
  template RCP<IntegratorBuilder< SCALAR > > \
  integratorBuilder(); \
  \
  template RCP<IntegratorBuilder< SCALAR > > \
  integratorBuilder(const RCP<ParameterList> &paraList); \
  \
  template RCP<IntegratorBase< SCALAR > > \
  createForwardSensitivityIntegrator( \
    const RCP<const Thyra::ModelEvaluator< SCALAR > >& model, \
    const int& p_index, \
    const Thyra::ModelEvaluatorBase::InArgs< SCALAR >& model_ic, \
    const RCP<Thyra::NonlinearSolverBase< SCALAR > >& nlSolver, \
    const RCP<ParameterList>& integratorBuilderPL \
    );


#endif //Rythmos_INTEGRATOR_BUILDER_DEF_H
