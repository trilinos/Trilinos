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


#ifndef Rythmos_INTEGRATOR_BUILDER_H
#define Rythmos_INTEGRATOR_BUILDER_H

// Rythmos classes:
#include "Rythmos_Types.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_IntegrationControlStrategyAcceptingIntegratorBase.hpp"
#include "Rythmos_StepperBuilder.hpp"
#include "Rythmos_StepControlStrategyBase.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Rythmos_InterpolationBufferBase.hpp"
#include "Rythmos_InterpolationBufferAppenderBase.hpp"
#include "Rythmos_ErrWtVecCalcBase.hpp"

// Specific objects to seed the builder:
#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_ImplicitBDFStepperStepControl.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_PointwiseInterpolationBufferAppender.hpp"
//#include "Rythmos_SmartInterpolationBufferAppender.hpp"
#include "Rythmos_ImplicitBDFStepperErrWtVecCalc.hpp"

// Teuchos:
#include "Teuchos_ObjectBuilder.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_as.hpp"

// Thyra:
#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_ModelEvaluator.hpp"

namespace {

  // Top level Parameter Sublist names:
  static std::string integrator_name = "Integrator";
  static std::string integrationControl_name = "Integration Control";
  static std::string stepper_name = "Stepper";
  static std::string stepControl_name = "Step Control";
  static std::string interpolationBuffer_name = "Trailing Interpolation Buffer";
  static std::string interpolationBufferAppender_name = "Interpolation Buffer Appender";
  static std::string errWtVecCalc_name = "Error Weight Vector Calculator";
  // Top level Parameter names:
  static std::string finalTime_name = "Final Time";
  static int finalTime_default = 1; // Should be Scalar(1.0)
  static std::string landOnFinalTime_name = "Land On Final Time";
  static bool landOnFinalTime_default = true;

  // Builder names:
  static std::string integratorBuilder_name = "Rythmos::Integrator";
  static std::string integratorBuilderType_name = "Integrator Type";
  static std::string integrationControlBuilder_name = "Rythmos::IntegrationControlStrategy";
  static std::string integrationControlBuilderType_name = "Integration Control Strategy Type";
  static std::string stepControlBuilder_name = "Rythmos::StepControlStrategy";
  static std::string stepControlBuilderType_name = "Step Control Strategy Type";
  static std::string interpolationBufferBuilder_name = "Rythmos::InterpolationBuffer";
  static std::string interpolationBufferBuilderType_name = "Interpolation Buffer Type";
  static std::string interpolationBufferAppenderBuilder_name = "Rythmos::InterpolationBufferAppender";
  static std::string interpolationBufferAppenderBuilderType_name = "Interpolation Buffer Appender Type";
  static std::string errWtVecCalcBuilder_name = "Rythmos::ErrWtVecCalc";
  static std::string errWtVecCalcBuilderType_name = "Error Weight Vector Calculator Type";

  // Specific object names:
  static std::string defaultIntegrator_name = "Default Integrator";
  static std::string simpleIntegrationControl_name = "Simple Integration Control Strategy";
  static std::string implicitBDFStepControl_name = "Implicit BDF Step Control Strategy";
  static std::string defaultInterpolationBuffer_name = "Interpolation Buffer";
  static std::string pointwiseInterpolationBufferAppender_name = "Pointwise Interpolation Buffer Appender";
//  static std::string smartInterpolationBufferAppender_name = "Smart Interpolation Buffer Appender";
  static std::string implicitBDFStepperErrWtVecCalc_name = "Implicit BDF Stepper Error Weight Vector Calculator";

} // namespace

namespace Rythmos {

template<class Scalar>
  class IntegratorBuilder : virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \brief . */
  IntegratorBuilder();

  /** \brief . */
  virtual ~IntegratorBuilder();

  /** \brief Set a new Integrator factory object. */
  void setIntegratorFactory(
    const RCP<const Teuchos::AbstractFactory<IntegratorBase<Scalar> > > &integratorFactory,
    const std::string &integratorFactoryName
    );

  /** \brief Set a new Integration Control Strategy factory object. */
  void setIntegrationControlFactory(
    const RCP<const Teuchos::AbstractFactory<IntegrationControlStrategyBase<Scalar> > > &integrationControlFactory,
    const std::string &integrationControlName
    );

  /** \brief Set the Stepper Builder object. */
  void setStepperBuilder(
    const RCP<StepperBuilder<Scalar> > &stepperBuilder
    );

  /** \brief Set a new Step Control Strategy factory object. */
  void setStepControlFactory(
    const RCP<const Teuchos::AbstractFactory<StepControlStrategyBase<Scalar> > > &stepControlStrategyFactory,
    const std::string &stepControlName
    );

  /** \brief Set an InterpolationBuffer factory object. */
  void setInterpolationBufferFactory(
    const RCP<const Teuchos::AbstractFactory<InterpolationBufferBase<Scalar> > > &interpolationBufferFactory,
    const std::string &interpolationBufferName
    );

  /** \brief Set an InterpolationBufferAppender factory object. */
  void setInterpolationBufferAppenderFactory(
    const RCP<const Teuchos::AbstractFactory<InterpolationBufferAppenderBase<Scalar> > > &interpolationBufferAppenderFactory,
    const std::string &interpolationBufferAppenderName
    );

  /** \brief Set an ErrWtVecCalc factory object. */
  void setErrWtVecCalcFactory(
    const RCP<const Teuchos::AbstractFactory<ErrWtVecCalcBase<Scalar> > > &errWtVecCalcFactory,
    const std::string &errWtVecCalcFactoryName
    );
  
  /** \brief . */
  RCP<IntegratorBase<Scalar> > create(
    const RCP<const Thyra::ModelEvaluator<Scalar> > model,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver
      ) const;
  
  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(const RCP<Teuchos::ParameterList> & paramList);
  
  /** \brief. */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  /** \brief. */
  RCP<ParameterList> getNonconstParameterList();

  /** \brief. */
  RCP<ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const ParameterList> getParameterList() const;
 
  //@}

private:

  // //////////////////////////////////////
  // Private data members

  RCP<Teuchos::ObjectBuilder<IntegratorBase<Scalar> > > integratorBuilder_;
  RCP<Teuchos::ObjectBuilder<IntegrationControlStrategyBase<Scalar> > > integrationControlBuilder_;
  RCP<StepperBuilder<Scalar> > stepperBuilder_;
  RCP<Teuchos::ObjectBuilder<StepControlStrategyBase<Scalar> > > stepControlBuilder_;
  RCP<Teuchos::ObjectBuilder<InterpolationBufferBase<Scalar> > > interpolationBufferBuilder_;
  RCP<Teuchos::ObjectBuilder<InterpolationBufferAppenderBase<Scalar> > > interpolationBufferAppenderBuilder_;
  RCP<Teuchos::ObjectBuilder<ErrWtVecCalcBase<Scalar> > > errWtVecCalcBuilder_;

  RCP<ParameterList> paramList_;
  mutable RCP<ParameterList> validPL_;

  // //////////////////////////////////////
  // Private member functions

  void initializeDefaults_();

};


// Nonmember constructor
template<class Scalar>
RCP<IntegratorBuilder<Scalar> > integratorBuilder()
{
  RCP<IntegratorBuilder<Scalar> > sb = rcp(new IntegratorBuilder<Scalar> );
  return sb;
}

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
  integrationControlBuilder_->setObjectFactory(integrationControlFactory, integrationControlName);
  validPL_ = Teuchos::null;
}

template<class Scalar>
void IntegratorBuilder<Scalar>::setStepperBuilder(
    const RCP<StepperBuilder<Scalar> > &stepperBuilder
    ) 
{
  TEST_FOR_EXCEPT(is_null(stepperBuilder));
  stepperBuilder_ = stepperBuilder;
  validPL_ = Teuchos::null;
}

template<class Scalar>
void IntegratorBuilder<Scalar>::setStepControlFactory(
  const RCP<const Teuchos::AbstractFactory<StepControlStrategyBase<Scalar> > > &stepControlStrategyFactory,
  const std::string &stepControlName
  )
{
  stepControlBuilder_->setObjectFactory(stepControlStrategyFactory, stepControlName);
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
  interpolationBufferAppenderBuilder_->setObjectFactory(interpolationBufferAppenderFactory, interpolationBufferAppenderName);
  validPL_ = Teuchos::null;
}

template<class Scalar>
void IntegratorBuilder<Scalar>::setErrWtVecCalcFactory(
    const RCP<const Teuchos::AbstractFactory<ErrWtVecCalcBase<Scalar> > > &errWtVecCalcFactory,
    const std::string &errWtVecCalcFactoryName
    )
{
  errWtVecCalcBuilder_->setObjectFactory(errWtVecCalcFactory,errWtVecCalcFactoryName);
  validPL_ = Teuchos::null;
}

template<class Scalar>
void IntegratorBuilder<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*this->getValidParameters());
  paramList_ = paramList;
}



template<class Scalar>
RCP<const Teuchos::ParameterList>
IntegratorBuilder<Scalar>::getValidParameters() const
{
  if (is_null(validPL_)) {
    RCP<ParameterList> pl = Teuchos::parameterList();

    // Integrator 
    ParameterList& integratorPL = 
      pl->sublist(integrator_name).disableRecursiveValidation();
    integratorPL.setParameters(*(integratorBuilder_->getValidParameters()));

    // Integration Control 
    ParameterList& integrationControlPL = 
      pl->sublist(integrationControl_name).disableRecursiveValidation();
    integrationControlPL.setParameters(*(integrationControlBuilder_->getValidParameters()));

    // Stepper
    ParameterList& stepperPL = 
      pl->sublist(stepper_name).disableRecursiveValidation();
    stepperPL.setParameters(*(stepperBuilder_->getValidParameters()));

    // Step Control
    ParameterList& stepControlPL = 
      pl->sublist(stepControl_name).disableRecursiveValidation();
    stepControlPL.setParameters(*(stepControlBuilder_->getValidParameters()));

    // Interpolation Buffer
    ParameterList& interpolationBufferPL = 
      pl->sublist(interpolationBuffer_name).disableRecursiveValidation();
    interpolationBufferPL.setParameters(*(interpolationBufferBuilder_->getValidParameters()));

    // Interpolation Buffer Appender
    ParameterList& interpolationBufferAppenderPL = 
      pl->sublist(interpolationBufferAppender_name).disableRecursiveValidation();
    interpolationBufferAppenderPL.setParameters(*(interpolationBufferAppenderBuilder_->getValidParameters()));

    // ErrWtVecCalc
    ParameterList& errWtVecCalcPL = 
      pl->sublist(errWtVecCalc_name).disableRecursiveValidation();
    errWtVecCalcPL.setParameters(*(errWtVecCalcBuilder_->getValidParameters()));

    // These parameters are necessary in order to set the stepper on the integrator
    pl->set(finalTime_name,Teuchos::as<Scalar>(finalTime_default));
    pl->set(landOnFinalTime_name,landOnFinalTime_default);

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
//
// a.  Its okay if the integration control comes back null, the 
//     IntegrationControlStrategyAcceptingIntegratorBase will deal with it
// b.  Its okay if the step control comes back null, the 
//     StepControlStrategyAcceptingStepperBase will deal with it
template<class Scalar>
RCP<IntegratorBase<Scalar> >
IntegratorBuilder<Scalar>::create(
    const RCP<const Thyra::ModelEvaluator<Scalar> > model,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver
    ) const
{
  TEST_FOR_EXCEPTION( is_null(model), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  The model passed in is null!"
      );
  TEST_FOR_EXCEPTION( is_null(paramList_), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  Please set a parameter list on this class before calling create."
      );
  // Create the integrator first
  RCP<ParameterList> integratorPL = sublist(paramList_,integrator_name);
  integratorBuilder_->setParameterList(integratorPL);
  RCP<IntegratorBase<Scalar> > integrator = integratorBuilder_->create();
  TEST_FOR_EXCEPTION( is_null(integrator), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  The integrator came back null from the ObjectBuilder!"
      );
  // Check for a trailing interpolation buffer and set it on the integrator
  RCP<TrailingInterpolationBufferAcceptingIntegratorBase<Scalar> > tibaIntegrator = 
    Teuchos::rcp_dynamic_cast<TrailingInterpolationBufferAcceptingIntegratorBase<Scalar> >(integrator,false);
  if (!is_null(tibaIntegrator)) {
    RCP<ParameterList> interpolationBufferPL = sublist(paramList_,interpolationBuffer_name);
    interpolationBufferBuilder_->setParameterList(interpolationBufferPL);
    RCP<InterpolationBufferBase<Scalar> > ib = interpolationBufferBuilder_->create();
    if (!is_null(ib)) {
      tibaIntegrator->setTrailingInterpolationBuffer(ib);
    }
  }
  // Check for an InterpolationBufferAppender and set it on the integrator
  RCP<InterpolationBufferAppenderAcceptingIntegratorBase<Scalar> > ibaaIntegrator = 
    Teuchos::rcp_dynamic_cast<InterpolationBufferAppenderAcceptingIntegratorBase<Scalar> >(integrator,false);
  if (!is_null(ibaaIntegrator)) {
    RCP<ParameterList> interpolationBufferAppenderPL = sublist(paramList_,interpolationBufferAppender_name);
    interpolationBufferAppenderBuilder_->setParameterList(interpolationBufferAppenderPL);
    RCP<InterpolationBufferAppenderBase<Scalar> > interpolationBufferAppender = interpolationBufferAppenderBuilder_->create();
    if (!is_null(interpolationBufferAppender)) {
      ibaaIntegrator->setInterpolationBufferAppender(interpolationBufferAppender);
    }
  }
  // Check for IntegrationControlStrategy and set it on the integrator
  RCP<IntegrationControlStrategyAcceptingIntegratorBase<Scalar> > icsaIntegrator = 
    Teuchos::rcp_dynamic_cast<IntegrationControlStrategyAcceptingIntegratorBase<Scalar> >(integrator,false);
  if (!is_null(icsaIntegrator)) {
    RCP<ParameterList> integrationControlPL = sublist(paramList_,integrationControl_name);
    integrationControlBuilder_->setParameterList(integrationControlPL);
    RCP<IntegrationControlStrategyBase<Scalar> > integrationControl = integrationControlBuilder_->create();
    if (!is_null(integrationControl)) {
      icsaIntegrator->setIntegrationControlStrategy(integrationControl);
    }
  }
  // Create the Stepper
  RCP<ParameterList> stepperPL = sublist(paramList_,stepper_name);
  stepperBuilder_->setParameterList(stepperPL);
  RCP<StepperBase<Scalar> > stepper = stepperBuilder_->create();
  TEST_FOR_EXCEPTION( is_null(stepper), std::logic_error,
      "Error!  IntegratorBuilder::create(...)  The stepper came back null from the StepperBuilder!"
      );
  // Create the Step Control
  RCP<StepControlStrategyAcceptingStepperBase<Scalar> > scsaStepper = 
    Teuchos::rcp_dynamic_cast<StepControlStrategyAcceptingStepperBase<Scalar> >(stepper,false);
  if (!is_null(scsaStepper)) {
    RCP<ParameterList> stepControlPL = sublist(paramList_,stepControl_name);
    stepControlBuilder_->setParameterList(stepControlPL);
    RCP<StepControlStrategyBase<Scalar> > stepControl = stepControlBuilder_->create();
    if (!is_null(stepControl)) {
      // Create the ErrWtVecCalc
      RCP<ErrWtVecCalcAcceptingStepControlStrategyBase<Scalar> > ewvcaStepControl = 
        Teuchos::rcp_dynamic_cast<ErrWtVecCalcAcceptingStepControlStrategyBase<Scalar> >(stepControl,false);
      if (!is_null(ewvcaStepControl)) {
        RCP<ParameterList> errWtVecCalcPL = sublist(paramList_,errWtVecCalc_name);
        errWtVecCalcBuilder_->setParameterList(errWtVecCalcPL);
        RCP<ErrWtVecCalcBase<Scalar> > errWtVecCalc = errWtVecCalcBuilder_->create();
        if (!is_null(errWtVecCalc)) {
          ewvcaStepControl->setErrWtVecCalc(errWtVecCalc);
        }
      }
      scsaStepper->setStepControlStrategy(stepControl);
    }
  }
  // Set model on stepper
  stepper->setModel(model);
  // Set initial condition on stepper
  stepper->setInitialCondition(initialCondition);
  // Set nonlinear solver on stepper
  RCP<SolverAcceptingStepperBase<Scalar> > saStepper = Teuchos::rcp_dynamic_cast<SolverAcceptingStepperBase<Scalar> >(stepper,false);
  if(!is_null(saStepper)) {
    TEST_FOR_EXCEPTION( is_null(nlSolver), std::logic_error,
        "Error!  IntegratorBuilder::create(...)  The nonlinear solver passed in is null and the stepper is implicit!"
        );
    saStepper->setSolver(nlSolver);
  }
  Scalar finalTime = paramList_->get<Scalar>(finalTime_name,Teuchos::as<Scalar>(finalTime_default));
  bool landOnFinalTime = paramList_->get<bool>(landOnFinalTime_name,landOnFinalTime_default);
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
      defaultIntegrator_name
      );

  // Integration Control Strategy
  integrationControlBuilder_ = Teuchos::objectBuilder<IntegrationControlStrategyBase<Scalar> >();
  integrationControlBuilder_->setObjectName(integrationControlBuilder_name);
  integrationControlBuilder_->setObjectTypeName(integrationControlBuilderType_name);
  integrationControlBuilder_->setObjectFactory(
      abstractFactoryStd< IntegrationControlStrategyBase<Scalar>, SimpleIntegrationControlStrategy<Scalar> >(),
      simpleIntegrationControl_name
      );

  // Stepper Builder
  stepperBuilder_ = stepperBuilder<Scalar>();

  // Step Control Strategy
  stepControlBuilder_ = Teuchos::objectBuilder<StepControlStrategyBase<Scalar> >();
  stepControlBuilder_->setObjectName(stepControlBuilder_name);
  stepControlBuilder_->setObjectTypeName(stepControlBuilderType_name);
  stepControlBuilder_->setObjectFactory(
      abstractFactoryStd< StepControlStrategyBase<Scalar>, ImplicitBDFStepperStepControl<Scalar> >(),
      implicitBDFStepControl_name
      );

  // Trailing Interpolation Buffer 
  interpolationBufferBuilder_ = Teuchos::objectBuilder<InterpolationBufferBase<Scalar> >();
  interpolationBufferBuilder_->setObjectName(interpolationBufferBuilder_name);
  interpolationBufferBuilder_->setObjectTypeName(interpolationBufferBuilderType_name);
  interpolationBufferBuilder_->setObjectFactory(
      abstractFactoryStd< InterpolationBufferBase<Scalar>, InterpolationBuffer<Scalar> >(),
      defaultInterpolationBuffer_name
      );

  // Interpolation Buffer Appender
  interpolationBufferAppenderBuilder_ = Teuchos::objectBuilder<InterpolationBufferAppenderBase<Scalar> >();
  interpolationBufferAppenderBuilder_->setObjectName(interpolationBufferAppenderBuilder_name);
  interpolationBufferAppenderBuilder_->setObjectTypeName(interpolationBufferAppenderBuilderType_name);
//  interpolationBufferAppenderBuilder_->setObjectFactory(
//      abstractFactoryStd< InterpolationBufferAppenderBase<Scalar>, SmartInterpolationBufferAppender<Scalar> >(),
//      smartInterpolationBufferAppender_name
//      );
  interpolationBufferAppenderBuilder_->setObjectFactory(
      abstractFactoryStd< InterpolationBufferAppenderBase<Scalar>, PointwiseInterpolationBufferAppender<Scalar> >(),
      pointwiseInterpolationBufferAppender_name
      );

  // ErrWtVecCalc
  errWtVecCalcBuilder_ = Teuchos::objectBuilder<ErrWtVecCalcBase<Scalar> >();
  errWtVecCalcBuilder_->setObjectName(errWtVecCalcBuilder_name);
  errWtVecCalcBuilder_->setObjectTypeName(errWtVecCalcBuilderType_name);
  errWtVecCalcBuilder_->setObjectFactory(
      abstractFactoryStd< ErrWtVecCalcBase<Scalar>, ImplicitBDFStepperErrWtVecCalc<Scalar> >(),
      implicitBDFStepperErrWtVecCalc_name
      );

}


} // namespace Rythmos

#endif //Rythmos_INTEGRATOR_BUILDER_H

//#endif // Rythmos_EXPERIMENTAL
