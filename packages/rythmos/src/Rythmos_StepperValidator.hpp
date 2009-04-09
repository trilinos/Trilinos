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

#ifndef Rythmos_STEPPER_VALIDATOR_H
#define Rythmos_STEPPER_VALIDATOR_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Rythmos_IntegratorBuilder.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp" 
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Thyra_ModelEvaluator.hpp" 

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_Types.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"



namespace Rythmos {


/** \brief Class for validating steppers
 *
 * There are a variety of requirements placed on the steppers by the higher
 * level objects like the integrator.  These requirements are documented and
 * tested in this class.
 * 
 */
template<class Scalar> 
class StepperValidator 
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<StepperValidator<Scalar> >
{
public:

  StepperValidator();

  virtual ~StepperValidator();

  /** \brief Set a Integrator Builder object. */
  void setIntegratorBuilder(
    const RCP<IntegratorBuilder<Scalar> > &integratorBuilder
    );

  /** \brief Validate the stepper by the StepperBuilder. */
  void validateStepper() const;

  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);

  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();

  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();

  /** \brief . */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

private:

  // Private member functions:
  RCP<StepperBase<Scalar> > getStepper_(
      const RCP<const Thyra::ModelEvaluator<Scalar> > model,
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition,
      const RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver
      ) const;

  void validateIC_() const;

  // Private member data:
  RCP<IntegratorBuilder<Scalar> > integratorBuilder_;
  RCP<ParameterList> paramList_;
  mutable std::string stepperName_;

};

//
// StepperValidator nonmember constructor:
//
template<class Scalar>
RCP<StepperValidator<Scalar> > stepperValidator()
{
  return rcp(new StepperValidator<Scalar>() );
}

//
//  Mock Model class for validating a stepper
//
template<class Scalar>
class StepperValidatorMockModel
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  StepperValidatorMockModel();

  virtual ~StepperValidatorMockModel();

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getPassedInArgs() const;

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;

  //@}
  
  /** \name Public functions overridden from ParameterListAcceptor. */
  //@{
  
  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}
  
private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar
    ) const;

  //@}
  
  void initialize_();

  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> passedInArgs_;
  bool isImplicit_;
  bool isInitialized_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  int Np_; // number of parameter vectors (1)
  int Ng_; // number of observation functions (0)
  int np_; // length of parameter vector (1)
  int dim_; // length of solution vector

};

//
// StepperValidatorMockModel nonmember constructors:
//
template<class Scalar>
RCP<StepperValidatorMockModel<Scalar> > stepperValidatorMockModel() 
{
  return(rcp(new StepperValidatorMockModel<Scalar>()));
}

template<class Scalar>
RCP<StepperValidatorMockModel<Scalar> > stepperValidatorMockModel(
    bool isImplicit
    ) 
{
  RCP<StepperValidatorMockModel<Scalar> > model = rcp(new StepperValidatorMockModel<Scalar>());
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Is Implicit",isImplicit);
  model->setParameterList(pl);
  return model;
}

//
// StepperValidatorMockModel Definitions:
//
template<class Scalar>
StepperValidatorMockModel<Scalar>::StepperValidatorMockModel()
{ 
  isImplicit_ = false;
  Np_ = 1;
  Ng_ = 0;
  np_ = 1;
  dim_ = 1;
  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  // Create p_space 
  p_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(np_);
  this->initialize_();
}

template<class Scalar>
StepperValidatorMockModel<Scalar>::~StepperValidatorMockModel()
{ }

template<class Scalar>
void StepperValidatorMockModel<Scalar>::initialize_()
{
  if (!isInitialized_) {
    if (isImplicit_) {
      TEST_FOR_EXCEPT(true);
    } 
    else {
      // Set up prototypical InArgs
      Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
      inArgs.setModelEvalDescription(this->description());
      inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_t );
      inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x );
      inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_beta );
      inArgs.set_Np(1); 
      inArgs_ = inArgs;
      // Set up prototypical OutArgs
      Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
      outArgs.setModelEvalDescription(this->description());
      outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_f );
      outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_W_op );
      outArgs.set_Np_Ng(Np_,Ng_);
      outArgs_ = outArgs;
      // Set up nominal values
      nominalValues_ = inArgs_;
      nominalValues_.set_t(Scalar(1.0));
      const RCP<VectorBase<Scalar> > x_ic = Thyra::createMember(x_space_);
      Thyra::V_S(Teuchos::outArg(*x_ic),Scalar(2.0));
      nominalValues_.set_x(x_ic);
      const RCP<VectorBase<Scalar> > p_ic = Thyra::createMember(p_space_);
      Thyra::V_S(Teuchos::outArg(*p_ic),Scalar(3.0));
      nominalValues_.set_p(0,p_ic);
    }
    isInitialized_ = true;
  }
}

template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> StepperValidatorMockModel<Scalar>::getPassedInArgs() const
{
  return passedInArgs_;
}
template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> > StepperValidatorMockModel<Scalar>::get_x_space() const
{
  TEST_FOR_EXCEPT(!isInitialized_);
  return x_space_;
}
template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> > StepperValidatorMockModel<Scalar>::get_f_space() const
{
  TEST_FOR_EXCEPT(!isInitialized_);
  return f_space_;
}
template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> StepperValidatorMockModel<Scalar>::getNominalValues() const
{
  TEST_FOR_EXCEPT(!isInitialized_);
  return nominalValues_;
}
template<class Scalar>
RCP<Thyra::LinearOpWithSolveBase<Scalar> > StepperValidatorMockModel<Scalar>::create_W() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}
template<class Scalar>
RCP<Thyra::LinearOpBase<Scalar> > StepperValidatorMockModel<Scalar>::create_W_op() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}
template<class Scalar>
RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > StepperValidatorMockModel<Scalar>::get_W_factory() const
{
  return Teuchos::null;
}
template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> StepperValidatorMockModel<Scalar>::createInArgs() const
{
  TEST_FOR_EXCEPT(!isInitialized_);
  return inArgs_;
}
template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> > StepperValidatorMockModel<Scalar>::get_p_space(int l) const
{
  TEST_FOR_EXCEPT(!isInitialized_);
  return p_space_;
}
template<class Scalar>
RCP<const Teuchos::Array<std::string> > StepperValidatorMockModel<Scalar>::get_p_names(int l) const
{
  RCP<Teuchos::Array<std::string> > p_strings = 
    Teuchos::rcp(new Teuchos::Array<std::string>());
  p_strings->push_back("Model Coefficient");
  return p_strings;
}
template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> > StepperValidatorMockModel<Scalar>::get_g_space(int j) const
{
  return Teuchos::null;
}
template<class Scalar>
void StepperValidatorMockModel<Scalar>::setParameterList(RCP<ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT( is_null(paramList) );
  paramList->validateParameters(*this->getValidParameters());
  Teuchos::readVerboseObjectSublist(&*paramList,this);
  bool test_isImplicit = paramList->get("Is Implicit",false);
  if (test_isImplicit != isImplicit_) {
    isImplicit_ = test_isImplicit;
    this->initialize_();
  }
  setMyParamList(paramList);
}
template<class Scalar>
RCP<const ParameterList> StepperValidatorMockModel<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Is Implicit",false);
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}
template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> StepperValidatorMockModel<Scalar>::createOutArgsImpl() const
{
  TEST_FOR_EXCEPT(!isInitialized_);
  return outArgs_;
}

template<class Scalar>
void StepperValidatorMockModel<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  passedInArgs_ = inArgs;
}


//
// StepperValidator Definitions:
//

template<class Scalar>
StepperValidator<Scalar>::StepperValidator()
{ }

template<class Scalar>
  StepperValidator<Scalar>::~StepperValidator()
{ }

template<class Scalar>
  void StepperValidator<Scalar>::setIntegratorBuilder(
    const RCP<IntegratorBuilder<Scalar> > &integratorBuilder
    )
{
  TEUCHOS_ASSERT( !is_null(integratorBuilder) );
  integratorBuilder_ = integratorBuilder;
}

template<class Scalar>
  void StepperValidator<Scalar>::validateStepper() const
{
  // Extract the name of the stepper for later
  {
    RCP<const ParameterList> pl = integratorBuilder_->getParameterList();
    if (is_null(pl)) {
      pl = integratorBuilder_->getValidParameters();
    }
    const Teuchos::ParameterList& stepperSelectionPL = pl->sublist("Stepper Settings").sublist("Stepper Selection");
    stepperName_ = stepperSelectionPL.get<std::string>("Stepper Type");
  }
  // Verify that the stepper passes parameters to the model in evalModel:
  this->validateIC_();
  // Verify that the stepper states are correct 
  //   uninitialized => getTimeRange == invalidTimeRange
  //   initialized, but no step => getTimeRange.length() == 0, [t_ic,t_ic]
  //   initialized, step taken => getTimeRange.length() > 0
  // Verify that getPoints(t_ic) returns the IC after initialization
  // Verify that getPoints(t_ic) returns the IC correctly after the first step
}

template<class Scalar>
  void StepperValidator<Scalar>::setParameterList(RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
  RCP<Teuchos::ParameterList> StepperValidator<Scalar>::getNonconstParameterList()
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
  RCP<Teuchos::ParameterList> StepperValidator<Scalar>::unsetParameterList()
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
  RCP<const Teuchos::ParameterList> StepperValidator<Scalar>::getValidParameters() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}

template<class Scalar>
  std::string StepperValidator<Scalar>::description() const
{
  TEST_FOR_EXCEPT(true);
  return("");
}

template<class Scalar>
  void StepperValidator<Scalar>::describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
RCP<StepperBase<Scalar> > StepperValidator<Scalar>::getStepper_(
    const RCP<const Thyra::ModelEvaluator<Scalar> > model,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition,
    const RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver
    ) const
{
  RCP<IntegratorBase<Scalar> > integrator = integratorBuilder_->create(model,initialCondition,nlSolver);
  RCP<StepperBase<Scalar> > stepper = integrator->getNonconstStepper();
  return stepper;
}


template<class Scalar>
void StepperValidator<Scalar>::validateIC_() const
{
  // Determine if the stepper is implicit or not:
  bool isImplicit;
  {
    RCP<const StepperBuilder<Scalar> > sb = integratorBuilder_->getStepperBuilder();
    RCP<StepperBase<Scalar> > stepper = sb->create(stepperName_);
    isImplicit = stepper->isImplicit();
  }
  RCP<StepperValidatorMockModel<Scalar> > model = 
    stepperValidatorMockModel<Scalar>(isImplicit);
  // Set up initial condition:
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stepper_ic = model->createInArgs();
  {
    // Explicit InArgs items:
    // time
    Scalar time = Scalar(0.1);
    stepper_ic.set_t(time);
    // x
    RCP<VectorBase<Scalar> > x_ic = Thyra::createMember(*(model->get_x_space()));
    {
      Thyra::DetachedVectorView<Scalar> x_ic_view( *x_ic );
      x_ic_view[0] = Scalar(10.0);
    }
    stepper_ic.set_x(x_ic);
    // parameter 
    RCP<VectorBase<Scalar> > p_ic = Thyra::createMember(*(model->get_p_space(0)));
    {
      Thyra::DetachedVectorView<Scalar> p_ic_view( *p_ic );
      p_ic_view[0] = Scalar(11.0);
    }
    stepper_ic.set_p(0,p_ic);
  }
  if (isImplicit) {
    TEST_FOR_EXCEPT(true);
    // Set x_dot, alpha, beta, etc. on stepper_ic
  }
  RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver;
  if (isImplicit) {
    TEST_FOR_EXCEPT(true);
    // Set up Nonlinear Solver
  }
  RCP<StepperBase<Scalar> > stepper = this->getStepper_(model,stepper_ic,nlSolver);
  Scalar dt = Scalar(0.1);
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Scalar dt_taken = ST::nan();
  dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
  // Verify that the parameter got set on the model by asking the model for the
  // inArgs passed to it through evalModel
  Thyra::ModelEvaluatorBase::InArgs<Scalar> passedInArgs = model->getPassedInArgs();
  TEST_FOR_EXCEPTION( passedInArgs.get_t() != stepper_ic.get_t(), std::logic_error,
      "Error!  StepperValidator::validateIC:  Time did not get correctly set on the model through StepperBase::setInitialCondition!"
      );
  RCP<const VectorBase<Scalar> > x_out = passedInArgs.get_x();
  {
    Thyra::ConstDetachedVectorView<Scalar> x_out_view( *x_out );
    TEST_FOR_EXCEPTION( x_out_view[0] != Scalar(10.0), std::logic_error,
        "Error!  StepperValidator::validateIC:  X did not get set correctly on the model through StepperBase::setInitialCondition!"
        );
  }
  RCP<const VectorBase<Scalar> > p_out = passedInArgs.get_p(0);
  {
    Thyra::ConstDetachedVectorView<Scalar> p_out_view( *p_out );
    TEST_FOR_EXCEPTION( p_out_view[0] != Scalar(11.0), std::logic_error,
        "Error!  StepperValidator::validateIC:  Parameter 0 did not get set correctly on the model through StepperBase::setInitialCondition!"
        );
  }
  if (isImplicit) {
    TEST_FOR_EXCEPT(true);
    // Test that x_dot, alpha, beta got passed to the model
  }
}

} // namespace Rythmos

#endif //Rythmos_STEPPER_VALIDATOR_H
