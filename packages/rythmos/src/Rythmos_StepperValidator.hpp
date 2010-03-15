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
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"

#include "Teuchos_StandardCatchMacros.hpp"


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
  void defaultInitializeAll_();
  RCP<StepperBase<Scalar> > getStepper_(
      const RCP<Thyra::ModelEvaluator<Scalar> >& model,
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& initialCondition,
      const RCP<Thyra::NonlinearSolverBase<Scalar> >& nlSolver
      ) const;
  bool isImplicitStepper_() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getSomeIC_(
    const Thyra::ModelEvaluator<Scalar>& model) const;

  // Validate that parameters set through setInitialCondition actually get set
  // on the model when evalModel is called.
  void validateIC_() const;
  // Validate that getTimeRange and getStepStatus return the correct states of
  // initialization.
  void validateStates_() const;
  // Validate that we can get the initial condition through getPoints after
  // setInitialCondition has been set and after the first step.
  void validateGetIC_() const;
  // Validate that we can get the initial condition through
  // getInitialCondition
  void validateGetIC2_() const;
  // Validate that the stepper supports getNodes, which is used by the
  // Trailing Interpolation Buffer feature of the Integrator.
  void validateGetNodes_() const;

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

  RCP<const std::vector<Thyra::ModelEvaluatorBase::InArgs<Scalar> > > getPassedInArgs() const;

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
  
  void defaultInitializeAll_();
  void initialize_();

  mutable RCP<std::vector<Thyra::ModelEvaluatorBase::InArgs<Scalar> > > passedInArgs_;
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
  this->defaultInitializeAll_();
  Np_ = 1;
  Ng_ = 0;
  np_ = 1;
  dim_ = 1;
  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  // Create p_space 
  p_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(np_);
  passedInArgs_ = rcp(new std::vector<Thyra::ModelEvaluatorBase::InArgs<Scalar> >);
  this->initialize_();
}

template<class Scalar>
void StepperValidatorMockModel<Scalar>::defaultInitializeAll_()
{
  passedInArgs_ = Teuchos::null;
  isImplicit_ = false;
  isInitialized_ = false;
  //inArgs_;
  //outArgs_;
  //nominalValues_;
  x_space_ = Teuchos::null;
  f_space_ = Teuchos::null;
  p_space_ = Teuchos::null;
  Np_ = -1;
  Ng_ = -1;
  np_ = -1;
  dim_ = -1;
}

template<class Scalar>
StepperValidatorMockModel<Scalar>::~StepperValidatorMockModel()
{ }

template<class Scalar>
void StepperValidatorMockModel<Scalar>::initialize_()
{
  if (!isInitialized_) {
    // Set up prototypical InArgs
    Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_t );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_beta );
#ifdef HAVE_THYRA_ME_POLYNOMIAL
    // For ExplicitTaylorPolynomialStepper
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x_poly );
#endif // HAVE_THYRA_ME_POLYNOMIAL
    inArgs.set_Np(1); 
    if (isImplicit_) {
      inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x_dot );
      inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_alpha );
    }
    inArgs_ = inArgs;
    // Set up prototypical OutArgs
    Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_f );
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_W_op );
#ifdef HAVE_THYRA_ME_POLYNOMIAL
    // For ExplicitTaylorPolynomialStepper
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_f_poly );
#endif // HAVE_THYRA_ME_POLYNOMIAL
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
    isInitialized_ = true;
  }
}


template<class Scalar>
RCP<const std::vector<Thyra::ModelEvaluatorBase::InArgs<Scalar> > > StepperValidatorMockModel<Scalar>::getPassedInArgs() const
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
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory = this->get_W_factory();
  RCP<Thyra::LinearOpBase<Scalar> > matrix = this->create_W_op();
  {
    // 01/20/09 tscoffe:  This is a total hack to provide a full rank matrix to
    // linearOpWithSolve because it ends up factoring the matrix during
    // initialization, which it really shouldn't do, or I'm doing something
    // wrong here.   The net effect is that I get exceptions thrown in
    // optimized mode due to the matrix being rank deficient unless I do this.
    RCP<Thyra::MultiVectorBase<Scalar> > multivec = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(matrix,true);
    {
      RCP<Thyra::VectorBase<Scalar> > vec = Thyra::createMember(x_space_);
      {
        Thyra::DetachedVectorView<Scalar> vec_view( *vec );
        vec_view[0] = 1.0;
      }
      V_V(&*(multivec->col(0)),*vec);
    }
  }
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > W = 
    Thyra::linearOpWithSolve<Scalar>(
      *W_factory,
      matrix
      );
  return W;
}


template<class Scalar>
RCP<Thyra::LinearOpBase<Scalar> > StepperValidatorMockModel<Scalar>::create_W_op() const
{
  RCP<Thyra::MultiVectorBase<Scalar> > matrix = Thyra::createMembers(x_space_, dim_);
  return(matrix);
}


template<class Scalar>
RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > StepperValidatorMockModel<Scalar>::get_W_factory() const
{
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
  return W_factory;
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
  if (isImplicit_ != test_isImplicit) {
    isImplicit_ = test_isImplicit;
    isInitialized_ = false;
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
  typedef Teuchos::ScalarTraits<Scalar> ST;
  passedInArgs_->push_back(inArgs);
  // Fill f with zeros.
  RCP<VectorBase<Scalar> > f_out = outArgs.get_f();
  if (!is_null(f_out)) {
    Thyra::V_S(Teuchos::outArg(*f_out),ST::zero());
  }
#ifdef HAVE_THYRA_ME_POLYNOMIAL
  if (outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_f_poly)) {
    RCP<Teuchos::Polynomial<VectorBase<Scalar> > > f_poly_out = outArgs.get_f_poly();
    if (!is_null(f_poly_out)) {
      //Thyra::V_S(Teuchos::outArg(*f_poly_out),ST::zero());
    }
  }
#endif // HAVE_THYRA_ME_POLYNOMIAL

}


//
// StepperValidator Definitions:
//

template<class Scalar>
StepperValidator<Scalar>::StepperValidator()
{ 
  this->defaultInitializeAll_();
}

template<class Scalar>
void StepperValidator<Scalar>::defaultInitializeAll_()
{
  integratorBuilder_ = Teuchos::null;
  paramList_ = Teuchos::null;
  stepperName_ = "";
}

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
  bool verbose = true;
  Array<bool> success_array;
  bool local_success = true;

  // Verify that the stepper passes parameters to the model in evalModel:
  local_success = true;
  try {
    this->validateIC_();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,local_success);
  success_array.push_back(local_success);

  // Verify that the stepper states are correct 
  //   uninitialized => getTimeRange == invalidTimeRange
  //   initialized, but no step => getTimeRange.length() == 0, [t_ic,t_ic]
  //   initialized, step taken => getTimeRange.length() > 0
  local_success = true;
  try {
    this->validateStates_();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,local_success);
  success_array.push_back(local_success);

  // Verify that getPoints(t_ic) returns the IC after initialization and after the first step
  local_success = true;
  try {
    this->validateGetIC_();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,local_success);
  success_array.push_back(local_success);

  // Verify that getInitialCondition() returns the IC after
  // setInitialCondition(...)
  local_success = true;
  try {
    this->validateGetIC2_();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,local_success);
  success_array.push_back(local_success);

  // Validate that the stepper supports getNodes, which is used by the Trailing Interpolation Buffer feature of the Integrator.
  local_success = true;
  try {
    this->validateGetNodes_();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,local_success);
  success_array.push_back(local_success);

  // Verify that getPoints(t) returns the same vectors as getStepStatus
  // TODO

  bool global_success = true;
  for (int i=0 ; i < Teuchos::as<int>(success_array.size()) ; ++i) {
    global_success = global_success && success_array[i];
  }

  TEST_FOR_EXCEPTION( !global_success, std::logic_error,
      "Error!  StepperValidator:  The stepper " << stepperName_ << " did not pass stepper validation."
      );
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
    const RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& initialCondition,
    const RCP<Thyra::NonlinearSolverBase<Scalar> >& nlSolver
    ) const
{
  RCP<IntegratorBase<Scalar> > integrator = integratorBuilder_->create(model,initialCondition,nlSolver);
  RCP<StepperBase<Scalar> > stepper = integrator->getNonconstStepper();
  return stepper;
}

template<class Scalar>
bool StepperValidator<Scalar>::isImplicitStepper_() const
{
  RCP<const StepperBuilder<Scalar> > sb = integratorBuilder_->getStepperBuilder();
  RCP<StepperBase<Scalar> > stepper = sb->create(stepperName_);
  return stepper->isImplicit();
}

template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> StepperValidator<Scalar>::getSomeIC_(
    const Thyra::ModelEvaluator<Scalar>& model 
    ) const
{
  // Set up some initial condition:
  Thyra::ModelEvaluatorBase::InArgs<Scalar> model_ic = model.createInArgs();
  if (model_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_t)) {
    Scalar time = Scalar(0.125);
    model_ic.set_t(time);
  }
  if (model_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_x)) {
    RCP<VectorBase<Scalar> > x_ic = Thyra::createMember(*(model.get_x_space()));
    {
      Thyra::DetachedVectorView<Scalar> x_ic_view( *x_ic );
      x_ic_view[0] = Scalar(10.0);
    }
    model_ic.set_x(x_ic);
  }
  for (int i=0 ; i<model_ic.Np() ; ++i) {
    RCP<VectorBase<Scalar> > p_ic = Thyra::createMember(*(model.get_p_space(i)));
    {
      Thyra::DetachedVectorView<Scalar> p_ic_view( *p_ic );
      p_ic_view[i] = Scalar(11.0+i);
    }
    model_ic.set_p(i,p_ic);
  }
  if (model_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot)) {
    RCP<VectorBase<Scalar> > xdot_ic = Thyra::createMember(*(model.get_x_space()));
    {
      Thyra::DetachedVectorView<Scalar> xdot_ic_view( *xdot_ic );
      xdot_ic_view[0] = Scalar(12.0);
    }
    model_ic.set_x_dot(xdot_ic);
  }
  return model_ic;
}


template<class Scalar>
void StepperValidator<Scalar>::validateIC_() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // Determine if the stepper is implicit or not:
  bool isImplicit = this->isImplicitStepper_();
  RCP<StepperValidatorMockModel<Scalar> > model = 
    stepperValidatorMockModel<Scalar>(isImplicit);
  // Set up some initial condition:
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stepper_ic = this->getSomeIC_(*model);
  // Create nonlinear solver (if needed)
  RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver;
  if (isImplicit) {
    nlSolver = Rythmos::timeStepNonlinearSolver<Scalar>();
  }
  RCP<StepperBase<Scalar> > stepper = this->getStepper_(model,stepper_ic,nlSolver);
  Scalar dt = Scalar(0.1);
  Scalar dt_taken = ST::nan();
  dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
  // Verify that the parameter got set on the model by asking the model for the
  // inArgs passed to it through evalModel
  RCP<const std::vector<Thyra::ModelEvaluatorBase::InArgs<Scalar> > >
    passedInArgs_ptr = model->getPassedInArgs();
  const std::vector<Thyra::ModelEvaluatorBase::InArgs<Scalar> >&
    passedInArgs = *passedInArgs_ptr;
  bool valid_t = false;
  // Technically this will fail for any Runge Kutta Butcher Tableau where:
  //  c(0) != 0 and c(s) != 1, where s = # of stages.
  if ( (passedInArgs[0].get_t() == stepper_ic.get_t()   )   || 
       (passedInArgs[0].get_t() == stepper_ic.get_t()+dt) )
  {
    valid_t = true;
  }
  TEST_FOR_EXCEPTION( !valid_t, std::logic_error,
    "Error!  StepperValidator::validateIC:  Time did not get correctly set on"
    " the model through StepperBase::setInitialCondition!"
    );
  // If a stepper uses a predictor, then the x passed into the model will not
  // be the same as the IC, so we can't check it.
  RCP<const VectorBase<Scalar> > p_out = passedInArgs[0].get_p(0);
  TEST_FOR_EXCEPTION( is_null(p_out), std::logic_error,
    "Error!  StepperValidator::validateIC:  Parameter 0 did not get set on the"
    " model through StepperBase::setInitialCondition!"
    );
  {
    Thyra::ConstDetachedVectorView<Scalar> p_out_view( *p_out );
    TEST_FOR_EXCEPTION( p_out_view[0] != Scalar(11.0), std::logic_error,
        "Error!  StepperValidator::validateIC:  Parameter 0 did not get set correctly on the model through StepperBase::setInitialCondition!"
        );
  }
  if (isImplicit) {
    // Each time integration method will approximate xdot according to its own
    // algorithm, so we can't really test it here.
  }
}

template<class Scalar>
void StepperValidator<Scalar>::validateStates_() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const StepperBuilder<Scalar> > sb = integratorBuilder_->getStepperBuilder();
  RCP<StepperBase<Scalar> > stepper = sb->create(stepperName_);
  {
    // Uninitialized Stepper:
    TimeRange<Scalar> tr = stepper->getTimeRange();
    TEST_FOR_EXCEPTION( tr.isValid(), std::logic_error,
        "Error!  StepperValidator::validateStates:  Uninitialized stepper returned a valid time range!"
        );
    const StepStatus<Scalar> ss = stepper->getStepStatus();
    TEST_FOR_EXCEPTION( ss.stepStatus != STEP_STATUS_UNINITIALIZED, std::logic_error,
        "Error!  StepperValidator::validateStates:  Uninitialized stepper returned a valid step status!"
        );
  }
  bool isImplicit = stepper->isImplicit();
  RCP<StepperValidatorMockModel<Scalar> > model = 
    stepperValidatorMockModel<Scalar>(isImplicit);
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stepper_ic = this->getSomeIC_(*model);
  stepper->setInitialCondition(stepper_ic);
  {
    // Has initial condition:
    TimeRange<Scalar> tr = stepper->getTimeRange();
    TEST_FOR_EXCEPTION( compareTimeValues(tr.lower(),tr.upper()) != 0, std::logic_error,
        "Error!  StepperValidator::validateStates:  Stepper with initial condition returned a non zero time range!"
        );
//    const StepStatus<Scalar> ss = stepper->getStepStatus();
//    TEST_FOR_EXCEPTION( ss.stepStatus != STEP_STATUS_UNKNOWN, std::logic_error,
//        "Error!  StepperValidator::validateStates:  Stepper with initial condition did not return STEP_STATUS_UNKNOWN!"
//        );
  }
  // 04/16/09 tscoffe:  Now we use the integratorBuilder to create a fully
  // initialized stepper which we can use for taking a step.  We can't just
  // continue setting them up because we don't know what they all require.
  RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver;
  if (isImplicit) {
    nlSolver = timeStepNonlinearSolver<Scalar>();
  }
  stepper = this->getStepper_(model,stepper_ic,nlSolver);
  {
    // Still has initial condition:
    TimeRange<Scalar> tr = stepper->getTimeRange();
    TEST_FOR_EXCEPTION( compareTimeValues(tr.lower(),tr.upper()) != 0, std::logic_error,
        "Error!  StepperValidator::validateStates:  Fully initialized stepper returned a non zero time range!"
        );
//    const StepStatus<Scalar> ss = stepper->getStepStatus();
//    TEST_FOR_EXCEPTION( ss.stepStatus != STEP_STATUS_UNKNOWN, std::logic_error,
//        "Error!  StepperValidator::validateStates:  Fully initialized stepper, prior to taking a step, did not return STEP_STATUS_UNKNOWN!"
//        );
  }
  Scalar dt = Scalar(0.1);
  Scalar dt_taken = ST::nan();
  dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
  {
    // Taken Step:
    TimeRange<Scalar> tr = stepper->getTimeRange();
    TEST_FOR_EXCEPTION( compareTimeValues(tr.lower(),tr.upper()) >= 0, std::logic_error,
        "Error!  StepperValidator::validateStates:  Stepper returned a zero (or invalid) time range after taking a step!"
        );
    const StepStatus<Scalar> ss = stepper->getStepStatus();
    TEST_FOR_EXCEPTION( ss.stepStatus != STEP_STATUS_CONVERGED, std::logic_error,
        "Error!  StepperValidator::validateStates:  Stepper did not return converged step status after taking a step!"
        );
  }
}

template<class Scalar>
void StepperValidator<Scalar>::validateGetIC_() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  // Determine if the stepper is implicit or not:
  bool isImplicit = this->isImplicitStepper_();
  RCP<StepperValidatorMockModel<Scalar> > model = 
    stepperValidatorMockModel<Scalar>(isImplicit);
  // Set up some initial condition:
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stepper_ic = this->getSomeIC_(*model);
  // Create nonlinear solver (if needed)
  RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver;
  if (isImplicit) {
    nlSolver = Rythmos::timeStepNonlinearSolver<Scalar>();
  }
  RCP<StepperBase<Scalar> > stepper = this->getStepper_(model,stepper_ic,nlSolver);
  // Verify we can get the IC through getPoints prior to taking a step:
  {
    Array<Scalar> time_vec;
    Array<RCP<const VectorBase<Scalar> > > x_vec;
    Array<RCP<const VectorBase<Scalar> > > xdot_vec;
    Array<ScalarMag> accuracy_vec;
    time_vec.push_back(stepper_ic.get_t());
    stepper->getPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
    {
      Thyra::ConstDetachedVectorView<Scalar> x_view( *x_vec[0] );
      TEST_FOR_EXCEPTION( compareTimeValues<Scalar>(x_view[0],Scalar(10.0)) != 0,
        std::logic_error,
        "Error!  StepperValidator::validateGetIC:  Stepper did not return the initial"
        " condition for X through getPoints prior to taking a step!"
        );
    }
    if (isImplicit && !is_null(xdot_vec[0])) {
      Thyra::ConstDetachedVectorView<Scalar> xdot_view( *xdot_vec[0] );
      TEST_FOR_EXCEPTION( compareTimeValues<Scalar>(xdot_view[0],Scalar(12.0)) != 0,
        std::logic_error,
        "Error!  StepperValidator::validateGetIC:  Stepper did not return the initial"
        " condition for XDOT through getPoints prior to taking a step!"
        );
    }
  }
  // Verify we can get the IC through getPoints after taking a step:
  {
    Scalar dt = Scalar(0.1);
    Scalar dt_taken = ST::nan();
    dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
    Array<Scalar> time_vec;
    Array<RCP<const VectorBase<Scalar> > > x_vec;
    Array<RCP<const VectorBase<Scalar> > > xdot_vec;
    Array<ScalarMag> accuracy_vec;
    time_vec.push_back(stepper_ic.get_t());
    stepper->getPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
    {
      Thyra::ConstDetachedVectorView<Scalar> x_view( *x_vec[0] );
      TEST_FOR_EXCEPTION( compareTimeValues<Scalar>(x_view[0],Scalar(10.0)) != 0,
        std::logic_error,
        "Error!  StepperValidator::validateGetIC:  Stepper did not return the initial"
        " condition for X through getPoints after taking a step!"
        );
    }
    if (isImplicit && !is_null(xdot_vec[0])) {
      Thyra::ConstDetachedVectorView<Scalar> xdot_view( *xdot_vec[0] );
      TEST_FOR_EXCEPTION( compareTimeValues<Scalar>(xdot_view[0],Scalar(12.0)) != 0,
        std::logic_error,
        "Error!  StepperValidator::validateGetIC:  Stepper did not return the initial"
        " condition for XDOT through getPoints after taking a step!"
        );
    }
  }
}


template<class Scalar>
void StepperValidator<Scalar>::validateGetIC2_() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  // Determine if the stepper is implicit or not:
  bool isImplicit = this->isImplicitStepper_();
  RCP<StepperValidatorMockModel<Scalar> > model = 
    stepperValidatorMockModel<Scalar>(isImplicit);
  // Set up some initial condition:
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stepper_ic = this->getSomeIC_(*model);
  // Create nonlinear solver (if needed)
  RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver;
  if (isImplicit) {
    nlSolver = Rythmos::timeStepNonlinearSolver<Scalar>();
  }
  RCP<StepperBase<Scalar> > stepper = this->getStepper_(model,stepper_ic,nlSolver);
  // Verify we can get the IC back through getInitialCondition:
  Thyra::ModelEvaluatorBase::InArgs<Scalar> new_ic = stepper->getInitialCondition();
  //TEUCHOS_ASSERT( new_ic == stepper_ic );
  // Verify new_ic == stepper_ic
  {
    TEUCHOS_ASSERT( new_ic.get_t() == stepper_ic.get_t() );
    TEUCHOS_ASSERT( new_ic.get_x() == stepper_ic.get_x() );
    for (int i=0 ; i<stepper_ic.Np() ; ++i) {
      TEUCHOS_ASSERT( new_ic.get_p(i) == stepper_ic.get_p(i) );
    }
    if (isImplicit) {
      TEUCHOS_ASSERT( new_ic.get_x_dot() == stepper_ic.get_x_dot() );
    }
  }
}


template<class Scalar>
void StepperValidator<Scalar>::validateGetNodes_() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // Create uninitialized stepper and verify we get no nodes back
  {
    RCP<const StepperBuilder<Scalar> > sb = integratorBuilder_->getStepperBuilder();
    RCP<StepperBase<Scalar> > stepper = sb->create(stepperName_);
    Array<Scalar> nodes;
    stepper->getNodes(&nodes);
    TEST_FOR_EXCEPTION( nodes.size() != 0, std::logic_error,
        "Error!  StepperValidator::validateGetNodes:  Uninitialized stepper returned non-empty node list!"
        );
  }
  // Create fully initialize stepper and verify we get back one node for IC
  bool isImplicit = this->isImplicitStepper_();
  RCP<StepperValidatorMockModel<Scalar> > model = 
    stepperValidatorMockModel<Scalar>(isImplicit);
  Thyra::ModelEvaluatorBase::InArgs<Scalar> stepper_ic = this->getSomeIC_(*model);
  RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver;
  if (isImplicit) {
    nlSolver = Rythmos::timeStepNonlinearSolver<Scalar>();
  }
  RCP<StepperBase<Scalar> > stepper = this->getStepper_(model,stepper_ic,nlSolver);
  {
    Array<Scalar> nodes;
    stepper->getNodes(&nodes);
    TEST_FOR_EXCEPTION( nodes.size() == 0, std::logic_error,
        "Error!  StepperValidator::validateGetNodes:  Initialized stepper returned empty node list!"
        );
    TEST_FOR_EXCEPTION( nodes.size() > 1, std::logic_error,
        "Error!  StepperValidator::validateGetNodes:  Initialized stepper returned node list with more than one node!"
        );
  }
  // Take a step with the stepper and verify we get back two nodes
  Scalar dt = Scalar(0.1);
  Scalar dt_taken = ST::nan();
  dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED);
  {
    Array<Scalar> nodes;
    stepper->getNodes(&nodes);
    TEST_FOR_EXCEPTION( nodes.size() == 0, std::logic_error,
        "Error!  StepperValidator::validateGetNodes:  After taking a step, stepper returned empty node list!"
        );
    TEST_FOR_EXCEPTION( nodes.size() == 1, std::logic_error,
        "Error!  StepperValidator::validateGetNodes:  After taking a step, stepper returned node list with only one node!"
        );
    TEST_FOR_EXCEPTION( nodes.size() > 2, std::logic_error,
        "Error!  StepperValidator::validateGetNodes:  After taking a step, stepper returned node list with more than two nodes!"
        );
  }
}

} // namespace Rythmos

#endif //Rythmos_STEPPER_VALIDATOR_H
