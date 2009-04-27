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

#ifndef Rythmos_FORWARDEULER_STEPPER_DEF_H
#define Rythmos_FORWARDEULER_STEPPER_DEF_H

#include "Rythmos_ForwardEulerStepper_decl.hpp"
#include "Rythmos_StepperHelpers.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"

namespace Rythmos {

// Nonmember constructor
template<class Scalar>
RCP<ForwardEulerStepper<Scalar> > forwardEulerStepper()
{
  RCP<ForwardEulerStepper<Scalar> > stepper = rcp(new ForwardEulerStepper<Scalar>());
  return stepper;
}

// Nonmember constructor
template<class Scalar>
RCP<ForwardEulerStepper<Scalar> > forwardEulerStepper(const RCP<const Thyra::ModelEvaluator<Scalar> > &model)
{
  RCP<ForwardEulerStepper<Scalar> > stepper = forwardEulerStepper<Scalar>();
  stepper->setModel(model);
  return stepper;
}

template<class Scalar>
ForwardEulerStepper<Scalar>::ForwardEulerStepper()
{
  this->defaultInitializAll_();
  typedef Teuchos::ScalarTraits<Scalar> ST;
  dt_ = ST::zero();
  numSteps_ = 0;
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::defaultInitializAll_()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = Teuchos::null;
  solution_vector_ = Teuchos::null;
  residual_vector_ = Teuchos::null;
  t_ = ST::nan();
  dt_ = ST::nan();
  t_old_ = ST::nan();
  solution_vector_old_ = Teuchos::null;
  //basePoint_;
  numSteps_ = -1;
  haveInitialCondition_ = false;
  parameterList_ = Teuchos::null;
  isInitialized_ = false;
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::initialize_()
{
  if (!isInitialized_) {
    TEST_FOR_EXCEPTION( is_null(model_), std::logic_error,
       "Error!  Please set a model on the stepper.\n" 
       );
    residual_vector_ = Thyra::createMember(model_->get_f_space());
    isInitialized_ = true;
  }
}

template<class Scalar>
ForwardEulerStepper<Scalar>::~ForwardEulerStepper()
{
}

template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> > ForwardEulerStepper<Scalar>::get_x_space() const
{
  TEST_FOR_EXCEPTION(!haveInitialCondition_,std::logic_error,"Error, attempting to call get_x_space before setting an initial condition!\n");
  return(solution_vector_->space());
}

template<class Scalar>
Scalar ForwardEulerStepper<Scalar>::takeStep(Scalar dt, StepSizeType flag)
{
  TEST_FOR_EXCEPTION( !haveInitialCondition_, std::logic_error,
     "Error!  Attempting to call takeStep before setting an initial condition!\n" 
     );
  this->initialize_();
  if (flag == STEP_TYPE_VARIABLE) { 
    // print something out about this method not supporting automatic variable step-size
    typedef Teuchos::ScalarTraits<Scalar> ST;
    return(-ST::one());
  }
  //Thyra::eval_f<Scalar>(*model_,*solution_vector_,t_+dt,&*residual_vector_);
  eval_model_explicit<Scalar>(*model_,basePoint_,*solution_vector_,t_+dt,Teuchos::outArg(*residual_vector_));

  // solution_vector_old_ = solution_vector_
  Thyra::V_V(Teuchos::outArg(*solution_vector_old_),*solution_vector_);
  // solution_vector = solution_vector + dt*residual_vector
  Thyra::Vp_StV(&*solution_vector_,dt,*residual_vector_); 
  t_old_ = t_;
  t_ += dt;
  dt_ = dt;
  numSteps_++;

  return(dt);
}

template<class Scalar>
const StepStatus<Scalar> ForwardEulerStepper<Scalar>::getStepStatus() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  StepStatus<Scalar> stepStatus;

  if (!haveInitialCondition_) {
    stepStatus.stepStatus = STEP_STATUS_UNINITIALIZED;
  } 
  else if (numSteps_ == 0) {
    stepStatus.stepStatus = STEP_STATUS_UNKNOWN;
    stepStatus.order = this->getOrder();
    stepStatus.time = t_;
    stepStatus.solution = solution_vector_;
  } 
  else {
    stepStatus.stepStatus = STEP_STATUS_CONVERGED;
    stepStatus.stepSize = dt_; 
    stepStatus.order = this->getOrder();
    stepStatus.time = t_;
    stepStatus.stepLETValue = Scalar(-ST::one()); 
    stepStatus.solution = solution_vector_;
    stepStatus.residual = residual_vector_;
  }

  return(stepStatus);
}

template<class Scalar>
std::string ForwardEulerStepper<Scalar>::description() const
{
  std::string name = "Rythmos::ForwardEulerStepper";
  return(name);
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::describe(
      Teuchos::FancyOStream            &out,
      const Teuchos::EVerbosityLevel   verbLevel
      ) const
{
  if ( (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_DEFAULT) ) ||
       (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)     )
     ) {
    out << description() << "::describe" << std::endl;
    out << "model = " << model_->description() << std::endl;
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)) {
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM)) {
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH)) {
    out << "model = " << std::endl;
    model_->describe(out,verbLevel);
    out << "solution_vector = " << std::endl;
    solution_vector_->describe(out,verbLevel);
    out << "residual_vector = " << std::endl;
    residual_vector_->describe(out,verbLevel);
  }
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::addPoints(
    const Array<Scalar>& time_vec
    ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec
    ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    )
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, addPoints is not implemented for ForwardEulerStepper.\n");
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::getPoints(
    const Array<Scalar>& time_vec
    ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec
    ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
    ,Array<ScalarMag>* accuracy_vec) const
{
  TEUCHOS_ASSERT( haveInitialCondition_ );
  using Teuchos::constOptInArg;
  using Teuchos::null;
  defaultGetPoints<Scalar>(
      t_old_,constOptInArg(*solution_vector_old_),null,
      t_,constOptInArg(*solution_vector_),null,
      time_vec,ptr(x_vec),ptr(xdot_vec),ptr(accuracy_vec),
      null
      );
}

template<class Scalar>
TimeRange<Scalar> ForwardEulerStepper<Scalar>::getTimeRange() const
{
  if (!haveInitialCondition_) {
    return(invalidTimeRange<Scalar>());
  } else {
    return(TimeRange<Scalar>(t_old_,t_));
  }
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  TEUCHOS_ASSERT( time_vec != NULL );
  time_vec->clear();
  if (!haveInitialCondition_) {
    return; 
  } else {
    time_vec->push_back(t_old_);
  }
  if (numSteps_ > 0) {
    time_vec->push_back(t_);
  }
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, removeNodes is not implemented for ForwardEulerStepper.\n");
}

template<class Scalar>
int ForwardEulerStepper<Scalar>::getOrder() const
{
  return(1);
}

template <class Scalar>
void ForwardEulerStepper<Scalar>::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> ForwardEulerStepper<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> ForwardEulerStepper<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
ForwardEulerStepper<Scalar>::getValidParameters() const
{
  using Teuchos::ParameterList;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model)
{
  TEST_FOR_EXCEPT( is_null(model) );
  assertValidModel( *this, *model );
  model_ = model;
}

template<class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
ForwardEulerStepper<Scalar>::getModel() const
{
  return model_;
}

template<class Scalar>
void ForwardEulerStepper<Scalar>::setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    )
{
  basePoint_ = initialCondition;
  t_ = initialCondition.get_t();
  t_old_ = t_;
  solution_vector_ = initialCondition.get_x()->clone_v();
  solution_vector_old_ = solution_vector_->clone_v();
  haveInitialCondition_ = true;
}


// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_FORWARD_EULER_STEPPER_INSTANT(SCALAR) \
  \
  template class ForwardEulerStepper< SCALAR >; \
  \
  template RCP<ForwardEulerStepper< SCALAR > > forwardEulerStepper(); \
  template RCP<ForwardEulerStepper< SCALAR > > forwardEulerStepper( \
      const RCP<const Thyra::ModelEvaluator< SCALAR > > &model \
      ); 


} // namespace Rythmos

#endif //Rythmos_FORWARDEULER_STEPPER_DEF_H
