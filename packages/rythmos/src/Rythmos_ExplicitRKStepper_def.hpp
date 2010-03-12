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

#ifndef Rythmos_ExplicitRK_STEPPER_DEF_H
#define Rythmos_ExplicitRK_STEPPER_DEF_H

#include "Rythmos_ExplicitRKStepper_decl.hpp"

#include "Rythmos_RKButcherTableau.hpp"
#include "Rythmos_RKButcherTableauHelpers.hpp"
#include "Rythmos_RKButcherTableauBuilder.hpp"
#include "Rythmos_StepperHelpers.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Rythmos {

// Non-member constructors
template<class Scalar>
RCP<ExplicitRKStepper<Scalar> > explicitRKStepper()
{
  RCP<ExplicitRKStepper<Scalar> > stepper = rcp(new ExplicitRKStepper<Scalar>());
  return stepper;
}

template<class Scalar>
RCP<ExplicitRKStepper<Scalar> > explicitRKStepper(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model 
    )
{
  RCP<RKButcherTableauBase<Scalar> > rkbt = createRKBT<Scalar>("Explicit 4 Stage");
  //RCP<RKButcherTableauBase<Scalar> > rkbt = rcp(new Explicit4Stage4thOrder_RKBT<Scalar>());
  RCP<ExplicitRKStepper<Scalar> > stepper = rcp(new ExplicitRKStepper<Scalar>());
  stepper->setModel(model);
  stepper->setRKButcherTableau(rkbt);
  return stepper;
}

template<class Scalar>
RCP<ExplicitRKStepper<Scalar> > explicitRKStepper(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const RCP<const RKButcherTableauBase<Scalar> >& rkbt 
    )
{
  RCP<ExplicitRKStepper<Scalar> > stepper = rcp(new ExplicitRKStepper<Scalar>());
  stepper->setModel(model);
  stepper->setRKButcherTableau(rkbt);
  return stepper;
}

template<class Scalar>
ExplicitRKStepper<Scalar>::ExplicitRKStepper()
{
  this->defaultInitializeAll_();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(15);
  erkButcherTableau_ = rKButcherTableau<Scalar>();
  numSteps_ = 0;
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::defaultInitializeAll_()
{
  model_ = Teuchos::null;
  solution_vector_ = Teuchos::null;
  solution_vector_old_ = Teuchos::null;
  //k_vector_;
  ktemp_vector_ = Teuchos::null;
  //basePoint_;
  erkButcherTableau_ = Teuchos::null;
  t_ = ST::nan();
  t_old_ = ST::nan();
  dt_ = ST::nan();
  numSteps_ = -1;
  parameterList_ = Teuchos::null;
  isInitialized_ = false;
  haveInitialCondition_ = false;
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::setRKButcherTableau(const RCP<const RKButcherTableauBase<Scalar> > &rkbt)
{
  TEUCHOS_ASSERT( !is_null(rkbt) );
  validateERKButcherTableau(*rkbt);
  int numStages_old = erkButcherTableau_->numStages();
  int numStages_new = rkbt->numStages();
  TEST_FOR_EXCEPTION( numStages_new == 0, std::logic_error,
      "Error!  The Runge-Kutta Butcher tableau has no stages!"
      );
  if (!is_null(model_)) {
    int numNewStages = numStages_new - numStages_old;
    if ( numNewStages > 0 ) {
      k_vector_.reserve(numStages_new);
      for (int i=0 ; i<numNewStages ; ++i) {
        k_vector_.push_back(Thyra::createMember(model_->get_f_space()));
      }
    }
  }
  erkButcherTableau_ = rkbt;
}

template<class Scalar>
RCP<const RKButcherTableauBase<Scalar> > ExplicitRKStepper<Scalar>::getRKButcherTableau() const
{
  return erkButcherTableau_;
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::initialize_()
{
  if (!isInitialized_) {
    TEUCHOS_ASSERT( !is_null(model_) );
    TEUCHOS_ASSERT( !is_null(erkButcherTableau_) );
    TEUCHOS_ASSERT( haveInitialCondition_ );
    TEST_FOR_EXCEPTION( erkButcherTableau_->numStages() == 0, std::logic_error,
        "Error!  The Runge-Kutta Butcher tableau has no stages!"
        );
    ktemp_vector_ = Thyra::createMember(model_->get_f_space());
    // Initialize the stage vectors
    int numStages = erkButcherTableau_->numStages();
    k_vector_.reserve(numStages);
    for (int i=0 ; i<numStages ; ++i) {
      k_vector_.push_back(Thyra::createMember(model_->get_f_space()));
    }
  }
#ifdef RYTHMOS_DEBUG
  THYRA_ASSERT_VEC_SPACES(
    "Rythmos::ExplicitRKStepper::initialize_(...)",
    *solution_vector_->space(), *model_->get_x_space() );
#endif // RYTHMOS_DEBUG
  isInitialized_ = true;
}


template<class Scalar>
ExplicitRKStepper<Scalar>::~ExplicitRKStepper()
{
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > ExplicitRKStepper<Scalar>::get_x_space() const
{
  TEUCHOS_ASSERT( !is_null(model_) );
  return(model_->get_x_space());
}

template<class Scalar>
Scalar ExplicitRKStepper<Scalar>::takeStep(Scalar dt, StepSizeType flag)
{
  typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag TScalarMag;
  this->initialize_();
#ifdef RYTHMOS_DEBUG
    TEST_FOR_EXCEPTION( flag == STEP_TYPE_VARIABLE, std::logic_error,
        "Error!  ExplicitRKStepper does not support variable time steps at this time."
        );
#endif // RYTHMOS_DEBUG
  if ((flag == STEP_TYPE_VARIABLE) || (dt == ST::zero())) {
    return(Scalar(-ST::one()));
  }
  // Store old solution & old time
  V_V(&*solution_vector_old_, *solution_vector_);
  t_old_ = t_;

  dt_ = dt;

  int stages = erkButcherTableau_->numStages();
  Teuchos::SerialDenseMatrix<int,Scalar> A = erkButcherTableau_->A();
  Teuchos::SerialDenseVector<int,Scalar> b = erkButcherTableau_->b();
  Teuchos::SerialDenseVector<int,Scalar> c = erkButcherTableau_->c();
  // Compute stage solutions
  for (int s=0 ; s < stages ; ++s) {
    Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
    for (int j=0 ; j < s ; ++j) { // assuming Butcher matix is strictly lower triangular
      if (A(s,j) != ST::zero()) {
        Thyra::Vp_StV(&*ktemp_vector_, A(s,j), *k_vector_[j]); // ktemp = ktemp + a_{s+1,j+1}*k_{j+1}
      }
    }
    TScalarMag ts = t_ + c(s)*dt;
    eval_model_explicit<Scalar>(*model_,basePoint_,*ktemp_vector_,ts,Teuchos::outArg(*k_vector_[s]));
    Thyra::Vt_S(&*k_vector_[s],dt); // k_s = k_s*dt
  } 
  // Sum for solution:
  for (int s=0 ; s < stages ; ++s) {
    if (b(s) != ST::zero()) {
      Thyra::Vp_StV(&*solution_vector_, b(s), *k_vector_[s]); // solution_vector += b_{s+1}*k_{s+1}
    }
  }

  // update current time:
  t_ = t_ + dt;

  numSteps_++;

  return(dt);
}

template<class Scalar>
const StepStatus<Scalar> ExplicitRKStepper<Scalar>::getStepStatus() const
{
  StepStatus<Scalar> stepStatus;

  if (!haveInitialCondition_) {
    stepStatus.stepStatus = STEP_STATUS_UNINITIALIZED;
  } else if (numSteps_ == 0) {
    stepStatus.stepStatus = STEP_STATUS_UNKNOWN;
    stepStatus.order = erkButcherTableau_->order();
    stepStatus.time = t_;
    stepStatus.solution = solution_vector_;
  } else {
    stepStatus.stepStatus = STEP_STATUS_CONVERGED;
    stepStatus.stepSize = dt_;
    stepStatus.order = erkButcherTableau_->order();
    stepStatus.time = t_;
    stepStatus.solution = solution_vector_;
  }

  return(stepStatus);
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const
{
  if ( (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_DEFAULT) ) ||
       (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)     )
     ) {
    out << this->description() << "::describe" << std::endl;
    out << "model = " << model_->description() << std::endl;
    out << erkButcherTableau_->numStages() << " stage Explicit RK method" << std::endl;
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)) {
    out << "solution_vector = " << std::endl;
    out << Teuchos::describe(*solution_vector_,verbLevel);
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM)) {
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH)) {
    out << "model = " << std::endl;
    out << Teuchos::describe(*model_,verbLevel);
    int stages = erkButcherTableau_->numStages();
    for (int i=0 ; i<stages ; ++i) {
      out << "k_vector[" << i << "] = " << std::endl;
      out << Teuchos::describe(*k_vector_[i],verbLevel);
    }
    out << "ktemp_vector = " << std::endl;
    out << Teuchos::describe(*ktemp_vector_,verbLevel);
    out << "ERK Butcher Tableau A matrix: " << erkButcherTableau_->A() << std::endl;
    out << "ERK Butcher Tableau b vector: " << erkButcherTableau_->b() << std::endl;
    out << "ERK Butcher Tableau c vector: " << erkButcherTableau_->c() << std::endl;
    out << "t = " << t_ << std::endl;
  }
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::addPoints(
    const Array<Scalar>& time_vec
    ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec
    ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    )
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, addPoints is not implemented for ExplicitRKStepper at this time.\n");
}

template<class Scalar>
TimeRange<Scalar> ExplicitRKStepper<Scalar>::getTimeRange() const
{
  if (!haveInitialCondition_) {
    return(invalidTimeRange<Scalar>());
  } else {
    return(TimeRange<Scalar>(t_old_,t_));
  }
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{
  TEUCHOS_ASSERT( haveInitialCondition_ );
  using Teuchos::constOptInArg;
  using Teuchos::null;
  defaultGetPoints<Scalar>(
      t_old_, constOptInArg(*solution_vector_old_),
      Ptr<const VectorBase<Scalar> >(null),
      t_, constOptInArg(*solution_vector_),
      Ptr<const VectorBase<Scalar> >(null),
      time_vec,ptr(x_vec), ptr(xdot_vec), ptr(accuracy_vec),
      Ptr<InterpolatorBase<Scalar> >(null)
      );
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  TEUCHOS_ASSERT( time_vec != NULL );
  time_vec->clear();
  if (!haveInitialCondition_) {
    return;
  }
  time_vec->push_back(t_old_);
  if (t_ != t_old_) {
    time_vec->push_back(t_);
  }
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, removeNodes is not implemented for ExplicitRKStepper at this time.\n");
}

template<class Scalar>
int ExplicitRKStepper<Scalar>::getOrder() const
{
  return(erkButcherTableau_->order());
}

template <class Scalar>
void ExplicitRKStepper<Scalar>::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> ExplicitRKStepper<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> ExplicitRKStepper<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
ExplicitRKStepper<Scalar>::getValidParameters() const
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
void ExplicitRKStepper<Scalar>::setModel(const RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  TEST_FOR_EXCEPT( is_null(model) );
  TEST_FOR_EXCEPT( !is_null(model_) ); // For now you can only call this once.
  assertValidModel( *this, *model );
  model_ = model;
}


template<class Scalar>
void ExplicitRKStepper<Scalar>::setNonconstModel(const RCP<Thyra::ModelEvaluator<Scalar> >& model)
{
  this->setModel(model); // TODO 09/09/09 tscoffe:  use ConstNonconstObjectContainer!
}


template<class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
ExplicitRKStepper<Scalar>::getModel() const
{
  return model_;
}


template<class Scalar>
Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >
ExplicitRKStepper<Scalar>::getNonconstModel() 
{
  return Teuchos::null;
}


template<class Scalar>
void ExplicitRKStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  basePoint_ = initialCondition;

  // x

  RCP<const Thyra::VectorBase<Scalar> >
    x_init = initialCondition.get_x();

#ifdef RYTHMOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(x_init), std::logic_error,
    "Error, if the client passes in an intial condition to setInitialCondition(...),\n"
    "then x can not be null!" );
#endif

  solution_vector_ = x_init->clone_v();
  solution_vector_old_ = x_init->clone_v();

  // t
  
  t_ = initialCondition.get_t();
  t_old_ = t_;

  haveInitialCondition_ = true;

}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> 
ExplicitRKStepper<Scalar>::getInitialCondition() const
{
  return basePoint_;
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::supportsCloning() const
{
  return true;
}

template<class Scalar>
RCP<StepperBase<Scalar> > ExplicitRKStepper<Scalar>::cloneStepperAlgorithm() const
{
  // Just use the interface to clone the algorithm in a basically
  // uninitialized state
  RCP<ExplicitRKStepper<Scalar> >
    stepper = Teuchos::rcp(new ExplicitRKStepper<Scalar>());

  if (!is_null(model_)) {
    stepper->setModel(model_); // Shallow copy is okay!
  }

  if (!is_null(erkButcherTableau_)) {
    // 06/16/09 tscoffe:  should we clone the RKBT here?
    stepper->setRKButcherTableau(erkButcherTableau_);
  }

  if (!is_null(parameterList_)) {
    stepper->setParameterList(Teuchos::parameterList(*parameterList_));
  }

  return stepper;

}
    


// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_EXPLICIT_RK_STEPPER_INSTANT(SCALAR) \
  \
  template class ExplicitRKStepper< SCALAR >; \
  \
  template RCP< ExplicitRKStepper< SCALAR > > \
  explicitRKStepper();  \
  \
  template RCP< ExplicitRKStepper< SCALAR > > \
  explicitRKStepper( \
    const RCP<Thyra::ModelEvaluator< SCALAR > >& model \
      ); \
  \
  template RCP< ExplicitRKStepper< SCALAR > > \
  explicitRKStepper( \
    const RCP<Thyra::ModelEvaluator< SCALAR > >& model, \
    const RCP<const RKButcherTableauBase< SCALAR > >& rkbt \
      ); \
   
} // namespace Rythmos

#endif //Rythmos_ExplicitRK_STEPPER_DEF_H

