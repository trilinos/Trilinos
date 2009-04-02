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

#include "Rythmos_StepperHelpers.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Thyra_ModelEvaluatorHelpers.hpp"

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
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model 
    )
{
  //RCP<RKButcherTableauBase<Scalar> > rkbt = createRKBT<Scalar>("Explicit 4 Stage");
  RCP<RKButcherTableauBase<Scalar> > rkbt = rcp(new Explicit4Stage4thOrder_RKBT<Scalar>());
  RCP<ExplicitRKStepper<Scalar> > stepper = rcp(new ExplicitRKStepper<Scalar>(model,rkbt));
  return stepper;
}

template<class Scalar>
RCP<ExplicitRKStepper<Scalar> > explicitRKStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    const RCP<const RKButcherTableauBase<Scalar> > &rkbt 
    )
{
  RCP<ExplicitRKStepper<Scalar> > stepper = rcp(new ExplicitRKStepper<Scalar>(model,rkbt));
  return stepper;
}

template<class Scalar>
ExplicitRKStepper<Scalar>::ExplicitRKStepper(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &model, 
    const RCP<const RKButcherTableauBase<Scalar> > &rkbt
    )
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(15);

  t_ = ST::nan();
  t_old_ = ST::nan();
  dt_ = ST::nan();
  erkButcherTableau_ = rKButcherTableau<Scalar>();
  numSteps_ = 0;
  haveInitialCondition_ = false;
  this->setModel(model);
  this->setRKButcherTableau(rkbt);
  initialize_();
}

template<class Scalar>
ExplicitRKStepper<Scalar>::ExplicitRKStepper()
  : isInitialized_(false)
{
  t_ = ST::nan();
  t_old_ = ST::nan();
  dt_ = ST::nan();
  erkButcherTableau_ = rKButcherTableau<Scalar>();
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
  int numNewStages = numStages_new - numStages_old;
  if ( numNewStages > 0 ) {
    k_vector_.reserve(numStages_new);
    for (int i=0 ; i<numNewStages ; ++i) {
      k_vector_.push_back(Thyra::createMember(model_->get_f_space()));
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
  haveInitialCondition_ = setDefaultInitialConditionFromNominalValues<Scalar>(
    *model_, Teuchos::ptr(this) );
  ktemp_vector_ = Thyra::createMember(model_->get_f_space());

  isInitialized_ = true;
}


template<class Scalar>
ExplicitRKStepper<Scalar>::~ExplicitRKStepper()
{
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > ExplicitRKStepper<Scalar>::get_x_space() const
{
  TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,"Error, attempting to call get_x_space before initialization!\n");
  return(solution_vector_->space());
}

template<class Scalar>
Scalar ExplicitRKStepper<Scalar>::takeStep(Scalar dt, StepSizeType flag)
{
  typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag ScalarMag;
  TEST_FOR_EXCEPTION( erkButcherTableau_->numStages() == 0, std::logic_error,
      "Error!  The Runge-Kutta Butcher Tableau has no stages!"
      );
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
    ScalarMag ts = t_ + c(s)*dt;
    Thyra::eval_f<Scalar>(*model_,*ktemp_vector_,ts,&*k_vector_[s]);
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
  typedef Teuchos::ScalarTraits<Scalar> ST;
  if (x_vec != NULL) {
    x_vec->clear();
  }
  if (xdot_vec != NULL) {
    xdot_vec->clear();
  }
  if (accuracy_vec != NULL) {
    accuracy_vec->clear();
  }
  typename Array<Scalar>::const_iterator time_it = time_vec.begin();
  RCP<Thyra::VectorBase<Scalar> > tmpVec;
  for (; time_it != time_vec.end() ; time_it++) {
    Scalar t = *time_it;
    if (compareTimeValues(t,t_old_)==0) {
      tmpVec = solution_vector_old_;
    } else if (compareTimeValues(t,t_)==0) {
      tmpVec = solution_vector_;
    } else {
      TEST_FOR_EXCEPTION(true,std::logic_error,"Error, ExplicitRKStepper::getPoints only supports time values on the boundaries!\n");
    }
    if (!Teuchos::is_null(tmpVec)) {
      if (x_vec != NULL) {
        x_vec->push_back(tmpVec->clone_v());
      }
      if (xdot_vec != NULL) {
        xdot_vec->push_back(Teuchos::null);
      }
      if (accuracy_vec != NULL) {
        accuracy_vec->push_back(ST::zero());
      }
      tmpVec = Teuchos::null;
    }
  }
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  TEST_FOR_EXCEPTION(!isInitialized_,std::logic_error,"Error, attempting to call getNodes before initialization!\n");
  if (time_vec != NULL) {
    time_vec->clear();
    time_vec->push_back(t_old_);
    if (t_ != t_old_) {
      time_vec->push_back(t_);
    }
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
void ExplicitRKStepper<Scalar>::setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model)
{
  TEST_FOR_EXCEPT( is_null(model) );
  TEST_FOR_EXCEPT( !is_null(model_) ); // For now you can only call this once.
  assertValidModel( *this, *model );
  model_ = model;
}

template<class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
ExplicitRKStepper<Scalar>::getModel() const
{
  return model_;
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  TEST_FOR_EXCEPT( is_null(model_) );

  basePoint_ = initialCondition;

  // x

  RCP<const Thyra::VectorBase<Scalar> >
    x_init = initialCondition.get_x();

#ifdef HAVE_RYTHMOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(x_init), std::logic_error,
    "Error, if the client passes in an intial condition to setInitialCondition(...),\n"
    "then x can not be null!" );
  THYRA_ASSERT_VEC_SPACES(
    "Rythmos::ExplicitRKStepper::setInitialCondition(...)",
    *x_init->space(), *model_->get_x_space() );
#endif

  solution_vector_ = x_init->clone_v();
  solution_vector_old_ = x_init->clone_v();

  // t
  
  t_ = initialCondition.get_t();
  t_old_ = t_;

  haveInitialCondition_ = true;

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
    const RCP<const Thyra::ModelEvaluator< SCALAR > > &model \
      ); \
  \
  template RCP< ExplicitRKStepper< SCALAR > > \
  explicitRKStepper( \
    const RCP<const Thyra::ModelEvaluator< SCALAR > > &model, \
    const RCP<const RKButcherTableauBase< SCALAR > > &rkbt \
      ); \
   
} // namespace Rythmos

#endif //Rythmos_ExplicitRK_STEPPER_DEF_H

