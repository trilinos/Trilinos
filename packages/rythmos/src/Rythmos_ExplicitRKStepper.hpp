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

#ifndef Rythmos_ExplicitRK_STEPPER_H
#define Rythmos_ExplicitRK_STEPPER_H

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_StepperHelpers.hpp"
#include "Rythmos_RKButcherTableau.hpp"
#include "Teuchos_RCP.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"

namespace Rythmos {

/** \brief . */
template<class Scalar>
class ExplicitRKStepper : virtual public StepperBase<Scalar>
{
  public:
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    
    /** \brief . */
    ExplicitRKStepper();
    ExplicitRKStepper(
        const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
        RKButcherTableau<Scalar> rkbt
        );

    /** \brief. */
    void setRKButcherTableau(RKButcherTableau<Scalar> rkbt);

    /** \brief. */
    RKButcherTableau<Scalar> getRKButcherTableau() const;

    /** \brief . */
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

    /** \brief . */
    void setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model);

    /** \brief . */
    Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
    getModel() const;
    
    /** \brief . */
    ~ExplicitRKStepper();

    /** \brief . */
    void setInitialCondition(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
      );

    /** \brief . */
    Scalar takeStep(Scalar dt, StepSizeType flag);

    /** \brief . */
    const StepStatus<Scalar> getStepStatus() const;

    /** \brief . */
    std::ostream& describe(
      std::ostream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ,const std::string          leadingIndent
      ,const std::string          indentSpacer
      ) const;

    /// Redefined from InterpolationBufferBase 
    /// Add points to buffer
    void addPoints(
      const Array<Scalar>& time_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
      );

    /// Get values from buffer
    void getPoints(
      const Array<Scalar>& time_vec
      ,Array<RCP<const VectorBase<Scalar> > >* x_vec
      ,Array<RCP<const VectorBase<Scalar> > >* xdot_vec
      ,Array<ScalarMag>* accuracy_vec) const;

    /** \brief . */
    TimeRange<Scalar> getTimeRange() const;

    /// Get interpolation nodes
    void getNodes(Array<Scalar>* time_vec) const;

    /// Remove interpolation nodes
    void removeNodes(Array<Scalar>& time_vec);

    /// Get order of interpolation
    int getOrder() const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    

  private:

    Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_;
    Teuchos::RCP<Thyra::VectorBase<Scalar> > solution_vector_;
    Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > k_vector_;
    Teuchos::RCP<Thyra::VectorBase<Scalar> > ktemp_vector_;

    Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;

    RKButcherTableau<Scalar> erkButcherTableau_;

    Scalar t_;
    Scalar dt_;
    int numSteps_;

    Teuchos::RCP<Teuchos::ParameterList> parameterList_;

    bool isInitialized_;

    bool haveInitialCondition_;

    // Private member functions:
    void initialize_();

};

// Non-member constructor
template<class Scalar>
RCP<ExplicitRKStepper<Scalar> > explicitRKStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model 
    )
{
  RKButcherTableau<Scalar> rkbt = createExplicit4Stage4thOrder_RKBT<Scalar>();
  RCP<ExplicitRKStepper<Scalar> > stepper = rcp(new ExplicitRKStepper<Scalar>(model,rkbt));
  return stepper;
}

template<class Scalar>
RCP<ExplicitRKStepper<Scalar> > explicitRKStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    const RKButcherTableau<Scalar> rkbt 
    )
{
  RCP<ExplicitRKStepper<Scalar> > stepper = rcp(new ExplicitRKStepper<Scalar>(model,rkbt));
  return stepper;
}

template<class Scalar>
ExplicitRKStepper<Scalar>::ExplicitRKStepper(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model, RKButcherTableau<Scalar> rkbt)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(15);

  numSteps_ = 0;
  haveInitialCondition_ = false;
  this->setModel(model);
  this->setRKButcherTableau(rkbt);
  initialize_();
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::setRKButcherTableau(RKButcherTableau<Scalar> rkbt)
{
  validateERKButcherTableau(rkbt);
  int numStages_old = erkButcherTableau_.numStages();
  int numStages_new = rkbt.numStages();
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
RKButcherTableau<Scalar> ExplicitRKStepper<Scalar>::getRKButcherTableau() const
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
ExplicitRKStepper<Scalar>::ExplicitRKStepper()
  : isInitialized_(false)
{
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
  if ((flag == STEP_TYPE_VARIABLE) || (dt == ST::zero())) {
    return(Scalar(-ST::one()));
  }
  dt_ = dt;

  int stages = erkButcherTableau_.numStages();
  Teuchos::SerialDenseMatrix<int,Scalar> A = erkButcherTableau_.A();
  Teuchos::SerialDenseVector<int,Scalar> b = erkButcherTableau_.b();
  Teuchos::SerialDenseVector<int,Scalar> c = erkButcherTableau_.c();
  // Compute stage solutions
  for (int s=0 ; s < stages ; ++s) {
    Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
    for (int j=0 ; j < s ; ++j) { // assuming Butcher matix is strictly lower triangular
      if (A(s,j) != ST::zero())
        Thyra::Vp_StV(&*ktemp_vector_, A(s,j), *k_vector_[j]); // ktemp = ktemp + a_{s+1,j+1}*k_{j+1}
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
    stepStatus.order = erkButcherTableau_.order();
    stepStatus.time = t_;
    stepStatus.solution = solution_vector_;
  } else {
    stepStatus.stepStatus = STEP_STATUS_CONVERGED;
    stepStatus.stepSize = dt_;
    stepStatus.order = erkButcherTableau_.order();
    stepStatus.time = t_;
    stepStatus.solution = solution_vector_;
  }

  return(stepStatus);
}

template<class Scalar>
std::ostream& ExplicitRKStepper<Scalar>::describe(
      std::ostream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ,const std::string          leadingIndent
      ,const std::string          indentSpacer
      ) const
{
  if ( (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_DEFAULT) ) ||
       (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)     )
     ) {
    out << this->description() << "::describe" << std::endl;
    out << "model = " << model_->description() << std::endl;
    out << erkButcherTableau_.numStages() << " stage Explicit RK method" << std::endl;
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)) {
    out << "solution_vector = " << std::endl;
    solution_vector_->describe(out,verbLevel,leadingIndent,indentSpacer); 
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM)) {
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH)) {
    out << "model = " << std::endl;
    model_->describe(out,verbLevel,leadingIndent,indentSpacer); 
    int stages = erkButcherTableau_.numStages();
    for (int i=0 ; i<stages ; ++i) {
      out << "k_vector[" << i << "] = " << std::endl;
      out << k_vector_[i]->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    }
    out << "ktemp_vector = " << std::endl;
    ktemp_vector_->describe(out,verbLevel,leadingIndent,indentSpacer); 
    out << "ERK Butcher Tableau A matrix: " << erkButcherTableau_.A() << std::endl;
    out << "ERK Butcher Tableau b vector: " << erkButcherTableau_.b() << std::endl;
    out << "ERK Butcher Tableau c vector: " << erkButcherTableau_.c() << std::endl;
    out << "t = " << t_ << std::endl;
  }
  return(out);
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
    return(TimeRange<Scalar>(t_,t_));
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
  if ((time_vec.size() == 1) && (time_vec[0] == t_)) {
    x_vec->clear();
    x_vec->push_back(solution_vector_->clone_v());
    xdot_vec->clear();
    xdot_vec->push_back(Teuchos::null);
    accuracy_vec->clear();
    accuracy_vec->push_back(ST::zero());
  } else {
    TEST_FOR_EXCEPTION(true,std::logic_error,"Error, ExplicitRKStepper::getPoints is only implemented for the last time point computed.\n");
  }
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, getNodes is not implemented for ExplicitRKStepper at this time.\n");
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::removeNodes(Array<Scalar>& time_vec) 
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, removeNodes is not implemented for ExplicitRKStepper at this time.\n");
}

template<class Scalar>
int ExplicitRKStepper<Scalar>::getOrder() const
{
  return(erkButcherTableau_.order());
}

template <class Scalar>
void ExplicitRKStepper<Scalar>::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  parameterList_ = paramList;
  int outputLevel = parameterList_->get( "outputLevel", int(-1) );
  outputLevel = std::min(std::max(outputLevel,-1),4);
  this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
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
void ExplicitRKStepper<Scalar>::setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model)
{
  TEST_FOR_EXCEPT(model == Teuchos::null);
  TEST_FOR_EXCEPT( !Teuchos::is_null(model_) ); // For now you can only call this once.
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

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    is_null(x_init), std::logic_error,
    "Error, if the client passes in an intial condition to setInitialCondition(...),\n"
    "then x can not be null!" );
  THYRA_ASSERT_VEC_SPACES(
    "Rythmos::BackwardEulerStepper::setInitialCondition(...)",
    *x_init->space(), *model_->get_x_space() );
#endif

  solution_vector_ = x_init->clone_v();

  // t
  
  t_ = initialCondition.get_t();

  haveInitialCondition_ = true;

}

} // namespace Rythmos

#endif //Rythmos_ExplicitRK_STEPPER_H

