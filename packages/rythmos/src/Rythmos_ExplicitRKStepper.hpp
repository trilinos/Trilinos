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
    ExplicitRKStepper(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model_);

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
      ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec
      ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
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

    int stages_; // Number of stages of RK
    Array<Array<Scalar> > b_A_; // Butcher tableau A matrix
    Array<Scalar> b_b_; // Butcher b vector
    Array<Scalar> b_c_; // Butcher c vector

    Scalar t_;
    Scalar dt_;

    Teuchos::RCP<Teuchos::ParameterList> parameterList_;

    bool isInitialized_;

    // Private member functions:
    void initialize_();

};

template<class Scalar>
ExplicitRKStepper<Scalar>::ExplicitRKStepper(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(15);

  this->setModel(model);
  initialize_();
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::initialize_()
{
  t_ = ST::zero();
  solution_vector_ = model_->getNominalValues().get_x()->clone_v();
  stages_ = 4; // 4 stage ERK
  k_vector_.reserve(stages_);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
    f_space = model_->get_f_space();
  for (int i=0 ; i<stages_ ; ++i) {
    k_vector_.push_back(Thyra::createMember(f_space));
  }
  ktemp_vector_ = Thyra::createMember(f_space);

  // initialize the Butcher tableau and its vectors
  Scalar zero = ST::zero();
  b_A_.reserve(stages_);
  b_b_.reserve(stages_);
  b_c_.reserve(stages_);

  // fill b with zeros
  b_b_.push_back(zero); 
  b_b_.push_back(zero); 
  b_b_.push_back(zero);
  b_b_.push_back(zero);
  // fill c with zeros
  b_c_.push_back(zero); 
  b_c_.push_back(zero); 
  b_c_.push_back(zero);
  b_c_.push_back(zero);
  // fill A with zeros
  for (int i=0 ; i<stages_ ; ++i) {
    b_A_.push_back(b_b_);
  }

  // Runge-Kutta methods in Butcher tableau form:
  //  c | A    A = sxs matrix
  //   -|---   b = s vector
  //      b'   c = s vector
  
/*
  // 3/8 Rule Runge-Kutta Method: (not implemented yet)
  // c = [  0  1/3 2/3  1  ]'
  // A = [  0              ]
  //     [ 1/3  0          ]
  //     [-1/3  1   0      ]
  //     [  1  -1   1   0  ]
  // b = [ 1/8 3/8 3/8 1/8 ]'
  Scalar one = ST::one();
  Scalar one_third    = Scalar(ST::one()/(3*ST::one()));
  Scalar two_third    = Scalar(2*ST::one()/(3*ST::one()));
  Scalar one_eighth   = Scalar(ST::one()/(8*ST::one()));
  Scalar three_eighth = Scalar(3*ST::one()/(8*ST::one()));

  // fill b_A_
  b_A_[0][0] = zero;
  b_A_[0][1] = zero;
  b_A_[0][2] = zero;
  b_A_[0][3] = zero;

  b_A_[1][0] = one_third;
  b_A_[1][1] = zero;
  b_A_[1][2] = zero;
  b_A_[1][3] = zero;

  b_A_[2][0] = Scalar(-one_third);
  b_A_[2][1] = one;
  b_A_[2][2] = zero;
  b_A_[2][3] = zero;

  b_A_[3][0] = one;
  b_A_[3][1] = Scalar(-one);
  b_A_[3][2] = one;
  b_A_[3][3] = zero;

  // fill b_b_
  b_b_[0] = one_eighth;
  b_b_[1] = three_eighth;
  b_b_[2] = three_eighth;
  b_b_[3] = one_eighth;
  
  // fill b_c_
  b_c_[0] = zero;
  b_c_[1] = one_third;
  b_c_[2] = two_third;
  b_c_[3] = one;
*/


  // "The" Runge-Kutta Method: (implemented below)
  // c = [  0  1/2 1/2  1  ]'
  // A = [  0              ] 
  //     [ 1/2  0          ]
  //     [  0  1/2  0      ]
  //     [  0   0   1   0  ]
  // b = [ 1/6 1/3 1/3 1/6 ]'
  Scalar one = ST::one();
  Scalar onehalf = ST::one()/(2*ST::one());
  Scalar onesixth = ST::one()/(6*ST::one());
  Scalar onethird = ST::one()/(3*ST::one());

  // fill b_A_
  b_A_[0][0] = zero;
  b_A_[0][1] = zero;
  b_A_[0][2] = zero;
  b_A_[0][3] = zero;

  b_A_[1][0] = onehalf;
  b_A_[1][1] = zero;
  b_A_[1][2] = zero;
  b_A_[1][3] = zero;

  b_A_[2][0] = zero;
  b_A_[2][1] = onehalf;
  b_A_[2][2] = zero;
  b_A_[2][3] = zero;

  b_A_[3][0] = zero;
  b_A_[3][1] = zero;
  b_A_[3][2] = one;
  b_A_[3][3] = zero;

  // fill b_b_
  b_b_[0] = onesixth;
  b_b_[1] = onethird;
  b_b_[2] = onethird;
  b_b_[3] = onesixth;
  
  // fill b_c_
  b_c_[0] = zero;
  b_c_[1] = onehalf;
  b_c_[2] = onehalf;
  b_c_[3] = one;

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

  // Compute stage solutions
  for (int s=0 ; s < stages_ ; ++s) {
    Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
    for (int j=0 ; j < s ; ++j) { // assuming Butcher matix is strictly lower triangular
      if (b_A_[s][j] != ST::zero())
        Thyra::Vp_StV(&*ktemp_vector_, b_A_[s][j], *k_vector_[j]); // ktemp = ktemp + a_{s+1,j+1}*k_{j+1}
    }
    ScalarMag ts = t_ + b_c_[s]*dt;
    Thyra::eval_f<Scalar>(*model_,*ktemp_vector_,ts,&*k_vector_[s]);
    Thyra::Vt_S(&*k_vector_[s],dt); // k_s = k_s*dt
  } 
  // Sum for solution:
  for (int s=0 ; s < stages_ ; ++s) {
    if (b_b_[s] != ST::zero()) {
      Thyra::Vp_StV(&*solution_vector_, b_b_[s], *k_vector_[s]); // solution_vector += b_{s+1}*k_{s+1}
    }
  }

  // update current time:
  t_ = t_ + dt;

  return(dt);
}

template<class Scalar>
const StepStatus<Scalar> ExplicitRKStepper<Scalar>::getStepStatus() const
{
  StepStatus<Scalar> stepStatus;

  stepStatus.stepSize = dt_;
  stepStatus.order = -1;
  stepStatus.time = t_;
  stepStatus.solution = solution_vector_;

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
    out << stages_ << " stage Explicit RK method" << std::endl;
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)) {
    out << "solution_vector = " << std::endl;
    solution_vector_->describe(out,verbLevel,leadingIndent,indentSpacer); 
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM)) {
  } else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH)) {
    out << "model = " << std::endl;
    model_->describe(out,verbLevel,leadingIndent,indentSpacer); 
    for (int i=0 ; i<stages_ ; ++i) {
      out << "k_vector[" << i << "] = " << std::endl;
      out << k_vector_[i]->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    }
    out << "ktemp_vector = " << std::endl;
    ktemp_vector_->describe(out,verbLevel,leadingIndent,indentSpacer); 
    for (int i=0 ; i<stages_ ; ++i) {
      for (int j=0 ; j<stages_ ; ++j) {
        out << "b_A_[" << i << "][" << j << "] = " << b_A_[i][j] << std::endl;
      }
    }
    for (int i=0 ; i<stages_ ; ++i) {
      out << "b_b_[" << i << "] = " << b_b_[i] << std::endl;
    }
    for (int i=0 ; i<stages_ ; ++i) {
      out << "b_c_[" << i << "] = " << b_c_[i] << std::endl;
    }
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
void ExplicitRKStepper<Scalar>::getPoints(
    const Array<Scalar>& time_vec
    ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* x_vec
    ,Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
    ,Array<ScalarMag>* accuracy_vec) const
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, getPoints is not implemented for ExplicitRKStepper at this time.\n");
}

template<class Scalar>
TimeRange<Scalar> ExplicitRKStepper<Scalar>::getTimeRange() const
{
  return invalidTimeRange<Scalar>();
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
  return(4);
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
  TEST_FOR_EXCEPT(model == Teuchos::null)
  model_ = model;
}

template<class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
ExplicitRKStepper<Scalar>::getModel() const
{
  return model_;
}

} // namespace Rythmos

#endif //Rythmos_ExplicitRK_STEPPER_H
