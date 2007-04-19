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
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"

namespace Rythmos {

/** \brief . */
template<class Scalar>
class ExplicitRKStepper : virtual public StepperBase<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    
    /** \brief . */
    ExplicitRKStepper();
    ExplicitRKStepper(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model_);

    /** \brief . */
    void setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model);
    
    /** \brief . */
    ~ExplicitRKStepper();

    /** \brief . */
    Scalar TakeStep(Scalar dt, StepSizeType flag);

    /** \brief . */
    const StepStatus<Scalar> getStepStatus();

    /** \brief . */
    std::string description() const;

    /** \brief . */
    std::ostream& describe(
      std::ostream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ,const std::string          leadingIndent
      ,const std::string          indentSpacer
      ) const;

    /// Redefined from InterpolationBufferBase 
    /// Add points to buffer
    bool SetPoints(
      const std::vector<Scalar>& time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
      ,const std::vector<ScalarMag> & accuracy_vec 
      );
    
    /// Get values from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
      ,std::vector<ScalarMag>* accuracy_vec) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar> & IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_vec) const;

    /// Remove interpolation nodes
    bool RemoveNodes(std::vector<Scalar>& time_vec);

    /// Get order of interpolation
    int GetOrder() const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
    

  private:

    Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > solution_vector_;
    std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > k_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > ktemp_vector_;

    int stages_; // Number of stages of RK
    std::vector<std::vector<Scalar> > b_A_; // Butcher tableau A matrix
    std::vector<Scalar> b_b_; // Butcher b vector
    std::vector<Scalar> b_c_; // Butcher c vector

    Scalar t_;
    Scalar dt_;

    Teuchos::RefCountPtr<Teuchos::ParameterList> parameterList_;

};

template<class Scalar>
ExplicitRKStepper<Scalar>::ExplicitRKStepper(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model)
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(15);
  out->setMaxLenLinePrefix(30);
  //out->pushLinePrefix("Rythmos::ExplicitRKStepper");
  //out->setShowLinePrefix(true);
  //out->setTabIndentStr("    ");

  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = model;
  t_ = ST::zero();
  solution_vector_ = model_->getNominalValues().get_x()->clone_v();
  stages_ = 4; // 4 stage ERK
  k_vector_.reserve(stages_);
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
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

}

template<class Scalar>
ExplicitRKStepper<Scalar>::ExplicitRKStepper()
{
}

template<class Scalar>
ExplicitRKStepper<Scalar>::~ExplicitRKStepper()
{
}

template<class Scalar>
Scalar ExplicitRKStepper<Scalar>::TakeStep(Scalar dt, StepSizeType flag)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag ScalarMag;
  if ((flag == VARIABLE_STEP) || (dt == ST::zero())) {
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
const StepStatus<Scalar> ExplicitRKStepper<Scalar>::getStepStatus()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  StepStatus<Scalar> stepStatus;

  stepStatus.stepSize = dt_;
  stepStatus.order = -1;
  stepStatus.time = t_;
  stepStatus.solution = solution_vector_;

  return(stepStatus);
}

template<class Scalar>
std::string ExplicitRKStepper<Scalar>::description() const
{
  std::string name = "Rythmos::ExplicitRKStepper";
  return(name);
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
    out << description() << "::describe" << std::endl;
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
bool ExplicitRKStepper<Scalar>::SetPoints(
    const std::vector<Scalar>& time_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
    ,const std::vector<ScalarMag> & accuracy_vec 
    )
{
  return(false);
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::GetPoints(
    const std::vector<Scalar>& time_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
    ,std::vector<ScalarMag>* accuracy_vec) const
{
  return(false);
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::SetRange(
    const Scalar& time_lower
    ,const Scalar& time_upper
    ,const InterpolationBufferBase<Scalar>& IB)
{
  return(false);
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::GetNodes(std::vector<Scalar>* time_vec) const
{
  return(false);
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::RemoveNodes(std::vector<Scalar>& time_vec) 
{
  return(false);
}

template<class Scalar>
int ExplicitRKStepper<Scalar>::GetOrder() const
{
  return(4);
}

template <class Scalar>
void ExplicitRKStepper<Scalar>::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  parameterList_ = paramList;
  int outputLevel = parameterList_->get( "outputLevel", int(-1) );
  outputLevel = min(max(outputLevel,-1),4);
  this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> ExplicitRKStepper<Scalar>::getParameterList()
{
  return(parameterList_);
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> ExplicitRKStepper<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
void ExplicitRKStepper<Scalar>::setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model)
{
  TEST_FOR_EXCEPT(model == Teuchos::null)
  model_ = model;
}

} // namespace Rythmos

#endif //Rythmos_ExplicitRK_STEPPER_H
