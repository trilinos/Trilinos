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

#ifndef Rythmos_STEPPER_ExplicitRK_H
#define Rythmos_STEPPER_ExplicitRK_H

#include "Rythmos_Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"

namespace Rythmos {

/** \brief . */
template<class Scalar>
class ExplicitRKStepper : virtual public Stepper<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    
    /** \brief . */
    ExplicitRKStepper();
    ExplicitRKStepper(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model_);
    
    /** \brief . */
    ~ExplicitRKStepper();

    /** \brief . */
    Scalar TakeStep(Scalar dt);
   
    /** \brief . */
    Scalar TakeStep();

    /** \brief . */
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_solution() const;
    
    /** \brief . */
    std::string description() const;

    /** \brief . */
    std::ostream& describe(
      std::ostream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ,const std::string          leadingIndent
      ,const std::string          indentSpacer
      ) const;

    /// Redefined from InterpolationBuffer 
    /// Add points to buffer
    bool SetPoints(
      const std::vector<Scalar>& time_list
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_list
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_list);
    
    /// Get values from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_list
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_list
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_list
      ,std::vector<ScalarMag>* accuracy_list) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBuffer<Scalar> & IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_list) const;

    /// Remove interpolation nodes
    bool RemoveNodes(std::vector<Scalar>& time_list) const;

    /// Get order of interpolation
    int GetOrder() const;

  private:

    Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > solution_vector_;
    std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > k_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > ktemp_vector_;

    int stages_; // Number of stages of RK
    std::vector<std::vector<Scalar> > b_A; // Butcher tableau A matrix
    std::vector<Scalar> b_b; // Butcher b vector
    std::vector<Scalar> b_c; // Butcher c vector

    Scalar t_;

};

template<class Scalar>
ExplicitRKStepper<Scalar>::ExplicitRKStepper(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = model;
  t_ = ST::zero();
  solution_vector_ = model_->getNominalValues().get_x()->clone_v();
  stages_ = 4; // 4 stage ERK
  k_vector_.reserve(stages_);
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
    f_space = model_->get_f_space();
  for (int i=0 ; i<stages_ ; ++i)
  {
    k_vector_.push_back(Thyra::createMember(f_space));
  }
  ktemp_vector_ = Thyra::createMember(f_space);

  // initialize the Butcher tableau and its vectors
  Scalar zero = ST::zero();
  b_A.reserve(stages_);
  b_b.reserve(stages_);
  b_c.reserve(stages_);

  // fill b with zeros
  b_b.push_back(zero); 
  b_b.push_back(zero); 
  b_b.push_back(zero);
  b_b.push_back(zero);
  // fill c with zeros
  b_c.push_back(zero); 
  b_c.push_back(zero); 
  b_c.push_back(zero);
  b_c.push_back(zero);
  // fill A with zeros
  for (int i=0 ; i<stages_ ; ++i)
  {
    b_A.push_back(b_b);
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

  // fill b_A
  b_A[0][0] = zero;
  b_A[0][1] = zero;
  b_A[0][2] = zero;
  b_A[0][3] = zero;

  b_A[1][0] = one_third;
  b_A[1][1] = zero;
  b_A[1][2] = zero;
  b_A[1][3] = zero;

  b_A[2][0] = Scalar(-one_third);
  b_A[2][1] = one;
  b_A[2][2] = zero;
  b_A[2][3] = zero;

  b_A[3][0] = one;
  b_A[3][1] = Scalar(-one);
  b_A[3][2] = one;
  b_A[3][3] = zero;

  // fill b_b
  b_b[0] = one_eighth;
  b_b[1] = three_eighth;
  b_b[2] = three_eighth;
  b_b[3] = one_eighth;
  
  // fill b_c
  b_c[0] = zero;
  b_c[1] = one_third;
  b_c[2] = two_third;
  b_c[3] = one;
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

  // fill b_A
  b_A[0][0] = zero;
  b_A[0][1] = zero;
  b_A[0][2] = zero;
  b_A[0][3] = zero;

  b_A[1][0] = onehalf;
  b_A[1][1] = zero;
  b_A[1][2] = zero;
  b_A[1][3] = zero;

  b_A[2][0] = zero;
  b_A[2][1] = onehalf;
  b_A[2][2] = zero;
  b_A[2][3] = zero;

  b_A[3][0] = zero;
  b_A[3][1] = zero;
  b_A[3][2] = one;
  b_A[3][3] = zero;

  // fill b_b
  b_b[0] = onesixth;
  b_b[1] = onethird;
  b_b[2] = onethird;
  b_b[3] = onesixth;
  
  // fill b_c
  b_c[0] = zero;
  b_c[1] = onehalf;
  b_c[2] = onehalf;
  b_c[3] = one;

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
Scalar ExplicitRKStepper<Scalar>::TakeStep()
{
  // print something out about this method not supporting automatic variable step-size
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return(-ST::one());
}

template<class Scalar>
Scalar ExplicitRKStepper<Scalar>::TakeStep(Scalar dt)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag ScalarMag;

  // Compute stage solutions
  for (int s=0 ; s < stages_ ; ++s)
  {
    Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
    for (int j=0 ; j < s ; ++j) // assuming Butcher matix is strictly lower triangular
    {
      if (b_A[s][j] != ST::zero())
        Thyra::Vp_StV(&*ktemp_vector_, b_A[s][j], *k_vector_[j]); // ktemp = ktemp + a_{s+1,j+1}*k_{j+1}
    }
    ScalarMag ts = t_ + b_c[s]*dt;
    Thyra::eval_f<Scalar>(*model_,*ktemp_vector_,ts,&*k_vector_[s]);
    Thyra::Vt_S(&*k_vector_[s],dt); // k_s = k_s*dt
  
  } 
  // Sum for solution:
  for (int s=0 ; s < stages_ ; ++s)
  {
    if (b_b[s] != ST::zero())
      Thyra::Vp_StV(&*solution_vector_, b_b[s], *k_vector_[s]); // solution_vector += b_{s+1}*k_{s+1}
  }

  // update current time:
  t_ = t_ + dt;

  return(dt);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ExplicitRKStepper<Scalar>::get_solution() const
{
  return(solution_vector_);
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
  if (verbLevel == Teuchos::VERB_EXTREME)
  {
    out << description() << "::describe" << std::endl;
    out << "model_ = " << std::endl;
    out << model_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "solution_vector_ = " << std::endl;
    out << solution_vector_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    for (int i=0 ; i<stages_ ; ++i)
    {
      out << "k_vector_[" << i << "] = " << std::endl;
      out << k_vector_[i]->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    }
    out << "ktemp_vector_ = " << std::endl;
    out << ktemp_vector_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    for (int i=0 ; i<stages_ ; ++i)
      for (int j=0 ; j<stages_ ; ++j)
        out << "b_A[" << i << "][" << j << "] = " << b_A[i][j] << std::endl;
    for (int i=0 ; i<stages_ ; ++i)
      out << "b_b[" << i << "] = " << b_b[i] << std::endl;
    for (int i=0 ; i<stages_ ; ++i)
      out << "b_c[" << i << "] = " << b_c[i] << std::endl;
    out << "t_ = " << t_ << std::endl;
  }
  return(out);
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::SetPoints(
    const std::vector<Scalar>& time_list
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_list
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_list)
{
  return(false);
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::GetPoints(
    const std::vector<Scalar>& time_list
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_list
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_list
    ,std::vector<ScalarMag>* accuracy_list) const
{
  return(false);
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::SetRange(
    const Scalar& time_lower
    ,const Scalar& time_upper
    ,const InterpolationBuffer<Scalar>& IB)
{
  return(false);
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::GetNodes(std::vector<Scalar>* time_list) const
{
  return(false);
}

template<class Scalar>
bool ExplicitRKStepper<Scalar>::RemoveNodes(std::vector<Scalar>& time_list) const
{
  return(false);
}

template<class Scalar>
int ExplicitRKStepper<Scalar>::GetOrder() const
{
  return(4);
}

} // namespace Rythmos

#endif //Rythmos_STEPPER_ExplicitRK_H
