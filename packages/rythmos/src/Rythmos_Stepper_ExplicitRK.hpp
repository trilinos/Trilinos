//
// @HEADER
// ***********************************************************************
// 
//                           Rythmos Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef Rythmos_STEPPER_ExplicitRK_H
#define Rythmos_STEPPER_ExplicitRK_H

#include "Rythmos_Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "Rythmos_ModelEvaluator.hpp"
#include <vector>

namespace Rythmos {

//-----------------------------------------------------------------------------
// Class         : ExplicitRK
// Purpose       : Define explicit Runge-Kutta methods
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 06/09/05
//-----------------------------------------------------------------------------
template<class Scalar>
class ExplicitRK : public Stepper<Scalar>
{
  public:
    
    // Constructor
    ExplicitRK();
    ExplicitRK(const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > &model_);
    
    // Destructor
    ~ExplicitRK();

    // Take a step _no larger_ than dt 
    Scalar TakeStep(Scalar dt);
   
    // Take a step 
    Scalar TakeStep();

    // Get solution vector
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_solution() const;

    // Get residual vector
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_residual() const;

  protected:

    Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > model_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > solution_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > residual_vector_;
    std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > k_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > ktemp_vector_;

    std::vector<std::vector<Scalar> > b_A; // Butcher tableau A matrix
    std::vector<Scalar> b_b; // Butcher b vector
    std::vector<Scalar> b_c; // Butcher c vector

    Scalar t_;

};

//-----------------------------------------------------------------------------
// Function      : ExplicitRK::ExplicitRK
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 06/09/05
//-----------------------------------------------------------------------------
template<class Scalar>
ExplicitRK<Scalar>::ExplicitRK(const Teuchos::RefCountPtr<const Rythmos::ModelEvaluator<Scalar> > &model)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = model;
  t_ = ST::zero();
  solution_vector_ = (*model_).get_vector();
  residual_vector_ = (*model_).get_vector();
  int stages = 4; // 4 stage ERK
  k_vector_.reserve(stages);
  for (int i=0 ; i<stages ; ++i)
  {
    k_vector_.push_back(model_->get_vector());
  }
  ktemp_vector_ = model_->get_vector();

  // Runge-Kutta methods in Butcher tableau form:
  //  c | A    A = sxs matrix
  //   -|---   b = s vector
  //      b'   c = s vector
  
  // 3/8 Rule Runge-Kutta Method: (not implemented yet)
  // c = [  0  1/3 2/3  1  ]'
  // A = [  0              ]
  //     [ 1/3  0          ]
  //     [-1/3  1   0      ]
  //     [  1  -1   1   0  ]
  // b = [ 1/8 3/8 3/8 1/8 ]'

  
  // "The" Runge-Kutta Method: (implemented below)
  // c = [  0  1/2 1/2  1  ]'
  // A = [  0              ] 
  //     [ 1/2  0          ]
  //     [  0  1/2  0      ]
  //     [  0   0   1   0  ]
  // b = [ 1/6 1/3 1/3 1/6 ]'
  Scalar zero = ST::zero();
  Scalar one = ST::one();
  Scalar onehalf = ST::one()/(2*ST::one());
  Scalar onesixth = ST::one()/(6*ST::one());
  Scalar onethird = ST::one()/(3*ST::one());
  b_A.reserve(stages);
  b_b.reserve(stages);
  b_c.reserve(stages);

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
  for (int i=0 ; i<stages ; ++i)
  {
    b_A.push_back(b_b);
  }

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
//
//-----------------------------------------------------------------------------
// Function      : ExplicitRK::ExplicitRK
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
template<class Scalar>
ExplicitRK<Scalar>::ExplicitRK()
{
}

//-----------------------------------------------------------------------------
// Function      : ~ExplicitRK::ExplicitRK
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
template<class Scalar>
ExplicitRK<Scalar>::~ExplicitRK()
{
}

//-----------------------------------------------------------------------------
// Function      : ExplicitRK::TakeStep
// Purpose       : Take a step 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
template<class Scalar>
Scalar ExplicitRK<Scalar>::TakeStep()
{
  // print something out about this method not supporting automatic variable step-size
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return(-ST::one());
}

//-----------------------------------------------------------------------------
// Function      : ExplicitRK::TakeStep
// Purpose       : Take a step no larger than dt
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
template<class Scalar>
Scalar ExplicitRK<Scalar>::TakeStep(Scalar dt)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  InArgs<Scalar> inargs;
  OutArgs<Scalar> outargs;

  // k_1:
  inargs.set_x(solution_vector_);
  double t1 = t_ + b_c[1-1]*dt;
  inargs.set_t(t1);
  outargs.request_F(k_vector_[1-1]);
  model_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k_vector_[1-1],dt); // k_1 = k_1*dt
  // k_2:
  Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
  Thyra::Vp_StV(&*ktemp_vector_, b_A[2-1][1-1], *k_vector_[1-1]); // ktemp = ktemp + a_{2,1}*k_1
  inargs.set_x(ktemp_vector_);
  double t2 = t_ + b_c[2-1]*dt;
  inargs.set_t(t2);
  outargs.request_F(k_vector_[2-1]);
  model_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k_vector_[2-1],dt); // k_2 = k_2*dt
  // k_3:
  Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
  Thyra::Vp_StV(&*ktemp_vector_, b_A[3-1][1-1], *k_vector_[1-1]); // ktemp = ktemp + a_{3,1}*k_1
  Thyra::Vp_StV(&*ktemp_vector_, b_A[3-1][2-1], *k_vector_[2-1]); // ktemp = ktemp + a_{3,2}*k_2
  inargs.set_x(ktemp_vector_);
  double t3 = t_ + b_c[3-1]*dt;
  inargs.set_t(t3);
  outargs.request_F(k_vector_[3-1]);
  model_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k_vector_[3-1],dt); // k_3 = k_3*dt
  // k_4:
  Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
  Thyra::Vp_StV(&*ktemp_vector_, b_A[4-1][1-1], *k_vector_[1-1]); // ktemp = ktemp + a_{4,1}*k_1
  Thyra::Vp_StV(&*ktemp_vector_, b_A[4-1][2-1], *k_vector_[2-1]); // ktemp = ktemp + a_{4,2}*k_2
  Thyra::Vp_StV(&*ktemp_vector_, b_A[4-1][3-1], *k_vector_[3-1]); // ktemp = ktemp + a_{4,3}*k_3
  inargs.set_x(ktemp_vector_);
  double t4 = t_ + b_c[4-1]*dt;
  inargs.set_t(t4);
  outargs.request_F(k_vector_[4-1]);
  model_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k_vector_[4-1],dt); // k_4 = k_4*dt
  // Sum for solution:
  Thyra::Vp_StV(&*solution_vector_, b_b[1-1], *k_vector_[1-1]); // solution_vector += b_1*k_1
  Thyra::Vp_StV(&*solution_vector_, b_b[2-1], *k_vector_[2-1]); // solution_vector += b_2*k_2
  Thyra::Vp_StV(&*solution_vector_, b_b[3-1], *k_vector_[3-1]); // solution_vector += b_3*k_3
  Thyra::Vp_StV(&*solution_vector_, b_b[4-1], *k_vector_[4-1]); // solution_vector += b_4*k_4

  // update current time:
  t_ = t_ + dt;

  return(dt);
}

//-----------------------------------------------------------------------------
// Function      : ExplicitRK::get_solution
// Purpose       : return current solution
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 06/07/05
//-----------------------------------------------------------------------------
template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ExplicitRK<Scalar>::get_solution() const
{
  return(solution_vector_);
}

//-----------------------------------------------------------------------------
// Function      : ExplicitRK::get_residual
// Purpose       : return current residual
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 06/07/05
//-----------------------------------------------------------------------------
template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ExplicitRK<Scalar>::get_residual() const
{
  // Since this RK method doesn't actually compute the residual, if you want
  // it, then I need to compute it on the fly.  Ideally, we'd keep track of
  // whether we've already evaluated it so we don't keep hitting the nonlinear
  // problem.
  InArgs<Scalar> inargs;
  OutArgs<Scalar> outargs;
  
  inargs.set_x(solution_vector_);
  inargs.set_t(t_);
  outargs.request_F(residual_vector_);
  model_->evalModel(inargs,outargs);
  return(residual_vector_);
}

} // namespace Rythmos

#endif //Rythmos_STEPPER_ExplicitRK_H
