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

#include "Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "NonlinearModel.hpp"

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
    ExplicitRK(const Teuchos::RefCountPtr<const NonlinearModel<Scalar> > &model_);
    
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

    Teuchos::RefCountPtr<const NonlinearModel<Scalar> > model_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > solution_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > residual_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > k1_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > k2_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > k3_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > k4_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > ktemp_vector_;

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
ExplicitRK<Scalar>::ExplicitRK(const Teuchos::RefCountPtr<const Rythmos::NonlinearModel<Scalar> > &model)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = model;
  t_ = ST::zero();
  solution_vector_ = (*model_).get_vector();
  residual_vector_ = (*model_).get_vector();
  k1_vector_ = model_->get_vector();
  k2_vector_ = model_->get_vector();
  k3_vector_ = model_->get_vector();
  k4_vector_ = model_->get_vector();
  ktemp_vector_ = model_->get_vector();
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

//#include<vector>
//  std::vector<std::vector<Scalar> > butcher_tableau;
//  std::vector<Scalar> butcher_b;
//  std::vector<Scalar> butcher_c;
  
  // Runge-Kutta methods in Butcher tableau form:
  //  c | A    A = sxs matrix
  //   -|---   b = s vector
  //      b'   c = s vector
  
  // "The" Runge-Kutta Method:
  // c = [  0  1/2 1/2  1  ]'
  // A = [  0              ] 
  //     [ 1/2  0          ]
  //     [  0  1/2  0      ]
  //     [  0   0   1   0  ]
  // b = [ 1/6 1/3 1/3 1/6 ]'
//  int s = 4;
//  butcher_tableau.reserve(4);
//  butcher_b.reserve(s);
//  butcher_c.reserve(s);
//  Scalar onesixth = ST::one()/(6*ST::one());
//  Scalar onethird = ST::one()/(3*ST::one());
//  butcher_b.push_back(onesixth);
//  butcher_b.push_back(onethird);
//  butcher_b.push_back(onethird);
//  butcher_b.push_back(onesixth);
//  Scalar onehalf = ST::one()/(2*ST::one());
//  butcher_c.push_back(ST::zero());
//  butcher_c.push_back(onehalf);
//  butcher_c.push_back(onehalf);
//  butcher_c.push_back(ST::one());

  // 3/8 Rule Runge-Kutta Method:
  // c = [  0  1/3 2/3  1  ]'
  // A = [  0              ]
  //     [ 1/3  0          ]
  //     [-1/3  1   0      ]
  //     [  1  -1   1   0  ]
  // b = [ 1/8 3/8 3/8 1/8 ]'


  // k1:
  inargs.set_x(solution_vector_);
  double t1 = t_ + ST::zero();
  inargs.set_t(t1);
  outargs.request_F(k1_vector_);
  model_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k1_vector_,dt); // k1 = k1*dt
  // k2:
  Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
  Thyra::Vp_StV(&*ktemp_vector_, 0.5*dt, *k1_vector_); // ktemp = ktemp + k1*0.5
  inargs.set_x(ktemp_vector_);
  double t2 = t_ + 0.5*dt;
  inargs.set_t(t2);
  outargs.request_F(k2_vector_);
  model_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k2_vector_,dt); // k2 = k2*dt
  // k3:
  Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
  Thyra::Vp_StV(&*ktemp_vector_, 0.5*dt, *k2_vector_); // ktemp = ktemp + k2*0.5
  inargs.set_x(ktemp_vector_);
  double t3 = t_ + 0.5*dt;
  inargs.set_t(t3);
  outargs.request_F(k3_vector_);
  model_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k3_vector_,dt); // k3 = k3*dt
  // k4:
  Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
  Thyra::Vp_StV(&*ktemp_vector_, 1.0*dt, *k3_vector_); // ktemp = ktemp + k3*1.0
  inargs.set_x(ktemp_vector_);
  double t4 = t_ + dt;
  inargs.set_t(t4);
  outargs.request_F(k4_vector_);
  model_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k4_vector_,dt); // k4 = k4*dt
  // Sum for solution:
  Thyra::Vp_StV(&*solution_vector_, 1.0/6.0, *k1_vector_); // solution_vector += (1/6)*k1
  Thyra::Vp_StV(&*solution_vector_, 1.0/3.0, *k2_vector_); // solution_vector += (1/3)*k2
  Thyra::Vp_StV(&*solution_vector_, 1.0/3.0, *k3_vector_); // solution_vector += (1/3)*k3
  Thyra::Vp_StV(&*solution_vector_, 1.0/6.0, *k4_vector_); // solution_vector += (1/6)*k4

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
