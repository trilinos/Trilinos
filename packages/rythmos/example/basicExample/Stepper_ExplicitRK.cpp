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

namespace Rythmos {

//-----------------------------------------------------------------------------
// Function      : ExplicitRK::ExplicitRK
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 06/09/05
//-----------------------------------------------------------------------------
ExplicitRK::ExplicitRK(const Teuchos::RefCountPtr<Rythmos::NonlinearModel<Scalar> > &model)
{
  model_ = model;
  t_ = 0.0; // this should come from Scalar class
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
ExplicitRK::ExplicitRK()
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
ExplicitRK::~ExplicitRK()
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
Scalar ExplicitRK::TakeStep()
{
  // print something out about this method not supporting automatic variable step-size
  return(-1);
}

//-----------------------------------------------------------------------------
// Function      : ExplicitRK::TakeStep
// Purpose       : Take a step no larger than dt
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
Scalar ExplicitRK::TakeStep(Scalar dt)
{
  InArgs<Scalar> inargs;
  OutArgs<Scalar> outargs;

  // Fourth order Runge-Kutta 
  // k1:
  inargs.set_x(solution_vector_);
  double t1 = t_ + 0.0;
  inargs.set_t(t1);
  outargs.set_F(k1_vector_);
  problem_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k1_vector_,dt); // k1 = k1*dt
  // k2:
  Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
  Thyra::Vp_StV(&*ktemp_vector_, 0.5, *k1_vector_); // ktemp = ktemp + k1*0.5
  inargs.set_x(ktemp_vector_);
  double t2 = t_ + 0.5*dt;
  inargs.set_t(t2);
  outargs.set_F(k2_vector_);
  problem_->evalMoedl(inargs,outargs);
  Thyra::Vt_S(&*k2_vector_,dt); // k2 = k2*dt
  // k3:
  Thyra::assign(&*ktemp_vector_, *solution_vector_); // ktemp = solution_vector
  Thyra::Vp_StV(&*ktemp_vector_, 0.5, *k2_vector_); // ktemp = ktemp + k1*0.5
  inargs.set_x(ktemp_vector_);
  double t3 = t_ + 0.5*dt;
  inargs.set_t(t3);
  outargs.set_F(k3_vector_);
  problem_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k3_vector_,dt); // k3 = k3*dt
  // k4:
  Thyra::assign(&*ktemp_vector_, *solution_vector); // ktemp = solution_vector
  Thyra::Vp_StV(&*ktemp_vector_, 1.0, k3_vector_); // ktemp = ktemp + k3
  inargs.set_x(ktemp_vector_);
  double t4 = t_ + dt;
  inargs.set_t(t4);
  outargs.set_F(k4_vector_);
  problem_->evalModel(inargs,outargs);
  Thyra::Vt_S(&*k4_vector_,dt); // k4 = k4*dt
  // Sum for solution:
  Thyra::Vp_StV(&*solution_vector_, 1/6, k1_vector_); // solution_vector += (1/6)*k1
  Thyra::Vp_StV(&*solution_vector_, 1/3, k2_vector_); // solution_vector += (1/3)*k2
  Thyra::Vp_StV(&*solution_vector_, 1/3, k3_vector_); // solution_vector += (1/3)*k3
  Thyra::Vp_StV(&*solution_vector_, 1/6, k4_vector_); // solution_vector += (1/6)*k4

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
const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &Forward_Euler::get_solution()
{
  return(solution_vector_);
};

//-----------------------------------------------------------------------------
// Function      : ExplicitRK::get_residual
// Purpose       : return current residual
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 06/07/05
//-----------------------------------------------------------------------------
const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &Forward_Euler::get_residual()
{
  // Since this RK method doesn't actually compute the residual, if you want
  // it, then I need to compute it on the fly.  Ideally, we'd keep track of
  // whether we've already evaluated it so we don't keep hitting the nonlinear
  // problem.
  InArgs<Scalar> inargs;
  OutArgs<Scalar> outargs;
  
  inargs.set_x(solution_vector_);
  inargs.set_t(t_);
  outargs.set_F(residual_vector_);
  problem_->evalModel(inargs,outargs);
  return(residual_vector_);
};

} // namespace Rythmos
