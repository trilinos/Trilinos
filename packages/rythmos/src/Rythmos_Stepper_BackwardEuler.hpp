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

#ifndef Rythmos_STEPPER_BACKWARDEULER_H
#define Rythmos_STEPPER_BACKWARDEULER_H

#include "Rythmos_Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "Rythmos_ModelEvaluator.hpp"

namespace Rythmos {

template<class Scalar>
class BackwardEuler : public Stepper<Scalar>
{
  public:
    
    // Constructor
    BackwardEuler() {};
    BackwardEuler(const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > &model);
    
    // Destructor
    ~BackwardEuler() {};

    // Take a step _no larger_ than dt 
    Scalar TakeStep(Scalar dt);
   
    // Take a step 
    Scalar TakeStep();

    // Get solution vector
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_solution() const;

  private:

    Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > model_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > solution_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > residual_vector_;
    Scalar t_;

};


template<class Scalar>
BackwardEuler<Scalar>::BackwardEuler(const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> > &model)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = model;
  t_ = ST::zero();
  solution_vector_ = model_->get_vector();
  residual_vector_ = model_->get_vector();
}

template<class Scalar>
Scalar BackwardEuler<Scalar>::TakeStep()
{
  // print something out about this method not supporting automatic variable step-size
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return(-ST::one());
}

template<class Scalar>
Scalar BackwardEuler<Scalar>::TakeStep(Scalar dt)
{
  InArgs<Scalar> inargs;
  OutArgs<Scalar> outargs;

  // Currently this is exactly the same as ForwardEuler
  // Basically we need to write a new residual function for the nonlinear solver of the form:
  // f(x) = x(t+dt)-x(t)-dt*model(x(t+dt),t+dt) = 0
  // with Jacobian:
  // f'(x) = I - dt*model'(x(t+dt),t+dt)
  //
  inargs.set_x(solution_vector_);
  inargs.set_t(t_+dt);

  outargs.request_F(residual_vector_);

  model_->evalModel(inargs,outargs);

  // solution_vector = solution_vector + dt*residual_vector
  Thyra::Vp_StV(&*solution_vector_,dt,*residual_vector_); 
  t_ += dt;

  return(dt);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > BackwardEuler<Scalar>::get_solution() const
{
  return(solution_vector_);
}



} // namespace Rythmos

#endif //Rythmos_STEPPER_BACKWARDEULER_H
