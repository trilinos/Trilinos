//@HEADER

// ***********************************************************************
//
//                     Rythmos Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_STEPPER_FORWARDEULER
#define RYTHMOS_STEPPER_FORWARDEULER

#include "Stepper.hpp"
#include "ModelEvaluator.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Rythmos {

template<class Scalar>
class ForwardEulerStepper : public Stepper<Scalar>
{
  public:
    ForwardEulerStepper();
    ForwardEulerStepper(ModelEvaluator<Scalar> *model);
    ~ForwardEulerStepper();
    Scalar takeStep(Scalar dt);
    Scalar get_solution();
  protected:
    Scalar t_;
    Scalar x_;
    Scalar f_;
    ModelEvaluator<Scalar> *model_;
};


template<class Scalar>
ForwardEulerStepper<Scalar>::ForwardEulerStepper() 
{
}
template<class Scalar>
ForwardEulerStepper<Scalar>::ForwardEulerStepper(ModelEvaluator<Scalar> *model)
{
  model_ = model;
  t_ = Teuchos::ScalarTraits<Scalar>::zero();
  x_ = model_->get_vector();
}
template<class Scalar>
ForwardEulerStepper<Scalar>::~ForwardEulerStepper() 
{
}
template<class Scalar>
Scalar ForwardEulerStepper<Scalar>::takeStep(Scalar dt)
{
  f_ = model_->evalModel(x_,t_);
  x_ = x_ + dt*f_;
  t_ = t_ + dt;
  return(dt);
}
template<class Scalar>
Scalar ForwardEulerStepper<Scalar>::get_solution()
{
  return(x_);
}

} // namespace Rythmos


#endif // RYTHMOS_STEPPER_FORWARDEULER
