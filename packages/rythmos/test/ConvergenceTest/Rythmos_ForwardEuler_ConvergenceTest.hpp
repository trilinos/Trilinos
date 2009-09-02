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
#ifndef Rythmos_FORWARD_EULER_CONVERGENCETEST_H
#define Rythmos_FORWARD_EULER_CONVERGENCETEST_H

#include "Rythmos_Types.hpp"
#include "Rythmos_ConvergenceTestHelpers.hpp"
#include "Rythmos_ForwardEulerStepper.hpp"

namespace Rythmos {

template<class Scalar>
class ForwardEulerStepperFactory : public virtual StepperFactoryBase<Scalar>
{
  public:
    ForwardEulerStepperFactory(RCP<ModelFactoryBase<Scalar> > modelFactory) 
    { 
      modelFactory_ = modelFactory;
    }
    virtual ~ForwardEulerStepperFactory() {}
    RCP<StepperBase<Scalar> > getStepper() const 
    { 
      RCP<ModelEvaluator<Scalar> > model = modelFactory_->getModel();
      RCP<ForwardEulerStepper<Scalar> > stepper = forwardEulerStepper<Scalar>(model);
      Thyra::ModelEvaluatorBase::InArgs<Scalar> ic = model->getNominalValues();
      stepper->setInitialCondition(ic);
      return stepper;
    }
  private:
    RCP<ModelFactoryBase<Scalar> > modelFactory_;
};
// non-member constructor
template<class Scalar>
RCP<ForwardEulerStepperFactory<Scalar> > forwardEulerStepperFactory(
    RCP<ModelFactoryBase<Scalar> > modelFactory)
{
  RCP<ForwardEulerStepperFactory<Scalar> > feFactory = Teuchos::rcp(
      new ForwardEulerStepperFactory<Scalar>(modelFactory)
      );
  return feFactory;
}

} // namespace Rythmos 

#endif // Rythmos_FORWARD_EULER_CONVERGENCETEST_H

