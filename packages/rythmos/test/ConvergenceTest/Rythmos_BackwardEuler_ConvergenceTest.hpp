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
#ifndef Rythmos_BACKWARD_EULER_CONVERGENCETEST_H
#define Rythmos_BACKWARD_EULER_CONVERGENCETEST_H

#include "Rythmos_Types.hpp"
#include "Rythmos_ConvergenceTestHelpers.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"

namespace Rythmos {

template<class Scalar>
class BackwardEulerStepperFactory : public virtual StepperFactoryBase<Scalar>
{
  public:
    BackwardEulerStepperFactory(RCP<ModelFactoryBase<Scalar> > modelFactory) 
    { 
      modelFactory_ = modelFactory;
    }
    virtual ~BackwardEulerStepperFactory() {}
    RCP<StepperBase<Scalar> > getStepper() const 
    { 
      RCP<ModelEvaluator<Scalar> > model = modelFactory_->getModel();
      Thyra::ModelEvaluatorBase::InArgs<double> model_ic = model->getNominalValues();
      RCP<Rythmos::TimeStepNonlinearSolver<Scalar> >
        nonlinearSolver = Rythmos::timeStepNonlinearSolver<Scalar>();
      RCP<ParameterList> nonlinearSolverPL = Teuchos::parameterList();
      nonlinearSolverPL->get("Default Tol",1.0e-9); // Set default if not set
      nonlinearSolver->setParameterList(nonlinearSolverPL);
      RCP<BackwardEulerStepper<Scalar> > stepper = Teuchos::rcp(new BackwardEulerStepper<Scalar>(model,nonlinearSolver));
      stepper->setInitialCondition(model_ic);
      return stepper;
    }
  private:
    RCP<ModelFactoryBase<Scalar> > modelFactory_;
};
// non-member constructor
template<class Scalar>
RCP<BackwardEulerStepperFactory<Scalar> > backwardEulerStepperFactory(
    RCP<ModelFactoryBase<Scalar> > modelFactory)
{
  RCP<BackwardEulerStepperFactory<Scalar> > beFactory = Teuchos::rcp(
      new BackwardEulerStepperFactory<Scalar>(modelFactory)
      );
  return beFactory;
}

} // namespace Rythmos 

#endif // Rythmos_BACKWARD_EULER_CONVERGENCETEST_H

