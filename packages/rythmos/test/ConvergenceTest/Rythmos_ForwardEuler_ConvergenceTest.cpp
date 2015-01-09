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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_ForwardEuler_ConvergenceTest.hpp"

namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;

TEUCHOS_UNIT_TEST( Rythmos_ForwardEulerStepper, GlobalErrorConvergenceStudy ) {
  RCP<SinCosModelFactory> modelFactory = sinCosModelFactory(false);
  RCP<SinCosModelExactSolutionObject> exactSolution = sinCosModelExactSolutionObject(modelFactory);
  RCP<ForwardEulerStepperFactory<double> > stepperFactory = forwardEulerStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(stepperFactory,exactSolution);

  double slope = computeOrderByGlobalErrorConvergenceStudy(stepperFactoryAndExactSolution);
  int order = stepperFactoryAndExactSolution.getStepper()->getOrder();
  double tol = 2.0e-2;
  TEST_FLOATING_EQUALITY( slope, 1.0*order, tol );
}

TEUCHOS_UNIT_TEST( Rythmos_ForwardEulerStepper, LocalErrorConvergenceStudy ) {
  RCP<SinCosModelFactory> modelFactory = sinCosModelFactory(false);
  RCP<SinCosModelExactSolutionObject> exactSolution = sinCosModelExactSolutionObject(modelFactory);
  RCP<ForwardEulerStepperFactory<double> > stepperFactory = forwardEulerStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(stepperFactory,exactSolution);

  double slope = computeOrderByLocalErrorConvergenceStudy(stepperFactoryAndExactSolution);
  int order = stepperFactoryAndExactSolution.getStepper()->getOrder();
  double tol = 1.0e-1;
  int localOrder = order+1;
  TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol );
}



} // namespace Rythmos

