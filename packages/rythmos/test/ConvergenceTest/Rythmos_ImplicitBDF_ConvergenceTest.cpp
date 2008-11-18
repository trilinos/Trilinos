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


#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_ImplicitBDF_ConvergenceTest.hpp"

namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;

// 08/19/08 tscoffe:  These tests won't work very well until I can handle the
// fact that the first few steps will be at a lower order than what I specified by
// minorder.  Currently, I'm imagining adding a few ramp-up steps of a smaller
// size proportionate to the disparity between minOrder and currentOrder.
// Basically, if you're step-size is h (h<1), your desired order is p, and
// you're current order is q, then instead of taking one step of size h,
// take h^(q/p) steps of size h^(p/q).  That way your error, which is normally
// O(h^q), will be O((h^(p/q))^q) = O(h^p).
//
// So for now, just do order 1 and 2, as they seem to pass.
TEUCHOS_UNIT_TEST( Rythmos_ImplicitBDFStepper, GlobalErrorConvergenceStudy ) {
  RCP<SinCosModelFactory> modelFactory = sinCosModelFactory(true);
  RCP<SinCosModelExactSolutionObject> exactSolution = sinCosModelExactSolutionObject(modelFactory);
  RCP<ImplicitBDFStepperFactory<double> > stepperFactory = implicitBDFStepperFactory<double>(modelFactory);
  StepperFactoryAndExactSolutionObject<double> stepperFactoryAndExactSolution(stepperFactory,exactSolution);

  int N = stepperFactory->maxOrder();
  for (int order=1 ; order<=N ; ++order) {
    stepperFactory->setOrder(order);
    double slope = computeOrderByGlobalErrorConvergenceStudy(stepperFactoryAndExactSolution);
    double tol = 1.0e-1;
    if (order == 1) { tol = 2.0e-2; }
    TEST_FLOATING_EQUALITY( slope, 1.0*order, tol );
  }
}

/*
// 08/19/08 tscoffe: We can't do a local error convergence study on a
// multi-step method as it will always take its first step at first order
TEUCHOS_UNIT_TEST( Rythmos_ImplicitBDFStepper, LocalErrorConvergenceStudy ) {
  SinCosModelIBDFStepperFactory ibdfFactory;
  for (int order=1 ; order<=5 ; ++order) {
    ibdfFactory.setOrder(order);
    double slope = computeOrderByLocalErrorConvergenceStudy(ibdfFactory);
    double tol = 1.0e-1;
    int localOrder = order+1;
    TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol );
  }
}
*/



} // namespace Rythmos

