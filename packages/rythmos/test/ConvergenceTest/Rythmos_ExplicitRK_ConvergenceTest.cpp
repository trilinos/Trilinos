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

#include "Rythmos_ExplicitRK_ConvergenceTest.hpp"
#include "Rythmos_RKButcherTableau.hpp"



namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;

RCP<ExplicitRKStepper<double> > getERKStepperByOrder(int order) {
  RCP<SinCosModel> model = sinCosModel(false);
  DefaultRKButcherTableauFactory<double> rkbtFactory;
  Teuchos::ParameterList pl;
  pl.set("Selection Type", "Explicit method by order");
  pl.set("Explicit method by order", order);
  RKButcherTableau<double> rkbt = rkbtFactory.create(pl);
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model,rkbt);
  return stepper;
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, GlobalErrorConvergenceStudy ) {

  SinCosModelERKStepperFactory erkFactory;
  for (int order=1 ; order<=4 ; ++order) {
    erkFactory.setOrder(order);
    double slope = computeOrderByGlobalErrorConvergenceStudy(erkFactory);

    double tol = 1.0e-2;
    TEST_FLOATING_EQUALITY( slope, 1.0*order, tol ); // is slope close to order?
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, LocalErrorConvergenceStudy ) {
  SinCosModelERKStepperFactory erkFactory;
  for (int order=1 ; order<=4 ; ++order) {
    erkFactory.setOrder(order);
    double slope = computeOrderByLocalErrorConvergenceStudy(erkFactory);

    int localOrder = order+1; // I don't know why the order is coming out one higher than it should!?!
    double tol = 1.0e-3;
    if (order == 1) { tol = 1.0e-2; }
    else if (order == 4) { tol = 1.0e-1; }
    TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol ); // is slope close to order?
  }
}



} // namespace Rythmos

