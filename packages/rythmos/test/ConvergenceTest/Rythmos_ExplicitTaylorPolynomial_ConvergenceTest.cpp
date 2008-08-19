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

#include "Rythmos_ExplicitTaylorPolynomial_ConvergenceTest.hpp"

namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;

/*
TEUCHOS_UNIT_TEST( Rythmos_ExplicitTaylorPolynomialStepper, GlobalErrorConvergenceStudy ) {
  int order = 1;
  SinCosModelETPStepperFactory etpFactory;
  double slope = computeOrderByGlobalErrorConvergenceStudy(etpFactory);
  double tol = 1.0e-10;
  TEST_FLOATING_EQUALITY( slope, 1.0*order, tol );
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitTaylorPolynomialStepper, LocalErrorConvergenceStudy ) {
  int order = 1;
  SinCosModelETPStepperFactory etpFactory;
  double slope = computeOrderByLocalErrorConvergenceStudy(etpFactory);
  double tol = 1.0e-10;
  int localOrder = order+1;
  TEST_FLOATING_EQUALITY( slope, 1.0*localOrder, tol );
}
*/



} // namespace Rythmos

