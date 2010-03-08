/*
// @HEADER
// ***********************************************************************
// 
//    GlobiPack: Collection of Scalar 1D globalizaton utilities
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "GlobiPack_TestLagrPolyMeritFunc1D.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Tuple.hpp"


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::outArg;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::tuple;
using GlobiPack::TestLagrPolyMeritFunc1D;
using GlobiPack::computeValue;


template<class Scalar>
inline Scalar sqr(const Scalar &x) { return x*x; }


double g_tol = Teuchos::ScalarTraits<double>::eps()*100.0;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "tol", &g_tol, "Floating point tolerance" );
}


//
// Unit tests
//


//
// Here we represent the simple quadratic merit function:
//
//   phi(alpha) = 0.5 * (alpha - 2.0)^2 + 2.0
//
// with the derivative:
//
//   Dphi(alpha) = alpha - 2.0
//

template<class Scalar>
Scalar phi_quad1(const Scalar &alpha)
{
  return as<Scalar>(0.5)*sqr(alpha-as<Scalar>(2.0)) + as<Scalar>(2.0);
}

template<class Scalar>
Scalar Dphi_quad1(const Scalar &alpha)
{
  return alpha-as<Scalar>(2.0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TestLagrPolyMeritFunc1D, basic, Scalar )
{
  
  typedef Teuchos::ScalarTraits<Scalar> ST;
  
  ECHO(Array<Scalar> alphaPoints = tuple<Scalar>(0.0, 2.0, 4.0));
  ECHO(Array<Scalar> phiPoints = tuple<Scalar>(4.0, 2.0, 4.0));
  ECHO(TestLagrPolyMeritFunc1D<Scalar> meritFunc(alphaPoints, phiPoints));
  TEST_ASSERT(meritFunc.supportsDerivEvals());

  Array<Scalar> alphaTestPoints = tuple<Scalar>(0.0, 1.0, 2.0, 3.0, 4.0);
  for (int test_i = 0; test_i < as<int>(alphaTestPoints.size()); ++test_i) {
    out << "\ntest_i="<<test_i<<"\n\n";
    Teuchos::OSTab tab(out);
    ECHO(const Scalar alpha = alphaTestPoints[test_i]);
    out << "alpha="<<alpha<<"\n";
    ECHO(Scalar phi = as<Scalar>(-1.0));
    ECHO(Scalar Dphi = as<Scalar>(-1.0));
    ECHO(meritFunc.eval(alpha, outArg(phi), outArg(Dphi)));
    TEST_FLOATING_EQUALITY(phi, phi_quad1(alpha), g_tol);
    TEST_FLOATING_EQUALITY(Dphi, Dphi_quad1(alpha), g_tol);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( TestLagrPolyMeritFunc1D, basic )


} // namespace


