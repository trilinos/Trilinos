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


#include "GlobiPack_GoldenQuadInterpBracket.hpp"
#include "GlobiPack_TestLagrPolyMeritFunc1D.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


//
// Helper code and declarations
//


using GlobiPack::TestLagrPolyMeritFunc1D;
using GlobiPack::testLagrPolyMeritFunc1D;
using GlobiPack::GoldenQuadInterpBracket;
using GlobiPack::goldenQuadInterpBracket;
using GlobiPack::PointEval1D;
using GlobiPack::computePoint;
using Teuchos::as;
using Teuchos::inOutArg;
using Teuchos::outArg;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Teuchos::Array;
using Teuchos::tuple;
using Teuchos::ParameterList;
using Teuchos::parameterList;


template<class Scalar>
inline Scalar sqr(const Scalar &x) { return x*x; }


// Set up a quadratic merit function with minimizer at alpha=2.0, phi=3.0;
template<class Scalar>
const RCP<TestLagrPolyMeritFunc1D<Scalar> > quadPhi()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  Array<Scalar> alphaPoints = tuple<Scalar>(0.0, 2.0, 4.0);
  Array<ScalarMag> phiPoints = tuple<ScalarMag>(6.0, 3.0, 6.0);
  return testLagrPolyMeritFunc1D<Scalar>(alphaPoints, phiPoints);
}


double g_tol = Teuchos::ScalarTraits<double>::eps()*100.0;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "tol", &g_tol, "Floating point tolerance" );
}


//
// Unit tests for GoldenQuadInterpBracket
//


//
// Check that object can exactly interplate a quadratic merit function at the
// very first iteration
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( GoldenQuadInterpBracket, bracket, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  using std::min;

  const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>();
  
  RCP<GoldenQuadInterpBracket<Scalar> > bracket = goldenQuadInterpBracket<Scalar>();

  bracket->setOStream(rcpFromRef(out));

  const Array<Scalar> alpha =
    tuple<Scalar>(
      1e-14, 1e-10, 1e-7, 1e-4, 0.1, 1.0, 1.1, 1.5,
      1.9, 2.0, 2.1, 4.0, 8.0, 30.0);

  for (int i = 0; i < as<int>(alpha.size()); ++i ) {

    PointEval1D<Scalar> p_l = computePoint(*phi, ST::zero());
    PointEval1D<Scalar> p_m = computePoint(*phi, alpha[i]);
    PointEval1D<Scalar> p_u;
    int numIters = -1;

    if (alpha[i] > ST::eps()) {
    
      const bool bracketResult = bracket->bracketMinimum(
        *phi, inOutArg(p_l), inOutArg(p_m), outArg(p_u), outArg(numIters) );
      
      TEST_ASSERT(bracketResult);
      TEST_COMPARE(p_l.alpha, <, p_m.alpha);
      TEST_COMPARE(p_m.alpha, <, p_u.alpha);
      TEST_COMPARE(p_l.phi, >, p_m.phi);
      TEST_COMPARE(p_m.phi, <, p_u.phi);

    }
    else {

      // If alpha[i] < eps, then there will not be enough accuracy in the
      // merit function to allow for a successfull line search.

    }

  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( GoldenQuadInterpBracket, bracket )


} // namespace
