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
#include "GlobiPack_BrentsLineSearch.hpp"
#include "Teuchos_Tuple.hpp"

#include "meritFuncsHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


//
// Helper code and declarations
//


using GlobiPack::BrentsLineSearch;
using GlobiPack::brentsLineSearch;
using GlobiPack::computeValue;
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


double g_tol_scale = 100.0;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "tol", &g_tol_scale, "Floating point tolerance scaling of eps." );
}


//
// Unit tests for BrentsLineSearch
//


//
// Check that object can exactly interplate a quadratic merit function. This
// takes more than one iteration because of the golden search stuff.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BrentsLineSearch, quadExact, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>();
  
  RCP<BrentsLineSearch<Scalar> > linesearch = brentsLineSearch<Scalar>();

  linesearch->setOStream(rcpFromRef(out));

  const PointEval1D<Scalar> point_k = computePoint(*phi, ST::zero());
  PointEval1D<Scalar> point_kp1 = computePoint(*phi, as<Scalar>(8.0));

  int numIters = -1;
  const bool linesearchResult = linesearch->doLineSearch(
    *phi, point_k, inOutArg(point_kp1), outArg(numIters) );

  TEST_ASSERT(linesearchResult);
  TEST_FLOATING_EQUALITY(point_kp1.alpha, as<Scalar>(2.0),
    as<Scalar>(g_tol_scale*ST::squareroot(ST::eps())));
  TEST_FLOATING_EQUALITY(point_kp1.phi, as<ScalarMag>(3.0),
    as<Scalar>(g_tol_scale)*ST::eps());
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( BrentsLineSearch, quadExact )


//
// Check that object can approximately mimimize a cubic function.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( BrentsLineSearch, cubicApprox, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>();
  
  RCP<BrentsLineSearch<Scalar> > linesearch = brentsLineSearch<Scalar>();

  const RCP<ParameterList> pl = parameterList();
  pl->sublist("Minimize").set("Relative Tol", as<double>(g_tol_scale*ST::eps()));
  pl->sublist("Minimize").set("Bracket Tol", as<double>(ST::eps()));
  linesearch->setParameterList(pl);
  
  linesearch->setOStream(rcpFromRef(out));

  const PointEval1D<Scalar> point_k = computePoint(*phi, ST::zero());
  PointEval1D<Scalar> point_kp1 = computePoint(*phi, as<Scalar>(8.0));

  int numIters = -1;
  const bool linesearchResult = linesearch->doLineSearch(
    *phi, point_k, inOutArg(point_kp1), outArg(numIters) );

  TEST_ASSERT(linesearchResult);
  TEST_FLOATING_EQUALITY(point_kp1.alpha, as<Scalar>(2.0),
    as<Scalar>(g_tol_scale*ST::squareroot(ST::eps())));
  TEST_FLOATING_EQUALITY(point_kp1.phi, as<ScalarMag>(3.0),
    as<Scalar>(g_tol_scale)*ST::eps());

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( BrentsLineSearch, cubicApprox )


} // namespace
