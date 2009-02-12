
#include "GlobiPack_TestLagrPolyMeritFunc1D.hpp"
#include "GlobiPack_ArmijoPolyInterpLineSearch.hpp"
#include "Teuchos_Tuple.hpp"

#include "meritFuncsHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


//
// Helper code and declarations
//


using GlobiPack::TestLagrPolyMeritFunc1D;
using GlobiPack::testLagrPolyMeritFunc1D;
using GlobiPack::ArmijoPolyInterpLineSearch;
using GlobiPack::armijoQuadraticLineSearch;
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


double g_tol = Teuchos::ScalarTraits<double>::eps()*100.0;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "tol", &g_tol, "Floating point tolerance" );
}


//
// Unit tests for ArmijoPolyInterpLineSearch
//


//
// Check that internal default parameters are set correctly
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, defaultParams, Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  namespace AQLSU = GlobiPack::ArmijoPolyInterpLineSearchUtils;
  RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>();
  TEST_EQUALITY(linesearch->eta(), as<ScalarMag>(AQLSU::eta_default));
  TEST_EQUALITY(linesearch->minFrac(), as<ScalarMag>(AQLSU::minFrac_default));
  TEST_EQUALITY(linesearch->maxFrac(), as<ScalarMag>(AQLSU::maxFrac_default));
  TEST_EQUALITY(linesearch->doMaxIters(), as<ScalarMag>(AQLSU::doMaxIters_default));
  TEST_EQUALITY(linesearch->maxIters(), as<ScalarMag>(AQLSU::maxIters_default));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, defaultParams )


//
// Check that parameter list is parsed correctly
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, parseParams, Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  namespace AQLSU = GlobiPack::ArmijoPolyInterpLineSearchUtils;
  ECHO(RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>());
  const double eta = 0.99999;
  const double minFrac = 4.0;
  const double maxFrac = 5.0;
  const int minIters = 5;
  const int maxIters = 100;
  const bool doMaxIters = true;
  ECHO(const RCP<ParameterList> pl = parameterList());
  ECHO(pl->set("Armijo Slope Fraction", eta));
  ECHO(pl->set("Min Backtrack Fraction", minFrac));
  ECHO(pl->set("Max Backtrack Fraction", maxFrac));
  ECHO(pl->set("Min Num Iterations", minIters));
  ECHO(pl->set("Max Num Iterations", maxIters));
  ECHO(pl->set("Do Max Iterations", doMaxIters));
  ECHO(linesearch->setParameterList(pl));
  const Scalar tol = ST::eps();
  TEST_FLOATING_EQUALITY(linesearch->eta(), as<Scalar>(eta), tol);
  TEST_FLOATING_EQUALITY(linesearch->minFrac(), as<Scalar>(minFrac), tol);
  TEST_FLOATING_EQUALITY(linesearch->maxFrac(), as<Scalar>(maxFrac), tol);
  TEST_EQUALITY(linesearch->minIters(), minIters);
  TEST_EQUALITY(linesearch->maxIters(), maxIters);
  TEST_EQUALITY(linesearch->doMaxIters(), doMaxIters);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, parseParams )


//
// Check that the ArmijoPolyInterpLineSearch object validates its parameters
// and their values correctly.
//

/*

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, validateParams, Scalar )
{
  TEST_FOR_EXCEPT(true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, validateParams )

*/


//
// Check that object can exactly interplate a quadratic merit function at the
// very first iteration
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, quadExact, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>();
  
  RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>();

  linesearch->setOStream(rcpFromRef(out));
 
  const PointEval1D<Scalar> point_k = computePoint(*phi, ST::zero(), true, true);
  PointEval1D<Scalar> point_kp1 = computePoint(*phi, as<Scalar>(5.0));

  int numIters = -1;
  const bool linesearchResult = linesearch->doLineSearch(
    *phi, point_k, inOutArg(point_kp1), outArg(numIters) );

  TEST_ASSERT(linesearchResult);
  TEST_EQUALITY(numIters, 1);
  TEST_FLOATING_EQUALITY(point_kp1.alpha, as<Scalar>(2.0), g_tol);
  TEST_FLOATING_EQUALITY(point_kp1.phi, as<ScalarMag>(3.0), g_tol);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, quadExact )


//
// Check that object will accept the inital point passed in without doing any
// evaluations.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, noEval, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>();
  
  RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>();

  linesearch->setOStream(rcpFromRef(out));

  const Scalar alpha_k_init = as<Scalar>(2.1);
  const PointEval1D<Scalar> point_k = computePoint(*phi, ST::zero(), true, true);
  PointEval1D<Scalar> point_kp1 = computePoint(*phi, alpha_k_init);

  int numIters = -1;
  const bool linesearchResult = linesearch->doLineSearch(
    *phi, point_k, inOutArg(point_kp1), outArg(numIters) );

  TEST_ASSERT(linesearchResult);
  TEST_EQUALITY(numIters, 0);
  TEST_FLOATING_EQUALITY(point_kp1.alpha, alpha_k_init, g_tol);
  TEST_FLOATING_EQUALITY(point_kp1.phi, computeValue(*phi, alpha_k_init), g_tol);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, noEval )


//
// Check that object will force a minimum number of iterations if asked.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, minIters, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>();
  
  RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>();

  const RCP<ParameterList> pl = parameterList();
  pl->set("Max Backtrack Fraction", 1.0);
  pl->set("Min Num Iterations", 1);
  linesearch->setParameterList(pl);

  linesearch->setOStream(rcpFromRef(out));

  const Scalar alpha_k_init = as<Scalar>(2.1);
  const PointEval1D<Scalar> point_k = computePoint(*phi, ST::zero(), true, true);
  PointEval1D<Scalar> point_kp1 = computePoint(*phi, alpha_k_init);

  int numIters = -1;
  const bool linesearchResult = linesearch->doLineSearch(
    *phi, point_k, inOutArg(point_kp1), outArg(numIters) );

  TEST_ASSERT(linesearchResult);
  TEST_EQUALITY(numIters, 1);
  TEST_FLOATING_EQUALITY(point_kp1.alpha, as<Scalar>(2.0), g_tol);
  TEST_FLOATING_EQUALITY(point_kp1.phi, as<ScalarMag>(3.0), g_tol);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, minIters )


//
// ToDo:
//
// (*) Check that the object performs min bracketing correctly.
//
// (*) Check that the object performs max bracketing correctly.
//
// (*) Check that the ArmijoPolyInterpLineSearch object deals with NaN returns
// currectly.
//
// (*) Check that the ArmijoPolyInterpLineSearch object will throw the right
// exception if a positive initial derivative is presented.
//
// (*) Check that the ArmijoPolyInterpLineSearch object will deal with line search
// failure correctly.
//
// (*) Check that the ArmijoPolyInterpLineSearch object deals with invalid input
// correctly (i.e. check preconditions and throws).
//


} // namespace
