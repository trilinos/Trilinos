
#include "Teuchos_UnitTestHarness.hpp"

#include "GlobiPack_TestLagrPolyMeritFunc1D.hpp"
#include "GlobiPack_ArmijoPolyInterpLineSearch.hpp"
#include "Teuchos_Tuple.hpp"


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
  ECHO(RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>());
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
  ECHO(pl->set("Min Number of Iterations", minIters));
  ECHO(pl->set("Max Number of Iterations", maxIters));
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

  ECHO(const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>());
  
  ECHO(RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>());

  ECHO(linesearch->setOStream(rcpFromRef(out)));

  ECHO(const Scalar phi_k = computeValue<Scalar>(*phi, ST::zero()));
  ECHO(Scalar alpha_k = as<Scalar>(5.0)); // Initial point should not be accepted!
  ECHO(Scalar phi_kp1 = computeValue<Scalar>(*phi, alpha_k));

  ECHO(const bool linesearchResult =
    linesearch->doLineSearch(*phi, phi_k, outArg(alpha_k), outArg(phi_kp1), null));

  TEST_ASSERT(linesearchResult);
  TEST_FLOATING_EQUALITY(alpha_k, as<Scalar>(2.0), g_tol);
  TEST_FLOATING_EQUALITY(phi_kp1, as<ScalarMag>(3.0), g_tol);
  TEST_EQUALITY(linesearch->numIterations(), 1);

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

  ECHO(const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>());
  
  ECHO(RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>());

  ECHO(linesearch->setOStream(rcpFromRef(out)));

  ECHO(const Scalar phi_k = computeValue<Scalar>(*phi, ST::zero()));
  ECHO(const Scalar alpha_k_init = as<Scalar>(2.1));
  ECHO(Scalar alpha_k = alpha_k_init); // Initial point should be accepted!
  ECHO(Scalar phi_kp1 = computeValue<Scalar>(*phi, alpha_k));

  ECHO(const bool linesearchResult =
    linesearch->doLineSearch(*phi, phi_k, outArg(alpha_k), outArg(phi_kp1), null));

  TEST_ASSERT(linesearchResult);
  TEST_FLOATING_EQUALITY(alpha_k, alpha_k_init, g_tol);
  TEST_FLOATING_EQUALITY(phi_kp1, computeValue<Scalar>(*phi, alpha_k_init), g_tol);
  TEST_EQUALITY(linesearch->numIterations(), 0);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, noEval )


//
// Check that object will force a minimum number of iterations if asked.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, minIters, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  ECHO(const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>());
  
  ECHO(RCP<ArmijoPolyInterpLineSearch<Scalar> > linesearch =
    armijoQuadraticLineSearch<Scalar>());

  ECHO(const RCP<ParameterList> pl = parameterList());
  ECHO(pl->set("Max Backtrack Fraction", 1.0));
  ECHO(pl->set("Min Number of Iterations", 1));
  ECHO(linesearch->setParameterList(pl));

  ECHO(linesearch->setOStream(rcpFromRef(out)));

  ECHO(const Scalar phi_k = computeValue<Scalar>(*phi, ST::zero()));
  ECHO(const Scalar alpha_k_init = as<Scalar>(2.1));
  ECHO(Scalar alpha_k = alpha_k_init); // Initial point should be accepted!
  ECHO(Scalar phi_kp1 = computeValue<Scalar>(*phi, alpha_k));

  ECHO(const bool linesearchResult =
    linesearch->doLineSearch(*phi, phi_k, outArg(alpha_k), outArg(phi_kp1), null));

  TEST_ASSERT(linesearchResult);
  TEST_FLOATING_EQUALITY(alpha_k, as<Scalar>(2.0), g_tol);
  TEST_FLOATING_EQUALITY(phi_kp1, as<ScalarMag>(3.0), g_tol);
  TEST_EQUALITY(linesearch->numIterations(), 1);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, minIters )


//
// Check that the object performs min bracketing correctly.
//

/*

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, minBracketing, Scalar )
{
  TEST_FOR_EXCEPT(true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, minBracketing )

*/


//
// Check that the object performs max bracketing correctly.
//

/*

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, maxBracketing, Scalar )
{
  TEST_FOR_EXCEPT(true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, maxBracketing )

*/


//
// Check that the ArmijoPolyInterpLineSearch object deals with NaN returns
// currectly.
//

/*

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, nanBehavior, Scalar )
{
  TEST_FOR_EXCEPT(true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, nanBehavior )

*/

//
// Check that the ArmijoPolyInterpLineSearch object will throw the right
// exception if a negative initial derivative is presented.
//

/*

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, noDescentDirec, Scalar )
{
  TEST_FOR_EXCEPT(true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, noDescentDirec )

*/


//
// Check that the ArmijoPolyInterpLineSearch object will deal with line search
// failure correctly.
//

/*

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, linesearchFailure, Scalar )
{
  TEST_FOR_EXCEPT(true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, linesearchFailure )

*/


//
// Check that the ArmijoPolyInterpLineSearch object deals with invalid input
// correctly (i.e. check preconditions and throws).
//

/*

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArmijoPolyInterpLineSearch, invalidArgs, Scalar )
{
  TEST_FOR_EXCEPT(true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( ArmijoPolyInterpLineSearch, invalidArgs )

*/


} // namespace
