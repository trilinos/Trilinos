
#include "GlobiPack_Brents1DMinimization.hpp"
#include "GlobiPack_TestLagrPolyMeritFunc1D.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


//
// Helper code and declarations
//


using GlobiPack::TestLagrPolyMeritFunc1D;
using GlobiPack::testLagrPolyMeritFunc1D;
using GlobiPack::Brents1DMinimization;
using GlobiPack::brents1DMinimization;
using GlobiPack::PointEval1D;
using GlobiPack::computePoint;
using Teuchos::as;
using Teuchos::inOutArg;
using Teuchos::outArg;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using Teuchos::ParameterList;
using Teuchos::parameterList;
using Teuchos::OSTab;


template<class Scalar>
inline Scalar sqr(const Scalar &x) { return x*x; }


template<class Scalar>
inline Scalar cube(const Scalar &x) { return x*x*x; }


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

//
// Set up a cubic merit function with minimizer at alpha=2.0, phi=3.0;
//
// The function being represented approximated is:
//
//   phi(alpha) = (alpha - 2.0)^2 + 1e-3 * (alpha - 2.0)^3 + 3.0
//
// This function has the first and second derivatives derivatives:
//
//   Dphi(alpha) = 2.0 * (alpha - 2.0) + 3e-3 * (alpha - 2.0)^2
//
//   D2phi(alpha) = 2.0 + 6e-3 * (alpha - 2.0)
//
// At alpha=2.0, the function has Dphi=0.0 and D2phi = 2.0 and therefore, this
// is a local minimum.
//

const double cubicMut = 1e-3;

template<class Scalar>
inline Scalar cubicPhiVal(const Scalar &alpha)
{ return sqr(alpha - 2.0) + cubicMut * cube(alpha - 2.0) + 3.0; }


template<class Scalar>
const RCP<TestLagrPolyMeritFunc1D<Scalar> > cubicPhi()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  Array<Scalar> alphaPoints =
    tuple<Scalar>(0.0, 1.0, 3.0, 4.0);
  Array<ScalarMag> phiPoints =
    tuple<ScalarMag>(
      cubicPhiVal(alphaPoints[0]),
      cubicPhiVal(alphaPoints[1]),
      cubicPhiVal(alphaPoints[2]),
      cubicPhiVal(alphaPoints[3])
      );
  return testLagrPolyMeritFunc1D<Scalar>(alphaPoints, phiPoints);
}


double g_tol_scale = 100.0;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "tol", &g_tol_scale, "Floating point tolerance scaling of eps." );
}


//
// Unit tests for Brents1DMinimization
//


//
// Check that object can exactly interplate a quadratic merit function.  This
// takes more than one iteration due to the golden search bracketing algorithm
// which takes time to terminate.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Brents1DMinimization, quadExact, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = quadPhi<Scalar>();
  
  RCP<Brents1DMinimization<Scalar> > brentsMin = brents1DMinimization<Scalar>();

  // Set up to do one iteration and that is it.  With the quadratic
  // interplation that should be enough.
  const RCP<ParameterList> pl = parameterList();
  //pl->set("Relative Tol", 1.0);
  //pl->set("Bracket Tol", 1.0);
  //pl->set("Max Iterations", 3);
  brentsMin->setParameterList(pl);

  brentsMin->setOStream(rcpFromRef(out));

  const Array<Array<double> > brackets =
    tuple<Array<double> >(
      tuple(0.0, 2.0, 4.0), // This is the midpoint and the solution!
      tuple(0.5, 2.5, 4.5), // This is the midpoint but *not* the solution!
      tuple(0.0, 1.0, 3.0),
      tuple(1.0, 3.0, 4.0),
      tuple(1.9, 2.0, 4.0),
      tuple(1.9, 3.9, 4.0),
      tuple(0.0, 2.0, 2.1),
      tuple(0.0, 0.1, 2.1)
      );

  for (int i = 0; i < as<int>(brackets.size()); ++i ) {

    const ArrayView<const double> bracket = brackets[i]();

    out << "\ni = "<<i<<": bracket = "<<bracket()<<"\n";

    OSTab tab(out);
    
    PointEval1D<Scalar> p_l = computePoint<Scalar>(*phi, bracket[0]);
    PointEval1D<Scalar> p_m = computePoint<Scalar>(*phi, bracket[1]);
    PointEval1D<Scalar> p_u = computePoint<Scalar>(*phi, bracket[2]);
    int numIters = -1;
    
    const bool mimimized = brentsMin->approxMinimize(
      *phi, p_l, inOutArg(p_m), p_u, outArg(numIters) );
    
    TEST_ASSERT(mimimized);
    //TEST_EQUALITY_CONST(numIters, 1);
    TEST_FLOATING_EQUALITY(p_m.alpha, as<Scalar>(2.0),
      as<Scalar>(g_tol_scale*ST::squareroot(ST::eps())));
    TEST_FLOATING_EQUALITY(p_m.phi, as<Scalar>(3.0),
      as<Scalar>(g_tol_scale)*ST::eps());
    TEST_COMPARE(p_l.alpha, <=, p_m.alpha);
    TEST_COMPARE(p_m.alpha, <=, p_u.alpha);
    TEST_COMPARE(p_m.phi, <=, p_l.phi);
    TEST_COMPARE(p_m.phi, <=, p_u.phi);
    
  }
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( Brents1DMinimization, quadExact )


//
// Check that object can approximately mimimize a cubic function.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Brents1DMinimization, cubicApprox, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  const RCP<TestLagrPolyMeritFunc1D<Scalar> > phi = cubicPhi<Scalar>();
  
  RCP<Brents1DMinimization<Scalar> > brentsMin = brents1DMinimization<Scalar>();

  // Set up to do one iteration and that is it.  With the quadratic
  // interplation that should be enough.
  const RCP<ParameterList> pl = parameterList();
  pl->set("Relative Tol", as<double>(g_tol_scale*ST::eps()));
  pl->set("Bracket Tol", as<double>(ST::eps()));
  brentsMin->setParameterList(pl);

  brentsMin->setOStream(rcpFromRef(out));

  const Array<Array<double> > brackets =
    tuple<Array<double> >(
      tuple(0.0, 2.0, 4.0), // This is the midpoint and the solution!
      tuple(0.5, 2.5, 4.5), // This is the midpoint but *not* the solution!
      tuple(0.0, 1.0, 3.0),
      tuple(1.0, 3.0, 4.0),
      tuple(1.9, 2.0, 4.0),
      tuple(1.9, 3.9, 4.0),
      tuple(0.0, 2.0, 2.1),
      tuple(0.0, 0.1, 2.1)
      );

  for (int i = 0; i < as<int>(brackets.size()); ++i ) {

    const ArrayView<const double> bracket = brackets[i]();

    out << "\ni = "<<i<<": bracket = "<<bracket()<<"\n";

    OSTab tab(out);
    
    PointEval1D<Scalar> p_l = computePoint<Scalar>(*phi, bracket[0]);
    PointEval1D<Scalar> p_m = computePoint<Scalar>(*phi, bracket[1]);
    PointEval1D<Scalar> p_u = computePoint<Scalar>(*phi, bracket[2]);
    int numIters = -1;
    
    const bool mimimized = brentsMin->approxMinimize(
      *phi, p_l, inOutArg(p_m), p_u, outArg(numIters) );
    
    TEST_ASSERT(mimimized);
    //TEST_EQUALITY_CONST(numIters, 1);
    TEST_FLOATING_EQUALITY(p_m.alpha, as<Scalar>(2.0),
      as<Scalar>(g_tol_scale*ST::squareroot(ST::eps())));
    TEST_FLOATING_EQUALITY(p_m.phi, as<Scalar>(3.0),
      as<Scalar>(g_tol_scale*ST::eps()));
    TEST_COMPARE(p_l.alpha, <=, p_m.alpha);
    TEST_COMPARE(p_m.alpha, <=, p_u.alpha);
    TEST_COMPARE(p_m.phi, <=, p_l.phi);
    TEST_COMPARE(p_m.phi, <=, p_u.phi);
    
  }
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES( Brents1DMinimization, cubicApprox )


} // namespace
