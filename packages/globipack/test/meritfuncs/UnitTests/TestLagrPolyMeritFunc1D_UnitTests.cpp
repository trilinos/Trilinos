
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


