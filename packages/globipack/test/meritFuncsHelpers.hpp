
#include "GlobiPack_TestLagrPolyMeritFunc1D.hpp"
#include "Teuchos_Tuple.hpp"


namespace {


using GlobiPack::TestLagrPolyMeritFunc1D;
using GlobiPack::testLagrPolyMeritFunc1D;
using GlobiPack::PointEval1D;
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::tuple;


template<class Scalar>
inline Scalar sqr(const Scalar &x) { return x*x; }


template<class Scalar>
inline Scalar cube(const Scalar &x) { return x*x*x; }


//
// Set up a quadratic merit function with minimizer at alpha=2.0, phi=3.0.
//

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


} // namespace
