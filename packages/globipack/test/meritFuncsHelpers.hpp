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
