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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
