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

#ifndef GLOBIPACK_TEST_LAGR_POLY_MERIT_FUNC_1D_DECL_HPP
#define GLOBIPACK_TEST_LAGR_POLY_MERIT_FUNC_1D_DECL_HPP


#include "GlobiPack_MeritFunc1DBase.hpp"


namespace GlobiPack {


/** \brief Lagrange Polynomial Merit Function used in testing.
 
This test class implements an arbitrary order polynomial specified as a set
points.

Let the order-n polynomial approximation be:

\verbatim

  phi(alpha) =
    sum( phi_k * L(n,k)(alpha), k = 0...n-1 )

\endverbatim

where <tt>L(n,k)(alpha) are the nth order Lagrange polynomials:

\verbatim

  L(n,k)(alpha) =
    product( (alpha - alpha[i]) / (alpha[k] - alpha[i]), i=0...n-1, i!=k )

\endverbatim

The derivative of <tt>phi(alpha)</tt> with respect to <tt>alpha</tt>
<tt>Dphi</tt> is given by:

\verbatim

  Dphi(alpha) =
    sum( phi_k * DL(n,k)(alpha), k = 0...n-1 )

\endverbatim

where:

\verbatim

  DL(n,k)(alpha) = sum(
      1/(alpha-alpha[j])
        * product( (alpha-alpha[i])/(alpha[k]-alpha[i]), i=0...n-1, i!=k, i!=j ),
      j=0...n-1, j!=k
      )

\endverbatim

Above, <tt>DL(n,k)(alpha)</tt> is derived using the simple product rule.


*/
template<typename Scalar>
class TestLagrPolyMeritFunc1D : public MeritFunc1DBase<Scalar>
{
public:

  /** \brief Constructor. */
  TestLagrPolyMeritFunc1D(
    const ArrayView<const Scalar> &alpha,
    const ArrayView<const Scalar> &phi
    );

  /** \name Overridden from MeritFunc1DBase */
  //@{

  /** \brief . */
  virtual bool supportsDerivEvals() const;
  /** \brief . */
  virtual void eval( const Scalar &alpha, const Ptr<Scalar> &phi,
    const Ptr<Scalar> &Dphi ) const;

  //@}

private:

  Array<Scalar> alpha_;
  Array<Scalar> phi_;

};


template<typename Scalar>
const RCP<TestLagrPolyMeritFunc1D<Scalar> >
testLagrPolyMeritFunc1D(
  const ArrayView<const Scalar> &alpha,
  const ArrayView<const Scalar> &phi
  )
{
  return Teuchos::rcp(new TestLagrPolyMeritFunc1D<Scalar>(alpha, phi));
}


} // namespace GlobiPack


#endif // GLOBIPACK_TEST_LAGR_POLY_MERIT_FUNC_1D_DECL_HPP
