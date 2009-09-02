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
