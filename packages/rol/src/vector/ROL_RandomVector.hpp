// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_RANDOMVECTOR_H
#define ROL_RANDOMVECTOR_H

#include "ROL_Vector.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Elementwise_Function.hpp"


namespace ROL {

/** @ingroup la_group
    \function RandomizeVector
    \brief Fill a ROL::Vector with uniformly-distributed random numbers
           in the interval [lower,upper]

    \f[ x_i \in \mathcal{U}(l,u)\;\forall\l i \f]

*/

template<class Real> 
void RandomizeVector( Vector<Real> &x, const Real &lower=0.0, const Real &upper=1.0 ) {

  Elementwise::UniformlyRandom<Real> ur(lower,upper);
  x.applyUnary(ur);
}

/** @ingroup la_group
    \function RandomizeFeasibleVector
    \brief Fill a ROL::Vector with uniformly-distributed random numbers
           which satisfy the supplied bound constraint

    \f[ x_i \in \mathcal{U}(l_i,u_i)\;\forall\l i \f]
*/

template<class Real>
void RandomizeFeasibleVector( Vector<Real> &x, BoundConstraint<Real> &bnd ) {
  const ROL::Ptr<const Vector<Real> > u = bnd.getUpperBound();
  const ROL::Ptr<const Vector<Real> > l = bnd.getLowerBound();

  Elementwise::UniformlyRandomMultiply<Real> urm;
  
  x.set(*u);
  x.axpy(-1.0,*l);
  x.applyUnary(urm);
  x.plus(*l);
}


} // namespace ROL

#endif // ROL_RANDOMVECTOR_H
